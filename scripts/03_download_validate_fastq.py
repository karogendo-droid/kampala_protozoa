#!/usr/bin/env python3
"""Download one audited ENA paired FASTQ run; verify bytes, MD5 and all read pairs.

Run from the project, with this file in scripts/:
    python scripts/03_download_validate_fastq.py --audit metadata/PRJNA1091428_audit_TIMESTAMP --run SRR28434174

Requires Python >=3.10 and requests. No classifier, filtering or upload is run.
Downloads use HTTPS with certificate verification and are resumable when supported.
Input audit records are never edited. Existing final FASTQs are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time
from datetime import datetime, timezone
from urllib.parse import urlsplit

import requests

PROJECT = "PRJNA1091428"
VERSION = "1.0.0"
CHUNK = 1024 * 1024
ALLOWED_HOSTS = {"ftp.sra.ebi.ac.uk", "ftp.ebi.ac.uk"}
BASES = frozenset(b"ACGTRYSWKMBDHVNacgtryswkmbdhvn.")


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header in {path.name}")
        return [{str(k): (v or "").strip() for k, v in row.items() if k is not None}
                for row in reader]


def checksum(path: Path, algorithm: str) -> str:
    # MD5 is used to compare with ENA's transfer checksum, not for security.
    digest = hashlib.md5(usedforsecurity=False) if algorithm == "md5" else hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def load_audit(audit: Path, run: str) -> tuple[dict, list[dict]]:
    summary = (audit / "audit_summary.txt").read_text(encoding="utf-8")
    if f"BioProject: {PROJECT}" not in summary:
        raise ValueError("The audit does not identify the expected BioProject.")
    samples = [r for r in read_tsv(audit / "sample_manifest.tsv")
               if r.get("run_accession") == run]
    if len(samples) != 1:
        raise ValueError("Expected exactly one sample-manifest row for the selected run.")
    sample = samples[0]
    for key, expected in {"library_layout": "PAIRED", "instrument_platform": "ILLUMINA",
                          "library_strategy": "WGS", "library_source": "METAGENOMIC"}.items():
        if sample.get(key, "").upper() != expected:
            raise ValueError(f"Review required: {key} is not {expected} for {run}.")
    files = [r for r in read_tsv(audit / "fastq_file_manifest.tsv")
             if r.get("run_accession") == run]
    expected_names = {f"{run}_1.fastq.gz", f"{run}_2.fastq.gz"}
    if len(files) != 2 or {r.get("filename") for r in files} != expected_names:
        raise ValueError("Expected exactly the two paired FASTQ files. Extra/missing files require review.")
    for row in files:
        url = urlsplit(row.get("https_url", ""))
        if (url.scheme != "https" or url.hostname not in ALLOWED_HOSTS or url.port not in (None, 443)
                or url.username or url.password or url.query or url.fragment
                or Path(url.path).name != row["filename"]):
            raise ValueError("Unexpected URL in FASTQ manifest; no download attempted.")
        if not re.fullmatch(r"[0-9a-fA-F]{32}", row.get("expected_md5", "")):
            raise ValueError("Missing/invalid archive MD5 checksum.")
        if not row.get("expected_bytes", "").isdigit() or int(row["expected_bytes"]) <= 0:
            raise ValueError("Missing/invalid expected file size.")
        row["expected_md5"] = row["expected_md5"].lower()
    return sample, sorted(files, key=lambda x: x["filename"])


def verify_file(path: Path, row: dict) -> dict:
    actual_bytes = path.stat().st_size
    if actual_bytes != int(row["expected_bytes"]):
        raise ValueError(f"Size mismatch for {path.name}: {actual_bytes} bytes; expected {row['expected_bytes']}.")
    actual_md5 = checksum(path, "md5")
    if actual_md5 != row["expected_md5"]:
        raise ValueError(f"MD5 mismatch for {path.name}. Do not use this file for analysis.")
    return {**row, "actual_bytes": actual_bytes, "actual_md5": actual_md5,
            "byte_count_matches": True, "md5_matches": True, "checked_utc": utc()}


def download_file(session: requests.Session, row: dict, folder: Path, emit=print) -> tuple[Path, dict]:
    final = folder / row["filename"]
    part = folder / (row["filename"] + ".part")
    expected = int(row["expected_bytes"])
    if final.is_symlink() or part.is_symlink():
        raise ValueError("Refusing to follow a symbolic link in the download directory.")
    if final.exists():
        emit(f"Existing file: {final.name}; verifying without overwriting.")
        result = verify_file(final, row)
        final.chmod(0o444)
        return final, {**result, "transfer": "reused_existing_verified_file"}
    if part.exists() and part.stat().st_size > expected:
        raise ValueError(f"Partial file {part.name} is larger than expected; manual review required.")

    for attempt in range(1, 4):
        offset = part.stat().st_size if part.exists() else 0
        if offset == expected:
            break
        headers = {"Accept-Encoding": "identity"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
            emit(f"Resuming {final.name} from {offset:,} bytes (attempt {attempt}/3).")
        else:
            emit(f"Downloading {final.name} ({expected / 1e6:.1f} MB; attempt {attempt}/3).")
        try:
            with session.get(row["https_url"], headers=headers, stream=True,
                             timeout=(30, 120), allow_redirects=True) as response:
                response.raise_for_status()
                dest = urlsplit(response.url)
                if dest.scheme != "https" or dest.hostname not in ALLOWED_HOSTS:
                    raise ValueError("Unexpected download redirect; stopped without accepting the file.")
                if response.headers.get("Content-Encoding", "identity").lower() not in ("", "identity"):
                    raise ValueError("Unexpected HTTP compression; cannot verify exact archived bytes.")
                if response.status_code == 206:
                    match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", response.headers.get("Content-Range", ""))
                    if not match or int(match[1]) != offset or int(match[3]) != expected:
                        raise ValueError("Unexpected partial-response range; stopped to protect file integrity.")
                    mode = "ab" if offset else "wb"
                elif response.status_code == 200:
                    if offset:
                        emit("Server did not resume; restarting this partial file from byte zero.")
                    offset, mode = 0, "wb"
                else:
                    raise ValueError(f"Unexpected HTTP response: {response.status_code}")
                length = response.headers.get("Content-Length")
                if length and int(length) != expected - offset:
                    raise ValueError("Server length differs from the archived file size; review required.")
                count, last_mark = offset, int(100 * offset / expected) // 10
                with part.open(mode) as handle:
                    for block in response.iter_content(chunk_size=CHUNK):
                        if not block:
                            continue
                        if count + len(block) > expected:
                            raise ValueError("Transfer exceeded the expected file size.")
                        handle.write(block)
                        count += len(block)
                        mark = int(100 * count / expected) // 10
                        if mark > last_mark:
                            emit(f"  {final.name}: {100 * count / expected:.0f}%")
                            last_mark = mark
            if part.stat().st_size != expected:
                raise requests.ConnectionError("Transfer ended before all expected bytes arrived.")
            break
        except requests.RequestException as exc:
            if attempt == 3:
                raise RuntimeError("Download interrupted. Partial file retained. Rerun the SAME command to retry.") from exc
            emit(f"Network attempt failed ({type(exc).__name__}); retrying with any partial file retained.")
            time.sleep(2 * attempt)

    emit(f"Checking file size and MD5: {final.name}")
    try:
        result = verify_file(part, row)
    except ValueError:
        rejected = folder / f"{part.name}.rejected_{stamp()}"
        part.rename(rejected)
        raise ValueError(f"Downloaded bytes failed verification; retained as {rejected.name}. "
                         "No final FASTQ was created. Rerun the command to try a fresh copy.")
    if final.exists():
        raise ValueError("A final FASTQ appeared during this download. No file was overwritten.")
    part.rename(final)
    final.chmod(0o444)
    emit(f"PASS: {final.name} matches the archive checksum and is now read-only.")
    return final, {**result, "transfer": "downloaded_and_verified"}


def read_record(handle, file_label: str, number: int):
    header = handle.readline()
    if not header:
        return None
    seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
    if not seq or not plus or not qual:
        raise ValueError(f"Truncated FASTQ record {number:,} in {file_label}.")
    header, seq, plus, qual = [x.rstrip(b"\r\n") for x in (header, seq, plus, qual)]
    if not header.startswith(b"@") or len(header) < 2 or not plus.startswith(b"+"):
        raise ValueError(f"Invalid FASTQ structure at record {number:,} in {file_label}.")
    if not seq or len(seq) != len(qual):
        raise ValueError(f"Sequence/quality length mismatch at record {number:,} in {file_label}.")
    if not set(seq).issubset(BASES) or min(qual) < 33 or max(qual) > 126:
        raise ValueError(f"Unexpected sequence/quality characters at record {number:,} in {file_label}.")
    identifier = header.split(None, 1)[0][1:]
    if identifier.endswith((b"/1", b"/2")):
        identifier = identifier[:-2]
    if not identifier:
        raise ValueError(f"Missing read identifier at record {number:,} in {file_label}.")
    return identifier, len(seq)


def validate_pairs(r1: Path, r2: Path, emit=print) -> dict:
    emit("Checking gzip integrity, four-line FASTQ structure and matching identifiers for EVERY read pair ...")
    count = 0
    totals, mins, maxs = [0, 0], [None, None], [0, 0]
    with gzip.open(r1, "rb") as first, gzip.open(r2, "rb") as second:
        while True:
            left = read_record(first, "R1", count + 1)
            right = read_record(second, "R2", count + 1)
            if left is None and right is None:
                break
            if left is None or right is None:
                raise ValueError("R1 and R2 contain unequal numbers of records; no pairing changes were made.")
            if left[0] != right[0]:
                raise ValueError(f"Read identifiers disagree at pair {count + 1:,}. "
                                 "Do not rename or repair reads manually; request review.")
            count += 1
            for i, item in enumerate((left, right)):
                n = item[1]
                totals[i] += n
                mins[i] = n if mins[i] is None else min(mins[i], n)
                maxs[i] = max(maxs[i], n)
            if count % 1_000_000 == 0:
                emit(f"  {count:,} read pairs checked.")
    if count == 0:
        raise ValueError("The FASTQ files contain no read pairs.")
    return {"gzip_integrity": "PASS", "fastq_structure": "PASS", "paired_identifiers": "PASS",
            "read_pairs": count, "total_reads": 2 * count,
            "R1": {"reads": count, "bases": totals[0], "min_length": mins[0], "max_length": maxs[0],
                   "mean_length": totals[0] / count},
            "R2": {"reads": count, "bases": totals[1], "min_length": mins[1], "max_length": maxs[1],
                   "mean_length": totals[1] / count}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path, help="Existing metadata audit folder")
    parser.add_argument("--run", default="SRR28434174", help="ONE run to download")
    args = parser.parse_args()
    if not re.fullmatch(r"[SED]RR\d+", args.run):
        parser.error("Invalid run accession.")
    script = Path(__file__).resolve()
    root = script.parent.parent
    if script.parent.name != "scripts" or not (root / "metadata").is_dir():
        print("STOP: Put this file in the existing project's scripts folder.")
        return 1
    audit = args.audit.expanduser()
    audit = (audit if audit.is_absolute() else root / audit).resolve()
    record = {"project": PROJECT, "script_version": VERSION, "script_sha256": checksum(script, "sha256"),
              "started_utc": utc(), "run_accession": args.run, "python": sys.version,
              "requests_version": requests.__version__, "status": "STARTED", "files": [],
              "classification_performed": False, "filtering_performed": False, "uploads_performed": False}
    output = None
    log = None
    session = requests.Session()
    session.headers.update({"User-Agent": f"KampalaProtozoaDownload/{VERSION}"})
    def emit(message: str):
        print(message, flush=True)
        if log is not None:
            log.write(message + "\n")
            log.flush()
    try:
        sample, file_rows = load_audit(audit, args.run)
        output = root / "metadata" / f"download_{args.run}_{stamp()}"
        output.mkdir(parents=True, exist_ok=False)
        log = (output / "download_log.txt").open("w", encoding="utf-8")
        record.update(sample_metadata=sample, source_audit=str(audit.relative_to(root)) if audit.is_relative_to(root) else str(audit))
        for name in ("fastq_file_manifest.tsv", "sample_manifest.tsv", "audit_summary.txt"):
            shutil.copy2(audit / name, output / name)
        record["input_sha256"] = {name: checksum(output / name, "sha256") for name in
                                  ("fastq_file_manifest.tsv", "sample_manifest.tsv", "audit_summary.txt")}
        raw = root / "data" / "raw" / args.run
        raw.mkdir(parents=True, exist_ok=True)
        expected_total = sum(int(r["expected_bytes"]) for r in file_rows)
        free = shutil.disk_usage(raw).free
        # A conservative buffer for THIS transfer only, not an analysis-space estimate.
        if free < expected_total + 2 * 1024**3:
            raise ValueError("Less than the compressed download size plus a 2-GiB safety buffer is free.")
        emit(f"DOWNLOADING ONE RUN: {args.run} | {sample.get('sample_alias', 'NOT_REPORTED')}")
        emit(f"Collection date as deposited: {sample.get('collection_date_as_reported', 'NOT_REPORTED')}")
        emit(f"Two compressed FASTQ files: {expected_total / 1e9:.3f} GB total.")
        emit("No filtering, taxonomic classification or uploads will be performed.")
        paths = []
        for row in file_rows:
            path, result = download_file(session, row, raw, emit)
            paths.append(path)
            record["files"].append({**result, "local_path": str(path.relative_to(root)), "permissions": oct(path.stat().st_mode & 0o777)})
        metrics = validate_pairs(paths[0], paths[1], emit)
        record.update(status="DOWNLOAD_AND_STRUCTURE_CHECKS_PASSED", validation=metrics)
        lines = ["DOWNLOAD AND VALIDATION COMPLETE", f"BioProject: {PROJECT}",
                 f"Run: {args.run}", f"Sample alias: {sample.get('sample_alias', 'NOT_REPORTED')}",
                 f"Collection date as deposited: {sample.get('collection_date_as_reported', 'NOT_REPORTED')}",
                 "Files verified against expected byte counts and MD5 checksums: 2/2",
                 "Gzip integrity: PASS", "FASTQ structure: PASS", "Paired read identifiers: PASS",
                 f"Read pairs: {metrics['read_pairs']:,}", f"Total reads (R1 + R2): {metrics['total_reads']:,}"]
        for mate in ("R1", "R2"):
            s = metrics[mate]
            lines.append(f"{mate} length (bases): min={s['min_length']}; max={s['max_length']}; mean={s['mean_length']:.2f}")
        lines.extend([f"Raw archive files (read-only): {raw.relative_to(root)}",
                      f"Reports: {output.relative_to(root)}",
                      "", "These are transfer/structure checks, NOT a full quality-control or pathogen result.",
                      "No reads were trimmed, classified or uploaded. No original audit files were edited.",
                      "Unresolved metadata and control availability still require review."])
        summary = "\n".join(lines) + "\n"
        (output / "download_report.txt").write_text(summary, encoding="utf-8")
        emit("\n" + summary)
        return 0
    except KeyboardInterrupt:
        record.update(status="INTERRUPTED", error="Stopped by user; rerun the same command to resume/reverify.")
        emit("\nSTOP: Interrupted. Partial files retained. Rerun the SAME command to resume/reverify.")
        return 130
    except Exception as exc:
        record.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        emit(f"\nSTOP: {type(exc).__name__}: {exc}")
        emit("Do not proceed to analysis. Share the error text or small log, NOT sequencing reads.")
        return 1
    finally:
        session.close()
        record["finished_utc"] = utc()
        if output is not None:
            (output / "download_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
            if record["status"] not in ("DOWNLOAD_AND_STRUCTURE_CHECKS_PASSED",):
                (output / "download_report.txt").write_text(
                    f"Status: {record['status']}\nRun: {args.run}\nError: {record.get('error', '')}\n", encoding="utf-8")
            if log is not None:
                log.close()


if __name__ == "__main__":
    raise SystemExit(main())
