#!/usr/bin/env python3
"""Audit stored quality scores and ambiguous bases without modifying FASTQ files.

Run from the existing Kampala project, using the successful step-04 QC folder:
  python scripts/05_audit_quality_encoding.py --qc results/qc/raw/RUN_TIMESTAMP

Python >=3.10; standard library only. No external commands or network requests.
Reads both complete gzip FASTQ files, verifies archived MD5 checksums during the
scan, and writes aggregate reports only. No sequence or read identifiers are
included in the shareable review bundle. This is NOT a trimming command.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import sys
import zipfile

PROJECT = "PRJNA1091428"
VERSION = "1.0.0"
QUALITY_OFFSET = 33
MIN_LENGTH = 75
MAX_N = 5
BUFFER_LIMIT = 2 * 1024 * 1024


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inside(root: Path, path: Path) -> Path:
    p = (path if path.is_absolute() else root / path).resolve()
    if not p.is_relative_to(root):
        raise ValueError("Input paths must remain within the project folder.")
    return p


class HashingReader:
    """Hash every compressed byte as gzip consumes the complete input stream."""
    def __init__(self, handle):
        self.handle = handle
        self.md5 = hashlib.md5(usedforsecurity=False)
        self.count = 0
        self.name = str(handle.name)

    def read(self, size: int = -1) -> bytes:
        data = self.handle.read(size)
        self.md5.update(data)
        self.count += len(data)
        return data

    def seekable(self) -> bool:
        return False


def read_record(handle, mate: str, number: int):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise ValueError(f"Truncated FASTQ record in {mate}, record {number:,}.")
    header, sequence, plus, quality = [x.rstrip(b"\r\n") for x in (header, sequence, plus, quality)]
    if not header.startswith(b"@") or not plus.startswith(b"+"):
        raise ValueError(f"Invalid FASTQ structure in {mate}, record {number:,}.")
    if not sequence or len(sequence) != len(quality):
        raise ValueError(f"Sequence/quality length discrepancy in {mate}, record {number:,}.")
    tokens = header[1:].split()
    if not tokens:
        raise ValueError(f"Missing read identifier in {mate}, record {number:,}.")
    name = tokens[0]
    if name.endswith((b"/1", b"/2")):
        wanted = b"/1" if mate == "R1" else b"/2"
        if not name.endswith(wanted):
            raise ValueError(f"Unexpected mate suffix in {mate}, record {number:,}.")
        name = name[:-2]
    return name, sequence.upper(), quality


class Counts:
    def __init__(self):
        self.reads = 0
        self.bases = 0
        self.lengths = Counter()
        self.n_counts = Counter()
        self.quality_counts = Counter()
        self.quality_buffer = bytearray()
        self.n_bases = 0
        self.all_n = 0
        self.non_acgtn = 0
        self.short = 0
        self.many_n = 0
        self.length_n_pass = 0
        self.uniform_reads = 0
        self.uniform_q30_reads = 0

    def update(self, sequence: bytes, quality: bytes) -> tuple[bool, bool]:
        size = len(sequence)
        ns = sequence.count(b"N")
        self.reads += 1
        self.bases += size
        self.lengths[size] += 1
        self.n_counts[ns] += 1
        self.n_bases += ns
        self.all_n += ns == size
        self.non_acgtn += bool(sequence.translate(None, b"ACGTN"))
        self.short += size < MIN_LENGTH
        self.many_n += ns > MAX_N
        ok = size >= MIN_LENGTH and ns <= MAX_N
        self.length_n_pass += ok
        uniform = quality.count(quality[:1]) == len(quality)
        self.uniform_reads += uniform
        self.uniform_q30_reads += uniform and quality[0] == 63
        self.quality_buffer.extend(quality)
        if len(self.quality_buffer) >= BUFFER_LIMIT:
            self.flush()
        return ok, ns == size

    def flush(self):
        self.quality_counts.update(self.quality_buffer)
        self.quality_buffer.clear()
        if any(c < 33 or c > 126 for c in self.quality_counts):
            raise ValueError("A quality byte falls outside printable Phred+33. Review the input encoding.")

    def state(self) -> str:
        if set(self.quality_counts) == {63}:
            return "CONSTANT_Q30_ALL_BASES"
        if len(self.quality_counts) == 1:
            return "ONE_CONSTANT_SCORE_ALL_BASES"
        if self.uniform_reads == self.reads:
            return "CONSTANT_WITHIN_EACH_READ"
        return "WITHIN_READ_VARIATION_PRESENT_OR_MIXED"

    def summary(self) -> dict:
        return {
            "reads": self.reads, "bases": self.bases,
            "minimum_length": min(self.lengths), "maximum_length": max(self.lengths),
            "mean_length": self.bases / self.reads,
            "stored_phred33_scores": sorted(c - QUALITY_OFFSET for c in self.quality_counts),
            "quality_state": self.state(),
            "reads_uniform_within_read": self.uniform_reads,
            "reads_all_quality_Q30": self.uniform_q30_reads,
            "N_bases": self.n_bases, "reads_entirely_N": self.all_n,
            "reads_with_other_than_ACGTN": self.non_acgtn,
            "reads_shorter_than_75": self.short, "reads_with_more_than_5_N": self.many_n,
            "reads_passing_length_and_N_only": self.length_n_pass,
        }


def parse_fastqc(path: Path, filename: str, count: int) -> list[dict]:
    with zipfile.ZipFile(path) as z:
        members = [n for n in z.namelist() if n.endswith("/fastqc_data.txt")]
        if len(members) != 1:
            raise ValueError(f"Cannot identify FastQC data table in {path.name}.")
        text = z.read(members[0]).decode("utf-8")
    modules = {}
    current = None
    for line in text.splitlines():
        if line == ">>END_MODULE":
            current = None
        elif line.startswith(">>"):
            current = line[2:].split("\t")[0]
            modules[current] = []
        elif current is not None:
            modules[current].append(line)
    basic = dict(line.split("\t", 1) for line in modules.get("Basic Statistics", [])
                 if "\t" in line and not line.startswith("#"))
    if basic.get("Filename") != filename or int(basic.get("Total Sequences", "-1")) != count:
        raise ValueError("FastQC report does not describe the selected input and read count.")
    if basic.get("Encoding") != "Sanger / Illumina 1.9":
        raise ValueError("This audit requires review of an encoding other than the reported Phred+33.")
    if "Overrepresented sequences" not in modules:
        raise ValueError("The FastQC overrepresented-sequence module is missing; review required.")
    rows = []
    for line in modules["Overrepresented sequences"]:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError("Unexpected overrepresented-sequence table format.")
        seq = fields[0].upper()
        if not seq:
            raise ValueError("Empty reported sequence in FastQC table.")
        rows.append({
            "row": len(rows) + 1, "reported_count": int(fields[1]),
            "reported_percent": float(fields[2]), "reported_segment_length": len(seq),
            "N_bases_in_segment": seq.count("N"),
            "segment_all_N": "YES" if set(seq) == {"N"} else "NO",
            "segment_is_single_repeated_base": "YES" if len(set(seq)) == 1 else "NO",
            "segment_non_ACGTN_bases": sum(b not in "ACGTN" for b in seq),
            "possible_source_as_reported": fields[3] if len(fields) > 3 else "NOT_REPORTED",
        })
    return rows


def tsv(path: Path, header: list[str], rows):
    with path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc", type=Path, required=True, help="Successful step-04 QC directory")
    parser.add_argument("--run", default="SRR28434174")
    args = parser.parse_args()
    if not re.fullmatch(r"[SED]RR\d+", args.run):
        parser.error("Invalid run accession.")
    script = Path(__file__).resolve()
    root = script.parent.parent
    if script.parent.name != "scripts" or not (root / "metadata").is_dir():
        print("STOP: Put the script in your existing project's scripts folder.")
        return 1
    output = None
    log = None
    record = {"project": PROJECT, "run_accession": args.run, "script_version": VERSION,
              "script_sha256": sha256(script), "started_utc": utc(), "status": "STARTED",
              "python": sys.version, "architecture": platform.machine(),
              "conda_environment": os.environ.get("CONDA_DEFAULT_ENV", "NOT_REPORTED"),
              "quality_offset": QUALITY_OFFSET, "filtering_performed": False,
              "classification_performed": False, "network_requests_performed": False}

    def emit(message):
        print(message, flush=True)
        if log is not None:
            log.write(message + "\n")
            log.flush()

    try:
        qc = inside(root, args.qc.expanduser())
        qpath = qc / "qc_record.json"
        qrecord = json.loads(qpath.read_text(encoding="utf-8"))
        if (qrecord.get("project") != PROJECT or qrecord.get("run_accession") != args.run
                or qrecord.get("status") != "QC_REPORTS_GENERATED_INTERPRETATION_PENDING"):
            raise ValueError("The specified step-04 record is not a successful QC run for this accession.")
        source = qc / "provenance" / "download_record.json"
        if sha256(source) != qrecord.get("source_download_sha256"):
            raise ValueError("The preserved download record no longer matches the QC provenance.")
        prior = json.loads(source.read_text(encoding="utf-8"))
        if (prior.get("project") != PROJECT or prior.get("run_accession") != args.run
                or prior.get("status") != "DOWNLOAD_AND_STRUCTURE_CHECKS_PASSED"):
            raise ValueError("The preserved download record is not valid for this run.")
        validation = prior.get("validation", {})
        if any(validation.get(k) != "PASS" for k in ("gzip_integrity", "fastq_structure", "paired_identifiers")):
            raise ValueError("The original download did not pass all structural checks.")
        expected = validation.get("read_pairs")
        if not isinstance(expected, int) or expected < 1:
            raise ValueError("The original paired-read count is unavailable.")
        rows = prior.get("files", [])
        names = {f"{args.run}_1.fastq.gz", f"{args.run}_2.fastq.gz"}
        if len(rows) != 2 or {r.get("filename") for r in rows} != names:
            raise ValueError("Expected exactly two previously verified FASTQ inputs.")
        rows = sorted(rows, key=lambda r: r["filename"])
        paths = []
        stats = []
        overrepresented = {}
        for mate, row in zip(("R1", "R2"), rows):
            path = root / row.get("local_path", "")
            expected_path = root / "data" / "raw" / args.run / row["filename"]
            if path != expected_path or path.is_symlink() or not path.is_file():
                raise ValueError(f"Unexpected or missing raw input for {mate}.")
            path = inside(root, path)
            if not re.fullmatch(r"[0-9a-fA-F]{32}", str(row.get("expected_md5", ""))):
                raise ValueError("The archive MD5 is missing or invalid.")
            st = path.stat()
            if st.st_size != int(row["expected_bytes"]):
                raise ValueError(f"Input size differs from the archived size for {mate}.")
            paths.append(path)
            stats.append((st.st_size, st.st_mtime_ns, st.st_ino))
            fqzip = qc / "fastqc" / (row["filename"].removesuffix(".fastq.gz") + "_fastqc.zip")
            overrepresented[mate] = parse_fastqc(fqzip, row["filename"], expected)
        output = root / "results" / "qc" / "quality_audit" / f"{args.run}_{stamp()}"
        output.mkdir(parents=True, exist_ok=False)
        log = (output / "quality_audit_log.txt").open("w", encoding="utf-8")
        record.update(source_qc=str(qc.relative_to(root)), source_qc_sha256=sha256(qpath),
                      source_download_sha256=sha256(source), output_directory=str(output.relative_to(root)))
        emit(f"STORED-QUALITY AND AMBIGUOUS-BASE AUDIT: {args.run}")
        emit("Scanning both complete FASTQ files. No trimming, filtering, classification or uploads.")
        emit("Quality values are decoded as Phred+33; their original provenance remains unverified.")
        counts = [Counts(), Counts()]
        pair_states = Counter({"both_pass": 0, "R1_only_pass": 0, "R2_only_pass": 0, "neither_pass": 0})
        n_pairs = Counter({"both_all_N": 0, "R1_only_all_N": 0, "R2_only_all_N": 0, "neither_all_N": 0})
        i = 0
        with paths[0].open("rb") as a, paths[1].open("rb") as b:
            hashed = [HashingReader(a), HashingReader(b)]
            with gzip.GzipFile(fileobj=hashed[0], mode="rb") as r1, gzip.GzipFile(fileobj=hashed[1], mode="rb") as r2:
                while True:
                    x = read_record(r1, "R1", i + 1)
                    y = read_record(r2, "R2", i + 1)
                    if x is None and y is None:
                        break
                    if x is None or y is None:
                        raise ValueError("The FASTQ files contain unequal record counts.")
                    if x[0] != y[0]:
                        raise ValueError(f"Read-pair identifiers differ at pair {i + 1:,}.")
                    i += 1
                    ok1, n1 = counts[0].update(x[1], x[2])
                    ok2, n2 = counts[1].update(y[1], y[2])
                    pair_states["both_pass" if ok1 and ok2 else "R1_only_pass" if ok1 else "R2_only_pass" if ok2 else "neither_pass"] += 1
                    n_pairs["both_all_N" if n1 and n2 else "R1_only_all_N" if n1 else "R2_only_all_N" if n2 else "neither_all_N"] += 1
                    if i % 1_000_000 == 0:
                        emit(f"  Examined {i:,} / {expected:,} read pairs.")
                    if i > expected:
                        raise ValueError("Observed more read pairs than in the validated download.")
        if i != expected:
            raise ValueError("Observed pair count differs from the validated download.")
        record["input_checks"] = []
        for mate, path, row, initial, h, c in zip(("R1", "R2"), paths, rows, stats, hashed, counts):
            c.flush()
            if sum(c.quality_counts.values()) != c.bases:
                raise ValueError("Quality histogram does not account for every base.")
            st = path.stat()
            if (st.st_size, st.st_mtime_ns, st.st_ino) != initial:
                raise ValueError(f"Input changed during audit: {mate}.")
            observed = h.md5.hexdigest()
            if h.count != int(row["expected_bytes"]) or observed != row["expected_md5"].lower():
                raise ValueError(f"The fully scanned {mate} file does not match its archive checksum.")
            record["input_checks"].append({"mate": mate, "filename": row["filename"],
                                           "compressed_bytes": h.count, "archive_md5_matches": True,
                                           "md5": observed, "file_metadata_unchanged": True})
        summaries = {m: c.summary() for m, c in zip(("R1", "R2"), counts)}
        record.update(status="COMPLETE_REVIEW_REQUIRED", mates=summaries,
                      pair_length_N_forecast=dict(pair_states), all_N_pair_counts=dict(n_pairs),
                      overrepresented_segment_composition=overrepresented,
                      forecast_definition={"minimum_length": MIN_LENGTH, "maximum_N": MAX_N,
                      "quality_adapter_and_other_filters_included": False})
        tsv(output / "quality_histogram.tsv", ["mate", "ascii_code", "character", "phred33", "bases", "percent_of_bases"],
            ([m, q, chr(q), q - QUALITY_OFFSET, n, f"{100*n/c.bases:.9f}"]
             for m, c in zip(("R1", "R2"), counts) for q, n in sorted(c.quality_counts.items())))
        tsv(output / "length_histogram.tsv", ["mate", "length_bases", "reads"],
            ([m, k, v] for m, c in zip(("R1", "R2"), counts) for k, v in sorted(c.lengths.items())))
        tsv(output / "N_count_histogram.tsv", ["mate", "N_bases_per_read", "reads"],
            ([m, k, v] for m, c in zip(("R1", "R2"), counts) for k, v in sorted(c.n_counts.items())))
        cols = ["row", "reported_count", "reported_percent", "reported_segment_length", "N_bases_in_segment",
                "segment_all_N", "segment_is_single_repeated_base", "segment_non_ACGTN_bases", "possible_source_as_reported"]
        tsv(output / "overrepresented_composition_no_sequences.tsv", ["mate"] + cols,
            ([m] + [r[k] for k in cols] for m in ("R1", "R2") for r in overrepresented[m]))
        sample = prior.get("sample_metadata", {}).get("sample_alias", "NOT_REPORTED")
        lines = ["QUALITY ENCODING AUDIT COMPLETE - REVIEW REQUIRED", f"BioProject: {PROJECT}",
                 f"Run: {args.run}", f"Sample alias: {sample}", f"Read pairs scanned: {i:,}",
                 "Archive file sizes and MD5 checksums: MATCH (2/2)",
                 "FASTQ/gzip structure and paired identifiers: PASS", "Input file metadata unchanged: YES", ""]
        for m in ("R1", "R2"):
            s = summaries[m]
            lines += [f"{m}:", f"  Reads: {s['reads']:,}; bases: {s['bases']:,}",
                      f"  Stored Phred+33 scores observed: {s['stored_phred33_scores']}",
                      f"  Quality pattern: {s['quality_state']}",
                      f"  Reads with a single uniform score within the read: {s['reads_uniform_within_read']:,}",
                      f"  Reads with Q30 at EVERY position: {s['reads_all_quality_Q30']:,}",
                      f"  Reads consisting entirely of N: {s['reads_entirely_N']:,}",
                      f"  Reads with more than 5 N bases: {s['reads_with_more_than_5_N']:,}",
                      f"  Reads shorter than 75 bases: {s['reads_shorter_than_75']:,} ({100*s['reads_shorter_than_75']/i:.4f}%)",
                      f"  Reads containing symbols other than A/C/G/T/N: {s['reads_with_other_than_ACGTN']:,}", ""]
        lines += ["FASTQC OVERREPRESENTED SEGMENTS (no DNA strings included):"]
        for m in ("R1", "R2"):
            if not overrepresented[m]:
                lines.append(f"  {m}: No reported overrepresented segment.")
            for r in overrepresented[m]:
                lines.append(f"  {m}, row {r['row']}: count={r['reported_count']:,}; percent={r['reported_percent']:.6f}; "
                             f"segment_length={r['reported_segment_length']}; N_bases={r['N_bases_in_segment']}; "
                             f"segment_all_N={r['segment_all_N']}; possible_source={r['possible_source_as_reported']}")
        lines += ["  FastQC's listed segment need not be the complete original read.", "",
                  "LENGTH/N-ONLY FORECAST - NO FILTERING HAS BEEN PERFORMED:",
                  "  Each mate is evaluated at length >=75 and N count <=5, before adapter/quality processing."]
        for key in ("both_pass", "R1_only_pass", "R2_only_pass", "neither_pass"):
            lines.append(f"  {key}: {pair_states[key]:,} pairs ({100*pair_states[key]/i:.4f}%)")
        lines += [f"  Pairs with both reads entirely N: {n_pairs['both_all_N']:,}", "",
                  "INTERPRETATION LIMITS:",
                  "- Stored quality values are not proof that the original instrument qualities were preserved.",
                  "- Constant scores cannot identify lower-quality positions within a read.",
                  "- A Q30-only pattern is compatible with simplified/assigned qualities, but does not prove SRA Lite use.",
                  "- Within-read variation, if present, does not by itself establish original-quality provenance.",
                  "- This forecast does not include adapter trimming, quality trimming or classification.",
                  "- No reads were modified, removed, uploaded or classified. No pathogen result was produced.",
                  f"Reports: {output.relative_to(root)}", ""]
        summary = "\n".join(lines)
        (output / "quality_audit_summary.txt").write_text(summary, encoding="utf-8")
        share = ["quality_audit_summary.txt", "quality_histogram.tsv", "length_histogram.tsv",
                 "N_count_histogram.tsv", "overrepresented_composition_no_sequences.tsv"]
        with zipfile.ZipFile(output / "quality_audit_review.zip", "w", zipfile.ZIP_DEFLATED) as z:
            for name in share:
                z.write(output / name, arcname=name)
            z.writestr("README.txt", "Aggregate metrics only; no DNA sequences, read identifiers or raw FASTQ files.\n")
        emit("\n" + summary)
        emit("For review: quality_audit_summary.txt and quality_audit_review.zip; do not upload FASTQ files.")
        return 0
    except KeyboardInterrupt:
        record.update(status="INTERRUPTED", error="Interrupted by user; no input files intentionally modified.")
        emit("\nSTOP: Interrupted. A rerun starts a new full audit in a separate output folder.")
        return 130
    except Exception as exc:
        record.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        emit(f"\nSTOP: {type(exc).__name__}: {exc}")
        emit("Do not trim or classify reads. Share this error or quality_audit_log.txt.")
        return 1
    finally:
        record["finished_utc"] = utc()
        if output is not None:
            (output / "quality_audit_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
            if record["status"] != "COMPLETE_REVIEW_REQUIRED":
                (output / "quality_audit_summary.txt").write_text(
                    f"Status: {record['status']}\nError: {record.get('error', '')}\n", encoding="utf-8")
            emit(f"Audit folder: {output.relative_to(root)}")
        if log is not None:
            log.close()


if __name__ == "__main__":
    raise SystemExit(main())
