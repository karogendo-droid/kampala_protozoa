#!/usr/bin/env python3
"""Run read-only FastQC and MultiQC on one previously validated paired FASTQ run.

Place in the project's scripts directory and run from the project, for example:
  python scripts/04_run_raw_qc.py --download metadata/download_SRR28434174_TIMESTAMP

Requires Python >=3.10, FastQC 0.12.1, Java, and MultiQC on PATH.
No Python packages beyond the standard library are required by this wrapper.
The script does not trim, remove, classify, or upload reads. It preserves earlier
results by creating a new timestamped output directory on every invocation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
import zipfile

PROJECT = "PRJNA1091428"
VERSION = "1.0.0"
FASTQC_VERSION = "FastQC v0.12.1"
# These standard FastQC modules contain aggregate statistics, not read sequences.
REVIEW_MODULES = {
    "Basic Statistics", "Per base sequence quality", "Per tile sequence quality",
    "Per sequence quality scores", "Per base sequence content",
    "Per sequence GC content", "Per base N content", "Sequence Length Distribution",
    "Sequence Duplication Levels", "Adapter Content",
}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def digest(path: Path, name: str = "sha256") -> str:
    h = hashlib.md5(usedforsecurity=False) if name == "md5" else hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def resolve_inside(root: Path, path: Path) -> Path:
    resolved = (path if path.is_absolute() else root / path).resolve()
    if not resolved.is_relative_to(root):
        raise ValueError("The supplied path must be within the project folder.")
    return resolved


def load_inputs(root: Path, folder: Path, run: str) -> tuple[dict, list[tuple[Path, dict]]]:
    record_path = folder / "download_record.json"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("project") != PROJECT or record.get("run_accession") != run:
        raise ValueError("The download record is for a different project or run.")
    if record.get("status") != "DOWNLOAD_AND_STRUCTURE_CHECKS_PASSED":
        raise ValueError("The previous download/structure check did not finish successfully.")
    checks = record.get("validation", {})
    for key in ("gzip_integrity", "fastq_structure", "paired_identifiers"):
        if checks.get(key) != "PASS":
            raise ValueError(f"Previous validation does not report PASS for {key}.")
    expected_pairs = checks.get("read_pairs")
    if not isinstance(expected_pairs, int) or expected_pairs < 1:
        raise ValueError("The previous record lacks a positive paired-read count.")
    rows = record.get("files", [])
    names = {f"{run}_1.fastq.gz", f"{run}_2.fastq.gz"}
    if len(rows) != 2 or {r.get("filename") for r in rows} != names:
        raise ValueError("Expected exactly the two previously validated paired FASTQ files.")
    result = []
    for row in sorted(rows, key=lambda r: r["filename"]):
        candidate = root / row.get("local_path", "")
        expected = root / "data" / "raw" / run / row["filename"]
        if candidate.is_symlink() or candidate != expected or not candidate.is_file():
            raise ValueError(f"Missing or unexpected raw-file path: {row['filename']}")
        path = resolve_inside(root, candidate)
        if not re.fullmatch(r"[0-9a-fA-F]{32}", str(row.get("expected_md5", ""))):
            raise ValueError("The source record lacks a valid archive MD5.")
        result.append((path, row))
    return record, result


def parse_fastqc_archive(path: Path) -> dict:
    with zipfile.ZipFile(path) as z:
        candidates = [n for n in z.namelist() if n.endswith("/fastqc_data.txt")]
        if len(candidates) != 1:
            raise ValueError(f"Expected one fastqc_data.txt in {path.name}.")
        text = z.read(candidates[0]).decode("utf-8")
    sections: dict[str, dict] = {}
    current = None
    version = "NOT_REPORTED"
    for line in text.splitlines():
        if line.startswith("##FastQC\t"):
            version = line.split("\t", 1)[1]
        elif line == ">>END_MODULE":
            current = None
        elif line.startswith(">>"):
            pieces = line[2:].split("\t", 1)
            if len(pieces) != 2:
                raise ValueError("Unexpected FastQC module header.")
            current = pieces[0]
            sections[current] = {"status": pieces[1].upper(), "lines": []}
        elif current is not None:
            sections[current]["lines"].append(line)
    basic = {}
    for line in sections.get("Basic Statistics", {}).get("lines", []):
        if line and not line.startswith("#") and "\t" in line:
            k, v = line.split("\t", 1)
            basic[k] = v
    if "Filename" not in basic or "Total Sequences" not in basic:
        raise ValueError("FastQC output lacks required basic statistics.")
    count = int(basic["Total Sequences"].replace(",", ""))
    return {"version": version, "basic": basic, "read_count": count, "sections": sections}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--download", type=Path, required=True, help="Previous validated-download report folder")
    ap.add_argument("--run", default="SRR28434174")
    args = ap.parse_args()
    if not re.fullmatch(r"[SED]RR\d+", args.run):
        ap.error("Invalid sequencing-run accession.")
    script = Path(__file__).resolve()
    root = script.parent.parent
    if script.parent.name != "scripts" or not (root / "metadata").is_dir():
        print("STOP: Put this script in the existing project's scripts folder.")
        return 1
    output = None
    log = None
    record = {"project": PROJECT, "run_accession": args.run, "script_version": VERSION,
              "script_sha256": digest(script), "started_utc": utc(), "status": "STARTED",
              "python": sys.version, "architecture": platform.machine(),
              "conda_environment": os.environ.get("CONDA_DEFAULT_ENV", "NOT_REPORTED"),
              "commands": [], "input_checks": [], "versions": {},
              "trimming_performed": False, "classification_performed": False,
              "upload_commands_executed": False}

    def emit(message: str) -> None:
        print(message, flush=True)
        if log is not None:
            log.write(message + "\n")
            log.flush()

    def run_command(command: list[str], name: str, *, echo: bool = True) -> str:
        assert output is not None
        command_record = {"name": name, "argv": command, "started_utc": utc()}
        record["commands"].append(command_record)
        emit("COMMAND: " + shlex.join(command))
        target = output / "logs" / f"{name}.log"
        chunks = []
        proc = subprocess.Popen(command, cwd=output, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        try:
            with target.open("w", encoding="utf-8") as handle:
                assert proc.stdout is not None
                for line in proc.stdout:
                    handle.write(line)
                    handle.flush()
                    chunks.append(line)
                    if echo:
                        emit(line.rstrip("\r\n"))
                code = proc.wait()
        except BaseException:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
            raise
        finally:
            if proc.stdout is not None:
                proc.stdout.close()
            command_record.update(finished_utc=utc(), returncode=proc.poll())
        if code != 0:
            raise RuntimeError(f"{name} failed (exit {code}). Read {target.relative_to(root)}.")
        return "".join(chunks).strip()

    try:
        folder = resolve_inside(root, args.download.expanduser())
        previous, inputs = load_inputs(root, folder, args.run)
        executables = {}
        for tool in ("fastqc", "multiqc", "java"):
            path = shutil.which(tool)
            if not path:
                raise ValueError(f"{tool} is not available. Activate the kampala-qc environment and check installation.")
            executables[tool] = path
        output = root / "results" / "qc" / "raw" / f"{args.run}_{stamp()}"
        output.mkdir(parents=True, exist_ok=False)
        for name in ("logs", "fastqc", "multiqc", "provenance"):
            (output / name).mkdir()
        log = (output / "qc_log.txt").open("w", encoding="utf-8")
        record["output_directory"] = str(output.relative_to(root))
        record["source_download"] = str(folder.relative_to(root))
        shutil.copy2(folder / "download_record.json", output / "provenance" / "download_record.json")
        record["source_download_sha256"] = digest(folder / "download_record.json")
        emit(f"RAW-ARCHIVE QUALITY ASSESSMENT: {args.run}")
        emit("Reading the existing paired files. No reads will be trimmed, discarded, merged or classified.")
        emit("Here 'raw' means unchanged since archive download, not proven unprocessed by the submitter.")
        for path, row in inputs:
            emit(f"Verifying archived input: {path.name}")
            observed = digest(path, "md5")
            size = path.stat().st_size
            if observed != row["expected_md5"].lower() or size != int(row["expected_bytes"]):
                raise ValueError(f"Input no longer matches the verified archive record: {path.name}")
            record["input_checks"].append({"filename": path.name, "bytes": size, "md5_before": observed})
        for tool, flag in (("fastqc", "--version"), ("multiqc", "--version"), ("java", "-version")):
            record["versions"][tool] = run_command([executables[tool], flag], f"{tool}_version", echo=False)
        if FASTQC_VERSION not in record["versions"]["fastqc"]:
            raise ValueError("This baseline step uses FastQC 0.12.1. Check the kampala-qc environment rather than silently changing versions.")
        conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
        if conda:
            explicit = run_command([conda, "list", "--explicit"], "conda_packages", echo=False)
            (output / "provenance" / "conda-explicit.txt").write_text(explicit + "\n", encoding="utf-8")
        run_command([executables["fastqc"], "--java", executables["java"], "--threads", "2",
                     "--noextract", "--format", "fastq", "--outdir", str(output / "fastqc"),
                     *[str(p) for p, _ in inputs]], "fastqc")
        parsed = {}
        for mate, (path, _) in zip(("R1", "R2"), inputs):
            stem = path.name.removesuffix(".fastq.gz")
            z = output / "fastqc" / f"{stem}_fastqc.zip"
            html = output / "fastqc" / f"{stem}_fastqc.html"
            if not z.is_file() or not html.is_file():
                raise ValueError(f"FastQC did not produce both expected reports for {mate}.")
            parsed[mate] = parse_fastqc_archive(z)
            if parsed[mate]["basic"]["Filename"] != path.name:
                raise ValueError("FastQC report filename does not match the intended input.")
            if parsed[mate]["read_count"] != previous["validation"]["read_pairs"]:
                raise ValueError(f"FastQC {mate} count differs from the previous full FASTQ validation.")
        cfg = output / "provenance" / "multiqc_config.yaml"
        cfg.write_text("no_version_check: true\nplots_force_interactive: true\n", encoding="utf-8")
        run_command([executables["multiqc"], str(output / "fastqc"), "--module", "fastqc",
                     "--config", str(cfg), "--outdir", str(output / "multiqc"),
                     "--filename", "multiqc_report.html"], "multiqc")
        if not (output / "multiqc" / "multiqc_report.html").is_file():
            raise ValueError("MultiQC did not create the expected HTML report.")
        emit("Rechecking input checksums after QC ...")
        for entry, (path, row) in zip(record["input_checks"], inputs):
            entry["md5_after"] = digest(path, "md5")
            entry["unchanged"] = entry["md5_after"] == entry["md5_before"] and path.stat().st_size == entry["bytes"]
            if not entry["unchanged"]:
                raise ValueError("A raw input changed during processing. Stop and review.")
        sample = previous.get("sample_metadata", {}).get("sample_alias", "NOT_REPORTED")
        lines = ["RAW QC REPORTS GENERATED - INTERPRETATION PENDING",
                 f"BioProject: {PROJECT}", f"Run: {args.run}", f"Sample alias: {sample}",
                 f"FastQC: {record['versions']['fastqc']}", f"MultiQC: {record['versions']['multiqc']}",
                 f"Read pairs in previous complete validation: {previous['validation']['read_pairs']:,}",
                 "Raw file checksums unchanged: YES", "", "BASIC STATISTICS"]
        for mate in ("R1", "R2"):
            b = parsed[mate]["basic"]
            lines.extend([f"{mate}: {b['Filename']}", f"  Total sequences: {parsed[mate]['read_count']:,}",
                          f"  Sequence length: {b.get('Sequence length', 'NOT_REPORTED')}",
                          f"  GC percent: {b.get('%GC', 'NOT_REPORTED')}",
                          f"  Encoding: {b.get('Encoding', 'NOT_REPORTED')}"])
        lines.extend(["", "FASTQC MODULE FLAGS (not automatic decisions to discard a library)", "Module\tR1\tR2"])
        modules = list(dict.fromkeys(list(parsed["R1"]["sections"]) + list(parsed["R2"]["sections"])))
        for name in modules:
            flags = [parsed[m]["sections"].get(name, {}).get("status", "NOT_REPORTED") for m in ("R1", "R2")]
            lines.append(name + "\t" + "\t".join(flags))
        lines.extend(["", "No reads were trimmed, removed, merged or classified.",
                      "Two FastQC entries represent the paired files of one library, not two specimens.",
                      "Module warnings/failures require interpretation in this metagenomic context.",
                      "No pathogen or specimen-level eligibility conclusion has been made.",
                      f"Reports: {output.relative_to(root)}", ""])
        summary = "\n".join(lines)
        (output / "qc_summary.txt").write_text(summary, encoding="utf-8")
        # Make a small review archive of numeric QC tables, without raw read strings.
        with zipfile.ZipFile(output / "qc_review_bundle.zip", "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("qc_summary.txt", summary)
            for mate in ("R1", "R2"):
                p = parsed[mate]
                aggregate_lines = [f"##FastQC\t{p['version']}"]
                for name, section in p["sections"].items():
                    if name in REVIEW_MODULES:
                        aggregate_lines.append(f">>{name}\t{section['status'].lower()}")
                        aggregate_lines.extend(section["lines"])
                        aggregate_lines.append(">>END_MODULE")
                z.writestr(f"{mate}_aggregate_qc_tables.txt", "\n".join(aggregate_lines) + "\n")
                # Keep counts for overrepresented signals, but do not put DNA strings in the review archive.
                over = p["sections"].get("Overrepresented sequences", {}).get("lines", [])
                buf = io.StringIO()
                writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
                writer.writerow(["row", "count", "percentage", "possible_source"])
                index = 0
                for line in over:
                    if not line or line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    if len(fields) >= 3:
                        index += 1
                        writer.writerow([index, fields[1], fields[2], fields[3] if len(fields) > 3 else "NOT_REPORTED"])
                z.writestr(f"{mate}_overrepresented_counts_no_sequences.tsv", buf.getvalue())
            z.writestr("README.txt", "Aggregate QC metrics for review. No FASTQ reads or overrepresented DNA strings are included.\n"
                       "The complete, unchanged FastQC reports remain in the local fastqc directory.\n")
        record["status"] = "QC_REPORTS_GENERATED_INTERPRETATION_PENDING"
        record["fastqc_basic_statistics"] = {m: parsed[m]["basic"] for m in parsed}
        record["fastqc_flags"] = {m: {k: v["status"] for k, v in parsed[m]["sections"].items()} for m in parsed}
        emit("\n" + summary)
        emit("For review, share qc_summary.txt and qc_review_bundle.zip, NOT the FASTQ files.")
        return 0
    except KeyboardInterrupt:
        record.update(status="INTERRUPTED", error="Stopped by user; original FASTQ files were not intentionally modified.")
        emit("\nSTOP: Interrupted. Retain this partial output; a rerun creates a separate folder.")
        return 130
    except Exception as exc:
        record.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        emit(f"\nSTOP: {type(exc).__name__}: {exc}")
        emit("Do not trim or classify reads. Share this error or qc_log.txt for review.")
        return 1
    finally:
        record["finished_utc"] = utc()
        if output is not None:
            (output / "qc_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
            if record["status"] != "QC_REPORTS_GENERATED_INTERPRETATION_PENDING":
                (output / "qc_summary.txt").write_text(f"Status: {record['status']}\nRun: {args.run}\nError: {record.get('error', '')}\n", encoding="utf-8")
                emit(f"Diagnostic folder: {output.relative_to(root)}")
        if log is not None:
            log.close()


if __name__ == "__main__":
    raise SystemExit(main())
