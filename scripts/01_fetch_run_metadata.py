#!/usr/bin/env python3
"""Fetch ENA run metadata for PRJNA1091428; never download sequence files.

Run from the study project after placing this file in its scripts/ folder:
    python scripts/01_fetch_run_metadata.py

Dependencies: pandas, requests.  Suitable for Python 3.12 or later.
"""
from __future__ import annotations

import hashlib
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

PROJECT = "PRJNA1091428"
ENDPOINT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
PARAMS = {"accession": PROJECT, "result": "read_run", "format": "tsv"}
# No fields parameter: ENA documents this as returning all file-report columns.


def parse_report(text: str) -> pd.DataFrame:
    """Validate a metadata response without inferring missing values."""
    if not text.strip():
        raise ValueError("ENA returned an empty response; access is unresolved.")
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    if "run_accession" not in df.columns:
        raise ValueError("The response is not an ENA run table (run_accession missing).")
    if df.empty:
        raise ValueError("ENA returned a header but no runs; this is not proof of absent data.")
    runs = df["run_accession"].str.strip()
    if not runs.str.fullmatch(r"[SED]RR\d+").all():
        raise ValueError("One or more run identifiers are missing or unexpected.")
    if runs.duplicated().any():
        raise ValueError("Duplicate run rows returned; review the response before proceeding.")
    return df


def reported_size(value: str) -> int | None:
    """Sum ENA's semicolon-separated byte counts; unknown is not zero."""
    fields = str(value).strip().split(";")
    if not fields or not all(item.strip().isdigit() for item in fields):
        return None
    return sum(int(item.strip()) for item in fields)


def make_summary(df: pd.DataFrame) -> pd.DataFrame:
    fields = [
        "run_accession", "sample_accession", "secondary_sample_accession",
        "sample_alias", "sample_title", "experiment_accession", "library_name",
        "library_strategy", "library_source", "library_selection", "library_layout",
        "instrument_platform", "instrument_model", "collection_date", "country",
        "read_count", "base_count",
    ]
    result = pd.DataFrame(index=df.index)
    for field in fields:
        if field in df.columns:
            values = df[field].astype(str).str.strip()
            result[field] = values.mask(values.eq(""), "NOT_REPORTED")
        else:
            result[field] = "NOT_IN_REPORT"
    result["fastq_compressed_GB"] = [
        f"{size / 1_000_000_000:.3f}" if size is not None else "NOT_REPORTED"
        for size in (reported_size(x) for x in df.get("fastq_bytes", pd.Series("", index=df.index)))
    ]
    return result


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    if Path(__file__).resolve().parent.name != "scripts" or not (root / "metadata").is_dir():
        print("STOP: Put this file in the project's scripts folder. The project must also have a metadata folder.", file=sys.stderr)
        return 1

    now = datetime.now(timezone.utc)
    folder = root / "metadata" / f"{PROJECT}_{now.strftime('%Y%m%dT%H%M%S%fZ')}"
    folder.mkdir(parents=True, exist_ok=False)
    print(f"Retrieving metadata only for {PROJECT} ...", flush=True)
    print("No FASTQ files or reference databases will be downloaded.", flush=True)

    session = requests.Session()
    retries = Retry(total=2, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers["User-Agent"] = "KampalaProtozoaMetadataAudit/1.0"
    provenance = {
        "project": PROJECT, "retrieved_at_utc": now.isoformat(),
        "endpoint": ENDPOINT, "parameters": PARAMS,
        "python": sys.version.split()[0], "pandas": pd.__version__,
        "requests": requests.__version__, "sequence_files_downloaded": False,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "status": "started",
    }
    try:
        response = session.get(ENDPOINT, params=PARAMS, timeout=(20, 120))
        provenance.update({"request_url": response.url, "http_status": response.status_code})
        response.raise_for_status()
        payload = response.content
        (folder / "ena_response.txt").write_bytes(payload)
        df = parse_report(payload.decode("utf-8-sig"))
        (folder / "ena_run_report.tsv").write_bytes(payload)
        provenance["report_sha256"] = hashlib.sha256(payload).hexdigest()
        summary = make_summary(df)
        summary.to_csv(folder / "run_summary.tsv", sep="\t", index=False)
        df["run_accession"].to_csv(folder / "run_accessions.txt", index=False, header=False)

        lines = ["METADATA RETRIEVAL COMPLETE", f"BioProject: {PROJECT}", f"Run records returned: {len(df)}"]
        for field in ("sample_accession", "secondary_sample_accession"):
            if field in df:
                nonempty = df[field].str.strip()
                nonempty = nonempty[nonempty.ne("")]
                lines.append(f"Distinct {field} values reported: {nonempty.nunique()}")
        lines.append("Archive sample identifiers are NOT a verified count of independent specimens.")
        for field in ("library_strategy", "library_layout", "instrument_platform", "instrument_model"):
            if field in df:
                values = sorted({v.strip() for v in df[field] if v.strip()})
                lines.append(f"{field}: {', '.join(values) if values else 'NOT_REPORTED'}")
        urls = df.get("fastq_ftp", pd.Series("", index=df.index))
        lines.append(f"Runs with FASTQ links listed: {urls.str.strip().ne('').sum()} / {len(df)}")
        sizes = [reported_size(v) for v in df.get("fastq_bytes", pd.Series("", index=df.index))]
        known = [n for n in sizes if n is not None]
        if known:
            lines.append(f"Sum of reported compressed FASTQ sizes: {sum(known) / 1_000_000_000:.3f} GB ({len(known)}/{len(df)} runs with sizes)")
            lines.append("This is a sum of report entries, not a workspace requirement; download links are not yet tested.")
        else:
            lines.append("Compressed FASTQ sizes: NOT_REPORTED")
        if "collection_date" not in df:
            lines.append("Collection dates are not in this report; full sample records still need review.")
        else:
            count = df["collection_date"].str.strip().ne("").sum()
            lines.append(f"Nonblank collection_date entries: {count}/{len(df)} (values remain unverified)")
        lines.extend(["", "RUN SUMMARY"])
        display_fields = ["run_accession", "sample_accession", "sample_alias", "library_layout", "fastq_compressed_GB"]
        lines.append(summary[display_fields].to_string(index=False))
        lines.extend(["", "NOT_IN_REPORT = field absent from this response; NOT_REPORTED = blank value in this response.",
                      "Neither label establishes that information is absent from the full sample records.",
                      "Saved under: " + str(folder.relative_to(root)),
                      "Stop here: metadata retrieval is not a completed eligibility audit."])
        text = "\n".join(lines) + "\n"
        (folder / "audit_summary.txt").write_text(text, encoding="utf-8")
        provenance.update({"status": "retrieved", "run_records": len(df)})
        print("\n" + text)
    except (requests.RequestException, ValueError, UnicodeError, OSError) as exc:
        provenance.update({"status": "failed", "error": str(exc)})
        print(f"\nSTOP: Metadata retrieval did not finish successfully.\n{exc}\nNo sequence files were downloaded.\nCopy this error for review; do not infer that the dataset is absent.", file=sys.stderr)
        return 1
    finally:
        session.close()
        (folder / "retrieval_record.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
