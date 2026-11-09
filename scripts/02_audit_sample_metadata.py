#!/usr/bin/env python3
"""Audit ENA metadata for the runs returned by 01_fetch_run_metadata.py.

Run from the project: python scripts/02_audit_sample_metadata.py
Dependency: requests (already installed in kampala-protozoa).
Only small metadata records are fetched. FASTQ files are NEVER downloaded.
Original reports are preserved. Missing values and unresolved conflicts are explicit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

PROJECT = "PRJNA1091428"
REPORT_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
XML_API = "https://www.ebi.ac.uk/ena/browser/api/xml/"
FIELDS = [
    "run_accession", "study_accession", "secondary_study_accession",
    "sample_accession", "secondary_sample_accession", "experiment_accession",
    "sample_alias", "sample_title", "sample_description", "library_name",
    "library_strategy", "library_source", "library_selection", "library_layout",
    "instrument_platform", "instrument_model", "collection_date", "country",
    "isolation_source", "sampling_site", "library_construction_protocol",
    "extraction_protocol", "first_public", "last_updated", "read_count", "base_count",
    "fastq_ftp", "fastq_md5", "fastq_bytes",
]
EMPTY = {"", "not_in_report", "not_reported", "not reported", "missing", "not provided",
         "not collected", "not applicable", "unknown", "na", "n/a", "null", "none"}


def clean(value: str | None) -> str:
    return " ".join((value or "").split())


def present(value: str | None) -> bool:
    return clean(value).lower() not in EMPTY


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_report(text: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(text.lstrip("\ufeff")), delimiter="\t")
    if not reader.fieldnames or "run_accession" not in reader.fieldnames:
        raise ValueError("Response is not a run table: run_accession is missing.")
    rows = list(reader)
    ids = [clean(row.get("run_accession")) for row in rows]
    if not ids or any(not re.fullmatch(r"[SED]RR\d+", x) for x in ids):
        raise ValueError("No valid run records, or an unexpected run identifier.")
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate run rows; manual review is required.")
    return rows


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_record(data: bytes, kind: str, accession: str) -> ET.Element:
    root = ET.fromstring(data)
    # Remove namespaces to handle namespace-declared and legacy ENA records alike.
    for element in root.iter():
        element.tag = element.tag.rsplit("}", 1)[-1]
    for record in root.iter(kind):
        ids = {clean(record.get("accession"))}
        ids.update(clean(x.text) for x in record.findall("./IDENTIFIERS/*"))
        if accession in ids:
            return record
    raise ValueError(f"No matching {kind} record for {accession}; no substitution was made.")


def text_at(record: ET.Element | None, path: str) -> str:
    return clean(record.findtext(path)) if record is not None else ""


def reference(record: ET.Element | None, path: str, pattern: str) -> str:
    if record is None:
        return ""
    element = record.find(path)
    if element is None:
        return ""
    possibilities = [element.get("accession", "")]
    possibilities += [x.text or "" for x in element.findall("./IDENTIFIERS/*")]
    return next((x.strip() for x in possibilities if re.fullmatch(pattern, x.strip())), "")


def attributes(record: ET.Element | None) -> list[tuple[str, str, str]]:
    if record is None:
        return []
    return [(text_at(x, "TAG"), text_at(x, "VALUE"), text_at(x, "UNITS"))
            for x in record.findall(".//SAMPLE_ATTRIBUTE")]


def attribute_values(items: list[tuple[str, str, str]], tags: set[str]) -> str:
    values = []
    for tag, value, unit in items:
        normalized = re.sub(r"[_\s]+", " ", tag.lower()).strip()
        if normalized in tags and present(value):
            v = value + (f" [{unit}]" if unit else "")
            if v not in values:
                values.append(v)
    return " | ".join(values)


def pick(report_value: str | None, xml_value: str | None, field: str,
         run: str, issues: list[dict]) -> str:
    a, b = clean(report_value), clean(xml_value)
    if present(a) and present(b) and a.casefold() != b.casefold():
        issues.append({"run_accession": run, "issue": f"{field} differs between report and XML",
                       "details": f"Report={a}; XML={b}"})
        return f"CONFLICT: report={a} | XML={b}"
    return a if present(a) else b if present(b) else "NOT_REPORTED"


class Fetcher:
    def __init__(self, folder: Path) -> None:
        self.folder = folder
        self.log: list[dict] = []
        self.session = requests.Session()
        retry = Retry(total=2, backoff_factor=1,
                      status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
        self.session.mount("https://", HTTPAdapter(max_retries=retry))
        self.session.headers["User-Agent"] = "KampalaProtozoaMetadataAudit/2.0"
        self.cache: dict[tuple[str, str], ET.Element] = {}

    def get(self, url: str, path: Path, params: dict | None = None) -> bytes:
        time.sleep(0.35)
        entry = {"url": url, "parameters": params,
                 "requested_at_utc": datetime.now(timezone.utc).isoformat()}
        self.log.append(entry)
        try:
            response = self.session.get(url, params=params, timeout=(15, 75))
            entry.update(url=response.url, http_status=response.status_code)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(response.content)
            entry.update(saved_as=str(path.relative_to(self.folder)), sha256=sha256(response.content))
            response.raise_for_status()
            return response.content
        except requests.RequestException as exc:
            entry["error"] = str(exc)
            raise

    def xml(self, accession: str, kind: str, issues: list[dict], run: str) -> ET.Element | None:
        if not accession:
            issues.append({"run_accession": run, "issue": f"No linked {kind} accession", "details": ""})
            return None
        key = (accession, kind)
        if key in self.cache:
            return self.cache[key]
        try:
            data = self.get(XML_API + accession, self.folder / "xml" / f"{accession}.xml")
            record = parse_record(data, kind, accession)
            self.cache[key] = record
            return record
        except (requests.RequestException, ET.ParseError, ValueError) as exc:
            issues.append({"run_accession": run, "issue": f"{kind} record could not be verified",
                           "details": f"{accession}: {exc}"})
            return None


def make_file_rows(report: dict, issues: list[dict]) -> list[dict]:
    run = report["run_accession"]
    values = [[clean(x) for x in (report.get(k) or "").split(";") if clean(x)]
              for k in ("fastq_ftp", "fastq_md5", "fastq_bytes")]
    links, checksums, sizes = values
    if not links or len({len(x) for x in values}) != 1:
        issues.append({"run_accession": run, "issue": "FASTQ manifest incomplete",
                       "details": "File links, MD5 checksums and byte counts are absent or have unequal counts."})
        return []
    result = []
    for link, md5, size in zip(links, checksums, sizes):
        parsed = urlsplit(link if "://" in link else "https://" + link)
        good = (parsed.hostname in {"ftp.sra.ebi.ac.uk", "ftp.ebi.ac.uk"}
                and parsed.scheme in {"https", "ftp"} and parsed.path.endswith(".fastq.gz")
                and re.fullmatch(r"[a-fA-F0-9]{32}", md5) and size.isdigit() and int(size) > 0)
        if not good:
            issues.append({"run_accession": run, "issue": "Unexpected FASTQ manifest entry", "details": link})
            continue
        result.append({"run_accession": run, "filename": Path(parsed.path).name,
                       "https_url": urlunsplit(("https", parsed.netloc, parsed.path, "", "")),
                       "expected_md5": md5.lower(), "expected_bytes": size,
                       "compressed_GB": f"{int(size)/1e9:.6f}", "link_tested": "NO"})
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Original ena_run_report.tsv; otherwise find the latest saved report.")
    args = parser.parse_args()
    script = Path(__file__).resolve()
    root = script.parent.parent
    if script.parent.name != "scripts" or not (root / "metadata").is_dir():
        print("STOP: Put the script in the project's scripts folder; metadata/ must already exist.")
        return 1
    fetch: Fetcher | None = None
    folder: Path | None = None
    provenance = {"project": PROJECT, "script_sha256": sha256(script.read_bytes()),
                  "python": sys.version, "requests_version": requests.__version__,
                  "sequence_files_downloaded": False, "status": "started"}
    try:
        if args.source:
            source = args.source.expanduser().resolve()
        else:
            sources = sorted((root / "metadata").glob(f"{PROJECT}_*/ena_run_report.tsv"))
            if not sources:
                raise ValueError("No saved ena_run_report.tsv found. Complete Step 21 first.")
            source = sources[-1]
        original = source.read_bytes()
        old_rows = parse_report(original.decode("utf-8-sig"))
        wanted = {row["run_accession"] for row in old_rows}
        now = datetime.now(timezone.utc)
        folder = root / "metadata" / f"{PROJECT}_audit_{now.strftime('%Y%m%dT%H%M%S%fZ')}"
        folder.mkdir(parents=True, exist_ok=False)
        (folder / "input_run_report.tsv").write_bytes(original)
        provenance.update(retrieved_at_utc=now.isoformat(), source_report=str(source),
                          source_report_sha256=sha256(original), original_run_count=len(old_rows))
        fetch = Fetcher(folder)
        issues: list[dict] = []
        print(f"Auditing metadata for {len(wanted)} runs in {PROJECT}.", flush=True)
        print("No FASTQ files, databases or sequence alignments will be downloaded.\n", flush=True)
        print("Requesting explicitly selected run-report fields ...", flush=True)
        indexed: dict[str, dict] = {}
        try:
            payload = fetch.get(REPORT_API, folder / "expanded_run_report.tsv", {
                "accession": PROJECT, "result": "read_run", "format": "tsv", "fields": ",".join(FIELDS)})
            expanded = parse_report(payload.decode("utf-8-sig"))
            indexed = {row["run_accession"]: row for row in expanded}
            if set(indexed) != wanted:
                issues.append({"run_accession": "PROJECT", "issue": "Run list differs from the original report",
                               "details": f"Added={sorted(set(indexed)-wanted)}; absent={sorted(wanted-set(indexed))}. Only original runs audited."})
        except (requests.RequestException, ValueError, UnicodeError) as exc:
            issues.append({"run_accession": "PROJECT", "issue": "Expanded report unavailable",
                           "details": str(exc) + "; attempting linked XML records instead."})

        manifest, evidence, file_rows, protocols = [], [], [], []
        for index, old in enumerate(sorted(old_rows, key=lambda x: x["run_accession"]), 1):
            run = old["run_accession"]
            report = dict(old)
            report.update(indexed.get(run, {}))
            print(f"[{index}/{len(old_rows)}] Retrieving experiment and sample records for {run} ...", flush=True)
            for field, value in report.items():
                if value:
                    evidence.append({"run_accession": run, "source_type": "ENA_RUN_REPORT", "source_accession": run,
                                     "field": field, "value": value, "unit": ""})
            experiment_id = clean(report.get("experiment_accession"))
            if not re.fullmatch(r"[SED]RX\d+", experiment_id):
                run_xml = fetch.xml(run, "RUN", issues, run)
                experiment_id = reference(run_xml, "./EXPERIMENT_REF", r"[SED]RX\d+")
            experiment = fetch.xml(experiment_id, "EXPERIMENT", issues, run)
            sample_id = reference(experiment, "./DESIGN/SAMPLE_DESCRIPTOR", r"(?:[SED]RS\d+|SAM[END][A-Z]?\d+)")
            if not sample_id:
                sample_id = next((clean(report.get(key)) for key in ("sample_accession", "secondary_sample_accession")
                                  if re.fullmatch(r"(?:[SED]RS\d+|SAM[END][A-Z]?\d+)", clean(report.get(key)))), "")
            sample = fetch.xml(sample_id, "SAMPLE", issues, run)
            attrs = attributes(sample)
            for tag, value, unit in attrs:
                evidence.append({"run_accession": run, "source_type": "ENA_SAMPLE_XML", "source_accession": sample_id,
                                 "field": tag, "value": value, "unit": unit})
            sample_identifiers = " | ".join(clean(x.text) for x in sample.findall("./IDENTIFIERS/*") if clean(x.text)) if sample is not None else ""
            sample_title = text_at(sample, "./TITLE")
            sample_description = text_at(sample, "./DESCRIPTION")
            sample_alias = clean(sample.get("alias")) if sample is not None else ""
            library = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR") if experiment is not None else None
            layout_node = library.find("LIBRARY_LAYOUT") if library is not None else None
            layout = list(layout_node)[0].tag if layout_node is not None and len(layout_node) else ""
            platform_node = experiment.find("PLATFORM") if experiment is not None else None
            platform = list(platform_node)[0].tag if platform_node is not None and len(platform_node) else ""
            model = text_at(platform_node, "./*/INSTRUMENT_MODEL")
            for key, value in {"sample_title": sample_title, "sample_alias": sample_alias,
                               "sample_description": sample_description, "identifiers": sample_identifiers}.items():
                if value:
                    evidence.append({"run_accession": run, "source_type": "ENA_SAMPLE_XML", "source_accession": sample_id,
                                     "field": key, "value": value, "unit": ""})
            row = {"run_accession": run, "experiment_accession": experiment_id or "NOT_RETRIEVED",
                   "sample_accession": sample_id or "NOT_RETRIEVED", "sample_identifiers": sample_identifiers or "NOT_RETRIEVED",
                   "sample_alias": pick(report.get("sample_alias"), sample_alias, "sample_alias", run, issues),
                   "sample_title": pick(report.get("sample_title"), sample_title, "sample_title", run, issues),
                   "sample_description": pick(report.get("sample_description"), sample_description, "sample_description", run, issues),
                   "sample_xml_retrieved": "YES" if sample is not None else "NO",
                   "experiment_xml_retrieved": "YES" if experiment is not None else "NO"}
            for field, xml_value in {
                "library_name": text_at(library, "LIBRARY_NAME"),
                "library_strategy": text_at(library, "LIBRARY_STRATEGY"),
                "library_source": text_at(library, "LIBRARY_SOURCE"),
                "library_selection": text_at(library, "LIBRARY_SELECTION"),
                "library_layout": layout, "instrument_platform": platform, "instrument_model": model,
                "library_construction_protocol": text_at(library, "LIBRARY_CONSTRUCTION_PROTOCOL"),
            }.items():
                row[field] = pick(report.get(field), xml_value, field, run, issues)
                if xml_value:
                    evidence.append({"run_accession": run, "source_type": "ENA_EXPERIMENT_XML", "source_accession": experiment_id,
                                     "field": field, "value": xml_value, "unit": ""})
            date = attribute_values(attrs, {"collection date", "sample collection date", "sampling date"})
            site = attribute_values(attrs, {"sampling site", "sample collection site", "collection site", "site", "facility", "treatment plant", "wastewater treatment plant"})
            row["collection_date_as_reported"] = pick(report.get("collection_date"), date, "collection_date", run, issues)
            row["site_as_reported"] = pick(report.get("sampling_site"), site, "sampling_site", run, issues)
            row["country_as_reported"] = pick(report.get("country"), attribute_values(attrs, {"geo loc name", "geographic location (country and/or sea)", "country"}), "country", run, issues)
            row["isolation_source_as_reported"] = pick(report.get("isolation_source"), attribute_values(attrs, {"isolation source"}), "isolation_source", run, issues)
            row["season_as_reported"] = attribute_values(attrs, {"season", "sampling season", "collection season"}) or "NOT_REPORTED"
            row["treatment_stage_as_reported"] = attribute_values(attrs, {"treatment stage", "treatment status", "sampling treatment stage"}) or "NOT_REPORTED"
            row["pooling_or_replicate_as_reported"] = attribute_values(attrs, {"pooling", "pooled", "replicate", "biological replicate", "technical replicate", "composite sample"}) or "NOT_REPORTED"
            row["read_count_reported"] = report.get("read_count") or "NOT_REPORTED"
            row["base_count_reported"] = report.get("base_count") or "NOT_REPORTED"
            row["first_public"] = report.get("first_public") or "NOT_REPORTED"
            row["last_updated"] = report.get("last_updated") or "NOT_REPORTED"
            row["specimen_independence"] = "REQUIRES_REVIEW"
            row["whole_library_eligibility"] = "REQUIRES_REVIEW"
            files = make_file_rows(report, issues)
            file_rows += files
            row["fastq_files_with_manifest_entries"] = str(len(files))
            row["fastq_compressed_GB"] = f"{sum(int(f['expected_bytes']) for f in files)/1e9:.3f}" if files else "NOT_REPORTED"
            manifest.append(row)
            protocols.append("\n".join([f"RUN: {run}", f"EXPERIMENT: {experiment_id}", f"SAMPLE: {sample_id}",
                                         f"Sample title: {row['sample_title']}", f"Sample description: {row['sample_description']}",
                                         f"Design description: {text_at(experiment, './DESIGN/DESIGN_DESCRIPTION') or 'NOT_REPORTED'}",
                                         f"Library construction: {row['library_construction_protocol']}",
                                         f"Extraction protocol (indexed): {report.get('extraction_protocol') or 'NOT_REPORTED'}"]))
        write_tsv(folder / "sample_manifest.tsv", manifest, list(manifest[0]))
        write_tsv(folder / "all_metadata_evidence.tsv", evidence, ["run_accession", "source_type", "source_accession", "field", "value", "unit"])
        write_tsv(folder / "fastq_file_manifest.tsv", file_rows, ["run_accession", "filename", "https_url", "expected_md5", "expected_bytes", "compressed_GB", "link_tested"])
        write_tsv(folder / "issues_to_review.tsv", issues, ["run_accession", "issue", "details"])
        (folder / "library_and_sample_descriptions.txt").write_text("\n\n".join(protocols)+"\n", encoding="utf-8")
        lines = ["METADATA AUDIT RECORDS SAVED", f"BioProject: {PROJECT}", f"Run records audited: {len(manifest)}",
                 f"Experiment XML records retrieved: {sum(r['experiment_xml_retrieved']=='YES' for r in manifest)}/{len(manifest)}",
                 f"Sample XML records retrieved: {sum(r['sample_xml_retrieved']=='YES' for r in manifest)}/{len(manifest)}",
                 f"Issues requiring review: {len(issues)}", ""]
        for row in manifest:
            lines += [row["run_accession"],
                      f"  Sample: {row['sample_accession']} | alias: {row['sample_alias']} | title: {row['sample_title']}",
                      f"  Sequencing: {row['instrument_platform']} / {row['instrument_model']} / {row['library_layout']}",
                      f"  Library: {row['library_strategy']} / {row['library_source']} / {row['library_selection']}",
                      f"  Collection date (reported): {row['collection_date_as_reported']}",
                      f"  Site (reported): {row['site_as_reported']}",
                      f"  Source material (reported): {row['isolation_source_as_reported']}",
                      f"  Season: {row['season_as_reported']} | treatment stage: {row['treatment_stage_as_reported']}",
                      f"  Manifest-listed FASTQ: {row['fastq_compressed_GB']} GB; {row['fastq_files_with_manifest_entries']} files", ""]
        if issues:
            lines += ["ISSUES TO REVIEW:"]
            for issue in issues:
                lines.append(f"  {issue['run_accession']}: {issue['issue']} ({issue['details']})")
            lines.append("")
        lines += ["IMPORTANT:", "- Reported dates are not independently validated collection dates.",
                  "- No site, season or pooling status was inferred from run/sample numbering.",
                  "- NOT_REPORTED means no value was found in the inspected fields; review the descriptions and original study too.",
                  "- REQUIRES_REVIEW is not an exclusion or a confirmed eligibility decision.",
                  "- File-manifest entries do not establish that file downloads work; links were not tested.",
                  "- No sequence files were downloaded or analysed.",
                  "Saved under: " + str(folder.relative_to(root))]
        summary = "\n".join(lines)+"\n"
        (folder / "audit_summary.txt").write_text(summary, encoding="utf-8")
        provenance.update(status="metadata_saved_pending_review", run_records=len(manifest), issues=len(issues))
        print("\n" + summary)
        return 0
    except (OSError, ValueError, UnicodeError, csv.Error) as exc:
        provenance.update(status="failed", error=str(exc))
        print(f"\nSTOP: Metadata audit did not finish.\n{exc}\nNo sequence files were downloaded.", file=sys.stderr)
        return 1
    finally:
        if fetch is not None:
            provenance["requests"] = fetch.log
            fetch.session.close()
        if folder is not None:
            (folder / "retrieval_record.json").write_text(json.dumps(provenance, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
