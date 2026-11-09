#!/usr/bin/env python3
"""Check public archive metadata for alternatives to constant-quality FASTQ files.

Run from the existing project. Python 3.10+; standard library only.
This script reads the successful step-05 record, NOT the sequencing files.
It requests metadata from ENA and NCBI, with bounded response sizes and retries.
It never downloads a sequence archive, opens a cloud account, accepts charges,
changes quality scores, or filters reads. Listed files are NOT verified downloads.

Example:
  python scripts/06_check_quality_sources.py --quality-audit results/qc/quality_audit/SRR28434174_TIMESTAMP --run SRR28434174

Official API descriptions:
  https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html
  https://www.ncbi.nlm.nih.gov/Traces/sdl/2/
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import socket
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit, urlunsplit
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
import zipfile

PROJECT = 'PRJNA1091428'
VERSION = '1.0.0'
MAX_RESPONSE_BYTES = 5 * 1024 * 1024
TIMEOUT_SECONDS = 30
COLUMNS = ['source', 'category', 'run_accession', 'filename', 'file_type',
           'bytes_reported', 'md5_reported', 'quality_metadata',
           'noqual_reported', 'rehydration_required_reported',
           'location_reported_not_tested', 'extra_metadata']
ENA_FIELDS = ['run_accession', 'study_accession', 'secondary_study_accession',
              'experiment_accession', 'sample_accession', 'fastq_ftp',
              'fastq_md5', 'fastq_bytes', 'submitted_ftp', 'submitted_format',
              'submitted_md5', 'submitted_bytes', 'sra_ftp', 'sra_md5', 'sra_bytes']


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def stamp() -> str:
    return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_location(value: str) -> str:
    """Do not put query strings, credentials, or signed tokens in review outputs."""
    if not value:
        return 'NOT_REPORTED'
    value = value.strip()
    if '://' not in value:
        value = 'ftp://' + value
    parts = urlsplit(value)
    if not parts.hostname:
        return 'UNPARSED_LOCATION'
    host = parts.hostname
    if parts.port:
        host += ':' + str(parts.port)
    clean = urlunsplit((parts.scheme, host, parts.path, '', ''))
    if parts.query or parts.fragment or parts.username or parts.password:
        clean += ' [credentials/query/fragment omitted]'
    return clean


def row(source: str, category: str, run: str, **values) -> dict:
    out = {k: 'NOT_REPORTED' for k in COLUMNS}
    out.update(source=source, category=category, run_accession=run,
               quality_metadata='NOT_DETERMINED_FROM_METADATA', **values)
    return out


def text_value(value) -> str:
    if value is None or value == '':
        return 'NOT_REPORTED'
    if isinstance(value, bool):
        return 'true' if value else 'false'
    return str(value)


def field_values(value: str) -> list[str]:
    return [s.strip() for s in value.split(';')] if value else []


def parse_ena_report(payload: bytes, run: str) -> tuple[list[dict], list[str]]:
    text = payload.decode('utf-8-sig')
    reader = csv.DictReader(io.StringIO(text), delimiter='\t')
    if not reader.fieldnames or 'run_accession' not in reader.fieldnames:
        raise ValueError('ENA response is not a run report with a run_accession column.')
    records = list(reader)
    if len(records) != 1 or records[0].get('run_accession') != run:
        raise ValueError('ENA report did not return exactly the requested run.')
    record = records[0]
    studies = {record.get('study_accession', ''), record.get('secondary_study_accession', '')}
    issues = []
    if PROJECT not in studies:
        issues.append('ENA report does not independently confirm the BioProject; retain the earlier provenance linkage.')
    rows = []
    for prefix, category, default_format in (
        ('fastq', 'ENA_GENERATED_FASTQ', 'fastq'),
        ('submitted', 'ENA_SUBMITTED_FILE', ''),
        ('sra', 'ENA_SRA_FILE', 'sra'),
    ):
        key = prefix + '_ftp'
        if key not in record:
            issues.append(f'ENA field {key} was absent; this is NOT evidence that the file class is unavailable.')
            continue
        locations = field_values(record.get(key, ''))
        if not locations:
            issues.append(f'No {key} locations listed in this response; this is NOT a complete availability assessment.')
            continue
        parallel = {}
        for suffix in ('bytes', 'md5', 'format'):
            values = field_values(record.get(prefix + '_' + suffix, ''))
            if values and len(values) != len(locations):
                # Do not arbitrarily associate one checksum/size with multiple files.
                issues.append(f'{prefix}_{suffix} has {len(values)} values for {len(locations)} locations; association left unresolved.')
                values = []
            parallel[suffix] = values
        for i, location in enumerate(locations):
            clean = safe_location(location)
            name = location.split('?', 1)[0].rstrip('/').rsplit('/', 1)[-1]
            values = {k: (parallel[k][i] if parallel[k] else '') for k in parallel}
            rows.append(row('ENA_FILE_REPORT', category, run, filename=name,
                            file_type=values['format'] or default_format or 'NOT_REPORTED',
                            bytes_reported=values['bytes'] or 'NOT_REPORTED',
                            md5_reported=values['md5'] or 'NOT_REPORTED',
                            location_reported_not_tested=clean))
    return rows, issues


def parse_ena_xml(payload: bytes, run: str) -> tuple[list[dict], list[str]]:
    root = ET.fromstring(payload)
    for elem in root.iter():
        elem.tag = elem.tag.split('}')[-1]
    runs = [e for e in root.iter('RUN') if e.get('accession') == run]
    if len(runs) != 1:
        raise ValueError('ENA XML did not contain exactly the requested RUN record.')
    rows = []
    for elem in runs[0].iter('FILE'):
        a = dict(elem.attrib)
        name = a.get('filename') or a.get('name') or 'NOT_REPORTED'
        checksum = a.get('checksum', '') if a.get('checksum_method', '').lower() == 'md5' else ''
        rows.append(row('ENA_RUN_XML', 'RUN_XML_FILE_METADATA', run,
                        filename=name, file_type=a.get('filetype', 'NOT_REPORTED'),
                        bytes_reported=a.get('size', 'NOT_REPORTED'),
                        md5_reported=checksum or 'NOT_REPORTED',
                        extra_metadata=json.dumps(a, sort_keys=True)))
    issues = [] if rows else ['No FILE elements found in the ENA RUN XML; this does not establish that original files are unavailable.']
    return rows, issues


def is_true(value) -> bool:
    return value is True or (isinstance(value, str) and value.lower() == 'true')


def parse_sdl(payload: bytes, run: str, source: str) -> tuple[list[dict], list[str]]:
    obj = json.loads(payload)
    if isinstance(obj, dict):
        bundles = obj.get('result', [obj] if ('bundle' in obj or 'accession' in obj) else [])
    elif isinstance(obj, list):
        bundles = obj
    else:
        raise ValueError('Unrecognized NCBI SDL response structure.')
    if not isinstance(bundles, list):
        raise ValueError('NCBI SDL result is not a list.')
    bundles = [b for b in bundles if isinstance(b, dict)
               and (b.get('bundle') == run or b.get('accession') == run)]
    if len(bundles) != 1:
        raise ValueError('NCBI SDL response did not identify exactly the requested accession.')
    bundle = bundles[0]
    if str(bundle.get('status', '')) != '200':
        # HTTP 200 can contain an application-level 403/404/500.
        raise ValueError('NCBI application status ' + str(bundle.get('status')) + ': '
                         + str(bundle.get('msg', bundle.get('message', 'No explanation'))))
    files = bundle.get('files', [])
    if not isinstance(files, list):
        raise ValueError('NCBI SDL files entry is not a list.')
    rows, issues = [], []
    for f in files:
        if not isinstance(f, dict):
            raise ValueError('NCBI SDL returned an invalid file entry.')
        if f.get('accession') and f['accession'] != run:
            issues.append('An SDL file linked to another accession was excluded from this inventory.')
            continue
        name = str(f.get('name', 'NOT_REPORTED'))
        ftype = str(f.get('type', 'NOT_REPORTED'))
        nq = f.get('noqual', bundle.get('noqual'))
        category, quality = 'NCBI_FILE_METADATA', 'NOT_DETERMINED_FROM_METADATA'
        if is_true(nq):
            quality = 'SIMPLIFIED_QUALITIES_EXPLICITLY_REPORTED'
        elif any(token in name.lower() for token in ('.sralite', '.lite.', '.noqual.', '.zq.')) or ftype.lower() in ('sralite', 'noqual'):
            quality = 'SIMPLIFIED_FORMAT_INDICATED_BY_NAME_OR_TYPE'
        elif (source == 'NCBI_NORMALIZED_ONLY' and ftype.lower() == 'sra'
              and not any(t in name.lower() for t in ('pileup', 'realign', 'vdbcache'))):
            category = 'NORMALIZED_ARCHIVE_CANDIDATE'
            quality = 'NORMALIZED_FORMAT_CANDIDATE_CONTENT_NOT_VERIFIED'
        elif nq is False or nq == 'false':
            quality = 'NOQUAL_FALSE_CONTENT_NOT_VERIFIED'
        # meta-only=yes deliberately avoids signed download URLs.
        rows.append(row(source, category, run, filename=name, file_type=ftype,
                        bytes_reported=text_value(f.get('size')),
                        md5_reported=text_value(f.get('md5')),
                        noqual_reported=text_value(nq),
                        rehydration_required_reported=text_value(f.get('rehydrationRequired', bundle.get('rehydrationRequired'))),
                        extra_metadata=json.dumps({k: f[k] for k in ('object', 'modificationDate') if k in f}, sort_keys=True)))
        rows[-1]['quality_metadata'] = quality
    if not rows:
        issues.append(source + ': no file metadata listed in the successful response; no global absence conclusion.')
    return rows, issues


def request_metadata(label: str, url: str, out: Path, events: list[dict]) -> bytes | None:
    """Metadata only. Bounded size, certificate checking, two attempts at most."""
    for attempt in (1, 2):
        event = {'source': label, 'attempt': attempt, 'requested_utc': utc(), 'url': url}
        events.append(event)
        try:
            request = Request(url, headers={'User-Agent': 'KampalaProtozoaQualitySourceCheck/' + VERSION,
                                           'Accept-Encoding': 'identity'})
            with urlopen(request, timeout=TIMEOUT_SECONDS) as response:
                payload = response.read(MAX_RESPONSE_BYTES + 1)
                event['http_status'] = response.status
                event['response_url'] = response.geturl()
            if len(payload) > MAX_RESPONSE_BYTES:
                raise ValueError('Response exceeds the metadata-only size limit; not treated as file metadata.')
            if not payload.strip():
                raise ValueError('Empty metadata response.')
            dest = out / 'raw_metadata' / (label + '.response.txt')
            dest.write_bytes(payload)
            event.update(status='RETRIEVED_NOT_YET_PARSED', bytes_received=len(payload),
                         sha256=hashlib.sha256(payload).hexdigest(), local_response=dest.name)
            return payload
        except HTTPError as exc:
            event.update(status='REQUEST_FAILED', http_status=exc.code, error=str(exc))
            if exc.code not in (408, 429, 500, 502, 503, 504):
                return None
        except (URLError, TimeoutError, socket.timeout, OSError, ValueError) as exc:
            event.update(status='REQUEST_FAILED', error=f'{type(exc).__name__}: {exc}')
        if attempt == 1:
            time.sleep(2)
    return None


def write_tsv(path: Path, columns: list[str], rows: list[dict]):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for item in rows:
            writer.writerow({k: str(item.get(k, '')).replace('\t', ' ').replace('\r', ' ').replace('\n', ' ') for k in columns})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quality-audit', type=Path, required=True,
                        help='Successful step-05 output folder (not the ZIP file)')
    parser.add_argument('--run', default='SRR28434174')
    args = parser.parse_args(argv)
    if not re.fullmatch(r'[SED]RR\d+', args.run):
        parser.error('Invalid run accession.')
    script = Path(__file__).resolve()
    root = script.parent.parent
    if script.parent.name != 'scripts' or not (root / 'metadata').is_dir():
        print('STOP: Put this script in the existing project scripts folder.')
        return 1
    output = None
    record = {'project': PROJECT, 'run_accession': args.run, 'version': VERSION,
              'started_utc': utc(), 'status': 'STARTED', 'python': sys.version,
              'script_sha256': digest(script), 'requests': [],
              'sequence_files_downloaded': False, 'raw_reads_opened': False,
              'quality_scores_modified': False, 'cloud_charges_accepted': False}
    try:
        src = args.quality_audit.expanduser()
        src = (src if src.is_absolute() else root / src).resolve()
        if not src.is_relative_to(root):
            raise ValueError('Quality-audit path must be within this project.')
        src_json = src / 'quality_audit_record.json'
        audited = json.loads(src_json.read_text(encoding='utf-8'))
        if (audited.get('project') != PROJECT or audited.get('run_accession') != args.run
                or audited.get('status') != 'COMPLETE_REVIEW_REQUIRED'):
            raise ValueError('The supplied quality-audit record is not a completed step-05 audit for this run.')
        mates = audited.get('mates', {})
        if not all(mates.get(m, {}).get('quality_state') == 'CONSTANT_Q30_ALL_BASES' for m in ('R1', 'R2')):
            raise ValueError('This step is intended for the verified Q30-only pilot; review a different pattern first.')
        if mates['R1'].get('reads') != mates['R2'].get('reads') or not mates['R1'].get('reads'):
            raise ValueError('The audit record does not establish equal, nonzero paired-read counts.')
        checks = audited.get('input_checks', [])
        if len(checks) != 2 or not all(x.get('archive_md5_matches') is True for x in checks):
            raise ValueError('The quality audit does not contain two successful archive checksum checks.')
        output = root / 'metadata' / f'quality_sources_{args.run}_{stamp()}'
        (output / 'raw_metadata').mkdir(parents=True, exist_ok=False)
        (output / 'provenance').mkdir()
        shutil.copy2(src_json, output / 'provenance' / src_json.name)
        if (src / 'quality_audit_summary.txt').is_file():
            shutil.copy2(src / 'quality_audit_summary.txt', output / 'provenance' / 'quality_audit_summary.txt')
        record.update(source_quality_audit=str(src.relative_to(root)), source_quality_audit_sha256=digest(src_json))
        print(f'CHECKING PUBLIC FILE-SOURCE METADATA: {args.run}', flush=True)
        print('No FASTQ/SRA downloads, no new full-file scans, no trimming, no cloud-account actions.', flush=True)
        print('A file listing is not proof of recoverable original quality scores.', flush=True)
        ena_url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urlencode({
            'accession': args.run, 'result': 'read_run', 'fields': ','.join(ENA_FIELDS), 'format': 'tsv'})
        sdl_base = 'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?'
        jobs = [
            ('ENA_FILE_REPORT', ena_url, parse_ena_report),
            ('ENA_RUN_XML', 'https://www.ebi.ac.uk/ena/browser/api/xml/' + args.run, parse_ena_xml),
            ('NCBI_DEFAULT_METADATA', sdl_base + urlencode({'acc': args.run, 'meta-only': 'yes'}),
             lambda payload, run: parse_sdl(payload, run, 'NCBI_DEFAULT_METADATA')),
            ('NCBI_NORMALIZED_ONLY', sdl_base + urlencode({'acc': args.run, 'meta-only': 'yes', 'capability': 'zqa:r'}),
             lambda payload, run: parse_sdl(payload, run, 'NCBI_NORMALIZED_ONLY')),
        ]
        rows, issues, statuses = [], [], []
        for label, url, parse in jobs:
            print(f'  Requesting {label} ...', flush=True)
            payload = request_metadata(label, url, output, record['requests'])
            if payload is None:
                statuses.append({'source': label, 'status': 'NOT_RETRIEVED', 'note': 'See source_check_record.json. Availability remains unresolved.'})
                issues.append(label + ': request failed; NOT a finding that files or quality scores are absent.')
                continue
            try:
                new_rows, notes = parse(payload, args.run)
                rows.extend(new_rows)
                issues.extend(notes)
                statuses.append({'source': label, 'status': 'METADATA_PARSED', 'note': f'{len(new_rows)} file-metadata rows; download access and qualities not tested.'})
            except (ValueError, TypeError, KeyError, ET.ParseError) as exc:
                statuses.append({'source': label, 'status': 'UNRESOLVED_RESPONSE', 'note': str(exc)})
                issues.append(label + ': ' + str(exc) + '; no absence conclusion.')
        write_tsv(output / 'source_inventory.tsv', COLUMNS, rows)
        write_tsv(output / 'retrieval_status.tsv', ['source', 'status', 'note'], statuses)
        n_ok = sum(s['status'] == 'METADATA_PARSED' for s in statuses)
        submitted = [r for r in rows if r['category'] == 'ENA_SUBMITTED_FILE']
        normalized = [r for r in rows if r['category'] == 'NORMALIZED_ARCHIVE_CANDIDATE']
        simplified = [r for r in rows if r['quality_metadata'].startswith('SIMPLIFIED_')]
        lines = ['QUALITY SOURCE CHECK - METADATA ONLY', f'BioProject: {PROJECT}', f'Run: {args.run}',
                 f'Previously audited paired reads: {mates["R1"]["reads"]:,}',
                 'Previously audited FASTQ qualities: Q30 at every position in both mates.',
                 f'Metadata responses parsed: {n_ok}/4',
                 f'ENA submitted-file locations listed: {len(submitted)}',
                 f'NCBI normalized-format candidates listed by normalized-only request: {len(normalized)}',
                 f'Rows explicitly flagging simplified qualities or a simplified format: {len(simplified)}',
                 'Counts are metadata rows, NOT independent specimens or confirmed accessible files.', '', 'REQUEST RESULTS:']
        lines.extend(f'  {s["source"]}: {s["status"]} - {s["note"]}' for s in statuses)
        lines += ['', 'REPORTED FILES:']
        if not rows:
            lines.append('  None could be inventoried. Availability remains unresolved.')
        for r in rows:
            lines.append(f'  {r["source"]} | {r["filename"]} | type={r["file_type"]} | bytes={r["bytes_reported"]} | quality={r["quality_metadata"]}')
        lines += ['', 'INTERPRETATION:']
        if submitted or normalized:
            lines.append('  An alternative file/format is listed. It must be retrieved and its contents checked before replacing the current reads.')
        else:
            lines.append('  No alternative submitted-file location or normalized candidate was established by these queries. This is not proof of global unavailability.')
        lines += ['  The origin of the current constant-Q30 values is NOT established by these listings.',
                  '  A normalized-format label does NOT restore qualities that were missing or replaced before submission.',
                  '  noqual=true identifies simplified qualities for THAT file, not necessarily the ENA FASTQ already downloaded.',
                  '  No file links were followed; no original-quality sequence data were recovered in this step.',
                  '  No FASTQ files were reread, changed, filtered, uploaded or classified.', '', 'OUTSTANDING ISSUES:']
        lines.extend('  - ' + issue for issue in issues)
        if not issues:
            lines.append('  Contents, access and quality provenance still require review even though metadata requests succeeded.')
        lines += ['', 'NEXT DECISION:',
                  '  Review any listed alternatives before another large download.',
                  '  If informative qualities cannot be obtained, record an explicit protocol amendment before using sequence-only filtering.',
                  f'Reports: {output.relative_to(root)}']
        summary = '\n'.join(lines) + '\n'
        (output / 'quality_source_summary.txt').write_text(summary, encoding='utf-8')
        (output / 'README.txt').write_text(
            'This review bundle contains public file metadata and aggregate audit information, not sequencing reads.\n'
            'No download locations were tested. A metadata failure is not proof of file absence.\n'
            'Raw server responses and the local audit JSON are preserved separately on the computer.\n', encoding='utf-8')
        with zipfile.ZipFile(output / 'quality_source_review.zip', 'w', zipfile.ZIP_DEFLATED) as z:
            for name in ('quality_source_summary.txt', 'source_inventory.tsv', 'retrieval_status.tsv', 'README.txt'):
                z.write(output / name, name)
        record.update(status='METADATA_SAVED_REVIEW_REQUIRED' if n_ok else 'INCOMPLETE_ALL_RESPONSES_UNRESOLVED',
                      parsed_responses=n_ok, issues=issues, source_statuses=statuses)
        print('\n' + summary, flush=True)
        print('For review: quality_source_summary.txt and quality_source_review.zip. Do not send FASTQ files.', flush=True)
        return 0 if n_ok else 2
    except KeyboardInterrupt:
        record.update(status='INTERRUPTED', error='Interrupted; no sequence files were opened or modified.')
        print('\nSTOP: Interrupted. Original sequencing files are untouched.')
        return 130
    except Exception as exc:
        record.update(status='FAILED', error=f'{type(exc).__name__}: {exc}')
        print(f'STOP: {type(exc).__name__}: {exc}')
        return 1
    finally:
        record['finished_utc'] = utc()
        if output is not None:
            (output / 'source_check_record.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
            if record['status'] in ('FAILED', 'INTERRUPTED'):
                (output / 'quality_source_summary.txt').write_text(
                    f'Status: {record["status"]}\n{record.get("error", "")}\nOriginal FASTQ files were not opened.\n', encoding='utf-8')


if __name__ == '__main__':
    raise SystemExit(main())
