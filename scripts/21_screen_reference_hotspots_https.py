#!/usr/bin/env python3
"""Stage 3C alternative transport: NCBI BLAST HTTPS URL API.

Uses ONLY public reference genome windows already prepared by script 20.
Reuses script 20's candidate selection, TSV interpretation, and ZIP packaging.
Does not access environmental FASTQ files or transmit wastewater reads.

NCBI guidance: no requests more often than every 10 s; poll a given RID
no more often than once per minute. This script polls >=65 s apart.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'scripts/20_screen_reference_hotspots.py'
URL = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
QUERY_IDS = {f'KAMPALA_REFLOCUS_{i:02d}' for i in range(1, 13)}


def fail(message: str) -> None:
    raise RuntimeError('STOP: ' + message)


def load_stage20():
    if not OLD.is_file():
        fail(f'Original Stage 3C script not found: {OLD}')
    spec = importlib.util.spec_from_file_location('stage20_hotspot', OLD)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def check_prepared(old):
    fasta = old.OUT / 'public_reference_hotspots.fna'
    manifest = old.OUT / 'query_manifest.tsv'
    if not fasta.is_file() or not manifest.is_file():
        print('Pre-existing reference windows not found; preparing using original script 20...', flush=True)
        old.prepare()
    if not fasta.is_file() or not manifest.is_file():
        fail('Missing reference window FASTA and manifest after preparation')
    with manifest.open(newline='', encoding='utf-8') as source:
        rows = list(csv.DictReader(source, delimiter='\t'))
    if {r.get('query_id') for r in rows} != QUERY_IDS or len(rows) != 12:
        fail('Reference manifest is not the expected set of 12 public reference windows')
    sequences = {}
    ident = None
    with fasta.open(encoding='ascii') as src:
        for line in src:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                ident = line[1:].split()[0]
                if ident in sequences:
                    fail(f'Duplicate query identifier: {ident}')
                sequences[ident] = ''
            else:
                if ident is None:
                    fail('Malformed reference-window FASTA')
                sequences[ident] += line.upper()
    if set(sequences) != QUERY_IDS:
        fail('Reference-window FASTA has missing or extra query IDs')
    for row in rows:
        sequence = sequences[row['query_id']]
        if len(sequence) != int(row['query_length']):
            fail(f"Length mismatch for {row['query_id']}")
        if hashlib.sha256(sequence.encode('ascii')).hexdigest() != row['sequence_sha256']:
            fail(f"SHA256 mismatch for {row['query_id']}")
        if not set(sequence) <= set('ACGTN'):
            fail(f"Invalid nucleotides in {row['query_id']}")
    return fasta, rows


def curl_request(method: str, fields: dict[str, str], *, timeout: int = 140) -> str:
    # Uses macOS curl, which has already worked against this service on this Mac.
    curl = shutil.which('curl')
    if not curl:
        fail('curl is not available')
    cmd = [curl, '--fail', '--silent', '--show-error', '--location',
           '--connect-timeout', '15', '--max-time', str(timeout)]
    if method == 'POST':
        cmd += ['--request', 'POST']
    elif method == 'GET':
        cmd += ['--get']
    else:
        raise ValueError(method)
    for key, value in fields.items():
        cmd += ['--data-urlencode', f'{key}={value}']
    cmd.append(URL)
    try:
        completed = subprocess.run(cmd, text=True, capture_output=True,
                                   timeout=timeout + 15, check=False)
    except subprocess.TimeoutExpired as exc:
        fail(f'NCBI {method} request timed out: {exc}')
    if completed.returncode != 0:
        fail(f'NCBI HTTPS {method} request failed: {completed.stderr.strip()[-600:]}')
    if not completed.stdout.strip():
        fail(f'NCBI {method} returned an empty response')
    return completed.stdout


def write_json(path: Path, obj) -> None:
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')
    temp.replace(path)


def submission(old, fasta: Path, *, email: str, restart: bool):
    statepath = old.OUT / 'ncbi_blast_https_job.json'
    query_hash = hashlib.sha256(fasta.read_bytes()).hexdigest()
    if statepath.exists() and not restart:
        state = json.loads(statepath.read_text(encoding='utf-8'))
        if state.get('fasta_sha256') != query_hash:
            fail('Saved RID belongs to a different FASTA. Do not reuse it; inspect the saved job.')
        print('Resuming existing NCBI RID:', state['rid'], flush=True)
        return state

    # Submit the same twelve PUBLIC REFERENCE windows in one request.
    data = {
        'CMD': 'Put', 'PROGRAM': 'blastn', 'DATABASE': 'core_nt',
        'QUERY': fasta.read_text(encoding='ascii'),
        'EXPECT': '1e-10', 'FILTER': 'mL', 'HITLIST_SIZE': '100',
        'TOOL': 'kampala_protozoa_reference_screen',
    }
    if email:
        data['EMAIL'] = email
    print('Submitting 12 public reference windows using the NCBI HTTPS BLAST API...', flush=True)
    response = curl_request('POST', data, timeout=180)
    rid_match = re.search(r'\bRID\s*=\s*([A-Za-z0-9-]+)', response)
    rtoe_match = re.search(r'\bRTOE\s*=\s*(\d+)', response)
    if not rid_match:
        fail('NCBI did not return a valid RID: ' + re.sub(r'<[^>]+>', ' ', response)[-700:])
    state = {
        'rid': rid_match.group(1),
        'fasta_sha256': query_hash,
        'submitted_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'submitted_epoch': time.time(),
        'database': 'core_nt', 'program': 'blastn',
        'rtoe_seconds': int(rtoe_match.group(1)) if rtoe_match else None,
        'endpoint': URL, 'format': 'JSON2_S',
        'status': 'SUBMITTED',
    }
    write_json(statepath, state)
    print('NCBI accepted job; RID:', state['rid'], flush=True)
    # Always honor conservative polling interval after submission.
    time.sleep(65)
    return state


def poll(old, state, max_minutes: int):
    if state.get('status') == 'READY':
        return
    started = time.monotonic()
    submitted_epoch = state.get('submitted_epoch')
    if submitted_epoch is not None:
        age = time.time() - float(submitted_epoch)
        if age < 65:
            time.sleep(65 - age)
    while True:
        response = curl_request('GET', {'CMD': 'Get', 'RID': state['rid'],
                                        'FORMAT_OBJECT': 'SearchInfo'}, timeout=90)
        match = re.search(r'\bStatus\s*=\s*(\w+)', response)
        if not match:
            fail('BLAST SearchInfo returned unexpected output: ' + re.sub(r'<[^>]+>', ' ', response)[-500:])
        status = match.group(1).upper()
        state['status'] = status
        write_json(old.OUT / 'ncbi_blast_https_job.json', state)
        print('BLAST status:', status, flush=True)
        if status == 'READY':
            # ThereAreHits=no is legitimate: retain JSON report with empty hits.
            return
        if status in {'FAILED', 'UNKNOWN'}:
            fail(f'BLAST RID {state["rid"]} status {status}; do not treat as negative findings')
        if status != 'WAITING':
            fail(f'Unrecognized BLAST status: {status}')
        if time.monotonic() - started >= max_minutes * 60:
            fail('Polling limit reached. Re-run the same command to resume the saved RID.')
        time.sleep(65)


def blast_searches(payload):
    # NCBI JSON2_S can contain one search, multiple per-query reports, or iterations.
    container = payload.get('BlastOutput2') if isinstance(payload, dict) else None
    if container is None:
        fail('JSON is not a valid BLAST JSON2_S response')
    reports = container if isinstance(container, list) else [container]
    searches = []
    for item in reports:
        if 'error' in item:
            fail('BLAST reported an error: ' + str(item['error']))
        report = item.get('report', {})
        result = report.get('results', {})
        if 'search' in result:
            searches.append(result['search'])
        for iteration in result.get('iterations', []):
            if 'search' in iteration:
                searches.append(iteration['search'])
    if not searches:
        fail('BLAST JSON2_S contains no per-query search records')
    return searches


def convert_json_to_tsv(old, payload, manifest):
    raw = old.OUT / 'reference_hotspots_core_nt_blast.tsv'
    expected = {r['query_id']: int(r['query_length']) for r in manifest}
    seen = set()
    rows = []
    for search in blast_searches(payload):
        title = str(search.get('query_title', '')).split()[0]
        if title not in expected:
            fail(f'Unexpected BLAST query title: {search.get("query_title", "")!r}')
        if title in seen:
            fail(f'Duplicate BLAST query result: {title}')
        seen.add(title)
        if int(search.get('query_len', 0)) != expected[title]:
            fail(f'BLAST reported incorrect query length for {title}')
        for hit in search.get('hits', []):
            descriptions = hit.get('description', [])
            if not descriptions:
                fail(f'BLAST hit missing taxonomic descriptions for {title}')
            hsps = hit.get('hsps', [])
            if not hsps:
                continue
            # One top HSP for each listed taxonomic description, matching prior -max_hsps 1.
            hsp = hsps[0]
            aligned = int(hsp['align_len'])
            if aligned <= 0:
                fail(f'Invalid HSP length in {title}')
            percent_identity = 100 * int(hsp['identity']) / aligned
            for desc in descriptions:
                tax_name = str(desc.get('sciname') or '').strip()
                if not tax_name:
                    fail('NCBI returned a hit without a scientific name; cannot assess taxonomic specificity')
                rows.append({
                    'qseqid': title,
                    'saccver': str(desc.get('accession') or desc.get('id') or ''),
                    'pident': f'{percent_identity:.3f}',
                    'length': str(aligned),
                    'qlen': str(expected[title]),
                    'qstart': str(hsp['query_from']),
                    'qend': str(hsp['query_to']),
                    'evalue': str(hsp['evalue']),
                    'bitscore': str(hsp['bit_score']),
                    'staxids': str(desc.get('taxid', '')),
                    'sscinames': tax_name,
                    'stitle': str(desc.get('title', '')).replace('\t', ' ').replace('\n', ' '),
                })
    if seen != set(expected):
        fail(f'Incomplete BLAST results; missing queries: {sorted(set(expected)-seen)}')
    # Write only after ALL queries have been verified.
    with raw.open('w', encoding='utf-8', newline='') as dest:
        writer = csv.DictWriter(dest, fieldnames=old.BLAST_COLUMNS,
                                delimiter='\t', extrasaction='ignore')
        writer.writerows(rows)
    print(f'Validated {len(seen)} queries; converted {len(rows)} BLAST hit descriptions.', flush=True)
    return raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true', help='Only verify the prepared reference windows')
    parser.add_argument('--run', action='store_true', help='Use NCBI HTTPS URL API on public reference windows')
    parser.add_argument('--email', default='', help='Optional contact email for NCBI requests (not required for local analysis)')
    parser.add_argument('--restart', action='store_true', help='Submit a NEW RID rather than resume; do this only if previous RID is invalid/expired')
    parser.add_argument('--max-minutes', type=int, default=30, help='Maximum number of minutes to poll this run (default: 30)')
    args = parser.parse_args()
    if args.check_only == args.run:
        parser.error('Choose exactly one of --check-only or --run')
    if args.max_minutes < 1:
        parser.error('--max-minutes must be >= 1')
    old = load_stage20()
    fasta, rows = check_prepared(old)
    print('PREFLIGHT PASS:', len(rows), 'previously prepared PUBLIC REFERENCE windows, hashes verified.', flush=True)
    if args.check_only:
        print('Nothing submitted to NCBI.', flush=True)
        return
    # Reuse validated result on retries; keep no duplicate submissions.
    output_json = old.OUT / 'reference_hotspots_core_nt_https.json'
    saved_job = old.OUT / 'ncbi_blast_https_job.json'
    if args.restart:
        output_json.unlink(missing_ok=True)
    if output_json.exists() and not args.restart:
        print('Reusing previously downloaded NCBI HTTPS JSON report; no new BLAST job.', flush=True)
    else:
        state = submission(old, fasta, email=args.email, restart=args.restart)
        if state.get('status') != 'READY':
            poll(old, state, args.max_minutes)
        print('Retrieving NCBI BLAST JSON2_S report...', flush=True)
        data = curl_request('GET', {'CMD':'Get', 'RID':state['rid'],
                                    'FORMAT_TYPE':'JSON2_S'}, timeout=180)
        try:
            payload = json.loads(data)
        except json.JSONDecodeError:
            fail('NCBI returned non-JSON results; no output summary was written. Response: ' + data[:300])
        blast_searches(payload)  # Validate before persisting the JSON.
        output_json.write_text(json.dumps(payload) + '\n', encoding='utf-8')
    payload = json.loads(output_json.read_text(encoding='utf-8'))
    raw = convert_json_to_tsv(old, payload, rows)
    old.interpret(rows, raw)
    methods = old.OUT / 'stage3c_methods.json'
    metadata = json.loads(methods.read_text(encoding='utf-8'))
    metadata['submission'] = {
        'protocol': 'NCBI BLAST Common HTTPS URL API',
        'format': 'JSON2_S',
        'reference_only': True,
        'details': 'Single multi-query BLASTN submission; low-complexity mask; hitlist size 100; expect 1e-10',
        'limitation': 'Top-hit matches on reference windows do not establish the origin of wastewater reads',
    }
    write_json(methods, metadata)
    old.make_review(True)
    import zipfile
    with zipfile.ZipFile(old.REVIEW_ZIP, 'a', compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(saved_job, arcname='ncbi_blast_https_job.json')
        archive.write(raw, arcname='reference_hotspots_core_nt_blast.tsv')
    print('STAGE 3C HTTPS SCREEN FINISHED. These are reference-locus matches, NOT confirmed wastewater detections.', flush=True)
    print('UPLOAD:', old.REVIEW_ZIP, flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, json.JSONDecodeError) as exc:
        print('ERROR:', str(exc), file=sys.stderr)
        sys.exit(1)
