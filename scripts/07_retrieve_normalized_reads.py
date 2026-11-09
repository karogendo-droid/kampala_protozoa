#!/usr/bin/env python3
"""Retrieve the explicitly selected SRA Normalized pilot, extract, compare and QC.

Python 3.10+ standard library; SRA Toolkit 3.4.1, FastQC and MultiQC are external.
Run from the existing Kampala project with this file in scripts/.
No Lite fallback, cloud-payment consent, uploads, filtering or taxon calls.
Only SRR28434174 is supported in this pilot-specific script.

The deposited filename 'host_removed' suggests prior processing; this script
cannot establish that history or restore reads removed before submission.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import shutil
import socket
import subprocess
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit, urlunsplit, parse_qsl
from urllib.request import Request, HTTPRedirectHandler, build_opener
import zipfile

PROJECT = 'PRJNA1091428'
RUN = 'SRR28434174'
VERSION = '1.0.0'
TOOLKIT_VERSION = '3.4.1'
CHUNK = 1024 * 1024
MAX_METADATA = 5 * CHUNK
VALID_BASES = b'ACGTRYSWKMBDHVNacgtryswkmbdhvn.'
VALID_QUAL = bytes(range(33, 127))
SAFE_MODULES = {
    'Basic Statistics', 'Per base sequence quality', 'Per sequence quality scores',
    'Per base sequence content', 'Per sequence GC content', 'Per base N content',
    'Sequence Length Distribution', 'Sequence Duplication Levels', 'Adapter Content',
}


def utc():
    return datetime.now(timezone.utc).isoformat()


def stamp():
    return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')


def checksum(path, algorithm='sha256'):
    h = hashlib.md5(usedforsecurity=False) if algorithm == 'md5' else hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4 * CHUNK), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, fields, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: str(row.get(k, '')).replace('\t', ' ').replace('\n', ' ') for k in fields})


def within(root, value):
    p = Path(value).expanduser()
    p = (p if p.is_absolute() else root / p).resolve()
    if not p.is_relative_to(root):
        raise ValueError('Input/output path must remain inside this project.')
    return p


def true(value):
    return value is True or str(value).lower() in ('true', 'yes', '1')


def nonfalse(value):
    return value is not None and str(value).lower() not in ('', 'false', 'no', '0', 'none')


def safe_url(url):
    p = urlsplit(url)
    return urlunsplit((p.scheme, p.hostname or '', p.path, '', ''))


def check_url(url):
    p = urlsplit(url)
    host = (p.hostname or '').lower()
    allowed = (host == 'ncbi.nlm.nih.gov' or host.endswith('.ncbi.nlm.nih.gov')
               or host.endswith('.amazonaws.com') or host == 'storage.googleapis.com'
               or host.endswith('.storage.googleapis.com'))
    if (p.scheme != 'https' or not allowed or p.port not in (None, 443)
            or p.username or p.password or p.fragment):
        raise ValueError('Unexpected/non-HTTPS archive host; no transfer performed.')
    # Never send a requester-pays or accept-charges flag, even via redirect.
    for k, v in parse_qsl(p.query):
        if 'request-payer' in k.lower() or 'accept-charges' in k.lower() or 'userproject' in k.lower():
            raise ValueError('Payment-related query parameter rejected.')
    return url


class SafeRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        check_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


OPENER = build_opener(SafeRedirect())


def open_https(url, headers=None, timeout=90):
    check_url(url)
    h = {'User-Agent': 'KampalaProtozoaNormalizedPilot/' + VERSION,
         'Accept-Encoding': 'identity'}
    h.update(headers or {})
    return OPENER.open(Request(url, headers=h), timeout=timeout)


def load_inventory(source, run):
    rec = json.loads((source / 'source_check_record.json').read_text())
    if rec.get('project') != PROJECT or rec.get('run_accession') != run:
        raise ValueError('Wrong project/run in the source-check record.')
    if rec.get('status') != 'METADATA_SAVED_REVIEW_REQUIRED':
        raise ValueError('Expected a completed step-06 metadata check.')
    rows = read_tsv(source / 'source_inventory.tsv')
    selected = [r for r in rows if r.get('category') == 'NORMALIZED_ARCHIVE_CANDIDATE'
                and r.get('run_accession') == run and r.get('filename') == run
                and r.get('file_type') == 'sra']
    if len(selected) != 1:
        raise ValueError('Expected one explicitly identified normalized SRA candidate.')
    selected = selected[0]
    if true(selected.get('noqual_reported')):
        raise ValueError('Refusing simplified qualities as a normalized candidate.')
    size = int(selected['bytes_reported'])
    md5 = selected['md5_reported'].lower()
    if size <= 0 or not re.fullmatch('[0-9a-f]{32}', md5):
        raise ValueError('Missing normalized-file size/checksum.')
    baseline = []
    for mate in (1, 2):
        name = f'{run}_{mate}.fastq.gz'
        found = [r for r in rows if r.get('category') == 'ENA_GENERATED_FASTQ'
                 and r.get('filename') == name and r.get('run_accession') == run]
        if len(found) != 1:
            raise ValueError('Missing/ambiguous baseline FASTQ metadata.')
        r = found[0]
        if not re.fullmatch('[0-9a-fA-F]{32}', r.get('md5_reported', '')):
            raise ValueError('Missing baseline checksum.')
        baseline.append({'name': name, 'size': int(r['bytes_reported']), 'md5': r['md5_reported'].lower()})
    audit = json.loads((source / 'provenance' / 'quality_audit_record.json').read_text())
    if (audit.get('project') != PROJECT or audit.get('run_accession') != run
            or audit.get('status') != 'COMPLETE_REVIEW_REQUIRED'):
        raise ValueError('Missing matching, completed quality audit.')
    mates = audit.get('mates', {})
    if not all(mates.get(m, {}).get('quality_state') == 'CONSTANT_Q30_ALL_BASES' for m in ('R1', 'R2')):
        raise ValueError('This pilot expects the previously established Q30-only pattern.')
    count = mates['R1']['reads']
    if count <= 0 or count != mates['R2']['reads']:
        raise ValueError('Invalid baseline paired-read count.')
    return {'name': run, 'size': size, 'md5': md5}, baseline, count


def select_normalized(payload, expected, run=RUN):
    obj = json.loads(payload)
    bundles = obj.get('result', []) if isinstance(obj, dict) else obj
    if not isinstance(bundles, list):
        raise ValueError('Unexpected SDL response structure.')
    bundles = [b for b in bundles if isinstance(b, dict)
               and (b.get('bundle') == run or b.get('accession') == run)]
    if len(bundles) != 1 or str(bundles[0].get('status')) != '200':
        raise ValueError('NCBI did not return a successful location response for this run.')
    b = bundles[0]
    if nonfalse(b.get('rehydrationRequired')) or nonfalse(b.get('payRequired')):
        raise ValueError('Restore/payment requirement reported; no action authorized.')
    found = [f for f in b.get('files', []) if isinstance(f, dict)
             and f.get('name') == run and f.get('type') == 'sra'
             and f.get('accession', run) == run]
    if len(found) != 1:
        raise ValueError('No unique exact normalized file returned; Lite fallback is disabled.')
    f = found[0]
    if nonfalse(f.get('noqual', b.get('noqual'))):
        raise ValueError('Normalized candidate flagged as simplified/uncertain quality format.')
    if int(f.get('size', 0)) != expected['size'] or str(f.get('md5', '')).lower() != expected['md5']:
        raise ValueError('Current normalized-file size/checksum differs from the saved inventory. Review required.')
    if nonfalse(f.get('rehydrationRequired')) or nonfalse(f.get('payRequired')) or nonfalse(f.get('ceRequired')):
        raise ValueError('File requires cloud environment, restore or payment; no action taken.')
    links = []
    for loc in f.get('locations', []):
        if not isinstance(loc, dict):
            continue
        if any(nonfalse(loc.get(k)) for k in ('ceRequired', 'payRequired', 'rehydrationRequired', 'encryptedForProjectId')):
            continue
        link = loc.get('link', '')
        try:
            check_url(link)
        except (ValueError, TypeError):
            continue
        if link not in links:
            links.append(link)
    if not links:
        raise ValueError('No direct HTTPS location without cloud/payment requirements was returned. No data downloaded.')
    links.sort(key=lambda s: 0 if (urlsplit(s).hostname or '').endswith('.ncbi.nlm.nih.gov') else 1)
    return links


def get_locations(expected, output):
    url = 'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?' + urlencode({
        'acc': RUN, 'meta-only': 'no', 'capability': 'zqa:r', 'accept-alternate-locations': 'yes'})
    error = None
    for attempt in range(2):
        try:
            with open_https(url, timeout=40) as r:
                data = r.read(MAX_METADATA + 1)
                if len(data) > MAX_METADATA:
                    raise ValueError('Metadata response exceeds safe size limit.')
            p = output / 'private_location_response.json'
            p.write_bytes(data)
            p.chmod(0o600)  # May contain expiring signed links; excluded from review bundle.
            return select_normalized(data, expected)
        except (URLError, TimeoutError, socket.timeout, ConnectionError) as exc:
            error = f'{type(exc).__name__}: location request failed'
            if attempt == 0:
                time.sleep(2)
    raise RuntimeError(error or 'Location request failed; no absence conclusion can be made.')


def verify(path, expected):
    if path.is_symlink():
        raise ValueError('Refusing a symbolic-link input.')
    if path.stat().st_size != expected['size']:
        raise ValueError(f'Size mismatch: {path.name}.')
    observed = checksum(path, 'md5')
    if observed != expected['md5']:
        raise ValueError(f'MD5 mismatch: {path.name}; do not use this file.')
    return {'name': path.name, 'bytes': path.stat().st_size, 'md5': observed}


def download_normalized(links, final, expected, emit=print):
    if final.is_symlink() or final.with_suffix(final.suffix + '.part').is_symlink():
        raise ValueError('Refusing a symbolic-link download path.')
    if final.exists():
        emit('Reusing the existing normalized archive only after size/checksum verification.')
        result = verify(final, expected)
        final.chmod(0o444)
        return result
    part = final.with_suffix(final.suffix + '.part')
    if part.exists() and part.stat().st_size > expected['size']:
        raise ValueError('Partial download exceeds expected size; manual review required.')
    for link in links:
        for attempt in range(1, 4):
            offset = part.stat().st_size if part.exists() else 0
            if offset == expected['size']:
                result = verify(part, expected)
                part.rename(final)
                final.chmod(0o444)
                return {**result, 'name': final.name}
            emit(f'Downloading normalized SRA: {offset / 1e6:.1f} of {expected["size"] / 1e6:.1f} MB present; attempt {attempt}/3.')
            h = {'Range': f'bytes={offset}-'} if offset else {}
            try:
                with open_https(link, h) as r:
                    check_url(r.geturl())
                    if r.headers.get('Content-Encoding', 'identity').lower() not in ('', 'identity'):
                        raise ValueError('Unexpected HTTP compression; exact transfer cannot be checked.')
                    if r.status == 206:
                        m = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', r.headers.get('Content-Range', ''))
                        if not m or int(m[1]) != offset or int(m[3]) != expected['size'] or int(m[2]) != expected['size'] - 1:
                            raise ValueError('Unexpected download range; refusing to append.')
                        mode = 'ab' if offset else 'wb'
                    elif r.status == 200:
                        if offset:
                            emit('Server did not resume; restarting this partial file only.')
                        mode, offset = 'wb', 0
                    else:
                        raise ValueError('Unexpected transfer response status.')
                    clen = r.headers.get('Content-Length')
                    if clen and int(clen) != expected['size'] - offset:
                        raise ValueError('Transfer size disagrees with saved archive metadata.')
                    done, mark = offset, int(100 * offset / expected['size']) // 10
                    with part.open(mode) as f:
                        while chunk := r.read(CHUNK):
                            if done + len(chunk) > expected['size']:
                                raise ValueError('Transfer exceeded expected file size.')
                            f.write(chunk)
                            done += len(chunk)
                            newmark = int(100 * done / expected['size']) // 10
                            if newmark > mark:
                                emit(f'  Archive download: {100 * done / expected["size"]:.0f}%')
                                mark = newmark
                    if done != expected['size']:
                        raise ConnectionError('Incomplete transfer.')
                result = verify(part, expected)
                part.rename(final)
                final.chmod(0o444)
                return {**result, 'name': final.name}
            except HTTPError as exc:
                emit(f'  Download HTTP {exc.code}; no payment or authentication will be attempted.')
                if exc.code in (401, 403, 404, 410):
                    break
                time.sleep(2)
            except (URLError, TimeoutError, socket.timeout, ConnectionError, OSError):
                emit('  Connection interrupted; keeping the partial file.')
                if attempt < 3:
                    time.sleep(2)
    raise RuntimeError('Transfer not completed. Partial file retained. Rerun the same command for fresh links/resume; persistent failures need review.')


def fastq_record(handle, label, number):
    h = handle.readline(1000001)
    if not h:
        return None
    s, plus, q = (handle.readline(1000001) for _ in range(3))
    if not q or any(len(x) > 1000000 for x in (h, s, plus, q)):
        raise ValueError(f'Truncated or oversized FASTQ record in {label}, record {number}.')
    h, s, plus, q = (x.rstrip(b'\r\n') for x in (h, s, plus, q))
    if not h.startswith(b'@') or not plus.startswith(b'+') or not s or len(s) != len(q):
        raise ValueError(f'Invalid FASTQ structure in {label}, record {number}.')
    if s.translate(None, VALID_BASES) or q.translate(None, VALID_QUAL):
        raise ValueError(f'Invalid sequence/quality character in {label}, record {number}.')
    fields = h[1:].split()
    if not fields:
        raise ValueError('Empty FASTQ identifier.')
    key = fields[0]
    suffix = None
    if key.endswith((b'/1', b'/2')):
        suffix, key = int(key[-1:]), key[:-2]
    if len(fields) > 1 and re.match(br'^[12]:[YN]:', fields[1]):
        illumina_mate = int(fields[1][:1])
        if suffix is not None and suffix != illumina_mate:
            raise ValueError('Conflicting mate labels in FASTQ header.')
        suffix = illumina_mate
    if not key:
        raise ValueError('Empty normalized read identifier.')
    return key, suffix, s, q


def check_pair(a, b, label, number):
    if a is None and b is None:
        return False
    if a is None or b is None or a[0] != b[0]:
        raise ValueError(f'Unequal/mismatched paired records in {label} at pair {number}.')
    if a[1] not in (None, 1) or b[1] not in (None, 2):
        raise ValueError(f'Unexpected mate direction in {label} at pair {number}.')
    return True


def scan_and_compare(new_paths, old_paths, expected_count, emit=print):
    stats = [{
        'reads': 0, 'bases': 0, 'within_read_variable': 0, 'all_q30_reads': 0,
        'all_N_reads': 0, 'over_5_N_reads': 0, 'below_75_reads': 0,
        'lengths': Counter(), 'quality_chars': Counter(), 'quality_buffer': bytearray(),
        'sequence_sha256': hashlib.sha256(), 'quality_sha256': hashlib.sha256(),
    } for _ in range(2)]
    mismatches = [0, 0]
    old_count = count = overlap = 0
    forecasts = Counter()
    old_hash = [hashlib.sha256(), hashlib.sha256()]
    # New uncompressed FASTQ; existing ENA input stays compressed and read-only.
    with ExitStack() as stack:
        nh = [stack.enter_context(p.open('rb', buffering=4 * CHUNK)) for p in new_paths]
        oh = [stack.enter_context(gzip.open(p, 'rb')) for p in old_paths]
        index = 0
        while True:
            index += 1
            a = [fastq_record(h, f'normalized R{i+1}', index) for i, h in enumerate(nh)]
            b = [fastq_record(h, f'ENA R{i+1}', index) for i, h in enumerate(oh)]
            hasnew = check_pair(a[0], a[1], 'normalized', index)
            hasold = check_pair(b[0], b[1], 'ENA baseline', index)
            if not hasnew and not hasold:
                break
            if hasold:
                old_count += 1
                for i in range(2):
                    old_hash[i].update(b[i][2] + b'\n')
            if hasnew:
                count += 1
                passes = []
                for i in range(2):
                    s, q = a[i][2], a[i][3]
                    t = stats[i]
                    n = s.count(b'N') + s.count(b'n')
                    t['reads'] += 1
                    t['bases'] += len(s)
                    t['within_read_variable'] += q.count(q[:1]) != len(q)
                    t['all_q30_reads'] += q.count(b'?') == len(q)
                    t['all_N_reads'] += n == len(s)
                    t['over_5_N_reads'] += n > 5
                    t['below_75_reads'] += len(s) < 75
                    t['lengths'][len(s)] += 1
                    t['quality_buffer'].extend(q)
                    if len(t['quality_buffer']) >= CHUNK:
                        t['quality_chars'].update(t['quality_buffer'])
                        t['quality_buffer'].clear()
                    t['sequence_sha256'].update(s + b'\n')
                    t['quality_sha256'].update(q + b'\n')
                    passes.append(len(s) >= 75 and n <= 5)
                label = 'both_pass' if all(passes) else 'R1_only_pass' if passes[0] else 'R2_only_pass' if passes[1] else 'neither_pass'
                forecasts[label] += 1
            if hasnew and hasold:
                overlap += 1
                for i in range(2):
                    mismatches[i] += a[i][2] != b[i][2]
            if index % 1000000 == 0:
                emit(f'  Compared {index:,} paired-record positions.')
    if not count or old_count != expected_count:
        raise ValueError('Empty new files or baseline read count differs from the completed audit.')
    for t in stats:
        t['quality_chars'].update(t.pop('quality_buffer'))
        t['phred33_histogram'] = {str(k - 33): v for k, v in sorted(t.pop('quality_chars').items())}
        t['sequence_sha256'] = t['sequence_sha256'].hexdigest()
        t['quality_sha256'] = t['quality_sha256'].hexdigest()
        t['length_histogram'] = {str(k): v for k, v in sorted(t.pop('lengths').items())}
        scores = sorted(map(int, t['phred33_histogram']))
        t['quality_state'] = ('CONSTANT_Q30_ALL_BASES' if scores == [30] else
                              'ONE_CONSTANT_SCORE_ALL_BASES' if len(scores) == 1 else
                              'WITHIN_READ_VARIATION_PRESENT' if t['within_read_variable'] else
                              'CONSTANT_WITHIN_EACH_READ_ONLY')
    same = count == old_count and not any(mismatches)
    return {
        'paired_reads_new': count, 'paired_reads_baseline': old_count,
        'paired_positions_compared': overlap, 'sequence_mismatches_same_position': dict(zip(('R1', 'R2'), mismatches)),
        'same_sequences_same_order': same,
        'interpretation': 'SAME_SEQUENCE_STRINGS_IN_SAME_ORDER' if same else 'DIFFERENT_COUNTS_OR_ORDERED_SEQUENCES_REQUIRES_REVIEW',
        'mates': dict(zip(('R1', 'R2'), stats)), 'length_N_forecast': dict(forecasts),
        'baseline_sequence_sha256': dict(zip(('R1', 'R2'), [x.hexdigest() for x in old_hash])),
        'headers_compared_between_archives': False,
        'phred33_provenance': 'NCBI FASTQ export convention; variation does not prove untouched instrument qualities',
    }


def parse_fastqc(path):
    with zipfile.ZipFile(path) as z:
        files = [n for n in z.namelist() if n.endswith('/fastqc_data.txt')]
        if len(files) != 1:
            raise ValueError('Unexpected FastQC report structure.')
        text = z.read(files[0]).decode('utf-8')
    sections, current = {}, None
    for line in text.splitlines():
        if line == '>>END_MODULE':
            current = None
        elif line.startswith('>>'):
            bits = line[2:].split('\t')
            if len(bits) != 2:
                raise ValueError('Unexpected FastQC module header.')
            current = bits[0]
            sections[current] = {'status': bits[1].upper(), 'lines': []}
        elif current:
            sections[current]['lines'].append(line)
    basic = dict(line.split('\t', 1) for line in sections.get('Basic Statistics', {}).get('lines', [])
                 if not line.startswith('#') and '\t' in line)
    return sections, basic


def run_command(argv, output, name, record, emit, timeout=None):
    logpath = output / 'logs' / (name + '.log')
    entry = {'argv': argv, 'started_utc': utc(), 'log': str(logpath.name)}
    record['commands'].append(entry)
    emit('Running ' + name + ' ...')
    # No shell interpolation: $si etc. in explicit deflines pass literally.
    with logpath.open('w', encoding='utf-8') as log:
        proc = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            proc.wait(timeout=timeout)
        except BaseException:
            try:
                os.killpg(proc.pid, 15)
                proc.wait(timeout=10)
            except (ProcessLookupError, subprocess.TimeoutExpired):
                try:
                    os.killpg(proc.pid, 9)
                except ProcessLookupError:
                    pass
                proc.wait()
            raise
    entry.update(finished_utc=utc(), returncode=proc.returncode)
    text = logpath.read_text(encoding='utf-8', errors='replace')
    if proc.returncode:
        # Tool log is local. Do not expose any unexpected sequence in the shareable summary.
        raise RuntimeError(f'{name} failed (exit {proc.returncode}); inspect logs/{name}.log.')
    return text.strip()


def execute(args, script):
    root = script.parent.parent.resolve()
    if script.parent.name != 'scripts' or not (root / 'metadata').is_dir():
        raise ValueError('Put this script in the existing project scripts folder.')
    src = within(root, args.sources)
    expected, baseline, expected_count = load_inventory(src, args.run)
    toolkit = args.toolkit.expanduser().resolve()
    bins = {n: toolkit / n for n in ('vdb-validate', 'fasterq-dump')}
    for name, p in bins.items():
        if not p.is_file() or not os.access(p, os.X_OK):
            raise ValueError(f'Missing executable {name}; check the --toolkit folder.')
    qcbins = {n: shutil.which(n) for n in ('fastqc', 'multiqc')}
    if not all(qcbins.values()):
        raise ValueError('Activate kampala-qc; FastQC and MultiQC must be available.')
    if shutil.disk_usage(root).free < 30 * 1024**3:
        raise ValueError('Less than the planning reserve of 30 GiB free. Review space before extraction.')
    oldpaths = [within(root, root / 'data' / 'raw' / args.run / r['name']) for r in baseline]
    if not all(p.is_file() for p in oldpaths):
        raise ValueError('Previously downloaded paired ENA FASTQ files not found in data/raw/SRR28434174.')
    output = root / 'metadata' / f'normalized_{args.run}_{stamp()}'
    (output / 'logs').mkdir(parents=True)
    (output / 'review_tables').mkdir()
    (output / 'provenance').mkdir()
    record = {'project': PROJECT, 'run_accession': args.run, 'sample_alias': 'RES001',
              'version': VERSION, 'script_sha256': checksum(script), 'started_utc': utc(),
              'status': 'STARTED', 'python': sys.version, 'architecture': platform.machine(),
              'source_folder': str(src.relative_to(root)), 'commands': [], 'versions': {},
              'cloud_charges_accepted': False, 'filtering_performed': False,
              'classification_performed': False, 'quality_scores_assigned': False}
    for name in ('source_inventory.tsv', 'quality_source_summary.txt', 'source_check_record.json'):
        shutil.copy2(src / name, output / 'provenance' / name)
    logfile = output / 'normalized_log.txt'
    def emit(s):
        print(s, flush=True)
        with logfile.open('a', encoding='utf-8') as h:
            h.write(s + '\n')
    lock = None
    summary_name = 'normalized_summary.txt'
    try:
        emit(f'RETRIEVING AND TESTING NORMALIZED SRA: {args.run}')
        emit('No Lite fallback, cloud-payment consent, filtering, classification or uploads.')
        emit(f'Expected archive: {expected["size"]:,} bytes ({expected["size"] / 1e9:.3f} GB).')
        for name, p in bins.items():
            ver = run_command([str(p), '--version'], output, name + '_version', record, emit, timeout=30)
            record['versions'][name] = ver
            if TOOLKIT_VERSION not in ver:
                raise ValueError(f'{name} is not the selected version {TOOLKIT_VERSION}.')
        for name, p in qcbins.items():
            record['versions'][name] = run_command([p, '--version'], output, name + '_version', record, emit, timeout=60)
        archdir = root / 'data' / 'raw' / 'ncbi_normalized' / args.run
        archdir.mkdir(parents=True, exist_ok=True)
        lock = (archdir / 'download.lock').open('a+')
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another instance is using this archive. Do not run two copies simultaneously.')
        archive = archdir / (args.run + '.sra')
        if archive.exists():
            record['archive'] = download_normalized([], archive, expected, emit)
        else:
            emit('Requesting current HTTPS locations for the exact normalized file ...')
            links = get_locations(expected, output)
            record['candidate_locations_without_query'] = [safe_url(x) for x in links]
            record['archive'] = download_normalized(links, archive, expected, emit)
        record['archive']['path'] = str(archive.relative_to(root))
        run_command([str(bins['vdb-validate']), str(archive)], output, 'vdb_validate', record, emit)
        # A fresh directory prevents extraction from overwriting any prior attempt.
        tag = output.name.replace('normalized_', '', 1)
        fqdir = root / 'data' / 'raw' / 'normalized_fastq' / tag
        fqdir.mkdir(parents=True, exist_ok=False)
        tmp = fqdir / 'scratch'
        tmp.mkdir()
        record['fastq_directory'] = str(fqdir.relative_to(root))
        argv = [str(bins['fasterq-dump']), str(archdir), '--split-files', '--threads', '4',
                '--outdir', str(fqdir), '--temp', str(tmp),
                '--seq-defline', '@$ac.$si/$ri', '--qual-defline', '+']
        run_command(argv + ['--size-check', 'only'], output, 'fasterq_size_check', record, emit)
        run_command(argv, output, 'fasterq_extract', record, emit)
        newpaths = [fqdir / f'{args.run}_{i}.fastq' for i in (1, 2)]
        if not all(p.is_file() and p.stat().st_size for p in newpaths):
            raise ValueError('Conversion did not produce two nonempty paired FASTQ files.')
        extra = [p for p in fqdir.glob('*.fastq') if p not in newpaths and p.stat().st_size]
        if extra:
            raise ValueError('Unexpected additional FASTQ files; preserve them and review before comparison.')
        emit('Verifying existing ENA baseline files before comparison ...')
        record['baseline_checks'] = [verify(p, spec) for p, spec in zip(oldpaths, baseline)]
        oldstat = [(p.stat().st_size, p.stat().st_mtime_ns) for p in oldpaths]
        emit('Comparing the extracted paired reads with the ENA reads and measuring stored qualities ...')
        comparison = scan_and_compare(newpaths, oldpaths, expected_count, emit)
        record['comparison'] = comparison
        if oldstat != [(p.stat().st_size, p.stat().st_mtime_ns) for p in oldpaths]:
            raise ValueError('Baseline file size/modification time changed during comparison.')
        record['baseline_metadata_unchanged'] = True
        record['normalized_fastq'] = []
        for i, p in enumerate(newpaths):
            record['normalized_fastq'].append({'mate': 'R' + str(i+1), 'path': str(p.relative_to(root)),
                                              'bytes': p.stat().st_size, 'sha256': checksum(p)})
            p.chmod(0o444)
        qcdir = root / 'results' / 'qc' / 'normalized' / tag
        fastqc = qcdir / 'fastqc'
        fastqc.mkdir(parents=True)
        multiqc = qcdir / 'multiqc'
        multiqc.mkdir()
        record['qc_directory'] = str(qcdir.relative_to(root))
        run_command([qcbins['fastqc'], '--threads', '2', '--outdir', str(fastqc), *map(str, newpaths)],
                    output, 'normalized_fastqc', record, emit)
        flags = []
        for i, path in enumerate(newpaths):
            mate = 'R' + str(i+1)
            sections, basic = parse_fastqc(fastqc / f'{args.run}_{i+1}_fastqc.zip')
            if int(basic.get('Total Sequences', '-1').replace(',', '')) != comparison['paired_reads_new']:
                raise ValueError('FastQC read count disagrees with the complete file scan.')
            for mod, section in sections.items():
                flags.append({'mate': mate, 'module': mod, 'flag': section['status']})
                if mod in SAFE_MODULES:
                    name = mate + '_' + re.sub('[^A-Za-z0-9]+', '_', mod).strip('_') + '.tsv'
                    (output / 'review_tables' / name).write_text('\n'.join(section['lines']) + '\n')
        write_tsv(output / 'review_tables' / 'fastqc_flags.tsv', ['mate', 'module', 'flag'], flags)
        run_command([qcbins['multiqc'], str(fastqc), '--outdir', str(multiqc),
                     '--filename', 'multiqc_report.html', '--module', 'fastqc'],
                    output, 'normalized_multiqc', record, emit)
        qrows = [{'mate': m, 'phred33_score': score, 'bases': n}
                 for m, s in comparison['mates'].items() for score, n in s['phred33_histogram'].items()]
        write_tsv(output / 'review_tables' / 'quality_histogram.tsv', ['mate', 'phred33_score', 'bases'], qrows)
        rows = [{'mate': m, 'read_length': k, 'reads': n}
                for m, s in comparison['mates'].items() for k, n in s['length_histogram'].items()]
        write_tsv(output / 'review_tables' / 'length_histogram.tsv', ['mate', 'read_length', 'reads'], rows)
        write_tsv(output / 'review_tables' / 'length_N_forecast.tsv', ['outcome', 'pairs'],
                  [{'outcome': k, 'pairs': comparison['length_N_forecast'].get(k, 0)}
                   for k in ('both_pass', 'R1_only_pass', 'R2_only_pass', 'neither_pass')])
        lines = ['NORMALIZED ARCHIVE ASSESSMENT COMPLETE - REVIEW REQUIRED',
                 f'BioProject: {PROJECT}', f'Run: {args.run}', 'Sample alias: RES001',
                 f'Normalized archive bytes: {expected["size"]:,}',
                 'Normalized archive MD5: MATCH', 'vdb-validate: PASS',
                 f'New read pairs: {comparison["paired_reads_new"]:,}',
                 f'Previously audited ENA pairs: {comparison["paired_reads_baseline"]:,}',
                 f'Ordered sequence comparison: {comparison["interpretation"]}',
                 'Cross-archive read names were not used to assert sequence identity.',
                 'Original ENA files preserved: YES (checksum checked before comparison; file metadata stable).', '']
        for m, s in comparison['mates'].items():
            lines += [m + ':', '  Quality pattern: ' + s['quality_state'],
                      '  Stored Phred+33 scores: ' + ', '.join(s['phred33_histogram']),
                      f'  Reads with within-read score variation: {s["within_read_variable"]:,}',
                      f'  Reads with Q30 at every base: {s["all_q30_reads"]:,}',
                      f'  Entirely-N reads: {s["all_N_reads"]:,}',
                      f'  Below 75 bases: {s["below_75_reads"]:,}',
                      f'  Sequence mismatches at matching positions: {comparison["sequence_mismatches_same_position"][m]:,}', '']
        variable = all(s['within_read_variable'] > 0 for s in comparison['mates'].values())
        lines += ['INTERPRETATION:',
                  ('  Same ordered sequence strings with variable quality information recovered in both mates.'
                   if comparison['same_sequences_same_order'] and variable else
                   '  Do not switch analysis inputs automatically; inspect sequence/quality differences.'),
                  '  Variable scores are not proof that all original instrument qualities survived submission.',
                  '  A normalized SRA file cannot restore pre-submission reads or qualities that were removed.',
                  '  Submitted filenames contain host_removed; the removal method and any other prior filtering remain unverified.',
                  '  FastQC flags are not pathogen detections or automatic filtering decisions.',
                  '  No reads were trimmed, taxonomically selected, or merged with the earlier version.',
                  '  The two versions are NOT independent specimens or extra sequencing depth.', '',
                  f'New uncompressed paired FASTQ: {fqdir.relative_to(root)}',
                  f'FastQC and MultiQC reports: {qcdir.relative_to(root)}',
                  f'Reports: {output.relative_to(root)}']
        (output / summary_name).write_text('\n'.join(lines) + '\n')
        record['status'] = 'COMPLETE_REVIEW_REQUIRED'
        emit('\n' + '\n'.join(lines))
    except KeyboardInterrupt:
        record.update(status='INTERRUPTED', error='Interrupted; partial files and incomplete outputs retained, not approved for analysis.')
        emit('\nSTOP: Interrupted. Existing ENA files have not been intentionally modified.')
    except Exception as exc:
        # Avoid publishing signed URLs embedded in low-level exception text.
        message = re.sub(r'https?://\S+', '[URL omitted]', str(exc))
        record.update(status='FAILED', error=f'{type(exc).__name__}: {message}')
        emit('STOP: ' + record['error'])
    finally:
        if lock:
            lock.close()
        record['finished_utc'] = utc()
        write_json(output / 'normalized_record.json', record)
        if record['status'] != 'COMPLETE_REVIEW_REQUIRED':
            (output / summary_name).write_text(
                f'Status: {record["status"]}\nRun: {args.run}\n{record.get("error", "")}\n'
                'This is not proof that alternative data are globally unavailable.\n'
                'No sequence filtering or pathogen analysis has been performed.\n'
                f'Reports: {output.relative_to(root)}\n')
        # Never share raw sequence files, full FastQC archives, or signed URLs.
        with zipfile.ZipFile(output / 'normalized_review.zip', 'w', zipfile.ZIP_DEFLATED) as z:
            z.write(output / summary_name, summary_name)
            for p in sorted((output / 'review_tables').glob('*.tsv')):
                z.write(p, 'review_tables/' + p.name)
        emit(f'For review: {output.relative_to(root)}/{summary_name} and normalized_review.zip')
    return 0 if record['status'] == 'COMPLETE_REVIEW_REQUIRED' else 2


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sources', type=Path, required=True, help='Completed step-06 folder, not the ZIP')
    ap.add_argument('--toolkit', type=Path, required=True, help='SRA Toolkit 3.4.1 bin directory')
    ap.add_argument('--run', default=RUN, choices=[RUN])
    args = ap.parse_args()
    try:
        return execute(args, Path(__file__).resolve())
    except Exception as exc:
        print('STOP: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
