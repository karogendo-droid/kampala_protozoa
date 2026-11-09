#!/usr/bin/env python3
"""Apply the agreed paired-end baseline filters to the verified normalized pilot.

Run from the existing project with this script in scripts/. Requires Python 3.10+,
fastp 1.4.0, FastQC 0.12.1 and MultiQC. No network calls, downloads, taxonomic
selection, read merging, base correction, deduplication or human-read analysis.

Uses ONLY the normalized paired FASTQs recorded by step 07. Original files are
read-only inputs. Reruns create new output directories. Single passing mates are
saved separately; they are NEVER paired with each other or added to the primary
paired denominator. The outputs are QC-filtered reads, NOT verified non-human reads.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import zipfile

PROJECT = 'PRJNA1091428'
RUN = 'SRR28434174'
VERSION = '1.0.0'
FASTP_VERSION = '1.4.0'
FASTQC_VERSION = '0.12.1'
ADAPTER = 'CTGTCTCTTATACACATCT'
CHUNK = 4 * 1024 * 1024
VALID_BASES = b'ACGTNacgtn'
VALID_QUAL = bytes(range(33, 127))
Q20_OR_HIGHER = bytes(range(53, 127))
Q30_OR_HIGHER = bytes(range(63, 127))
SAFE_MODULES = {
    'Basic Statistics', 'Per base sequence quality', 'Per sequence quality scores',
    'Per base sequence content', 'Per sequence GC content', 'Per base N content',
    'Sequence Length Distribution', 'Sequence Duplication Levels', 'Adapter Content',
}
PARAMETERS = {
    'profile': 'proposal_baseline_v1',
    'quality_encoding': 'Phred+33, as exported and audited; not proof of original qualities',
    'adapter_read1': ADAPTER, 'adapter_read2': ADAPTER,
    'adapter_method': 'paired-overlap trimming plus explicit Nextera adapter sequence',
    'adapter_sequence_source': 'Illumina adapter documentation; archived library description is Nextera Flex',
    'cut_front': True, 'cut_tail': True,
    'cut_window_size': 4, 'cut_mean_quality': 20,
    'qualified_quality_phred': 20, 'unqualified_percent_limit': 20,
    'n_base_limit': 5, 'length_required': 75, 'average_qual': 0,
    'polyG_trimming': False, 'polyX_trimming': False,
    'polyG_reason': 'Explicitly off for this baseline; no undeclared homopolymer filter or header-based auto-detection.',
    'fixed_head_crop': 0, 'fixed_tail_crop': 0,
    'overlap_len_require': 30, 'overlap_diff_limit': 5, 'overlap_diff_percent_limit': 20,
    'gap_overlap_trimming': False, 'merging': False, 'base_correction': False,
    'deduplication': False, 'fastp_duplication_estimation': False,
    'low_complexity_filter': False,
    'passing_unpaired_mates': 'Retained separately for sensitivity analysis only',
    'threads': 4, 'compression_level': 4,
}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def stamp() -> str:
    return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(CHUNK), b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def within(root: Path, value: str | Path) -> Path:
    p = Path(value).expanduser()
    p = (p if p.is_absolute() else root / p).resolve()
    if not p.is_relative_to(root):
        raise ValueError('All study input/output paths must remain inside this project.')
    return p


def load_inputs(root: Path, source: Path, run: str) -> tuple[dict, list[dict], int]:
    rec = json.loads((source / 'normalized_record.json').read_text())
    if rec.get('project') != PROJECT or rec.get('run_accession') != run:
        raise ValueError('The normalized record belongs to a different project/run.')
    if rec.get('status') != 'COMPLETE_REVIEW_REQUIRED':
        raise ValueError('Step 07 did not finish successfully; do not use incomplete outputs.')
    comp = rec.get('comparison', {})
    if comp.get('same_sequences_same_order') is not True:
        raise ValueError('The normalized/baseline sequence comparison was not identical.')
    n = comp.get('paired_reads_new')
    if not isinstance(n, int) or isinstance(n, bool) or n <= 0 or n != comp.get('paired_reads_baseline'):
        raise ValueError('Invalid or inconsistent paired-read counts in the normalized record.')
    rows = rec.get('normalized_fastq', [])
    if len(rows) != 2:
        raise ValueError('Expected exactly two normalized FASTQ records.')
    selected = []
    for mate in ('R1', 'R2'):
        matches = [r for r in rows if r.get('mate') == mate]
        if len(matches) != 1:
            raise ValueError('Missing or duplicate mate in the normalized record.')
        r = dict(matches[0])
        p = within(root, r['path'])
        if not p.is_relative_to(root / 'data' / 'raw' / 'normalized_fastq') or p.suffix != '.fastq':
            raise ValueError('Expected normalized uncompressed FASTQ inputs, not the old Q30-only files.')
        if not p.is_file() or p.stat().st_size != r.get('bytes'):
            raise ValueError(f'{mate}: missing input or size differs from the normalized record.')
        if not re.fullmatch(r'[0-9a-f]{64}', r.get('sha256', '')):
            raise ValueError('Missing normalized-input SHA-256 checksum.')
        if comp.get('mates', {}).get(mate, {}).get('within_read_variable', 0) <= 0:
            raise ValueError('Meaningful within-read quality variation was not established for both mates.')
        if p.stat().st_mode & 0o222:
            raise ValueError(f'{mate} input is no longer read-only. Review the file history before proceeding.')
        selected.append({**r, 'absolute_path': p})
    if selected[0]['absolute_path'] == selected[1]['absolute_path']:
        raise ValueError('Read 1 and Read 2 cannot be the same input file.')
    before_qc = within(root, rec['qc_directory']) / 'fastqc'
    if not before_qc.is_dir():
        raise ValueError('The completed normalized FastQC directory was not found.')
    return rec, selected, n


def command(cmd: list[str], name: str, reports: Path, record: dict, emit, timeout=None) -> str:
    log = reports / 'logs' / f'{name}.log'
    entry = {'name': name, 'argv': cmd, 'started_utc': utc(), 'log': str(log)}
    record['commands'].append(entry)
    emit(f'Running {name}:\n  {shlex.join(cmd)}')
    with log.open('w', encoding='utf-8') as f:
        try:
            proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True, timeout=timeout, check=False)
        except BaseException:
            entry['finished_utc'] = utc()
            entry['interrupted_or_failed_to_start'] = True
            raise
    text = log.read_text(encoding='utf-8', errors='replace')
    entry.update(returncode=proc.returncode, finished_utc=utc())
    if text:
        # fastp/FastQC summaries are displayed after the external process finishes.
        emit(text[-16000:].rstrip())
    if proc.returncode != 0:
        raise RuntimeError(f'{name} failed (exit {proc.returncode}). See {log}.')
    return text.strip()


def build_fastp(tool: str, inputs: list[dict], outputs: dict[str, Path], qcdir: Path) -> list[str]:
    return [
        tool, '--in1', str(inputs[0]['absolute_path']), '--in2', str(inputs[1]['absolute_path']),
        '--out1', str(outputs['R1_paired']), '--out2', str(outputs['R2_paired']),
        '--unpaired1', str(outputs['R1_unpaired']), '--unpaired2', str(outputs['R2_unpaired']),
        '--adapter_sequence', ADAPTER, '--adapter_sequence_r2', ADAPTER,
        '--cut_front', '--cut_tail', '--cut_window_size', '4', '--cut_mean_quality', '20',
        '--qualified_quality_phred', '20', '--unqualified_percent_limit', '20',
        '--n_base_limit', '5', '--length_required', '75', '--average_qual', '0',
        '--disable_trim_poly_g', '--dont_eval_duplication',
        '--overlap_len_require', '30', '--overlap_diff_limit', '5', '--overlap_diff_percent_limit', '20',
        '--thread', '4', '--compression', '4', '--dont_overwrite',
        '--json', str(qcdir / 'fastp' / 'fastp.json'),
        '--html', str(qcdir / 'fastp' / 'fastp.html'),
        '--report_title', 'SRR28434174 RES001 - proposal baseline preprocessing',
    ]


def read_fastq(f, label: str, number: int):
    h = f.readline()
    if not h:
        return None
    s, plus, q = f.readline(), f.readline(), f.readline()
    if not s or not plus or not q:
        raise ValueError(f'{label}: incomplete FASTQ record {number}.')
    h, s, plus, q = (x.rstrip(b'\r\n') for x in (h, s, plus, q))
    if not h.startswith(b'@') or not plus.startswith(b'+') or not s or len(s) != len(q):
        raise ValueError(f'{label}: invalid FASTQ structure at record {number}.')
    if s.translate(None, VALID_BASES) or q.translate(None, VALID_QUAL):
        raise ValueError(f'{label}: unexpected base/quality character at record {number}.')
    token = h[1:].split(None, 1)[0]
    if not token:
        raise ValueError(f'{label}: empty read identifier at record {number}.')
    suffix = None
    if token.endswith((b'/1', b'/2')):
        suffix, token = token[-1:], token[:-2]
    return token, suffix, s, q


def stats_init() -> dict:
    return {'reads': 0, 'bases': 0, 'q20_bases': 0, 'q30_bases': 0, 'N_bases': 0,
            'length_min': None, 'length_max': None}


def add_read(stats: dict, r, label: str, mate: int) -> None:
    _, suffix, s, q = r
    if suffix is not None and suffix != str(mate).encode():
        raise ValueError(f'{label}: unexpected mate suffix.')
    n = s.upper().count(b'N')
    low20 = len(q.translate(None, Q20_OR_HIGHER))
    if len(s) < 75 or n > 5 or 100 * low20 > 20 * len(s):
        raise ValueError(f'{label}: an output read fails the configured length/N/quality limits.')
    stats['reads'] += 1
    stats['bases'] += len(s)
    stats['q20_bases'] += len(q) - low20
    stats['q30_bases'] += len(q) - len(q.translate(None, Q30_OR_HIGHER))
    stats['N_bases'] += n
    stats['length_min'] = len(s) if stats['length_min'] is None else min(stats['length_min'], len(s))
    stats['length_max'] = len(s) if stats['length_max'] is None else max(stats['length_max'], len(s))


def finish_stats(stats: dict) -> dict:
    stats['q20_percent'] = 100 * stats['q20_bases'] / stats['bases'] if stats['bases'] else None
    stats['q30_percent'] = 100 * stats['q30_bases'] / stats['bases'] if stats['bases'] else None
    stats['length_mean'] = stats['bases'] / stats['reads'] if stats['reads'] else None
    return stats


def scan_outputs(outputs: dict[str, Path], emit) -> dict:
    stats = {k: stats_init() for k in outputs}
    # Empty unpaired files are permitted; paired files must be readable gzip streams.
    with gzip.open(outputs['R1_paired'], 'rb') as a, gzip.open(outputs['R2_paired'], 'rb') as b:
        i = 0
        while True:
            r1, r2 = read_fastq(a, 'R1_paired', i+1), read_fastq(b, 'R2_paired', i+1)
            if r1 is None and r2 is None:
                break
            if r1 is None or r2 is None or r1[0] != r2[0]:
                raise ValueError('Processed Read 1/Read 2 counts or identifiers do not match.')
            add_read(stats['R1_paired'], r1, 'R1_paired', 1)
            add_read(stats['R2_paired'], r2, 'R2_paired', 2)
            i += 1
            if i % 1000000 == 0:
                emit(f'  Validated {i:,} retained pairs.')
    for mate in (1, 2):
        key = f'R{mate}_unpaired'
        p = outputs[key]
        if not p.is_file():
            raise ValueError(f'Missing unpaired output: {p.name}')
        if p.stat().st_size:
            with gzip.open(p, 'rb') as f:
                while True:
                    r = read_fastq(f, key, stats[key]['reads'] + 1)
                    if r is None:
                        break
                    add_read(stats[key], r, key, mate)
    return {k: finish_stats(v) for k, v in stats.items()}


def check_counts(fastp: dict, counts: dict, original_pairs: int) -> dict:
    summary = fastp.get('summary', {})
    if summary.get('fastp_version') != FASTP_VERSION:
        raise ValueError('fastp JSON version differs from the selected version.')
    if summary.get('before_filtering', {}).get('total_reads') != 2 * original_pairs:
        raise ValueError('fastp input read count differs from the validated normalized files.')
    pairs = counts['R1_paired']['reads']
    if counts['R2_paired']['reads'] != pairs:
        raise ValueError('Paired output counts differ.')
    if summary.get('after_filtering', {}).get('total_reads') != 2 * pairs:
        raise ValueError('fastp paired-read summary differs from the actual retained paired files.')
    if fastp.get('filtering_result', {}).get('passed_filter_reads') != 2 * pairs:
        raise ValueError('fastp filtering result differs from the retained paired-read count.')
    u1, u2 = counts['R1_unpaired']['reads'], counts['R2_unpaired']['reads']
    neither = original_pairs - pairs - u1 - u2
    if neither < 0:
        raise ValueError('Output counts exceed the number of input fragments.')
    return {'input_pairs': original_pairs, 'retained_pairs': pairs,
            'R1_only_pass': u1, 'R2_only_pass': u2, 'neither_pass': neither,
            'retained_pair_percent': 100 * pairs / original_pairs,
            'all_passing_reads': 2 * pairs + u1 + u2,
            'paired_read_bases': counts['R1_paired']['bases'] + counts['R2_paired']['bases']}


def fastqc_tables(path: Path):
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith('/fastqc_data.txt')]
        if len(names) != 1:
            raise ValueError('Unexpected FastQC ZIP layout.')
        text = z.read(names[0]).decode('utf-8')
    sections = {}
    current = None
    for line in text.splitlines():
        if line == '>>END_MODULE':
            current = None
        elif line.startswith('>>'):
            name, flag = line[2:].split('\t', 1)
            sections[name] = {'flag': flag.upper(), 'lines': []}
            current = name
        elif current is not None:
            sections[current]['lines'].append(line)
    basic = {}
    for line in sections.get('Basic Statistics', {}).get('lines', []):
        if '\t' in line and not line.startswith('#'):
            k, v = line.split('\t', 1)
            basic[k] = v
    return sections, basic


def export_qc(qcdir: Path, before_qc: Path, reports: Path, counts: dict, original: int) -> None:
    flags = []
    for stage, folder, names in [
        ('before', before_qc, {f'R{i}': f'{RUN}_{i}_fastqc.zip' for i in (1, 2)}),
        ('after', qcdir / 'fastqc', {f'R{i}': f'{RUN}_R{i}.paired_fastqc.zip' for i in (1, 2)}),
    ]:
        for mate, name in names.items():
            sections, basic = fastqc_tables(folder / name)
            expected = original if stage == 'before' else counts[f'{mate}_paired']['reads']
            if int(basic.get('Total Sequences', '-1').replace(',', '')) != expected:
                raise ValueError('FastQC count does not match the corresponding FASTQ count.')
            for module, sec in sections.items():
                flags.append({'stage': stage, 'mate': mate, 'module': module, 'flag': sec['flag']})
                if module in SAFE_MODULES:
                    safe = re.sub('[^A-Za-z0-9]+', '_', module).strip('_')
                    (reports / 'review_tables' / f'{stage}_{mate}_{safe}.tsv').write_text('\n'.join(sec['lines']) + '\n')
    write_tsv(reports / 'review_tables' / 'fastqc_flags_before_after.tsv', ['stage', 'mate', 'module', 'flag'], flags)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--normalized', required=True, help='Completed step-07 metadata folder')
    ap.add_argument('--run', default=RUN, choices=[RUN])
    args = ap.parse_args()
    script = Path(__file__).resolve()
    root = script.parent.parent
    if script.parent.name != 'scripts' or not (root / 'metadata').is_dir():
        raise ValueError('Copy this script into your existing project scripts folder before running it.')
    source = within(root, args.normalized)
    prior, inputs, original = load_inputs(root, source, args.run)
    tools = {n: shutil.which(n) for n in ('fastp', 'fastqc', 'multiqc')}
    if not all(tools.values()):
        raise ValueError('Activate kampala-qc and install fastp 1.4.0; fastp/FastQC/MultiQC must all be available.')
    if platform.system() == 'Darwin' and platform.machine() != 'arm64':
        raise ValueError('Use the native arm64 environment, not an Intel-emulated Python session.')
    if shutil.disk_usage(root).free < 15 * 1024**3:
        raise ValueError('Less than the conservative 15 GiB working-space reserve. Review storage before proceeding.')
    tag = f'{args.run}_{stamp()}'
    reports = root / 'metadata' / ('preprocess_' + tag)
    qcdir = root / 'results' / 'qc' / 'processed' / tag
    outdir = root / 'data' / 'processed' / 'fastp' / tag
    for p in [reports / 'logs', reports / 'review_tables', reports / 'provenance',
              qcdir / 'fastp', qcdir / 'fastqc', qcdir / 'multiqc', outdir]:
        p.mkdir(parents=True, exist_ok=False)
    record = {
        'project': PROJECT, 'run_accession': args.run, 'sample_alias': prior.get('sample_alias', 'RES001'),
        'status': 'STARTED', 'script_version': VERSION, 'script_sha256': digest(script),
        'started_utc': utc(), 'python': sys.version, 'architecture': platform.machine(),
        'source_folder': str(source.relative_to(root)), 'parameters': PARAMETERS,
        'commands': [], 'versions': {}, 'qc_directory': str(qcdir.relative_to(root)),
        'processed_directory': str(outdir.relative_to(root)),
        'classification_performed': False, 'human_read_screening_performed': False,
        'quality_scores_reassigned': False,
    }
    logfile = reports / 'preprocess_log.txt'
    def emit(s: str) -> None:
        print(s, flush=True)
        with logfile.open('a', encoding='utf-8') as h:
            h.write(s + '\n')
    lock = None
    try:
        emit(f'PREPROCESSING VERIFIED NORMALIZED READS: {args.run} | RES001')
        emit('This stage writes filtered copies; original FASTQ files remain unchanged.')
        emit('No classification, human-read removal, merging, correction or deduplication.')
        lockpath = root / 'data' / 'processed' / 'fastp' / f'{args.run}.lock'
        lock = lockpath.open('a+')
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another preprocessing process is active for this run.')
        for name, tool in tools.items():
            record['versions'][name] = command([tool, '--version'], name + '_version', reports, record, emit, timeout=60)
        if not re.search(r'\bfastp\s+1\.4\.0\b', record['versions']['fastp']):
            raise ValueError('Install the selected fastp version 1.4.0 before running this script.')
        if not re.search(r'FastQC v0\.12\.1\b', record['versions']['fastqc']):
            raise ValueError('FastQC differs from the baseline version 0.12.1; review before comparison.')
        for filename in ('normalized_record.json', 'normalized_summary.txt'):
            if (source / filename).is_file():
                shutil.copy2(source / filename, reports / 'provenance' / filename)
        shutil.copy2(script, reports / 'provenance' / script.name)
        save_json(reports / 'provenance' / 'processing_parameters.json', PARAMETERS)
        note = (
            'BASELINE PROCESSING DECISION\n'
            'Select normalized FASTQ only: sequence identity to ENA was established, and variable stored qualities recovered.\n'
            'Do not combine the normalized and ENA representations; they are the same specimen and reads.\n'
            'Apply proposal thresholds: 4-base end windows at mean Q20; length >=75; <=20% bases below Q20; <=5 Ns.\n'
            'Adapter sequence is the documented Nextera read-through motif, combined with paired-overlap trimming.\n'
            'Poly-G/Poly-X, fixed cropping, low-complexity filtering, merging, base correction and deduplication are OFF.\n'
            'Turning off poly-G auto-detection is explicit, not inferred from archive-renamed headers.\n'
            'Passing single mates are retained separately for later sensitivity analysis.\n'
            'host_removed filenames indicate unresolved prior processing; current processing does not establish absence of human reads.\n'
            'These are operational research thresholds, not a validated parasite-detection assay.\n'
        )
        (reports / 'provenance' / 'processing_decision.txt').write_text(note)
        emit('Checking normalized input checksums ...')
        for r in inputs:
            if digest(r['absolute_path']) != r['sha256']:
                raise ValueError(f'{r["mate"]}: input checksum changed since step 07.')
        record['inputs'] = [{k: v for k, v in r.items() if k != 'absolute_path'} for r in inputs]
        outputs = {f'R{i}_{kind}': outdir / f'{args.run}_R{i}.{kind}.fastq.gz'
                   for kind in ('paired', 'unpaired') for i in (1, 2)}
        cmd = build_fastp(tools['fastp'], inputs, outputs, qcdir)
        (reports / 'provenance' / 'fastp_command.txt').write_text(shlex.join(cmd) + '\n')
        command(cmd, 'fastp', reports, record, emit)
        fpres = json.loads((qcdir / 'fastp' / 'fastp.json').read_text())
        emit('Checking retained paired and single-mate files and counting actual outputs ...')
        counts = scan_outputs(outputs, emit)
        totals = check_counts(fpres, counts, original)
        record.update(output_counts=counts, fragment_accounting=totals)
        rows = [{'file_group': k, **v} for k, v in counts.items()]
        write_tsv(reports / 'review_tables' / 'output_counts.tsv', ['file_group', *stats_init(), 'q20_percent', 'q30_percent', 'length_mean'], rows)
        write_tsv(reports / 'review_tables' / 'fragment_accounting.tsv', ['outcome', 'count'],
                  [{'outcome': k, 'count': v} for k, v in totals.items()])
        # Only scalar numeric values are shared from fastp, not adapter/read sequences.
        frows = []
        for section in ('before_filtering', 'after_filtering'):
            for key, val in fpres.get('summary', {}).get(section, {}).items():
                if isinstance(val, (int, float)) and not isinstance(val, bool):
                    frows.append({'section': section, 'metric': key, 'value': val})
        for key, val in fpres.get('filtering_result', {}).items():
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                frows.append({'section': 'pair_based_filtering_result', 'metric': key, 'value': val})
        for key, val in fpres.get('adapter_cutting', {}).items():
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                frows.append({'section': 'adapter_cutting', 'metric': key, 'value': val})
        write_tsv(reports / 'review_tables' / 'fastp_numeric_summary.tsv', ['section', 'metric', 'value'], frows)
        if totals['retained_pairs'] == 0:
            raise ValueError('No paired reads were retained. Counts were saved; do not proceed to downstream analysis.')
        for p in outputs.values():
            p.chmod(0o444)
        record['outputs'] = [{'group': k, 'path': str(p.relative_to(root)),
                              'bytes': p.stat().st_size, 'sha256': digest(p)} for k, p in outputs.items()]
        command([tools['fastqc'], '--threads', '2', '--outdir', str(qcdir / 'fastqc'),
                 str(outputs['R1_paired']), str(outputs['R2_paired'])], 'postfilter_fastqc', reports, record, emit)
        before_qc = within(root, prior['qc_directory']) / 'fastqc'
        export_qc(qcdir, before_qc, reports, counts, original)
        command([tools['multiqc'], str(before_qc), str(qcdir / 'fastqc'), str(qcdir / 'fastp'),
                 '--module', 'fastqc', '--module', 'fastp', '--outdir', str(qcdir / 'multiqc'),
                 '--filename', 'multiqc_report.html'], 'before_after_multiqc', reports, record, emit)
        if not (qcdir / 'multiqc' / 'multiqc_report.html').is_file():
            raise ValueError('MultiQC did not produce the expected report.')
        emit('Verifying that the normalized input files are unchanged ...')
        for r in inputs:
            if digest(r['absolute_path']) != r['sha256']:
                raise ValueError('A normalized input checksum changed during processing.')
        record['normalized_inputs_unchanged_sha256'] = True
        record['status'] = 'COMPLETE_REVIEW_REQUIRED'
        lines = [
            'PREPROCESSING COMPLETE - POST-FILTER REVIEW REQUIRED', f'BioProject: {PROJECT}',
            f'Run: {args.run}', f'Sample alias: {record["sample_alias"]}',
            f'fastp: {record["versions"]["fastp"]}', 'Normalized inputs unchanged: YES (SHA-256 checked)',
            f'Input read pairs: {original:,}', f'Retained read pairs: {totals["retained_pairs"]:,}',
            f'Paired-fragment retention: {totals["retained_pair_percent"]:.4f}%',
            f'R1-only passing reads (separate): {totals["R1_only_pass"]:,}',
            f'R2-only passing reads (separate): {totals["R2_only_pass"]:,}',
            f'Input pairs with neither mate retained: {totals["neither_pass"]:,}',
            'Output gzip/FASTQ structure and paired identifiers: PASS',
            'Every retained read meets length >=75, N <=5 and <=20% stored bases below Q20: YES', '',
            'RETAINED PAIRED FILES:',
        ]
        for key in ('R1_paired', 'R2_paired'):
            s = counts[key]
            lines += [f'{key}: {s["reads"]:,} reads; length {s["length_min"]}-{s["length_max"]}; '
                      f'Q20 bases {s["q20_percent"]:.3f}%; Q30 bases {s["q30_percent"]:.3f}%']
        lines += ['', 'INTERPRETATION:',
            'Counts describe this preprocessing step, not numbers of organisms or independent specimens.',
            'Read quality summaries reflect stored scores, not externally validated sequencing accuracy.',
            'fastp paired-end failure categories are pair-based; do not interpret them as independent single-read losses.',
            'Single mates are separate sensitivity-analysis material, not two files that can be paired with each other.',
            'No parasite detection, new human-read removal, deduplication, correction or merging was performed.',
            'The prior length/N-only forecast is not a predicted result of adapter and quality processing.',
            'GC/content warnings need context; passing files are not certified human-free or contamination-free.',
            'Baseline poly-G/Poly-X trimming is off; assess any residual sequence artefacts during verification.',
            '', f'Processed files: {outdir.relative_to(root)}',
            f'Before/after QC: {qcdir.relative_to(root)}', f'Reports: {reports.relative_to(root)}']
        (reports / 'preprocess_summary.txt').write_text('\n'.join(lines) + '\n')
        emit('\n' + '\n'.join(lines))
    except KeyboardInterrupt:
        record.update(status='INTERRUPTED', error='Interrupted. Incomplete output files must not be used downstream.')
        emit('\nSTOP: ' + record['error'])
    except Exception as exc:
        record.update(status='FAILED', error=f'{type(exc).__name__}: {exc}')
        emit('STOP: ' + record['error'])
    finally:
        record['finished_utc'] = utc()
        save_json(reports / 'preprocess_record.json', record)
        if record['status'] != 'COMPLETE_REVIEW_REQUIRED':
            (reports / 'preprocess_summary.txt').write_text(
                f'Status: {record["status"]}\nRun: {args.run}\n{record.get("error", "")}\n'
                f'Original reads are not intentionally altered. Do not analyse incomplete outputs.\n'
                f'Reports: {reports.relative_to(root)}\n')
        with zipfile.ZipFile(reports / 'preprocess_review.zip', 'w', zipfile.ZIP_DEFLATED) as z:
            z.write(reports / 'preprocess_summary.txt', 'preprocess_summary.txt')
            for p in sorted((reports / 'review_tables').glob('*.tsv')):
                z.write(p, 'review_tables/' + p.name)
        if lock:
            fcntl.flock(lock, fcntl.LOCK_UN)
            lock.close()
        emit(f'For review: {reports.relative_to(root)}/preprocess_summary.txt and preprocess_review.zip')
    return 0 if record['status'] == 'COMPLETE_REVIEW_REQUIRED' else 2


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'STOP: {type(exc).__name__}: {exc}', file=sys.stderr)
        sys.exit(2)
