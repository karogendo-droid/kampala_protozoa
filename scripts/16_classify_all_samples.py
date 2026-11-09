#!/usr/bin/env python3
"""Stage 1 ONLY: classify six QC-filtered, human-screened Kampala wastewater read pairs.

This script does not calculate species abundances, confirm pathogens, or compare plants.
No raw FASTQs or read-level assignment files are uploaded.
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / 'references/kraken2/kampala_curated_v1'
BASE = ROOT / 'results/taxonomy/kraken2_curated_v1_baseline'
MANIFEST = ROOT / 'metadata/analysis_sample_manifest.csv'
REFERENCE_MANIFEST = ROOT / 'metadata/taxonomy/curated_v1_selected_assemblies.tsv'
REVIEW = ROOT / 'metadata/taxonomy/classification_v1_review.zip'
THREADS = 4
CONFIDENCE = '0.0'  # Sensitive discovery pass; every candidate requires later confirmation.
HIT_GROUPS = 2
RUN_IDS = ('SRR28434174', 'SRR28434173', 'SRR28434172',
           'SRR28434171', 'SRR28434170', 'SRR28434169')


def stop(message):
    raise RuntimeError('STOP: ' + message)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_table(path, delimiter):
    if not path.is_file():
        stop('Missing file ' + str(path))
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def in_project(relative):
    p = (ROOT / relative).resolve()
    if not p.is_relative_to(ROOT.resolve()):
        stop('Recorded file lies outside the project: ' + str(relative))
    return p


def load_samples():
    rows = read_table(MANIFEST, ',')
    by_id = {r['run_accession']: r for r in rows}
    if len(rows) != 6 or set(by_id) != set(RUN_IDS):
        stop('Sample manifest must contain exactly the six expected runs.')
    if len({r['sample_id'] for r in rows}) != 6:
        stop('Sample IDs must be distinct.')
    result = []
    for run in RUN_IDS:
        match = []
        for path in (ROOT / 'metadata').glob(f'human_screen_{run}_*/human_screen_record.json'):
            rec = json.loads(path.read_text(encoding='utf-8'))
            if rec.get('status') == 'COMPLETE_REVIEW_REQUIRED' and rec.get('run_accession') == run:
                match.append((path, rec))
        if len(match) != 1:
            stop(f'Expected ONE complete human-screen record for {run}; found {len(match)}. Review duplicates.')
        path, rec = match[0]
        if rec.get('sample_alias') != by_id[run]['sample_id']:
            stop(f'Sample alias mismatch for {run}.')
        out_list = rec.get('outputs', [])
        if len(out_list) != 2 or {r.get('mate') for r in out_list} != {'R1', 'R2'}:
            stop(f'Human-screen record has invalid paired output listing for {run}.')
        mates = {}
        expected = int(rec['fragment_accounting']['retained_pairs_both_unmapped'])
        if int(rec['output_validation']['pairs']) != expected or expected <= 0:
            stop(f'Invalid read-pair accounting for {run}.')
        for item in out_list:
            f = in_project(item['path'])
            if not f.is_file() or f.stat().st_size != int(item['bytes']):
                stop(f'Missing or changed FASTQ size for {run}: {f}')
            if f.suffix != '.gz' or not f.name.endswith('.human_screened.fastq.gz'):
                stop(f'Unexpected human-screen FASTQ name for {run}: {f}')
            if int(item['reads']) != expected:
                stop(f'Human-screen output record count mismatch for {run}.')
            mates[item['mate']] = (f, item)
        result.append({'run': run, 'sample_id': by_id[run]['sample_id'],
                       'plant': by_id[run]['treatment_plant'], 'pairs': expected,
                       'metadata': path, 'r1': mates['R1'], 'r2': mates['R2']})
    return result


def prepare_db():
    for filename in ('hash.k2d', 'opts.k2d', 'taxo.k2d', 'BUILD_COMPLETE'):
        p = DB / filename
        if not p.is_file() or (filename != 'BUILD_COMPLETE' and p.stat().st_size == 0):
            stop('Missing Kraken2 database component: ' + str(p))
    if not REFERENCE_MANIFEST.is_file():
        stop('Missing curated selection manifest.')
    if not shutil.which('kraken2'):
        stop('Kraken2 not available. Activate Conda environment kampala-taxanomy.')
    version = subprocess.run(['kraken2', '--version'], capture_output=True, text=True, check=True)
    print(version.stdout.strip() or version.stderr.strip(), flush=True)
    return {'hash.k2d': sha256(DB / 'hash.k2d'),
            'opts.k2d': sha256(DB / 'opts.k2d'),
            'taxo.k2d': sha256(DB / 'taxo.k2d'),
            'reference_manifest': sha256(REFERENCE_MANIFEST)}


def report_counts(report):
    taxa = {}
    with report.open(encoding='utf-8') as file:
        for line in file:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                stop('Malformed Kraken2 report line: ' + line[:150])
            taxid = parts[4].strip()
            if not taxid.isdigit():
                stop('Invalid taxonomy identifier in Kraken2 report.')
            taxa[int(taxid)] = {'clade': int(parts[1].strip()),
                               'direct': int(parts[2].strip()),
                               'rank': parts[3].strip(), 'name': parts[5].strip()}
    return taxa


def assignment_counts(path):
    total = classified = 0
    with gzip.open(path, 'rt', encoding='utf-8') as file:
        for line in file:
            parts = line.rstrip('\r\n').split('\t', 4)
            if len(parts) < 5 or parts[0] not in ('C', 'U') or not parts[2].isdigit():
                stop(f'Malformed Kraken2 read-assignment record near line {total + 1}.')
            total += 1
            if parts[0] == 'C':
                classified += 1
    return total, classified


def completed_ok(run_dir, item, db_hashes):
    marker = run_dir / 'completion.json'
    if not marker.exists():
        return False
    x = json.loads(marker.read_text(encoding='utf-8'))
    if x.get('run') != item['run'] or x.get('input_pairs') != item['pairs']:
        stop('Completion metadata does not match input: ' + item['run'])
    if x.get('database_sha256') != db_hashes or x.get('confidence') != CONFIDENCE or x.get('hit_groups') != HIT_GROUPS:
        stop('Database or classification parameters differ for previously completed run ' + item['run'])
    for name, digest in x['outputs_sha256'].items():
        p = run_dir / name
        if not p.is_file() or sha256(p) != digest:
            stop('Output integrity check failed: ' + str(p))
    return True


def classify_one(item, db_hashes):
    run = item['run']
    folder = BASE / run
    if folder.exists():
        if completed_ok(folder, item, db_hashes):
            print(f'REUSED VERIFIED: {run}', flush=True)
            return json.loads((folder / 'completion.json').read_text())
        stop(f'Incomplete previous output exists at {folder}. Do not overwrite it; share its contents for review.')

    for mate in ('r1', 'r2'):
        path, entry = item[mate]
        print(f'Checking stored SHA-256: {run} {mate.upper()}', flush=True)
        if sha256(path) != entry['sha256']:
            stop(f'Input FASTQ changed since human screening: {path}')

    folder.mkdir(parents=True, exist_ok=False)
    out_tmp = folder / 'assignments.tsv.gz.partial'
    report_tmp = folder / 'kraken2.report.txt.partial'
    log = folder / 'kraken2.stderr.log'
    args = ['kraken2', '--db', str(DB), '--threads', str(THREADS),
            '--paired', '--gzip-compressed', '--confidence', CONFIDENCE,
            '--minimum-hit-groups', str(HIT_GROUPS), '--report', str(report_tmp),
            str(item['r1'][0]), str(item['r2'][0])]
    print(f'CLASSIFY: {run} | {item["sample_id"]} | {item["pairs"]:,} pairs', flush=True)
    print('Running: kraken2 --db <curated_v1> --threads 4 --paired --gzip-compressed '
          '--confidence 0.0 --minimum-hit-groups 2 --report <report> <R1> <R2>', flush=True)
    try:
        with log.open('wb') as err, gzip.open(out_tmp, 'wb', compresslevel=3) as gz:
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=err)
            try:
                shutil.copyfileobj(proc.stdout, gz, length=1024 * 1024)
            finally:
                proc.stdout.close()
            code = proc.wait()
        if code:
            stop(f'Kraken2 returned exit code {code} for {run}. Read {log}. Partial files are preserved.')
        if not report_tmp.is_file():
            stop(f'Kraken2 did not create report for {run}.')
        total, classified = assignment_counts(out_tmp)
        taxa = report_counts(report_tmp)
        unclassified = total - classified
        root_classified = taxa.get(1, {}).get('clade', 0)
        report_unclassified = taxa.get(0, {}).get('clade', 0)
        if total != item['pairs'] or classified != root_classified or unclassified != report_unclassified:
            stop(f'Classification counts do not reconcile for {run}: assignments={total}, '
                 f'expected={item["pairs"]}, classified={classified}, root={root_classified}, '
                 f'unclassified={unclassified}, report_unclassified={report_unclassified}.')
        final_out = folder / 'assignments.tsv.gz'
        final_report = folder / 'kraken2.report.txt'
        os.replace(out_tmp, final_out)
        os.replace(report_tmp, final_report)
        record = {'run': run, 'sample_id': item['sample_id'], 'treatment_plant': item['plant'],
                  'input_pairs': total, 'classified_pairs': classified,
                  'unclassified_pairs': unclassified,
                  'classified_percent': round(100 * classified / total, 6),
                  'database_sha256': db_hashes, 'threads': THREADS,
                  'confidence': CONFIDENCE, 'hit_groups': HIT_GROUPS,
                  'input_record': str(item['metadata'].relative_to(ROOT)),
                  'input_fastqs': [str(item[m][0].relative_to(ROOT)) for m in ('r1', 'r2')],
                  'kraken_command': ['kraken2', '--db', '<project>/references/kraken2/kampala_curated_v1',
                                    '--threads', str(THREADS), '--paired', '--gzip-compressed',
                                    '--confidence', CONFIDENCE, '--minimum-hit-groups', str(HIT_GROUPS),
                                    '--report', 'kraken2.report.txt', '<R1>', '<R2>'],
                  'outputs_sha256': {f.name: sha256(f) for f in (final_out, final_report)},
                  'finished_utc': dt.datetime.now(dt.timezone.utc).isoformat()}
        temp_marker = folder / 'completion.json.partial'
        temp_marker.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
        os.replace(temp_marker, folder / 'completion.json')
        print(f'PASS: {run} | classified={classified:,} | unclassified={unclassified:,}', flush=True)
        return record
    except Exception:
        print(f'INCOMPLETE: {run} — files retained in {folder} for review.', file=sys.stderr)
        raise


def consolidate(records, db_hashes):
    summary = BASE / 'classification_summary.tsv'
    with summary.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['run_accession', 'sample_id', 'treatment_plant', 'input_pairs',
                    'classified_pairs', 'unclassified_pairs', 'classified_percent'])
        for rec in records:
            w.writerow([rec['run'], rec['sample_id'], rec['treatment_plant'], rec['input_pairs'],
                        rec['classified_pairs'], rec['unclassified_pairs'], rec['classified_percent']])

    targets = read_table(REFERENCE_MANIFEST, '\t')
    unique_targets = {(r['target_group'], r['taxid'], r['organism_name']) for r in targets if r['role'] == 'target'}
    signal = BASE / 'target_taxid_signals.tsv'
    with signal.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['run_accession', 'sample_id', 'treatment_plant', 'target_group',
                    'reference_organism', 'ncbi_taxid', 'clade_assigned_pairs', 'direct_assigned_pairs'])
        for rec in records:
            found = report_counts(BASE / rec['run'] / 'kraken2.report.txt')
            for group, taxid, organism in sorted(unique_targets):
                x = found.get(int(taxid), {'clade': 0, 'direct': 0})
                w.writerow([rec['run'], rec['sample_id'], rec['treatment_plant'], group,
                            organism, taxid, x['clade'], x['direct']])

    methods = BASE / 'classification_methods.json'
    methods.write_text(json.dumps({'project': 'PRJNA1091428',
                                   'analysis': 'Stage 1 Kraken2 exploratory taxonomic classification',
                                   'note': 'Preliminary assignments only; narrow reference background; '
                                           'do not infer confirmed detection, prevalence or abundance.',
                                   'db': str(DB.relative_to(ROOT)), 'database_sha256': db_hashes,
                                   'threads': THREADS, 'confidence': CONFIDENCE,
                                   'minimum_hit_groups': HIT_GROUPS,
                                   'unit': 'paired-end read fragments (one paired classification per pair)',
                                   'samples': RUN_IDS}, indent=2) + '\n')
    REVIEW.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(REVIEW, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for p in (summary, signal, methods):
            z.write(p, p.name)
        for rec in records:
            report = BASE / rec['run'] / 'kraken2.report.txt'
            z.write(report, f"reports/{rec['run']}_kraken2.report.txt")
    print('\nALL SIX CLASSIFIED AND RECONCILED. Review before stage 2.', flush=True)
    print('Summary:', summary.relative_to(ROOT), flush=True)
    print('Review ZIP:', REVIEW.relative_to(ROOT), flush=True)
    print('Do NOT upload raw FASTQs or read-level assignments.', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check-only', action='store_true', help='Verify inputs and database without running Kraken2.')
    a = p.parse_args()
    os.chdir(ROOT)
    db_hashes = prepare_db()
    items = load_samples()
    free = shutil.disk_usage(ROOT).free / (1024 ** 3)
    print(f'Free disk space: {free:.1f} GiB; minimum reserve: 25 GiB', flush=True)
    if free < 25:
        stop('Insufficient free storage; do not start classification.')
    for item in items:
        print(f"FOUND: {item['run']} {item['sample_id']} {item['plant']}: {item['pairs']:,} screened pairs", flush=True)
    if a.check_only:
        print('PREFLIGHT PASS. No reads classified.', flush=True)
        return
    BASE.mkdir(parents=True, exist_ok=True)
    records = [classify_one(x, db_hashes) for x in items]
    consolidate(records, db_hashes)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(2)
