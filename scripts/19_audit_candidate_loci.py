#!/usr/bin/env python3
"""Stage 3B: read-location concentration audit of Stage 3A alignments.

This is a distribution/quality audit, NOT independent taxonomic confirmation.
It examines local Stage-3A BAMs and does not download, upload or modify FASTQs.
All six samples are handled in one run. Outputs are summary tables only.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import shutil
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGE3A = ROOT / 'results/confirmation/stage3a_curated_v1'
REVIEW = ROOT / 'metadata/taxonomy/stage3b_locus_review.zip'
OUT = ROOT / 'results/confirmation/stage3b_locus_audit_v1'
RUNS = ('SRR28434174', 'SRR28434173', 'SRR28434172',
        'SRR28434171', 'SRR28434170', 'SRR28434169')
MAPQ_MIN = 20
ALIGNED_BP_MIN = 75
IDENTITY_MIN = 0.95
BIN_SIZE = 1000
REF_RE = re.compile(r'^(GC[AF]_\d+\.\d+)_\d+\|kraken:taxid\|\d+$')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def stop(msg: str):
    raise RuntimeError('STOP: ' + msg)


def rows(path: Path):
    if not path.is_file():
        stop('Missing file: ' + str(path))
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, records, columns):
    with path.open('w', encoding='utf-8', newline='') as f:
        out = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        out.writeheader()
        out.writerows(records)


def clean_id(qname):
    qname = qname.split()[0]
    if qname.endswith(('/1', '/2')):
        return qname[:-2]
    return qname


def aligned_bp(cigar):
    groups = CIGAR_RE.findall(cigar)
    if not groups or ''.join(n + op for n, op in groups) != cigar:
        return 0
    return sum(int(n) for n, op in groups if op in 'M=X')


def ref_blocks(cigar, pos0):
    groups = CIGAR_RE.findall(cigar)
    ref = pos0
    for n, op in groups:
        n = int(n)
        if op in 'M=X':
            yield ref, ref + n
            ref += n
        elif op in 'DN':
            ref += n


def decode_alignment(parts):
    """Use exactly the primary-alignment thresholds reported for Stage 3A."""
    if len(parts) < 11:
        stop('Malformed SAM record')
    flag = int(parts[1])
    if flag & 2308:  # Unmapped, secondary, supplementary
        return None
    ref = parts[2]
    m = REF_RE.fullmatch(ref)
    if not m:
        return None  # mitochondrial decoy or non-selected contig
    accession = m.group(1)
    mq = int(parts[4])
    cigar = parts[5]
    length = aligned_bp(cigar)
    if mq < MAPQ_MIN or length < ALIGNED_BP_MIN:
        return None
    nm_values = [part[5:] for part in parts[11:] if part.startswith('NM:i:')]
    if len(nm_values) != 1:
        return None
    nm = int(nm_values[0])
    if 1 - nm / length < IDENTITY_MIN:
        return None
    qname = clean_id(parts[0])
    mate = 'R1' if flag & 64 else 'R2' if flag & 128 else None
    if mate is None:
        stop('Passing alignment lacks mate flag')
    p0 = int(parts[3]) - 1
    if p0 < 0:
        stop('Negative alignment start')
    bins = set()
    for s, e in ref_blocks(cigar, p0):
        if e > s:
            bins.update(range(s // BIN_SIZE, (e - 1) // BIN_SIZE + 1))
    return accession, ref, qname, mate, p0, bins, bool(flag & 2)


def inspect_sample(run):
    directory = STAGE3A / run
    bam = directory / 'candidate_mapped.sorted.bam'
    metrics = directory / 'competitive_alignment_metrics.tsv'
    if not bam.is_file():
        stop('Missing Stage 3A BAM: ' + str(bam))
    expected_rows = rows(metrics)
    expected = {r['assembly_accession']: r for r in expected_rows if r['role'] == 'target'}
    if len(expected) != 23:
        stop('Expected 23 target reference assemblies in ' + str(metrics))
    if not shutil.which('samtools'):
        stop('SAMtools not installed on PATH. Run: conda activate kampala-align')
    subprocess.run(['samtools', 'quickcheck', '-v', str(bam)], check=True)
    pair_ids = defaultdict(set)
    concordant = defaultdict(lambda: defaultdict(set))
    mate_counts = Counter()
    bins = defaultdict(set)
    contigs = defaultdict(set)
    starts = defaultdict(set)
    proc = subprocess.Popen(['samtools', 'view', '-F', '2308', str(bam)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True)
    try:
        for line in proc.stdout:
            parsed = decode_alignment(line.rstrip('\n').split('\t'))
            if parsed is None:
                continue
            acc, ref, qname, mate, start, windows, proper = parsed
            if acc not in expected:
                continue
            pair_ids[acc].add(qname)
            mate_counts[(acc, mate)] += 1
            if proper:
                concordant[acc][qname].add(mate)
            contigs[(acc, ref)].add(qname)
            starts[(acc, ref, start)].add(qname)
            for window in windows:
                bins[(acc, ref, window)].add(qname)
    finally:
        if proc.stdout:
            proc.stdout.close()
    stderr = proc.stderr.read() if proc.stderr else ''
    code = proc.wait()
    if code:
        stop(f'samtools view failed for {run}: {stderr[:400]}')

    result_rows = []
    hotspot_rows = []
    for acc, old in expected.items():
        count = len(pair_ids[acc])
        confirmed_count = sum(1 for x in concordant[acc].values() if x == {'R1', 'R2'})
        if (count != int(old['high_quality_distinct_pairs']) or
            confirmed_count != int(old['high_quality_concordant_pairs']) or
            mate_counts[(acc, 'R1')] != int(old['high_quality_R1']) or
            mate_counts[(acc, 'R2')] != int(old['high_quality_R2'])):
            stop(f'Stage 3A count reconciliation failed for {run} {acc}')
        if count == 0:
            continue
        bitems = sorted(((key, s) for key, s in bins.items() if key[0] == acc),
                        key=lambda item: (-len(item[1]), item[0][1], item[0][2]))
        citems = sorted(((key, s) for key, s in contigs.items() if key[0] == acc),
                        key=lambda item: (-len(item[1]), item[0][1]))
        top_bins = bitems[:3]
        union_top3 = set().union(*(s for _, s in top_bins)) if top_bins else set()
        top_start = max((len(s) for (k, s) in starts.items() if k[0] == acc), default=0)
        r = {
            'run_accession': run, 'sample_id': old['sample_id'],
            'treatment_plant': old['treatment_plant'],
            'assembly_accession': acc, 'organism_name': old['organism_name'],
            'high_quality_distinct_pairs': count,
            'high_quality_concordant_pairs': confirmed_count,
            'reference_breadth_percent': old['reference_breadth_percent'],
            'total_1kb_bins': len(bitems),
            'bins_with_two_or_more_pairs': sum(len(s) >= 2 for _, s in bitems),
            'top_1kb_bin_pairs': len(bitems[0][1]) if bitems else 0,
            'top_1kb_bin_pair_fraction': f'{len(bitems[0][1])/count:.6f}' if bitems else '0',
            'top_3_bins_union_pairs': len(union_top3),
            'top_3_bins_union_pair_fraction': f'{len(union_top3)/count:.6f}',
            'top_contig_pairs': len(citems[0][1]) if citems else 0,
            'top_contig_pair_fraction': f'{len(citems[0][1])/count:.6f}' if citems else '0',
            'distinct_reference_start_sites': sum(1 for k in starts if k[0] == acc),
            'largest_exact_start_count': top_start,
        }
        result_rows.append(r)
        for (key, pset) in bitems[:10]:
            _, ref, binno = key
            hotspot_rows.append({
                'run_accession': run,
                'sample_id': old['sample_id'],
                'assembly_accession': acc,
                'organism_name': old['organism_name'],
                'reference_contig': ref,
                'bin_start_0based': binno * BIN_SIZE,
                'bin_end_0based_exclusive': (binno + 1) * BIN_SIZE,
                'distinct_pairs_in_bin': len(pset),
                'fraction_of_assembly_pairs': f'{len(pset)/count:.6f}',
            })
    return result_rows, hotspot_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check-only', action='store_true')
    mode.add_argument('--run', action='store_true')
    args = ap.parse_args()
    if not shutil.which('samtools'):
        stop('SAMtools not on PATH. Activate kampala-align.')
    if not (STAGE3A / 'stage3a_methods.json').is_file():
        stop('Stage 3A summary missing; classification not ready for this audit.')
    for run in RUNS:
        for name in ('candidate_mapped.sorted.bam', 'competitive_alignment_metrics.tsv'):
            if not (STAGE3A / run / name).is_file():
                stop(f'Missing input {run}/{name}')
    print('PREFLIGHT PASS: all six Stage 3A BAMs and metrics located.')
    if args.check_only:
        return
    if OUT.exists() or REVIEW.exists():
        stop('Stage 3B output already exists. Do not overwrite; review existing files.')
    OUT.mkdir(parents=True)
    records = []
    hotspots = []
    for run in RUNS:
        print(f'Auditing {run}...', flush=True)
        rr, hh = inspect_sample(run)
        records.extend(rr)
        hotspots.extend(hh)
        print(f'  Reconciled {len(rr)} target references with qualifying alignments', flush=True)
    fields = [
        'run_accession', 'sample_id', 'treatment_plant', 'assembly_accession',
        'organism_name', 'high_quality_distinct_pairs', 'high_quality_concordant_pairs',
        'reference_breadth_percent', 'total_1kb_bins', 'bins_with_two_or_more_pairs',
        'top_1kb_bin_pairs', 'top_1kb_bin_pair_fraction', 'top_3_bins_union_pairs',
        'top_3_bins_union_pair_fraction', 'top_contig_pairs',
        'top_contig_pair_fraction', 'distinct_reference_start_sites',
        'largest_exact_start_count',
    ]
    hotspot_fields = [
        'run_accession', 'sample_id', 'assembly_accession', 'organism_name',
        'reference_contig', 'bin_start_0based', 'bin_end_0based_exclusive',
        'distinct_pairs_in_bin', 'fraction_of_assembly_pairs',
    ]
    summary = OUT / 'locus_concentration_metrics.tsv'
    loci = OUT / 'top_reference_loci.tsv'
    methods = OUT / 'stage3b_locus_methods.json'
    write_tsv(summary, records, fields)
    write_tsv(loci, hotspots, hotspot_fields)
    methods.write_text(json.dumps({
        'stage': '3B locus concentration audit; NOT independent taxonomic confirmation',
        'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'input': 'Stage 3A target-reference BAM primary alignments, all six samples',
        'parameters': {
            'MAPQ_min': MAPQ_MIN,
            'aligned_reference_bases_min': ALIGNED_BP_MIN,
            'identity_min_approx': IDENTITY_MIN,
            'reference_bin_size_bp': BIN_SIZE,
        },
        'counts': {'sample_count': len(RUNS), 'target_assembly_rows_with_evidence': len(records)},
        'interpretation': [
            'Bins use aligned reference intervals; one read pair may contribute to multiple bins.',
            'Top-three-bin numerator is a union of read-pair identifiers (not a sum).',
            'Coverage patterns are not proof of species; high concentration warrants review.',
            'This step uses the same reference sequences as stages 1/3A.',
            'Independent taxonomic confirmation requires alternative marker/reference data and more comprehensive decoys.',
            'No infection/viability or plant-comparison inference.',
        ]
    }, indent=2) + '\n', encoding='utf-8')
    REVIEW.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(REVIEW, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in (summary, loci, methods):
            z.write(p, p.name)
    print('COMPLETE: locus audit reconciled with Stage 3A')
    print('UPLOAD: metadata/taxonomy/stage3b_locus_review.zip')


if __name__ == '__main__':
    try:
        main()
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as e:
        print(str(e), file=sys.stderr)
        raise SystemExit(2)
