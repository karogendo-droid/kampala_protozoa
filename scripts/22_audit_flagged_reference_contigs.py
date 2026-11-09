#!/usr/bin/env python3
"""Stage 3D: audit suspect cross-taxon reference contigs in existing Stage 3A BAMs.

Read-only with respect to existing references, FASTQs, BAMs, and previous outputs.
Does NOT confirm organisms; it estimates how dependent three prior signals are
on contigs with strong non-target matches in Stage 3C core_nt results.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
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
BLAST_REVIEW = ROOT / 'metadata/taxonomy/stage3c_reference_specificity_review.zip'
OUTPUT_DIR = ROOT / 'results/confirmation/stage3d_contamination_audit_v1'
OUTPUT_ZIP = ROOT / 'metadata/taxonomy/stage3d_contamination_audit_review.zip'
RUNS = ('SRR28434174', 'SRR28434173', 'SRR28434172',
        'SRR28434171', 'SRR28434170', 'SRR28434169')
TARGET = {
    'GCA_000209485.1': ('Enterocytozoon bieneusi H348', 'Enterocytozoon'),
    'GCA_000963415.1': ('Blastocystis sp. subtype 6', 'Blastocystis'),
    'GCA_002893445.1': ('Cyclospora cayetanensis', 'Cyclospora'),
}
MIN_MAPQ = 20
MIN_ALIGNED_BP = 75
MIN_IDENTITY = 0.95
REF_RE = re.compile(r'^(GC[AF]_\d+\.\d+)_\d+\|kraken:taxid\|\d+$')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def stop(message):
    raise RuntimeError('STOP: ' + message)


def read_tsv_path(path):
    with path.open('rt', encoding='utf-8', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def read_tsv_in_zip(z, name):
    if name not in z.namelist():
        stop(f'Stage 3C archive lacks {name}')
    with z.open(name) as stream, io.TextIOWrapper(stream, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, columns):
    with path.open('wt', encoding='utf-8', newline='') as fh:
        out = csv.DictWriter(fh, fieldnames=columns, extrasaction='ignore', delimiter='\t')
        out.writeheader()
        out.writerows(rows)


def load_stage3c():
    if not BLAST_REVIEW.is_file():
        stop('Missing Stage 3C archive: ' + str(BLAST_REVIEW))
    with zipfile.ZipFile(BLAST_REVIEW) as z:
        manifest = read_tsv_in_zip(z, 'query_manifest.tsv')
        summary = read_tsv_in_zip(z, 'reference_hotspot_blast_summary.tsv')
        job = json.loads(z.read('ncbi_blast_https_job.json'))
    if job.get('status') != 'READY':
        stop('Stage 3C remote BLAST job not recorded as READY')
    if len(manifest) != 12 or len(summary) != 12:
        stop('Expected 12 public reference loci and 12 BLAST result summaries')
    byquery = {r['query_id']: r for r in summary}
    if len(byquery) != 12 or len({r['query_id'] for r in manifest}) != 12:
        stop('Stage 3C query identifiers are duplicated')
    flagged = defaultdict(set)
    rows = []
    for item in manifest:
        query = item['query_id']
        if query not in byquery:
            stop('Missing BLAST results for ' + query)
        hit = byquery[query]
        acc = item['assembly_accession']
        if acc not in TARGET or item['genus'] != TARGET[acc][1]:
            stop('Unexpected target reference in Stage 3C: ' + str(acc))
        if hit['expected_genus'] != item['genus']:
            stop('Stage 3C reference / BLAST genus mismatch for ' + query)
        if not hit['best_scientific_name'].strip():
            stop('BLAST hit lacks scientific name for ' + query)
        pid = float(hit['best_percent_identity'])
        alen = int(hit['best_alignment_bp'])
        if not (pid >= 95 and alen >= 500):
            stop('Strong non-target match criterion not met for ' + query +
                 f' (identity={pid}, aligned={alen}). Reassess manually.')
        if hit['note'] != 'REVIEW_NONMATCHING_NEAR_BEST':
            stop('BLAST summary did not flag a near-best non-target match: ' + query)
        flagged[acc].add(item['reference_contig'])
        rows.append({
            'query_id': query, 'assembly_accession': acc,
            'organism_name': item['organism_name'],
            'reference_contig': item['reference_contig'],
            'queried_0based_start': item['bin_start_0based'],
            'query_length_bp': item['query_length'],
            'blast_best_organism': hit['best_scientific_name'],
            'blast_best_accession': hit['best_subject'],
            'blast_identity_percent': hit['best_percent_identity'],
            'blast_alignment_bp': hit['best_alignment_bp'],
            'action': 'FLAG_ENTIRE_CONTIG_FOR_CONSERVATIVE_SENSITIVITY_AUDIT',
        })
    if any(len(flagged[acc]) != 4 for acc in TARGET):
        stop('Expected four distinct suspect contigs per target assembly')
    return flagged, rows


def aligned_bp(cigar):
    ops = CIGAR_RE.findall(cigar)
    if not ops or ''.join(n + op for n, op in ops) != cigar:
        return 0
    return sum(int(n) for n, op in ops if op in 'M=X')


def reference_bins(cigar, pos0, bin_size=1000):
    refpos = pos0
    bins = set()
    for nstr, op in CIGAR_RE.findall(cigar):
        n = int(nstr)
        if op in 'M=X':
            if n:
                bins.update(range(refpos // bin_size,
                                  (refpos + n - 1) // bin_size + 1))
            refpos += n
        elif op in 'DN':
            refpos += n
    return bins


def decode(parts):
    if len(parts) < 11:
        stop('Malformed SAM alignment')
    flag = int(parts[1])
    if flag & 2308:
        return None
    ref = parts[2]
    m = REF_RE.fullmatch(ref)
    if not m or m.group(1) not in TARGET:
        return None
    mq = int(parts[4])
    length = aligned_bp(parts[5])
    if mq < MIN_MAPQ or length < MIN_ALIGNED_BP:
        return None
    nmvals = [p[5:] for p in parts[11:] if p.startswith('NM:i:')]
    if len(nmvals) != 1:
        return None
    if 1 - int(nmvals[0]) / length < MIN_IDENTITY:
        return None
    mate = 'R1' if flag & 64 else 'R2' if flag & 128 else None
    if mate is None:
        stop('Qualifying primary alignment lacks read mate flag')
    qname = parts[0].split()[0]
    if qname.endswith(('/1', '/2')):
        qname = qname[:-2]
    pos0 = int(parts[3]) - 1
    if pos0 < 0:
        stop('Negative alignment coordinate')
    return m.group(1), ref, qname, mate, bool(flag & 2), reference_bins(parts[5], pos0)


def preflight(flagged):
    if shutil.which('samtools') is None:
        stop('samtools is unavailable. Run: conda activate kampala-align')
    required = []
    for run in RUNS:
        folder = STAGE3A / run
        bam = folder / 'candidate_mapped.sorted.bam'
        metrics = folder / 'competitive_alignment_metrics.tsv'
        if not bam.is_file() or not metrics.is_file():
            stop('Stage 3A BAM or metrics missing for ' + run)
        subprocess.run(['samtools', 'quickcheck', '-v', str(bam)], check=True)
        records = read_tsv_path(metrics)
        byacc = {r['assembly_accession']: r for r in records}
        if not all(acc in byacc and byacc[acc]['role'] == 'target' for acc in flagged):
            stop('Stage 3A target rows missing for ' + run)
        required.append((run, bam, byacc))
    return required


def inspect(run, bam, metrics, flagged):
    all_pairs = defaultdict(set)
    flagged_pairs = defaultdict(set)
    clean_pairs = defaultdict(set)
    concordant = defaultdict(lambda: defaultdict(set))
    clean_concordant = defaultdict(lambda: defaultdict(set))
    mate_counts = Counter()
    clean_contigs = defaultdict(lambda: defaultdict(set))
    clean_bins = defaultdict(lambda: defaultdict(set))
    p = subprocess.Popen(['samtools', 'view', '-F', '2308', str(bam)],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        for line in p.stdout:
            rec = decode(line.rstrip('\n').split('\t'))
            if rec is None:
                continue
            acc, ref, rid, mate, proper, bins = rec
            all_pairs[acc].add(rid)
            mate_counts[(acc, mate)] += 1
            if proper:
                concordant[acc][rid].add(mate)
            if ref in flagged[acc]:
                flagged_pairs[acc].add(rid)
            else:
                clean_pairs[acc].add(rid)
                clean_contigs[acc][ref].add(rid)
                for b in bins:
                    clean_bins[acc][(ref, b)].add(rid)
                if proper:
                    clean_concordant[acc][rid].add(mate)
    finally:
        p.stdout.close()
    stderr = p.stderr.read() if p.stderr else ''
    code = p.wait()
    if code:
        stop('samtools view failed for ' + run + ': ' + stderr[:400])
    rows = []
    for acc in TARGET:
        previous = metrics[acc]
        count = len(all_pairs[acc])
        prevcount = int(previous['high_quality_distinct_pairs'])
        prevconcordant = int(previous['high_quality_concordant_pairs'])
        fullconcordant = sum(s == {'R1', 'R2'} for s in concordant[acc].values())
        if (count != prevcount or fullconcordant != prevconcordant or
            mate_counts[(acc, 'R1')] != int(previous['high_quality_R1']) or
            mate_counts[(acc, 'R2')] != int(previous['high_quality_R2'])):
            stop(f'Stage 3A count mismatch for {run} {acc}; no results written')
        contamination = flagged_pairs[acc]
        only_clean = clean_pairs[acc] - contamination
        both = clean_pairs[acc] & contamination
        clean_concordant_pairs = sum(
            clean_concordant[acc][rid] == {'R1', 'R2'} for rid in only_clean
        )
        retained_refs = {ref for ref, qids in clean_contigs[acc].items() if qids & only_clean}
        retained_bins = {k for k, qids in clean_bins[acc].items() if qids & only_clean}
        rows.append({
            'run_accession': run,
            'sample_id': previous['sample_id'],
            'treatment_plant': previous['treatment_plant'],
            'assembly_accession': acc,
            'organism_name': TARGET[acc][0],
            'stage3a_high_quality_distinct_pairs': count,
            'pairs_with_alignment_on_flagged_contig': len(contamination),
            'pct_pairs_with_flagged_alignment': f'{100*len(contamination)/count:.2f}' if count else '0.00',
            'pairs_with_both_flagged_and_unflagged_alignment': len(both),
            'remaining_pairs_not_touching_flagged_contigs': len(only_clean),
            'remaining_properly_paired_fragments': clean_concordant_pairs,
            'remaining_distinct_contigs': len(retained_refs),
            'remaining_distinct_1kb_bins': len(retained_bins),
            'evidence_status': 'UNCONFIRMED_REQUIRES_BROADER_REFERENCE_VALIDATION',
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check-only', action='store_true')
    mode.add_argument('--run', action='store_true')
    args = parser.parse_args()
    flagged, flags = load_stage3c()
    samples = preflight(flagged)
    print(f'PREFLIGHT PASS: {len(flags)} suspect contigs; {len(samples)} Stage 3A BAMs found.', flush=True)
    if args.check_only:
        print('No alignments scanned and no results written.', flush=True)
        return
    if OUTPUT_DIR.exists() or OUTPUT_ZIP.exists():
        stop('Stage 3D output already exists. Do not overwrite: ' + str(OUTPUT_DIR))
    audit = []
    for run, bam, metrics in samples:
        subset = inspect(run, bam, metrics, flagged)
        audit.extend(subset)
        print('AUDITED ' + run + ': ' + '; '.join(
            r['organism_name'] + ' ' + str(r['remaining_pairs_not_touching_flagged_contigs']) +
            '/' + str(r['stage3a_high_quality_distinct_pairs']) + ' remaining'
            for r in subset), flush=True)
    OUTPUT_DIR.mkdir(parents=True)
    a_path = OUTPUT_DIR / 'affected_species_by_sample.tsv'
    f_path = OUTPUT_DIR / 'flagged_reference_contigs.tsv'
    method_path = OUTPUT_DIR / 'stage3d_methods.json'
    summary_path = OUTPUT_DIR / 'stage3d_summary.txt'
    write_tsv(a_path, audit, list(audit[0].keys()))
    write_tsv(f_path, flags, list(flags[0].keys()))
    method = {
        'stage': '3D reference-contamination sensitivity audit',
        'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'source_blast_review': str(BLAST_REVIEW.relative_to(ROOT)),
        'source_blast_review_sha256': hashlib.sha256(BLAST_REVIEW.read_bytes()).hexdigest(),
        'source_stage3a': str(STAGE3A.relative_to(ROOT)),
        'mapping_filter': {'MAPQ_min': MIN_MAPQ, 'aligned_bp_min': MIN_ALIGNED_BP,
                           'approx_identity_min': MIN_IDENTITY, 'exclude_flag': 2308},
        'selection': 'Use 12 original contigs whose 1kb query windows have BLAST hits to non-target organisms with >=95% nucleotide identity over >=500 bp. Exclude whole contigs conservatively.',
        'meaning': 'Recompute qualifying read-pair support after excluding suspect contigs; not taxonomic confirmation and not a corrected whole-metagenome classification.',
        'limitations': [
            'Other target-reference contigs were not comprehensively screened for contaminants.',
            'Some removed contigs may contain a mixture of true and foreign sequence; whole-contig exclusion is conservative.',
            'The candidate reads were preselected by Stage 1 Kraken2, possibly hiding false negatives.',
            'Alignment uniqueness depends on the limited reference panel; remaining alignments are unconfirmed.',
            'Stage 3C core_nt is not exhaustive and top hits may be affected by database representation.',
            'No conclusion about absence, infection, viability or differences in plant populations can be drawn.'
        ]
    }
    method_path.write_text(json.dumps(method, indent=2) + '\n', encoding='utf-8')
    lines = ['STAGE 3D: impact of suspect reference contigs on three organism signals',
             'All read-pair counts below are Stage 3A high-quality alignment evidence; NOT confirmed detections.', '']
    for acc, (name, _) in TARGET.items():
        sub = [r for r in audit if r['assembly_accession'] == acc]
        original = sum(r['stage3a_high_quality_distinct_pairs'] for r in sub)
        affected = sum(r['pairs_with_alignment_on_flagged_contig'] for r in sub)
        remaining = sum(r['remaining_pairs_not_touching_flagged_contigs'] for r in sub)
        lines.append(f'{name}: original={original}, flagged={affected}, remaining={remaining}, flagged_fraction={affected/original:.1%}' if original else name + ': no qualifying pairs')
    lines += ['', 'IMPORTANT: These are conservative exclusion counts, not proof of reference contamination for every contig.',
              'Even residual alignments must not be interpreted as independent proof of protozoan DNA.',
              'Do not use uncorrected Kraken2 RPMs for Bugolobi-versus-Lubigi biological comparison.']
    summary_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    with zipfile.ZipFile(OUTPUT_ZIP, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for item in (a_path, f_path, method_path, summary_path):
            z.write(item, arcname=item.name)
    print('REVIEW PACKAGE: ' + str(OUTPUT_ZIP.relative_to(ROOT)), flush=True)
    print('Original reads, references, BAMs and earlier results were not changed.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError,
            zipfile.BadZipFile, json.JSONDecodeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)
