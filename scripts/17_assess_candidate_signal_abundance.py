#!/usr/bin/env python3
"""Stage 2: quantify exploratory Kraken2 candidate signals per million read pairs.

No sequence data are read or changed. This does NOT estimate confirmed parasite
abundance, biomass, viability or infection prevalence.
"""
import argparse
import csv
import hashlib
import io
import json
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ZIP = ROOT / 'metadata/taxonomy/classification_v1_review.zip'
DEFAULT_OUTPUT = ROOT / 'results/abundance/stage2_v1'
REQUIRED = {
    'classification_summary.tsv', 'target_taxid_signals.tsv',
    'classification_methods.json'
}
# NCBI genus taxids: include relatives, not just the named target species.
TARGET_GENERA = [
    ('Blastocystis', 12967),
    ('Cryptosporidium', 5806),
    ('Cyclospora', 44417),
    ('Entamoeba', 5758),
    ('Giardia', 5740),
    ('Enterocytozoon', 27971),
    ('Encephalitozoon', 6033),
]


def tsv_from_zip(zf, name):
    data = zf.read(name).decode('utf-8-sig')
    reader = csv.DictReader(io.StringIO(data), delimiter='\t')
    if reader.fieldnames is None:
        raise ValueError(f'Empty or invalid TSV: {name}')
    return list(reader)


def write_tsv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def numeric_count(raw, label):
    try:
        n = int(raw)
    except (TypeError, ValueError) as e:
        raise ValueError(f'Invalid count for {label}: {raw!r}') from e
    if n < 0:
        raise ValueError(f'Negative count for {label}')
    return n


def parse_report(raw_bytes, run):
    taxids = {}
    for line in raw_bytes.decode('utf-8').splitlines():
        cols = line.split('\t')
        if len(cols) != 6:
            raise ValueError(f'{run}: invalid Kraken report line: {line[:120]!r}')
        _, clade, direct, rank, taxid, label = cols
        taxid = int(taxid)
        if taxid in taxids:
            raise ValueError(f'{run}: duplicate taxonomy ID in report: {taxid}')
        clade = numeric_count(clade, f'{run} clade {taxid}')
        direct = numeric_count(direct, f'{run} direct {taxid}')
        if direct > clade:
            raise ValueError(f'{run}: direct count > clade count: {taxid}')
        taxids[taxid] = {'clade': clade, 'direct': direct, 'rank': rank, 'name': label.strip()}
    if 0 not in taxids or 1 not in taxids:
        raise ValueError(f'{run}: report lacks unclassified or root taxon')
    return taxids


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--review-zip', type=Path, default=DEFAULT_ZIP)
    ap.add_argument('--out-dir', type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()
    src = args.review_zip.expanduser().resolve()
    out = args.out_dir.expanduser().resolve()
    if not src.is_file():
        raise ValueError(f'Review ZIP not found: {src}. Locate it with: find metadata -name classification_v1_review.zip')
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Output folder is not empty; refusing to overwrite: {out}')

    with zipfile.ZipFile(src) as zf:
        bad = zf.testzip()
        if bad:
            raise ValueError(f'Damaged ZIP member: {bad}')
        existing = set(zf.namelist())
        missing = REQUIRED - existing
        if missing:
            raise ValueError(f'Missing ZIP members: {sorted(missing)}')
        summary = tsv_from_zip(zf, 'classification_summary.tsv')
        signals = tsv_from_zip(zf, 'target_taxid_signals.tsv')
        methods = json.loads(zf.read('classification_methods.json'))
        if methods.get('analysis') != 'Stage 1 Kraken2 exploratory taxonomic classification':
            raise ValueError('Unexpected or unsupported source classification methods')
        if len(summary) != 6:
            raise ValueError(f'Expected six samples, found {len(summary)}')
        method_runs = set(methods.get('samples', []))

        seen = set()
        sample_table = []
        reports = {}
        sample_lookup = {}
        for row in summary:
            run = row['run_accession']
            if run in seen:
                raise ValueError(f'Duplicate run: {run}')
            seen.add(run)
            name = f'reports/{run}_kraken2.report.txt'
            if name not in existing:
                raise ValueError(f'Missing Kraken2 report: {name}')
            inp = numeric_count(row['input_pairs'], f'{run} input_pairs')
            classified = numeric_count(row['classified_pairs'], f'{run} classified_pairs')
            unclassified = numeric_count(row['unclassified_pairs'], f'{run} unclassified_pairs')
            if inp == 0 or classified + unclassified != inp:
                raise ValueError(f'{run}: classified + unclassified != input pairs')
            percent = 100 * classified / inp
            if abs(percent - float(row['classified_percent'])) > 0.001:
                raise ValueError(f'{run}: inconsistent classified percentage')
            report = parse_report(zf.read(name), run)
            if report[0]['clade'] != unclassified or report[1]['clade'] != classified:
                raise ValueError(f'{run}: report root/unclassified counts do not match summary')
            reports[run] = report
            sample_lookup[run] = (row, inp)
            sample_table.append({
                'run_accession': run, 'sample_id': row['sample_id'],
                'treatment_plant': row['treatment_plant'],
                'input_pairs': inp, 'classified_pairs': classified,
                'unclassified_pairs': unclassified,
                'classified_percent': f'{percent:.6f}',
            })
        if method_runs != seen:
            raise ValueError('Run list does not match classification_methods.json')

        genus_table = []
        for sample in sample_table:
            run = sample['run_accession']
            inp = sample['input_pairs']
            for genus, taxid in TARGET_GENERA:
                record = reports[run].get(taxid)
                if not record:
                    # Not printed in a Kraken report means no reported assignments.
                    clade, direct = 0, 0
                else:
                    if record['rank'] != 'G':
                        raise ValueError(f'Unexpected taxonomy rank for {genus} ({taxid})')
                    clade, direct = record['clade'], record['direct']
                genus_table.append({
                    'run_accession': run, 'sample_id': sample['sample_id'],
                    'treatment_plant': sample['treatment_plant'],
                    'genus': genus, 'ncbi_genus_taxid': taxid,
                    'genus_clade_pairs': clade, 'genus_direct_pairs': direct,
                    'genus_clade_pairs_per_million': f'{clade / inp * 1_000_000:.4f}',
                    'genus_clade_percent_of_all_pairs': f'{clade / inp * 100:.7f}',
                })

        taxon_table = []
        seen_taxa = set()
        for row in signals:
            run = row['run_accession']
            if run not in sample_lookup:
                raise ValueError(f'Taxon entry for unexpected run {run}')
            sample, inp = sample_lookup[run]
            if row['sample_id'] != sample['sample_id'] or row['treatment_plant'] != sample['treatment_plant']:
                raise ValueError(f'{run}: taxon signal sample metadata mismatch')
            taxid = int(row['ncbi_taxid'])
            key = (run, taxid)
            if key in seen_taxa:
                raise ValueError(f'Duplicate candidate taxon entry: {key}')
            seen_taxa.add(key)
            clade = numeric_count(row['clade_assigned_pairs'], f'{run} clade taxid {taxid}')
            direct = numeric_count(row['direct_assigned_pairs'], f'{run} direct taxid {taxid}')
            report = reports[run].get(taxid)
            # Kraken2 reports omit taxa with zero assignments, while the
            # Stage 1 target table may include such taxa with zero counts.
            expected_clade = report['clade'] if report else 0
            expected_direct = report['direct'] if report else 0
            if clade != expected_clade or direct != expected_direct:
                raise ValueError(f'{run}: target taxid {taxid} disagrees with Kraken report')
            taxon_table.append({
                'run_accession': run, 'sample_id': sample['sample_id'],
                'treatment_plant': sample['treatment_plant'],
                'target_group': row['target_group'],
                'reference_organism': row['reference_organism'],
                'ncbi_taxid': taxid,
                'direct_assigned_pairs': direct,
                'clade_assigned_pairs': clade,
                'direct_pairs_per_million': f'{direct / inp * 1_000_000:.4f}',
                'clade_pairs_per_million': f'{clade / inp * 1_000_000:.4f}',
                'direct_percent_of_all_pairs': f'{direct / inp * 100:.7f}',
            })

    # All validations complete. Only now create outputs.
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'sample_classification_metrics.tsv', sample_table,
              ['run_accession','sample_id','treatment_plant','input_pairs','classified_pairs','unclassified_pairs','classified_percent'])
    write_tsv(out / 'genus_clade_signal_rpm.tsv', genus_table,
              ['run_accession','sample_id','treatment_plant','genus','ncbi_genus_taxid','genus_clade_pairs','genus_direct_pairs','genus_clade_pairs_per_million','genus_clade_percent_of_all_pairs'])
    write_tsv(out / 'candidate_taxon_direct_rpm.tsv', taxon_table,
              ['run_accession','sample_id','treatment_plant','target_group','reference_organism','ncbi_taxid','direct_assigned_pairs','clade_assigned_pairs','direct_pairs_per_million','clade_pairs_per_million','direct_percent_of_all_pairs'])
    provenance = {
        'stage': '2: exploratory read-pair-normalized candidate signal assessment',
        'source_zip': str(src),
        'source_zip_sha256': hashlib.sha256(src.read_bytes()).hexdigest(),
        'source_methods': methods,
        'formula': 'Read pairs per million (RPM) = number of Kraken2-assigned pairs / all input human-screened pairs x 1,000,000',
        'denominator': 'Input human-screened paired fragments in Stage 1, not classified pairs',
        'genus_count_type': 'Kraken2 clade counts include genus descendants and competing species in the curated database',
        'taxon_count_type': 'Direct assignments to exact named NCBI taxonomy ID; NOT clade counts',
        'limitations': [
            'Candidate classifier signals only. No confirmed parasite presence or viability.',
            'Not an estimate of biological abundance or prevalence.',
            'Reference database is deliberately small and competitive coverage is incomplete.',
            'Genome lengths and extraction recovery are not corrected by RPM.',
            'Closely related species and conserved/repetitive regions can cause false assignments.',
            'No plant comparison or hypothesis testing has been performed.',
            'Cystoisospora belli lacks genome coverage in the v1 curated database.',
        ],
    }
    (out / 'stage2_methods.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    total_pairs = sum(x['input_pairs'] for x in sample_table)
    classified_pairs = sum(x['classified_pairs'] for x in sample_table)
    lines = [
        '# Stage 2: Preliminary candidate-signal normalization',
        '',
        f'Samples: {len(sample_table)}',
        f'Total read pairs: {total_pairs:,}',
        f'Classified pairs: {classified_pairs:,} ({100 * classified_pairs / total_pairs:.3f}%)',
        '',
        'Every sample is normalized against its total human-screened read pairs.',
        'RPM is a sequencing-depth-adjusted read-assignment rate, **not organism abundance**.',
        'Do not report named species or subtypes as detected until independently confirmed.',
        'No inferential comparison between Bugolobi and Lubigi has been performed.',
        'Cystoisospora belli is not assessable using this database version.',
        '',
        'Output tables:',
        '- `sample_classification_metrics.tsv`: classified/unclassified counts.',
        '- `genus_clade_signal_rpm.tsv`: genus subtree counts (includes relatives).',
        '- `candidate_taxon_direct_rpm.tsv`: exactly assigned named taxids (may include nested taxa).',
    ]
    (out / 'stage2_review.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    review_zip = ROOT / 'metadata/taxonomy/abundance_v1_review.zip'
    if review_zip.exists():
        raise ValueError(f'Existing review ZIP would be overwritten; delete or rename only after checking: {review_zip}')
    review_zip.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(review_zip, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        for filename in ('sample_classification_metrics.tsv', 'genus_clade_signal_rpm.tsv',
                         'candidate_taxon_direct_rpm.tsv', 'stage2_methods.json','stage2_review.txt'):
            zf.write(out / filename, filename)
    print(f'PASS: {len(sample_table)} samples; {total_pairs:,} paired fragments; {classified_pairs:,} assigned to custom database.')
    print('Note: RPM values describe candidate sequence-assignment signals, not confirmed species abundance.')
    print(f'Outputs: {out}')
    print(f'For review: {review_zip}')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError, zipfile.BadZipFile) as exc:
        print(f'STOP: {exc}', file=sys.stderr)
        raise SystemExit(2)
