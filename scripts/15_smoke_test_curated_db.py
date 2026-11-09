#!/usr/bin/env python3
"""Limited operational self-check, NOT an independent taxonomic validation."""
import argparse
import csv
import gzip
import subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SELECTED = ROOT / 'metadata/taxonomy/curated_v1_selected_assemblies.tsv'

def make_reads(db):
    records = list(csv.DictReader(SELECTED.open(encoding='utf-8'), delimiter='\t'))
    reads = db / 'reference_selfcheck_150bp.fa'
    expectations = db / 'reference_selfcheck_expected.tsv'
    with reads.open('w') as out, expectations.open('w') as meta:
        w = csv.writer(meta, delimiter='\t')
        w.writerow(['read_id','reference_assembly','expected_taxid','organism','role'])
        for r in records:
            acc = r['assembly_accession']
            base = r['ftp_path'].rstrip('/').split('/')[-1]
            archive = db / 'downloads' / f'{base}_genomic.fna.gz'
            selected = None
            sequence = []
            with gzip.open(archive, 'rt', encoding='ascii') as fh:
                for line in fh:
                    if line.startswith('>'):
                        if len(''.join(sequence)) >= 600:
                            selected = ''.join(sequence)
                            break
                        sequence = []
                    else:
                        sequence.append(line.strip().upper())
                if selected is None:
                    selected = ''.join(sequence)
            # Find one window with no ambiguities to establish the classifier runs.
            fragment = None
            for pos in range(0, max(0, len(selected)-150), 20):
                frag = selected[pos:pos+150]
                if len(frag)==150 and all(c in 'ACGT' for c in frag):
                    fragment=frag
                    break
            if fragment:
                read_id = f'selfcheck_{acc}'
                out.write(f'>{read_id}\n{fragment}\n')
                w.writerow([read_id,acc,r['taxid'],r['organism_name'],r['role']])
    return reads,expectations

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db',required=True,type=Path)
    args=p.parse_args()
    db=args.db.resolve()
    for name in ('hash.k2d','opts.k2d','taxo.k2d'):
        if not (db/name).is_file() or (db/name).stat().st_size==0:
            raise RuntimeError('Missing database component: '+name)
    reads,expected=make_reads(db)
    out=db/'reference_selfcheck_kraken_output.tsv'
    report=db/'reference_selfcheck_kraken_report.txt'
    cmd=['kraken2','--db',str(db),'--threads','2','--output',str(out),'--report',str(report),str(reads)]
    subprocess.run(cmd,check=True)
    counts={}
    for line in out.open(encoding='utf-8'):
        fields=line.rstrip('\n').split('\t')
        if len(fields)>=3:
            counts[fields[1]]=(fields[0],fields[2])
    if len(counts) < 35:
        raise RuntimeError(f'Operational self-check generated too few Kraken2 assignments: {len(counts)}')
    with (ROOT/'metadata/taxonomy/curated_v1_reference_selfcheck.tsv').open('w',newline='') as dst:
        w=csv.writer(dst,delimiter='\t')
        w.writerow(['read_id','reference_assembly','expected_taxid','observed_taxid','kraken_assignment','interpretation'])
        for r in csv.DictReader(expected.open(encoding='utf-8'),delimiter='\t'):
            status,taxid=counts.get(r['read_id'],('MISSING','0'))
            check='EXACT_REFERENCE_TAXID' if status=='C' and taxid==r['expected_taxid'] else 'REVIEW_MATCH_OR_UNCLASSIFIED'
            w.writerow([r['read_id'],r['reference_assembly'],r['expected_taxid'],taxid,status,check])
    print(f'OPERATIONAL TEST COMPLETE: {len(counts)} reference-derived reads classified; inspect curated_v1_reference_selfcheck.tsv')
    print('WARNING: This is a same-reference smoke test only. It cannot estimate sensitivity/specificity or confirm wastewater detections.')
