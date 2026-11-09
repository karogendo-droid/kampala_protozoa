#!/usr/bin/env python3
"""Stage 3F: diagnose zero PR2 target-marker hits using existing BAMs ONLY.

No reads are remapped; existing results are not changed. Counts are READ
ALIGNMENTS (not read pairs), because both mates can align separately.
"""
from __future__ import annotations
import argparse
import collections
import csv
import json
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'confirmation'/'stage3f_pr2_zero_audit_v1'
REVIEW=ROOT/'metadata'/'taxonomy'/'stage3f_pr2_zero_audit_review.zip'
REF=ROOT/'references'/'markers'/'pr2_v5.1.1'
TAX=REF/'pr2_version_5.1.1_SSU_mothur.tax'
BAM_DIR=ROOT/'results'/'confirmation'/'stage3e_pr2_18s_v1'
RUNS=[('RES001','SRR28434174','Bugolobi'),('RES002','SRR28434173','Bugolobi'),
      ('RES003','SRR28434172','Lubigi'),('RES004','SRR28434171','Lubigi'),
      ('RES005','SRR28434170','Bugolobi'),('RES013','SRR28434169','Lubigi')]
TARGETS=['Giardia','Cryptosporidium','Entamoeba','Cyclospora','Cystoisospora',
         'Blastocystis','Enterocytozoon','Encephalitozoon']
CIGAR=re.compile(r'(\d+)([MIDNSHP=X])')

def stop(message): raise SystemExit('STOP: '+message)

def target_genus(taxonomy):
    levels=[v.strip().split('(')[0].strip() for v in taxonomy.rstrip(';').split(';')]
    for target in TARGETS:
        if any(re.match(r'^'+re.escape(target)+r'(?:[0-9]*)(?:$|[_\s])',
                        level,re.I) for level in levels):
            return target
    return None

def read_taxonomy():
    table={}
    with TAX.open(encoding='utf-8',errors='replace') as handle:
        for line in handle:
            cols=line.rstrip('\n\r').split('\t',1)
            if len(cols)==2:
                table[cols[0].split()[0]]=target_genus(cols[1])
    if len(table)<100: stop('PR2 taxonomy file has too few records')
    return table

def alignment_stats(cigar,nm):
    ops=CIGAR.findall(cigar)
    length=sum(int(n) for n,op in ops if op in ('M','=','X','I'))
    return length, 1.0 - (nm/max(length,1))

def audit_bam(bam,taxa):
    counters=collections.defaultdict(collections.Counter)
    total=collections.Counter()
    cmd=['samtools','view','-F','2308',str(bam)]
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,bufsize=1024*1024)
    try:
        for raw in proc.stdout:
            fields=raw.rstrip('\n').split('\t')
            if len(fields)<11: continue
            total['primary_mapped_alignments']+=1
            refname=fields[2]
            if refname not in taxa:
                total['no_taxonomy_record']+=1
                continue
            mapq=int(fields[4]); genus=taxa[refname]
            if mapq>=20: total['primary_mapq20']+=1
            if genus is None:
                if mapq>=20: total['nontarget_mapq20']+=1
                continue
            c=counters[genus]
            c['primary_target_alns_any_mapq']+=1
            if mapq>=20: c['target_alns_mapq20_any_identity']+=1
            nm=None
            for tag in fields[11:]:
                if tag.startswith('NM:i:'):
                    try: nm=int(tag[5:])
                    except ValueError: nm=None
                    break
            if nm is None:
                c['target_alns_missing_NM']+=1
                continue
            alnlen,identity=alignment_stats(fields[5],nm)
            if alnlen>=75:
                c['target_alns_length75']+=1
                if identity>=0.90: c['target_alns_length75_identity90']+=1
                if identity>=0.97:
                    c['target_alns_length75_identity97']+=1
                    if mapq>=20: c['target_alns_full_threshold']+=1
    finally:
        proc.stdout.close()
    stderr=proc.stderr.read()
    code=proc.wait()
    proc.stderr.close()
    if code: stop(f'samtools failed on {bam.name}: {stderr[:600]}')
    return counters,total

def write_tsv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--check-only',action='store_true')
    actions.add_argument('--run',action='store_true')
    args=parser.parse_args()
    if not shutil.which('samtools'): stop('samtools not installed; activate kampala-align')
    if not TAX.is_file(): stop(f'Missing PR2 taxonomy: {TAX}')
    for sample,run,plant in RUNS:
        bam=BAM_DIR/f'{sample}_{run}.pr2_aligned.bam'
        if not bam.is_file() or bam.stat().st_size==0: stop(f'Missing BAM: {bam}')
    print('PREFLIGHT PASS: PR2 taxonomy and six existing BAM files found.')
    if args.check_only: return
    taxa=read_taxonomy()
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]; summary=[]
    for sample,run,plant in RUNS:
        print('AUDITING',sample,run,flush=True)
        bam=BAM_DIR/f'{sample}_{run}.pr2_aligned.bam'
        hits,tot=audit_bam(bam,taxa)
        summary.append(dict(sample_id=sample,run_accession=run,treatment_plant=plant,
           primary_mapped_alignments=tot['primary_mapped_alignments'],
           primary_mapped_alignments_mapq20=tot['primary_mapq20'],
           nontarget_alignments_mapq20=tot['nontarget_mapq20'],
           missing_taxonomy=tot['no_taxonomy_record']))
        for genus in TARGETS:
            c=hits[genus]
            rows.append(dict(sample_id=sample,run_accession=run,treatment_plant=plant,
                genus=genus,
                target_alignments_any_mapq=c['primary_target_alns_any_mapq'],
                target_alignments_mapq20_any_identity=c['target_alns_mapq20_any_identity'],
                target_alignments_length75=c['target_alns_length75'],
                target_alignments_length75_identity90=c['target_alns_length75_identity90'],
                target_alignments_length75_identity97=c['target_alns_length75_identity97'],
                target_alignments_full_threshold=c['target_alns_full_threshold'],
                target_alignments_missing_NM=c['target_alns_missing_NM']))
        print('  Target alignments at any MAPQ:',sum(c['primary_target_alns_any_mapq'] for c in hits.values()),
              '| full threshold:',sum(c['target_alns_full_threshold'] for c in hits.values()),flush=True)
    table=OUT/'pr2_threshold_diagnostic_by_sample_genus.tsv'
    audit=OUT/'pr2_alignment_totals_by_sample.tsv'
    write_tsv(table,rows);write_tsv(audit,summary)
    metadata=OUT/'stage3f_methods.json'
    metadata.write_text(json.dumps(dict(stage='3F PR2 null-result diagnostic',
      created_utc=datetime.now(timezone.utc).isoformat(),
      original_bam_directory=str(BAM_DIR.relative_to(ROOT)),
      source_taxonomy=str(TAX.relative_to(ROOT)),
      unit='primary mapped read alignments, not read pairs',
      note='Counts are nested sensitivity checks, not distinct read-pair totals. This audit does not confirm taxonomic identity. Existing BAMs were not modified.'),indent=2)+'\n')
    REVIEW.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(REVIEW,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for path in (table,audit,metadata): z.write(path,arcname=path.name)
    print('DONE:',REVIEW)

if __name__=='__main__': main()
