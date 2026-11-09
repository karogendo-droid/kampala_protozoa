#!/usr/bin/env python3
"""Stage 3E: independent (reference-wise) PR2 SSU-marker screen of all six samples.

This is a marker-level corroboration screen, NOT biological/species confirmation.
Works with plain Python 3 and Bowtie2/SAMtools. No original FASTQs are modified.
"""
from __future__ import annotations
import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = '5.1.1'
REF = ROOT / 'references' / 'markers' / f'pr2_v{VERSION}'
OUT = ROOT / 'results' / 'confirmation' / 'stage3e_pr2_18s_v1'
REVIEW = ROOT / 'metadata' / 'taxonomy' / 'stage3e_pr2_marker_review.zip'
BASE = f'https://github.com/pr2database/pr2database/releases/download/v{VERSION}'
FILES = {
  'fasta': f'pr2_version_{VERSION}_SSU_mothur.fasta.gz',
  'tax': f'pr2_version_{VERSION}_SSU_mothur.tax.gz',
}
RUNS = [
  ('RES001','SRR28434174','Bugolobi'),
  ('RES002','SRR28434173','Bugolobi'),
  ('RES003','SRR28434172','Lubigi'),
  ('RES004','SRR28434171','Lubigi'),
  ('RES005','SRR28434170','Bugolobi'),
  ('RES013','SRR28434169','Lubigi'),
]
TARGETS = ['Giardia','Cryptosporidium','Entamoeba','Cyclospora',
           'Cystoisospora','Blastocystis','Enterocytozoon','Encephalitozoon']
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def stop(msg): raise SystemExit('STOP: ' + msg)


def command_required(name):
  if not shutil.which(name): stop(f'{name} not found in PATH. Activate kampala-align (and install if necessary).')


def inputs():
  root = ROOT/'data'/'processed'/'human_screen'
  if not root.is_dir(): stop(f'input directory missing: {root}')
  found=[]
  for sample,run,plant in RUNS:
    mates=[]
    for side in (1,2):
      pattern=f'{run}_R{side}.human_screened.fastq.gz'
      paths = sorted(root.rglob(pattern))
      if len(paths) != 1:
        stop(f'{sample}: expected exactly one {pattern}, found {len(paths)}. Check duplicate folders before running.')
      if paths[0].stat().st_size == 0: stop(f'empty FASTQ: {paths[0]}')
      mates.append(paths[0])
    found.append((sample,run,plant,*mates))
  return found


def check():
  for p in ('bowtie2','bowtie2-build','samtools','curl'): command_required(p)
  ins=inputs()
  print('PREFLIGHT PASS: six paired screened datasets located; tools available.')
  for sample,run,plant,r1,r2 in ins:
    print(f'  {sample} {plant}: {r1.name}, {r2.name}')
  print('No downloads or alignments performed.')
  return ins


def downloaded(path):
  if not path.exists() or path.stat().st_size < 1000: return False
  try:
    with gzip.open(path,'rb') as f: return bool(f.read(1024))
  except (OSError,EOFError): return False


def download_ref():
  REF.mkdir(parents=True,exist_ok=True)
  for key,name in FILES.items():
    path=REF/name
    if downloaded(path):
      print('REUSING', path)
      continue
    print('DOWNLOADING PR2', key, name, flush=True)
    temp=path.with_suffix(path.suffix+'.partial')
    url=f'{BASE}/{name}'
    p=subprocess.run(['curl','-fL','--retry','3','--connect-timeout','30',
                      '--max-time','1800','-o',str(temp),url], check=False)
    if p.returncode: stop(f'Could not download {url}; retry when connected (HTTP/network error).')
    if not downloaded(temp): stop(f'{name} was not a valid gzip. Check the PR2 release URL.')
    temp.replace(path)
  return REF/FILES['fasta'],REF/FILES['tax']


def unpack(gz_path):
  dest = gz_path.with_suffix('')
  if dest.exists() and dest.stat().st_size > 1000:
    return dest
  print('DECOMPRESSING',gz_path.name,flush=True)
  tmp=dest.with_suffix(dest.suffix+'.partial')
  with gzip.open(gz_path,'rb') as src, tmp.open('wb') as dst:
    shutil.copyfileobj(src,dst,1024*1024)
  tmp.replace(dest)
  return dest


def index_ready(prefix):
  return all((Path(str(prefix)+f'.{i}.bt2')).exists() for i in (1,2,3,4)) or all(
      (Path(str(prefix)+f'.{i}.bt2l')).exists() for i in (1,2,3,4))


def build_index(fasta, threads):
  index=REF/'pr2_ssu'
  if not index_ready(index):
    print('BUILDING BOWTIE2 PR2 INDEX',flush=True)
    p=subprocess.run(['bowtie2-build','--threads',str(threads),str(fasta),str(index)],
                     stdout=subprocess.DEVNULL,check=False)
    if p.returncode: stop('Bowtie2 index build failed')
  return index


def read_taxonomies(taxfile):
  table={}
  with taxfile.open(encoding='utf-8', errors='replace') as handle:
    for line in handle:
      cols=line.rstrip('\n\r').split('\t')
      if len(cols) >= 2:
        table[cols[0].split()[0]]=cols[1]
  if len(table)<100: stop('PR2 taxonomy table unexpectedly small or malformed')
  return table


def target_taxon(taxonomy):
  # Use genus-level names, not species/subtype calls.
  levels=[x.strip().split('(')[0] for x in taxonomy.rstrip(';').split(';')]
  for target in TARGETS:
    for level in levels:
      if re.match(r'^'+re.escape(target)+r'(?:[0-9]*)(?:$|[_\s])',level,re.I):
        return target
  return None


def run_alignment(row,index,threads):
  sample,run,plant,r1,r2=row
  bam=OUT/f'{sample}_{run}.pr2_aligned.bam'
  log=OUT/f'{sample}_{run}.bowtie2.log'
  if bam.is_file() and bam.stat().st_size > 1000:
    print('REUSING BAM for',sample)
    return bam,log
  print(f'ALIGNING {sample} ({plant}) against PR2 SSU...',flush=True)
  tmp=bam.with_suffix('.bam.partial')
  cmd=['bowtie2','--very-sensitive-local','--threads',str(threads),
       '--no-unal','-x',str(index),'-1',str(r1),'-2',str(r2)]
  with log.open('w') as err:
    bow=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err)
    sam=subprocess.run(['samtools','view','-b','-F','4','-o',str(tmp),'-'],
                       stdin=bow.stdout,stderr=err,check=False)
    bow.stdout.close()
    bowrc=bow.wait()
  if bowrc or sam.returncode:
    tmp.unlink(missing_ok=True)
    stop(f'{sample}: Bowtie2={bowrc}, samtools={sam.returncode}. See {log}')
  if not tmp.exists(): stop(f'{sample}: no BAM produced')
  tmp.rename(bam)
  return bam,log


def aln_stats(cigar,nm):
  ops=CIGAR.findall(cigar)
  if not ops: return 0,0
  alen=sum(int(n) for n,op in ops if op in ('M','=','X','I'))
  return alen, 1-(nm/max(alen,1))


def read_alignments(bam,taxa):
  # High-quality read hits are collapsed by fragment / marker genus.
  # Regions of 100bp count bins with at least one qualifying alignment.
  per_pair=collections.defaultdict(lambda:collections.defaultdict(set))
  bins=collections.defaultdict(set)
  accessions=collections.defaultdict(set)
  good_alns=collections.Counter()
  all_primary=0
  unknown_tax=0
  proc=subprocess.Popen(['samtools','view','-F','2308','-q','20',str(bam)],
                        stdout=subprocess.PIPE,text=True)
  for line in proc.stdout:
    row=line.rstrip('\n').split('\t')
    if len(row)<11: continue
    all_primary+=1
    qname,refname,pos,cigar=row[0],row[2],int(row[3]),row[5]
    lineage=taxa.get(refname)
    if lineage is None:
      unknown_tax+=1
      continue
    genus=target_taxon(lineage)
    if not genus: continue
    nm=None
    for field in row[11:]:
      if field.startswith('NM:i:'):
        try: nm=int(field[5:])
        except ValueError: pass
        break
    if nm is None: continue
    aligned,identity=aln_stats(cigar,nm)
    if aligned < 75 or identity < 0.97: continue
    good_alns[genus]+=1
    per_pair[qname][genus].add(refname)
    accessions[genus].add(refname)
    for b in range(max(0,(pos-1)//100),max(0,(pos-1)//100)+max(1,(aligned+99)//100)):
      bins[genus].add((refname,b))
  proc.stdout.close()
  if proc.wait(): stop(f'samtools view failed: {bam}')
  unique=collections.Counter()
  ambiguous=0
  for pair,taxa_found in per_pair.items():
    if len(taxa_found)==1: unique[next(iter(taxa_found))]+=1
    else: ambiguous+=1
  return {
    'unique':unique,'genus_alignments':good_alns,
    'marker_count':{g:len(v) for g,v in accessions.items()},
    'bins_count':{g:len(v) for g,v in bins.items()},
    'ambiguous_pairs':ambiguous,'primary_mapped_mapq20':all_primary,
    'taxids_without_lineage':unknown_tax,
  }


def write_tsv(path,fields,records):
  with path.open('w',newline='') as handle:
    out=csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
    out.writeheader()
    for rec in records: out.writerow(rec)


def main():
  a=argparse.ArgumentParser(description=__doc__)
  action=a.add_mutually_exclusive_group(required=True)
  action.add_argument('--check-only',action='store_true')
  action.add_argument('--run',action='store_true')
  a.add_argument('--threads',type=int,default=4)
  args=a.parse_args()
  if not 1<=args.threads<=16: stop('--threads must be between 1 and 16')
  rows=check()
  if args.check_only: return
  fasta_gz,tax_gz=download_ref()
  fasta=unpack(fasta_gz)
  tax=unpack(tax_gz)
  taxa=read_taxonomies(tax)
  genus_ref_counts=collections.Counter()
  for lineage in taxa.values():
    g=target_taxon(lineage)
    if g: genus_ref_counts[g]+=1
  print('PR2 TARGET GENUS REFERENCE COVERAGE:',dict(genus_ref_counts),flush=True)
  OUT.mkdir(parents=True,exist_ok=True)
  index=build_index(fasta,args.threads)
  all_records=[]; audit=[]
  for row in rows:
    sample,run,plant,*_=row
    bam,log=run_alignment(row,index,args.threads)
    data=read_alignments(bam,taxa)
    for genus in TARGETS:
      all_records.append(dict(sample_id=sample,run_accession=run,treatment_plant=plant,
        genus=genus,pr2_reference_sequences_for_genus=genus_ref_counts[genus],
        qualifying_unique_fragment_pairs=data['unique'][genus],
        qualifying_read_alignments=data['genus_alignments'][genus],
        distinct_pr2_accessions_hit=data['marker_count'].get(genus,0),
        distinct_100bp_reference_bins=data['bins_count'].get(genus,0),
        status='MARKER_LEVEL_CORROBORATION_ONLY_NOT_SPECIES_CONFIRMED'))
    audit.append(dict(sample_id=sample,run_accession=run,treatment_plant=plant,
      primary_mapped_alignments_mapq_ge20=data['primary_mapped_mapq20'],
      unknown_reference_lineages=data['taxids_without_lineage'],
      fragments_with_conflicting_target_genera=data['ambiguous_pairs'],
      bam_file=str(bam.relative_to(ROOT))))
    print(f'FINISHED {sample}: '+', '.join(f'{g}={data["unique"][g]}' for g in TARGETS),flush=True)
  tsv=OUT/'pr2_marker_genus_signal_by_sample.tsv'
  write_tsv(tsv,list(all_records[0]),all_records)
  logtable=OUT/'sample_alignment_audit.tsv'
  write_tsv(logtable,list(audit[0]),audit)
  cov=OUT/'pr2_target_genus_reference_coverage.tsv'
  write_tsv(cov,['genus','available_ssu_references'],
            [dict(genus=g,available_ssu_references=genus_ref_counts[g]) for g in TARGETS])
  methods=OUT/'stage3e_methods.json'
  methods.write_text(json.dumps({
    'stage':'3E independent PR2 marker-reference corroboration',
    'created_utc':datetime.now(timezone.utc).isoformat(),
    'reference':f'PR2 v{VERSION} SSU Mothur FASTA and taxonomy',
    'download_base':BASE,
    'read_source':'All six human-screened FASTQ read pairs; NOT preselected using Kraken2',
    'mapper':'bowtie2 --very-sensitive-local --no-unal; samtools',
    'high_quality_evidence':{'MAPQ_min':20,'aligned_query_bases_min':75,'approx_identity_min':0.97,
       'exclude_flags':2308,'counting_unit':'unique read-pair identifier within a target genus',
       'ambiguous_genera':'excluded from genus-specific counts'},
    'warning':'Not independent diagnostic confirmation or unbiased genome abundance. SSU may not distinguish species or subtypes. Target taxa absent from PR2 cannot be judged absent from samples. Masked gene regions and reference incompleteness remain limitations.'
  },indent=2)+'\n')
  REVIEW.parent.mkdir(parents=True,exist_ok=True)
  with zipfile.ZipFile(REVIEW,'w',zipfile.ZIP_DEFLATED) as z:
    for p in (tsv,logtable,cov,methods): z.write(p,p.name)
  print('SUCCESS: review archive:',REVIEW)
  print('CAVEAT: Marker-level alignments are independent of our suspect genome contigs, but not species/subtype confirmation.')

if __name__=='__main__': main()
