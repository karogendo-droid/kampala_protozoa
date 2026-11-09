#!/usr/bin/env python3
"""Stage 3C: screen public REFERENCE hotspot sequences against NCBI core_nt.

This is a specificity SCREEN of reference loci, NOT confirmation of wastewater organisms.
Only reference-derived FASTA sequences are submitted to NCBI; no wastewater reads
are uploaded. The reference panels in prior stages are not exhaustive.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
import re
import shutil
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / 'metadata/taxonomy/stage3b_locus_review.zip'
REFERENCE = ROOT / 'references/kraken2/kampala_curated_v1/staged_reference_library.fna'
OUT = ROOT / 'results/confirmation/stage3c_reference_specificity_v1'
REVIEW_ZIP = ROOT / 'metadata/taxonomy/stage3c_reference_specificity_review.zip'
FOCAL = ('Enterocytozoon bieneusi H348', 'Blastocystis sp. subtype 6', 'Cyclospora cayetanensis')
GENUS = {
    'Enterocytozoon bieneusi H348': 'Enterocytozoon',
    'Blastocystis sp. subtype 6': 'Blastocystis',
    'Cyclospora cayetanensis': 'Cyclospora',
}
MAX_CONTIGS_PER_ORGANISM = 4
WINDOW_BP = 1000
BLAST_COLUMNS = ['qseqid','saccver','pident','length','qlen','qstart','qend','evalue','bitscore','staxids','sscinames','stitle']


def stop(s):
    raise RuntimeError('STOP: ' + str(s))


def read_zip_tsv(z, name):
    try:
        text = z.read(name).decode('utf-8')
    except KeyError:
        stop(f'Missing {name} in {AUDIT}')
    return list(csv.DictReader(io.StringIO(text), delimiter='\t'))


def write_tsv(path, rows, columns):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def get_candidates():
    if not AUDIT.is_file(): stop(f'Missing audit ZIP: {AUDIT}')
    with zipfile.ZipFile(AUDIT) as archive:
        loci = read_zip_tsv(archive, 'top_reference_loci.tsv')
        metrics = read_zip_tsv(archive, 'locus_concentration_metrics.tsv')
    found = {(r['sample_id'], r['organism_name']) for r in metrics}
    for organism in FOCAL:
        samples = {sample for sample, name in found if name == organism}
        if len(samples) != 6: stop(f'Expected evidence in six samples for {organism}, found {sorted(samples)}')
    grouped = defaultdict(lambda: defaultdict(dict))
    for row in loci:
        org = row['organism_name']
        if org not in FOCAL: continue
        ref = row['reference_contig']
        start = int(row['bin_start_0based'])
        count = int(row['distinct_pairs_in_bin'])
        sample = row['sample_id']
        if start < 0 or count < 0: stop('Invalid locus coordinates or counts')
        # Keep per-sample results distinct; for a given bin, one pair is counted once.
        grouped[org][(ref, start)][sample] = count
    out=[]
    for org in FOCAL:
        if org not in grouped: stop(f'No loci for {org}')
        ranked = sorted(grouped[org].items(), key=lambda it: (-sum(it[1].values()), -len(it[1]), it[0]))
        contigs=set()
        selected=[]
        for (ref,start), values in ranked:
            if ref in contigs: continue  # deliberately sample distinct reference contigs
            contigs.add(ref)
            selected.append((ref,start,values))
            if len(selected)>=MAX_CONTIGS_PER_ORGANISM: break
        for ref,start,values in selected:
            acc_match=re.match(r'^(GC[AF]_\d+\.\d+)_',ref)
            if not acc_match: stop('Unexpected reference ID: '+ref)
            out.append(dict(organism_name=org, genus=GENUS[org],
                            assembly_accession=acc_match.group(1), reference_contig=ref,
                            bin_start_0based=start, sample_count=len(values),
                            aggregate_sample_bin_pair_touches=sum(values.values())))
    if len(out)<9: stop(f'Only {len(out)} candidate loci; expected at least 9')
    for ix,r in enumerate(out,1): r['query_id']=f'KAMPALA_REFLOCUS_{ix:02d}'
    return out


def extract_windows(records):
    if not REFERENCE.is_file(): stop(f'Genome reference FASTA missing: {REFERENCE}')
    by_ref = defaultdict(list)
    for r in records: by_ref[r['reference_contig']].append(r)
    found={}
    ref=None; offset=0; current=[]
    # Stream the large FASTA without loading the ~786 Mb reference into RAM.
    with REFERENCE.open('r', encoding='ascii') as stream:
        for line in stream:
            if line.startswith('>'):
                ref=line[1:].strip().split()[0]
                offset=0
                current=by_ref.get(ref, [])
                if current:
                    for r in current: found[r['query_id']]=[]
                continue
            if not current: continue
            chunk=line.strip().upper()
            end=offset+len(chunk)
            for r in current:
                start=r['bin_start_0based']
                a=max(start,offset)
                b=min(start+WINDOW_BP,end)
                if b>a: found[r['query_id']].append(chunk[a-offset:b-offset])
            offset=end
    result={k:''.join(v) for k,v in found.items()}
    for r in records:
        seq=result.get(r['query_id'],'')
        if len(seq)<75: stop(f'No adequate reference sequence for {r["query_id"]}: {r["reference_contig"]} {r["bin_start_0based"]}')
        if any(ch not in 'ACGTN' for ch in seq):
            stop(f'Unexpected FASTA base in {r["query_id"]}')
    return result


def shannon_entropy(seq):
    n=len(seq)
    c=Counter(seq)
    return round(-sum((v/n)*math.log2(v/n) for v in c.values()),4) if n else 0


def homopolymer_run(seq):
    last=None; size=0; high=0
    for x in seq:
        size = size+1 if x==last else 1
        high=max(high,size)
        last=x
    return high


def prepare():
    candidates=get_candidates()
    seqs=extract_windows(candidates)
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'public_reference_hotspots.fna').open('w',encoding='ascii') as fasta:
        for r in candidates:
            seq=seqs[r['query_id']]
            r.update(query_length=len(seq),gc_percent=round(100*(seq.count('C')+seq.count('G'))/len(seq),2),
                     N_percent=round(100*seq.count('N')/len(seq),2),
                     base_entropy_bits=shannon_entropy(seq),
                     max_homopolymer=homopolymer_run(seq),
                     sequence_sha256=hashlib.sha256(seq.encode('ascii')).hexdigest())
            fasta.write('>'+r['query_id']+'\n')
            for i in range(0,len(seq),80):fasta.write(seq[i:i+80]+'\n')
    cols=['query_id','organism_name','genus','assembly_accession','reference_contig',
          'bin_start_0based','sample_count','aggregate_sample_bin_pair_touches',
          'query_length','gc_percent','N_percent','base_entropy_bits','max_homopolymer',
          'sequence_sha256']
    write_tsv(OUT/'query_manifest.tsv',candidates,cols)
    meta=dict(stage='3C public-reference hotspot screening',
              created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
              inputs=dict(stage3b_zip=str(AUDIT.relative_to(ROOT)),
                          reference_fasta=str(REFERENCE.relative_to(ROOT))),
              method=('Select up to four most-supported 1-kb windows on distinct reference contigs per focal organism; '
                      'query only PUBLIC REFERENCE sequence against independent NCBI core_nt by BLASTN.'),
              privacy=('No wastewater reads, read identifiers, FASTQs, or BAMs are submitted to NCBI.'),
              caveats=['Top-locus bin-touch tallies may count a read in more than one bin.',
                       'Reference-locus specificity is not direct proof of the identity of environmental reads.',
                       'core_nt is a non-exhaustive, partially collapsed reference collection; incomplete metadata can mislead.',
                       'Remote BLAST results are an exploratory specificity screen, not confirmed detections.'],
              parameters=dict(max_contigs_per_organism=MAX_CONTIGS_PER_ORGANISM,
                              window_bp=WINDOW_BP,blast_database='core_nt'))
    (OUT/'stage3c_methods.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(f'PREPARED {len(candidates)} public reference windows across {len(FOCAL)} organism groups',flush=True)
    for r in candidates:
        print(f'  {r["query_id"]} {r["organism_name"]}: {r["reference_contig"]} position {r["bin_start_0based"]} length {r["query_length"]}',flush=True)
    return candidates


def run_blast():
    path=OUT/'public_reference_hotspots.fna'
    blast=shutil.which('blastn')
    if not blast: stop('blastn not found. Activate kampala-taxanomy; if needed install BLAST+ there.')
    raw=OUT/'reference_hotspots_core_nt_blast.tsv'
    command=[blast,'-remote','-db','core_nt','-task','blastn','-query',str(path),
             '-dust','yes','-evalue','1e-10','-max_target_seqs','100','-max_hsps','1',
             '-outfmt','6 '+' '.join(BLAST_COLUMNS),'-out',str(raw)]
    print('Submitting ONLY public assembly reference windows to NCBI BLAST core_nt...',flush=True)
    print('NCBI remote BLAST can be slow or unavailable; no local genome database download is needed.',flush=True)
    p=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    if p.returncode!=0:
        if raw.exists():raw.unlink()
        stop('Remote BLAST failed (no valid results recorded). NCBI response: '+(p.stderr or p.stdout)[-1200:])
    if not raw.is_file():stop('BLAST exited successfully but produced no output file')
    print('BLAST results saved:',raw,flush=True)
    return raw


def interpret(candidates, raw):
    byquery=defaultdict(list)
    with raw.open(encoding='utf-8') as handle:
        for lineno,line in enumerate(handle,1):
            if not line.strip() or line.startswith('#'):continue
            fields=line.rstrip('\r\n').split('\t')
            if len(fields)!=len(BLAST_COLUMNS):stop(f'Malformed BLAST output at line {lineno}; got {len(fields)} columns')
            r=dict(zip(BLAST_COLUMNS,fields))
            try:r['_score']=float(r['bitscore']);r['_pident']=float(r['pident']);r['_len']=int(r['length'])
            except ValueError:stop(f'Invalid numeric BLAST field at line {lineno}')
            byquery[r['qseqid']].append(r)
    summaries=[]
    top_rows=[]
    for row in candidates:
        q=row['query_id']
        matches=sorted(byquery.get(q,[]),key=lambda r:(-r['_score'],-r['_len'],r['saccver']))
        best=matches[0] if matches else None
        near=[]
        if best:
            near=[m for m in matches if m['_score']>=0.98*best['_score'] and m['_len']>=0.8*best['_len']]
        other=[m for m in near if row['genus'].lower() not in m['sscinames'].lower()]
        summaries.append(dict(query_id=q,expected_genus=row['genus'],
                              total_reported_hits=len(matches),
                              best_subject=best['saccver'] if best else '',
                              best_scientific_name=best['sscinames'] if best else '',
                              best_percent_identity=best['pident'] if best else '',
                              best_alignment_bp=best['length'] if best else '',
                              best_query_bp=row['query_length'],
                              near_best_hit_count=len(near),
                              near_best_nonmatching_genus_hits=len(other),
                              example_nonmatching_name=other[0]['sscinames'] if other else '',
                              note=('NO_HITS' if not best else
                                    'REVIEW_NONMATCHING_NEAR_BEST' if other else
                                    'NO_NONMATCHING_HIT_IN_REPORTED_NEAR_BEST_SET_NOT_PROOF_OF_SPECIFICITY')))
        for hit in matches[:20]:
            top_rows.append({k:hit[k] for k in BLAST_COLUMNS})
    write_tsv(OUT/'reference_hotspot_blast_summary.tsv',summaries,
              list(summaries[0].keys()))
    write_tsv(OUT/'top20_reference_blast_hits.tsv',top_rows,BLAST_COLUMNS)
    print('\nSCREEN SUMMARY (NOT direct read confirmation):',flush=True)
    for r in summaries:
        print(f'  {r["query_id"]} {r["expected_genus"]}: {r["best_scientific_name"] or "NO HIT"}; {r["note"]}',flush=True)
    return summaries


def make_review(include_hits):
    paths=['query_manifest.tsv','stage3c_methods.json']
    if include_hits:
        paths+=['reference_hotspot_blast_summary.tsv','top20_reference_blast_hits.tsv']
    REVIEW_ZIP.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(REVIEW_ZIP,'w',compression=zipfile.ZIP_DEFLATED) as output:
        for name in paths:output.write(OUT/name,arcname=name)
    print('Review ZIP created:',REVIEW_ZIP,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check-only',action='store_true',help='Check local stage 3B ZIP and public reference FASTA')
    mode.add_argument('--prepare-only',action='store_true',help='Extract selected PUBLIC REFERENCE windows; no internet')
    mode.add_argument('--run',action='store_true',help='Extract windows and submit only PUBLIC REFERENCE windows to NCBI BLAST')
    args=parser.parse_args()
    if args.check_only:
        c=get_candidates()
        if not REFERENCE.is_file():stop('Missing reference: '+str(REFERENCE))
        print(f'PREFLIGHT PASS: {len(c)} reference windows selectable; FASTA found. Nothing submitted.',flush=True)
        return
    candidates=prepare()
    if args.prepare_only:
        make_review(False)
        print('Preparation completed. No BLAST search or wastewater analysis was performed.')
        return
    raw=run_blast()
    interpret(candidates,raw)
    make_review(True)
    print('STAGE 3C REFERENCE SCREEN FINISHED. These results do NOT establish organism detection.',flush=True)


if __name__=='__main__':
    try: main()
    except (OSError,ValueError,RuntimeError,zipfile.BadZipFile) as error:
        print('ERROR:',error,file=sys.stderr)
        sys.exit(1)
