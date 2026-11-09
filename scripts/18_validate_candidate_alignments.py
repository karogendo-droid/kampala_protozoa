#!/usr/bin/env python3
"""Stage 3A ONLY: orthogonal competitive alignment of Stage-1 candidate read pairs.

This is *conditional* validation: reads were preselected by Kraken2 and are
aligned to the same curated, incomplete reference collection using Bowtie2.
Genome breadth, mapping specificity within this panel and paired-mate support
are examined. No species is automatically declared detected. Stage 3B must
check selected findings against broader independent references/markers.

Read-level files remain on the researcher's Mac, never in the review ZIP.
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / 'references/kraken2/kampala_curated_v1'
REF = DB / 'staged_reference_library.fna'
SELECTED = ROOT / 'metadata/taxonomy/curated_v1_selected_assemblies.tsv'
STAGE1 = ROOT / 'results/taxonomy/kraken2_curated_v1_baseline'
STAGE2 = ROOT / 'results/abundance/stage2_v1'
INDEX_DIR = ROOT / 'references/validation_bowtie2/curated_v1'
PREFIX = INDEX_DIR / 'competitive'
OUT = ROOT / 'results/confirmation/stage3a_curated_v1'
REVIEW_ZIP = ROOT / 'metadata/taxonomy/stage3a_alignment_review.zip'
RUNS = ('SRR28434174','SRR28434173','SRR28434172',
        'SRR28434171','SRR28434170','SRR28434169')
GENERA = {'Giardia':5740, 'Cryptosporidium':5806, 'Entamoeba':5758,
          'Cyclospora':44417, 'Blastocystis':12967,
          'Enterocytozoon':27971, 'Encephalitozoon':6033}
MIN_MAPQ = 20
MIN_ALIGNED_BP = 75
MIN_IDENTITY = 0.95
BIN_BP = 1000
THREADS = 4
ACCESSION_RE = re.compile(r'^(GC[AF]_\d+\.\d+)_\d+\|kraken:taxid\|(\d+)$')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def fail(s):
    raise RuntimeError('STOP: '+s)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    if not path.is_file():
        fail('Required file missing: '+str(path))
    with path.open(newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, columns):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer=csv.DictWriter(fh, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def required_files():
    if not REF.is_file() or not SELECTED.is_file():
        fail('Curated reference FASTA or selected-assembly manifest missing.')
    if not (DB/'BUILD_COMPLETE').is_file():
        fail('Curated Kraken2 database has not passed its operational check.')
    for exe in ('bowtie2-build','bowtie2','samtools'):
        if not shutil.which(exe):
            fail(f'{exe} not on PATH; activate environment: conda activate kampala-align')
    check = subprocess.run(['samtools','--version'],capture_output=True,text=True,check=True)
    print(check.stdout.splitlines()[0],flush=True)
    check = subprocess.run(['bowtie2','--version'],capture_output=True,text=True,check=True)
    print(check.stdout.splitlines()[0],flush=True)
    manifest=read_tsv(SELECTED)
    byacc={r['assembly_accession']:r for r in manifest}
    if len(manifest)!=44 or len(byacc)!=44:
        fail('Expected exactly 44 distinct selected reference assemblies.')
    if not (DB/'taxonomy/nodes.dmp').is_file():
        fail('NCBI taxonomy nodes.dmp missing from curated database.')
    if not STAGE2.joinpath('candidate_taxon_direct_rpm.tsv').is_file():
        fail('Stage 2 output not found. Do not run Stage 3 before Stage 2.')
    methods_path=STAGE1/'classification_methods.json'
    if not methods_path.is_file():
        fail('Missing Stage 1 methods record.')
    methods=json.loads(methods_path.read_text())
    if methods.get('analysis')!='Stage 1 Kraken2 exploratory taxonomic classification':
        fail('Unexpected Stage 1 methods record.')
    for filename, known in methods.get('database_sha256', {}).items():
        check_path=(SELECTED if filename=='reference_manifest' else DB/filename)
        if not check_path.is_file(): fail('Database fingerprint input absent: '+str(check_path))
        if filename=='reference_manifest' and sha256(check_path)!=known:
            fail('Stage 1 selection manifest has changed since classification.')
    sample_rows=read_tsv(STAGE2/'sample_classification_metrics.tsv')
    if len(sample_rows)!=6 or {r['run_accession'] for r in sample_rows}!=set(RUNS):
        fail('Stage 2 metrics must contain six expected runs.')
    samples={r['run_accession']:r for r in sample_rows}
    genus_rows=read_tsv(STAGE2/'genus_clade_signal_rpm.tsv')
    for g_row in genus_rows:
        run=g_row['run_accession']
        if run not in samples or g_row['genus'] not in GENERA:
            fail('Unexpected Stage 2 genus table entry')
        samples[run].setdefault('expected_genus_counts', {})[g_row['genus']]=int(g_row['genus_clade_pairs'])
    for run in RUNS:
        if set(samples[run].get('expected_genus_counts',{}))!=set(GENERA):
            fail('Missing Stage 2 genus counts for '+run)
    for run in RUNS:
        folder=STAGE1/run
        record_path=folder/'completion.json'
        if not record_path.is_file(): fail('Missing Stage 1 completion record for '+run)
        record=json.loads(record_path.read_text())
        if record.get('run')!=run or int(record.get('input_pairs',-1))!=int(samples[run]['input_pairs']):
            fail('Stage 1 / Stage 2 counts differ: '+run)
        if record.get('database_sha256')!=methods['database_sha256']:
            fail('Stage 1 database hashes mismatch: '+run)
        input_paths=record.get('input_fastqs',[])
        if len(input_paths)!=2: fail('Expected two screened FASTQs for '+run)
        paths=[]
        for rel in input_paths:
            p=(ROOT/rel).resolve()
            if not p.is_relative_to(ROOT.resolve()) or not p.is_file() or not p.name.endswith('.human_screened.fastq.gz'):
                fail('Unexpected screened FASTQ path for '+run+': '+str(p))
            paths.append(p)
        assignments=folder/'assignments.tsv.gz'
        if not assignments.is_file(): fail('Missing Kraken assignment output for '+run)
        if not (folder/'kraken2.report.txt').is_file(): fail('Missing Kraken report for '+run)
        samples[run]['r1'], samples[run]['r2'], samples[run]['assignments']=(*paths, assignments)
        print(f'FOUND {run}: {samples[run]["input_pairs"]} screened pairs',flush=True)
    if shutil.disk_usage(ROOT).free < 20*(1024**3):
        fail('Less than 20 GiB free disk space. Resolve storage before building index.')
    return manifest, samples


def read_taxonomy_parents():
    parents={}
    with (DB/'taxonomy/nodes.dmp').open(encoding='utf-8') as fh:
        for line in fh:
            fields=line.split('|')
            if len(fields)>=2:
                parents[int(fields[0].strip())]=int(fields[1].strip())
    for t in GENERA.values():
        if t not in parents: fail('Genus taxid missing from taxonomy: '+str(t))
    return parents


def taxid_genus(taxid, parents, cache):
    if taxid in cache: return cache[taxid]
    seen=set(); cur=taxid; hits=set()
    while cur in parents and cur not in seen:
        seen.add(cur)
        for name,gid in GENERA.items():
            if cur==gid: hits.add(name)
        if cur==parents[cur]: break
        cur=parents[cur]
    cache[taxid]=hits
    return hits


def clean_read_id(raw):
    text=raw.strip().split()[0]
    if text.startswith('@'):text=text[1:]
    if text.endswith('/1') or text.endswith('/2'):text=text[:-2]
    return text


def read_candidate_ids(path, parents):
    ids=set(); genuscounts=Counter(); cache={}
    total=0
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        for line in f:
            cols=line.rstrip('\r\n').split('\t',4)
            if len(cols)<4 or cols[0] not in ('C','U'):
                fail('Malformed Stage 1 Kraken assignment at line '+str(total+1))
            total+=1
            if cols[0]!='C': continue
            try: tid=int(cols[2])
            except ValueError:fail('Invalid Kraken taxid at line '+str(total))
            groups=taxid_genus(tid,parents,cache)
            if groups:
                rid=clean_read_id(cols[1])
                if rid in ids: fail('Duplicate selected Kraken read identifier: '+rid)
                ids.add(rid)
                for g in groups:genuscounts[g]+=1
    return ids,genuscounts,total


def extract_reads(source1,source2,ids,dest1,dest2,expected):
    extracted=set(); total=0
    with gzip.open(source1,'rt',encoding='ascii') as x, gzip.open(source2,'rt',encoding='ascii') as y, \
         gzip.open(dest1,'wt',encoding='ascii',compresslevel=4) as o1, \
         gzip.open(dest2,'wt',encoding='ascii',compresslevel=4) as o2:
        while True:
            r1=[x.readline() for _ in range(4)]
            r2=[y.readline() for _ in range(4)]
            if not r1[0] and not r2[0]:break
            if not all(r1) or not all(r2):fail('Broken or unequal FASTQ record near pair '+str(total+1))
            id1,id2=clean_read_id(r1[0]),clean_read_id(r2[0])
            if id1!=id2 or not r1[0].startswith('@') or not r2[0].startswith('@'):
                fail('Mismatched FASTQ pair names at pair '+str(total+1))
            if r1[2][0]!='+' or r2[2][0]!='+' or len(r1[1].strip())!=len(r1[3].strip()) or len(r2[1].strip())!=len(r2[3].strip()):
                fail('Malformed FASTQ at pair '+str(total+1))
            total+=1
            if id1 in ids:
                if id1 in extracted:fail('Duplicate selected FASTQ identifier: '+id1)
                extracted.add(id1)
                o1.writelines(r1);o2.writelines(r2)
            if total%10_000_000==0: print(f'  Scanned {total:,} read pairs',flush=True)
    if total!=expected or extracted!=ids:
        fail(f'Input or selected read ID mismatch: FASTQ pairs={total}, expected={expected}; selected={len(ids)}, found={len(extracted)}, missing={len(ids-extracted)}')
    return len(extracted)


def fasta_catalog(path):
    lengths=defaultdict(int); contigs=defaultdict(dict)
    cur_acc=cur_name=None; size=0
    def flush():
        if cur_acc:
            lengths[cur_acc]+=size
            contigs[cur_acc][cur_name]=size
    with path.open('rt',encoding='ascii') as fh:
        for line in fh:
            if line.startswith('>'):
                flush()
                token=line[1:].split()[0]
                m=ACCESSION_RE.match(token)
                if token=='GRCh38_chrM|kraken:taxid|9606':
                    cur_acc='GRCh38_chrM'
                elif m:cur_acc=m.group(1)
                else:fail('Unexpected FASTA ID: '+token)
                cur_name=token;size=0
            else:
                size+=len(line.strip())
        flush()
    return lengths,contigs


def index_files():
    for suffix in ('.bt2','.bt2l'):
        paths=[Path(str(PREFIX)+x+suffix) for x in ('.1','.2','.3','.4','.rev.1','.rev.2')]
        if all(p.is_file() and p.stat().st_size>0 for p in paths):return paths
    return []


def ensure_index(sha):
    marker=INDEX_DIR/'index_complete.json'
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    existing=list(INDEX_DIR.glob('competitive*.bt2*'))
    if marker.is_file():
        m=json.loads(marker.read_text())
        paths=index_files()
        if m.get('reference_sha256')!=sha or not paths:
            fail('Existing index marker disagrees with FASTA or index files. Do not overwrite.')
        print('REUSED VERIFIED INDEX:',str(PREFIX),flush=True)
        return
    if existing:
        fail('Unmarked Bowtie2 index files exist. Inspect before overwriting: '+str(INDEX_DIR))
    print('Building competitive Bowtie2 index from all 44 reference assemblies + mitochondrial decoy',flush=True)
    with (INDEX_DIR/'build.log').open('w') as log:
        cmd=['bowtie2-build','--threads',str(THREADS),str(REF),str(PREFIX)]
        code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
    if code!=0 or not index_files():fail('Bowtie2 index build failed: see '+str(INDEX_DIR/'build.log'))
    marker.write_text(json.dumps({'reference_sha256':sha,'files':[p.name for p in index_files()],
                                 'created_utc':dt.datetime.now(dt.timezone.utc).isoformat()},indent=2)+'\n')
    print('INDEX COMPLETE',flush=True)


def run_align(r1,r2,bam,log):
    cmd=['bowtie2','--very-sensitive-local','--phred33','--no-unal',
         '-k','3','--seed','42','--threads',str(THREADS),
         '-x',str(PREFIX),'-1',str(r1),'-2',str(r2)]
    view=['samtools','view','-b','-F','4','-']
    sort=['samtools','sort','-@','2','-m','512M','-o',str(bam),'-']
    with log.open('wb') as stderr:
        a=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=stderr)
        b=subprocess.Popen(view,stdin=a.stdout,stdout=subprocess.PIPE,stderr=stderr)
        a.stdout.close()
        c=subprocess.Popen(sort,stdin=b.stdout,stdout=stderr,stderr=stderr)
        b.stdout.close()
        cr=c.wait();br=b.wait();ar=a.wait()
    if any(x!=0 for x in (ar,br,cr)):
        fail(f'Bowtie2/samtools returned exit status {(ar,br,cr)}; inspect {log}')
    subprocess.run(['samtools','quickcheck','-v',str(bam)],check=True)
    subprocess.run(['samtools','index','-@','2',str(bam)],check=True)


def cigar_blocks(cigar, pos0):
    operations=CIGAR_RE.findall(cigar)
    if not operations or ''.join(n+code for n,code in operations)!=cigar:
        return [],0
    ref=pos0;blocks=[]; aligned=0
    for n,op in operations:
        k=int(n)
        if op in ('M','=','X'):
            blocks.append((ref,ref+k));ref+=k;aligned+=k
        elif op in ('D','N'):ref+=k
    return blocks,aligned


def union_bp(intervals):
    if not intervals:return 0
    intervals.sort();s,e=intervals[0];covered=0
    for ns,ne in intervals[1:]:
        if ns<=e:e=max(e,ne)
        else:covered+=e-s;s,e=ns,ne
    return covered+(e-s)


def evaluate_bam(bam, manifest, genome_lengths, contigs):
    roles={r['assembly_accession']:r for r in manifest}
    intervals=defaultdict(list); bins=defaultdict(set); pairs=defaultdict(set)
    proper=defaultdict(lambda:defaultdict(set)); counts=Counter(); all_aligned=Counter()
    p=subprocess.Popen(['samtools','view','-F','2308',str(bam)],stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,text=True)
    for line in p.stdout:
        parts=line.rstrip('\n').split('\t')
        if len(parts)<11:fail('Malformed BAM alignment')
        qname=clean_read_id(parts[0]);flag=int(parts[1]);ref=parts[2];pos0=int(parts[3])-1
        mapq=int(parts[4]);cigar=parts[5]
        m=ACCESSION_RE.match(ref)
        acc=m.group(1) if m else ('GRCh38_chrM' if ref.startswith('GRCh38_chrM|') else None)
        if acc is None:fail('Unknown alignment reference: '+ref)
        if acc=='GRCh38_chrM':
            all_aligned[acc]+=1;continue
        if acc not in roles or ref not in contigs[acc]:fail('Alignment contig not in curated FASTA')
        all_aligned[acc]+=1
        blocks,alen=cigar_blocks(cigar,pos0)
        if alen<MIN_ALIGNED_BP or mapq<MIN_MAPQ:continue
        nm=None
        for field in parts[11:]:
            if field.startswith('NM:i:'):
                try:nm=int(field[5:])
                except ValueError:fail('Invalid NM tag')
                break
        if nm is None:continue
        identity=1-nm/max(1,alen)
        if identity<MIN_IDENTITY:continue
        mate='R1' if flag&64 else 'R2' if flag&128 else '?'
        if mate=='?':fail('Alignment not marked as R1/R2')
        counts[(acc,mate)]+=1
        pairs[acc].add(qname)
        if flag&2:proper[acc][qname].add(mate)
        for s,e in blocks:
            if s<0 or e>contigs[acc][ref]:fail('Alignment coordinates beyond contig')
            intervals[(acc,ref)].append((s,e))
            for b in range(s//BIN_BP, (e-1)//BIN_BP+1):bins[acc].add((ref,b))
    err=p.stderr.read();returncode=p.wait()
    if returncode:fail('samtools view failed: '+err[:400])
    rows=[]
    for acc,r in sorted(roles.items()):
        breadth=sum(union_bp(v) for (a,ref),v in intervals.items() if a==acc)
        paired=sum(1 for mate_set in proper[acc].values() if mate_set=={'R1','R2'})
        if breadth>genome_lengths[acc]:fail('Covered bases exceed genome length for '+acc)
        rows.append({'assembly_accession':acc,'organism_name':r['organism_name'],
                     'role':r['role'],'target_group':r['target_group'],'genome_bp':genome_lengths[acc],
                     'primary_mapped_alignments_before_filters':all_aligned[acc],
                     'high_quality_R1':counts[(acc,'R1')], 'high_quality_R2':counts[(acc,'R2')],
                     'high_quality_distinct_pairs':len(pairs[acc]),
                     'high_quality_concordant_pairs':paired,
                     'covered_reference_bases':breadth,
                     'reference_breadth_percent':f'{100*breadth/genome_lengths[acc]:.6f}',
                     'distinct_1kb_reference_bins':len(bins[acc]),
                     'reference_contigs_with_coverage':sum(1 for (a,rname),v in intervals.items() if a==acc and v)})
    return rows


def completion_ok(run_dir,run,ref_sha):
    marker=run_dir/'completion.json'
    if not marker.exists():return False
    record=json.loads(marker.read_text())
    if record.get('run')!=run or record.get('reference_sha256')!=ref_sha:
        fail('Prior completion record changed for '+run)
    if record.get('selected_pairs',0)>0:
        for name in ('candidate_R1.fastq.gz','candidate_R2.fastq.gz',
                     'candidate_mapped.sorted.bam','candidate_mapped.sorted.bam.bai'):
            p=run_dir/name
            if not p.is_file() or p.stat().st_size==0:
                fail('Prior completed sequence/alignment file missing: '+str(p))
    for relative,digest in record.get('checksums',{}).items():
        path=run_dir/relative
        if not path.is_file() or sha256(path)!=digest:
            fail('Prior completed output altered: '+str(path))
    return True


def do_run(run,row,parents,manifest,genome_lengths,contigs,ref_sha):
    d=OUT/run
    if d.exists():
        if completion_ok(d,run,ref_sha):
            print('REUSED VERIFIED:',run,flush=True)
            return
        fail('Incomplete prior result exists; do not overwrite: '+str(d))
    d.mkdir(parents=True)
    print('\nSTAGE 3A:',run,row['sample_id'],row['treatment_plant'],flush=True)
    try:
        ids,groupcounts,assignment_total=read_candidate_ids(row['assignments'],parents)
        expected=int(row['input_pairs'])
        if assignment_total!=expected:fail('Kraken assignment count mismatch for '+run)
        if any(groupcounts[g]!=row['expected_genus_counts'][g] for g in GENERA):
            fail('Selected genus counts do not reconcile with Stage 2 for '+run+': '+str(groupcounts))
        print(f'  Selected {len(ids):,} candidate pairs by genus-lineage assignment',flush=True)
        (d/'candidate_selection.json').write_text(json.dumps(
            {'run':run,'selected_pairs':len(ids),'stage1_total_pairs':expected,'by_genus':groupcounts},indent=2)+'\n')
        r1=d/'candidate_R1.fastq.gz';r2=d/'candidate_R2.fastq.gz'
        if ids:
            n=extract_reads(row['r1'],row['r2'],ids,r1,r2,expected)
            if n!=len(ids):fail('FASTQ extraction count mismatch')
            bam=d/'candidate_mapped.sorted.bam'
            run_align(r1,r2,bam,d/'bowtie2_samtools.log')
            mapped=evaluate_bam(bam,manifest,genome_lengths,contigs)
        else:
            mapped=[{'assembly_accession':r['assembly_accession'], 'organism_name':r['organism_name'],
                     'role':r['role'], 'target_group':r['target_group'],'genome_bp':genome_lengths[r['assembly_accession']],
                     'primary_mapped_alignments_before_filters':0,'high_quality_R1':0,'high_quality_R2':0,
                     'high_quality_distinct_pairs':0,'high_quality_concordant_pairs':0,
                     'covered_reference_bases':0,'reference_breadth_percent':'0.000000',
                     'distinct_1kb_reference_bins':0,'reference_contigs_with_coverage':0} for r in manifest]
        for x in mapped:
            x['run_accession']=run;x['sample_id']=row['sample_id'];x['treatment_plant']=row['treatment_plant'];
            x['kraken_candidate_pairs_all_target_genera']=len(ids)
        cols=['run_accession','sample_id','treatment_plant','assembly_accession','organism_name','role',
              'target_group','genome_bp','kraken_candidate_pairs_all_target_genera',
              'primary_mapped_alignments_before_filters','high_quality_R1','high_quality_R2',
              'high_quality_distinct_pairs','high_quality_concordant_pairs',
              'covered_reference_bases','reference_breadth_percent',
              'distinct_1kb_reference_bins','reference_contigs_with_coverage']
        write_tsv(d/'competitive_alignment_metrics.tsv',mapped,cols)
        checksums={'candidate_selection.json':sha256(d/'candidate_selection.json'),
                   'competitive_alignment_metrics.tsv':sha256(d/'competitive_alignment_metrics.tsv')}
        (d/'completion.json').write_text(json.dumps({'run':run,'reference_sha256':ref_sha,
           'selected_pairs':len(ids),'checksums':checksums,
           'finished_utc':dt.datetime.now(dt.timezone.utc).isoformat()},indent=2)+'\n')
        print(f'PASS: {run} | candidate pairs {len(ids):,}',flush=True)
    except Exception:
        print('INCOMPLETE — stop and inspect '+str(d),file=sys.stderr)
        raise


def summary_zip(manifest,ref_sha):
    if REVIEW_ZIP.exists():fail('Review ZIP already exists: '+str(REVIEW_ZIP))
    rows=[]; selections=[]
    for run in RUNS:
        d=OUT/run
        x=json.loads((d/'candidate_selection.json').read_text())
        selections.append({'run_accession':run,'candidate_pairs':x['selected_pairs'],
                           'stage1_total_pairs':x['stage1_total_pairs'],
                           **{g:x['by_genus'].get(g,0) for g in GENERA}})
        rows.extend(read_tsv(d/'competitive_alignment_metrics.tsv'))
    cols=['run_accession','sample_id','treatment_plant','assembly_accession','organism_name','role',
          'target_group','genome_bp','kraken_candidate_pairs_all_target_genera',
          'primary_mapped_alignments_before_filters','high_quality_R1','high_quality_R2',
          'high_quality_distinct_pairs','high_quality_concordant_pairs',
          'covered_reference_bases','reference_breadth_percent',
          'distinct_1kb_reference_bins','reference_contigs_with_coverage']
    write_tsv(OUT/'all_assembly_alignment_metrics.tsv',rows,cols)
    selcols=['run_accession','stage1_total_pairs','candidate_pairs',*GENERA]
    write_tsv(OUT/'candidate_selection_summary.tsv',selections,selcols)
    note={'stage':'3A: orthogonal competitive alignment of Kraken2-selected candidate pairs',
          'reference':'44 curated genomes and human chrM, SAME reference set used for Kraken2 (not independent sequence source)',
          'reference_fasta_sha256':ref_sha, 'aligner':'Bowtie2 --very-sensitive-local -k 3',
          'filter':{'minimum_MAPQ':MIN_MAPQ,'minimum_aligned_bp':MIN_ALIGNED_BP,
                    'minimum_approx_identity':MIN_IDENTITY,'secondary_and_supplementary_excluded':True},
          'read_selection':'all Kraken2 paired assignments falling under the seven target genus taxids',
          'limitations':['Kraken-preselected reads: cannot reveal Kraken false negatives.',
                         'Same reference library as initial classifier, so validation is not independent of reference choice.',
                         'MAPQ measures uniqueness only among the limited included genomes.',
                         'Genome-breadth estimates use primary high-quality alignments, not biological abundance.',
                         'Repeated reference regions, sequence contaminants or shared loci may still produce misleading evidence.',
                         'No species automatically confirmed; further BLAST/independent marker and broader decoy checks required.',
                         'No clinical or viability conclusion, no treatment-plant comparison.']}
    (OUT/'stage3a_methods.json').write_text(json.dumps(note,indent=2)+'\n')
    REVIEW_ZIP.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(REVIEW_ZIP,'w',zipfile.ZIP_DEFLATED) as zipout:
        for name in ['all_assembly_alignment_metrics.tsv','candidate_selection_summary.tsv','stage3a_methods.json']:
            zipout.write(OUT/name,name)
    print('\nSTAGE 3A COMPLETE — competitive mapping reports ready; NO taxa automatically confirmed.',flush=True)
    print('Upload only: '+str(REVIEW_ZIP.relative_to(ROOT)),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check-only',action='store_true')
    mode.add_argument('--build-index-only',action='store_true')
    mode.add_argument('--run',action='store_true')
    args=ap.parse_args()
    os.chdir(ROOT)
    manifest,samples=required_files()
    if args.check_only:
        print('PREFLIGHT PASS. No index built, no reads extracted or aligned.',flush=True)
        return
    ref_sha=sha256(REF)
    ensure_index(ref_sha)
    if args.build_index_only:
        print('INDEX READY. No reads extracted or aligned.',flush=True)
        return
    parents=read_taxonomy_parents()
    genome_lengths,contigs=fasta_catalog(REF)
    for r in manifest:
        if r['assembly_accession'] not in genome_lengths:
            fail('Missing reference assembly in FASTA: '+r['assembly_accession'])
    OUT.mkdir(parents=True,exist_ok=True)
    for run in RUNS:
        do_run(run,samples[run],parents,manifest,genome_lengths,contigs,ref_sha)
    summary_zip(manifest,ref_sha)


if __name__=='__main__':
    try: main()
    except (OSError,RuntimeError,ValueError,subprocess.SubprocessError) as e:
        print('ERROR: '+str(e),file=sys.stderr)
        raise SystemExit(2)
