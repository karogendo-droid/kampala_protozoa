#!/usr/bin/env python3
"""PRJNA1091428: restartable, sequential workflow using the validated RES001 settings.

Five nonpilot samples: normalized SRA -> paired FASTQ -> full-file validation
-> before FastQC -> same fastp 1.4.0 baseline -> Bowtie2 human screen.
RES001 is checked and reused (never reprocessed here).

Run this file from PROJECT/scripts/ (root inferred from path). Python >= 3.10.
No specimen data are uploaded; only SRA public accessions are downloaded.
No original or intermediate user files are deleted by this orchestrator.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import fcntl
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

PROJECT = "PRJNA1091428"
PILOT = "SRR28434174"
RUN_TO_SAMPLE = {
    "SRR28434174": ("RES001", "Bugolobi", "Emptier_truck", "Dry_morning"),
    "SRR28434173": ("RES002", "Bugolobi", "Environmental", "Dry_morning"),
    "SRR28434172": ("RES003", "Lubigi", "Environmental", "Dry_morning"),
    "SRR28434171": ("RES004", "Lubigi", "Emptier_truck", "Dry_morning"),
    "SRR28434170": ("RES005", "Bugolobi", "Environmental", "Dry_evening"),
    "SRR28434169": ("RES013", "Lubigi", "Environmental", "Dry_evening"),
}
RUNS = list(RUN_TO_SAMPLE)
GIB = 1024**3
BLOCK = 4 * 1024 * 1024
EXPECTED_HUMAN_SCRIPT_SHA256 = "95125d099162be8ee18a53c9ee62ea7e4f0b918027cffc3978f88c1c18ddca71"
VALID_BASES = frozenset(b'ACGTNacgtn')


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def tag() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(BLOCK), b''):
            h.update(chunk)
    return h.hexdigest()


def rel(root: Path, p: Path) -> str:
    return str(p.relative_to(root))


def need_space(root: Path, minimum: int, context: str) -> None:
    available = shutil.disk_usage(root).free
    print(f'Disk space before {context}: {available / GIB:.1f} GiB free (required reserve: {minimum / GIB:.1f} GiB)', flush=True)
    if available < minimum:
        raise RuntimeError(f'Insufficient disk space for {context}. Stop and review existing files; no automatic deletion.')


def run_cmd(cmd: list[str], log: Path, *, env=None) -> str:
    """Stream command output to log and Terminal without buffering large outputs."""
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f'COMMAND: {shlex.join(cmd)}', flush=True)
    lines = []
    with log.open('a', encoding='utf-8') as out:
        out.write('\n' + now() + '\nCOMMAND: ' + shlex.join(cmd) + '\n')
        out.flush()
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              encoding='utf-8', errors='replace', bufsize=1, env=env) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                out.write(line)
                out.flush()
                print(line, end='', flush=True)
                if len(lines) < 300:
                    lines.append(line.rstrip('\n'))
            code = proc.wait()
        out.write(f'EXIT CODE: {code}\n')
    if code != 0:
        raise RuntimeError(f'Command failed with exit code {code}: {shlex.join(cmd)}; see {log}')
    return '\n'.join(lines)


def conda_cmd(envname: str, *args: str) -> list[str]:
    return ['conda', 'run', '--no-capture-output', '-n', envname, *args]


def check_runtime(root: Path, args) -> None:
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        raise RuntimeError('This release is configured for your Apple Silicon macOS setup (arm64).')
    if not (root / 'metadata').is_dir() or not (root / 'scripts').is_dir():
        raise RuntimeError('Run from the established Kampala project; required folders are missing.')
    if not (root / 'metadata' / 'analysis_sample_manifest.csv').is_file():
        raise RuntimeError('Missing metadata/analysis_sample_manifest.csv')
    if not shutil.which('conda'):
        raise RuntimeError('Conda is not on PATH. Activate a Conda shell first.')
    for tool in ('prefetch', 'vdb-validate', 'fasterq-dump'):
        if not (args.sra_bin / tool).is_file() or not os.access(args.sra_bin / tool, os.X_OK):
            raise RuntimeError(f'Missing SRA Toolkit executable: {args.sra_bin / tool}')
    need = [root/'scripts'/'10_batch_preprocess_reads.py', root/'scripts'/'09_screen_human_reads.py']
    for p in need:
        if not p.is_file():
            raise RuntimeError(f'Missing required local script: {p}')
    if sha256(need[1]) != EXPECTED_HUMAN_SCRIPT_SHA256:
        raise RuntimeError('The human-screen script has changed since review; stop to prevent a silent methods change.')
    with (root/'metadata'/'analysis_sample_manifest.csv').open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    if len(rows) != 6 or len({r.get('run_accession') for r in rows}) != 6:
        raise RuntimeError('Manifest must contain exactly the six unique reviewed runs.')
    for row in rows:
        run=row.get('run_accession')
        if run not in RUN_TO_SAMPLE:
            raise RuntimeError(f'Unexpected run in manifest: {run}')
        actual=tuple(row.get(k) for k in ('sample_id','treatment_plant','collection_source','sampling_condition'))
        if actual != RUN_TO_SAMPLE[run]:
            raise RuntimeError(f'Manifest mismatch for {run}: {actual} (expected {RUN_TO_SAMPLE[run]})')
    for envname, tools in [(args.qc_env, ['python','fastp','fastqc','multiqc']),
                           (args.align_env, ['python','bowtie2','bowtie2-build','samtools'])]:
        probe = ('import shutil,sys; tools=' + repr(tools) + ';'
                 'missing=[n for n in tools if shutil.which(n) is None];'
                 'print("ENV",sys.version.split()[0],"TOOLS",", ".join(tools));'
                 'sys.exit(0 if not missing else 7)')
        proc=subprocess.run(conda_cmd(envname,'python','-c',probe),capture_output=True,text=True)
        if proc.returncode:
            raise RuntimeError(f'Conda environment {envname} does not supply required tools. '
                               f'{proc.stdout} {proc.stderr}')
        print(f'Environment {envname}: tools found',flush=True)
    checks=[
        (args.qc_env,'fastp',['fastp','--version'],r'\bfastp\s+1\.4\.0\b'),
        (args.qc_env,'fastqc',['fastqc','--version'],r'FastQC v0\.12\.1\b'),
        (args.align_env,'bowtie2',['bowtie2','--version'],r'\bversion\s+2\.5\.4\b'),
        (args.align_env,'bowtie2-build',['bowtie2-build','--version'],r'\bversion\s+2\.5\.4\b'),
        (args.align_env,'samtools',['samtools','--version'],r'^samtools 1\.22\.1\b'),
    ]
    for envname,name,cmd,pattern in checks:
        proc=subprocess.run(conda_cmd(envname,*cmd),capture_output=True,text=True)
        output=(proc.stdout or '')+'\n'+(proc.stderr or '')
        if proc.returncode or not re.search(pattern,output,re.M):
            raise RuntimeError(f'{envname}: {name} version differs from the validated pilot: {output[:500]}')
        print(f'{name} version: approved baseline',flush=True)
    print('All preflight checks PASS. The original scripts are untouched.', flush=True)


def existing_complete(root: Path, pattern: str, record_name: str, run: str, outputs: str) -> tuple[Path,dict]|None:
    """Only trust completed records with intact checksummed paired outputs."""
    for directory in sorted((root/'metadata').glob(pattern),reverse=True):
        path=directory/record_name
        if not path.is_file():
            continue
        rec=json.loads(path.read_text())
        if rec.get('run_accession') != run or rec.get('project') != PROJECT:
            continue
        if rec.get('status') != 'COMPLETE_REVIEW_REQUIRED':
            continue
        rows=rec.get('outputs',[])
        if outputs=='screen':
            selected=[r for r in rows if r.get('mate') in ('R1','R2')]
            keys=[r.get('mate') for r in selected]
        else:
            selected=[r for r in rows if r.get('group') in ('R1_paired','R2_paired')]
            keys=[r.get('group') for r in selected]
        expected={'R1','R2'} if outputs=='screen' else {'R1_paired','R2_paired'}
        if len(selected)!=2 or set(keys)!=expected:
            raise RuntimeError(f'Completed {run} record in {directory} has invalid paired-output entries.')
        for r in selected:
            p=(root / r['path']).resolve()
            if not p.is_relative_to(root) or not p.is_file() or p.stat().st_size!=r.get('bytes'):
                raise RuntimeError(f'Completed {run} output missing or size changed: {p}')
            print(f'Checking existing completed {run} output SHA-256: {p.name}',flush=True)
            if sha256(p)!=r.get('sha256'):
                raise RuntimeError(f'Completed {run} output checksum changed: {p}. Do not overwrite.')
        return directory,rec
    return None


def locate_archive(root: Path, run: str) -> Path|None:
    top=root/'data'/'raw'/'ncbi_normalized'/run
    possibilities=[top/run/f'{run}.sra',top/f'{run}.sra']
    found=[p for p in possibilities if p.is_file()]
    if len(found)>1:
        raise RuntimeError(f'Found multiple SRA archives for {run}; review manually.')
    if any(p.is_file() for p in (top/run/f'{run}.sralite',top/f'{run}.sralite')):
        raise RuntimeError(f'SRA Lite candidate detected for {run}; review before proceeding.')
    return found[0] if found else None


def fetch_archive(root: Path, run: str, sra_bin: Path, logs: Path) -> Path:
    archive=locate_archive(root,run)
    if archive is None:
        need_space(root,25*GIB,f'{run} SRA download')
        output=root/'data'/'raw'/'ncbi_normalized'/run
        output.mkdir(parents=True,exist_ok=True)
        run_cmd([str(sra_bin/'prefetch'),'--type','sra','--output-directory',str(output),run],
                logs/f'{run}_prefetch.log')
        archive=locate_archive(root,run)
        if archive is None:
            raise RuntimeError(f'Prefetch completed but {run}.sra was not found at the expected path.')
    else:
        print(f'Reusing existing SRA archive: {rel(root,archive)}',flush=True)
    validate_log=logs/f'{run}_vdb_validate.log'
    run_cmd([str(sra_bin/'vdb-validate'),str(archive)],validate_log)
    text=validate_log.read_text(errors='replace')
    if f"Database '{run}.sra' is consistent" not in text:
        raise RuntimeError(f'No explicit consistent-database message for {run}; review validation log.')
    return archive


def fq_record(handle, mate: int, number: int):
    h=handle.readline()
    if not h:
        return None
    seq=handle.readline()
    plus=handle.readline()
    quality=handle.readline()
    if not seq or not plus or not quality:
        raise RuntimeError(f'Truncated mate {mate} FASTQ record {number}')
    h=h.rstrip(b'\r\n'); seq=seq.rstrip(b'\r\n'); plus=plus.rstrip(b'\r\n'); quality=quality.rstrip(b'\r\n')
    if not h.startswith(b'@') or not plus.startswith(b'+') or not seq or len(seq)!=len(quality):
        raise RuntimeError(f'Malformed mate {mate} FASTQ record {number}')
    if not set(seq).issubset(VALID_BASES) or any(q < 33 or q>126 for q in quality):
        raise RuntimeError(f'Unexpected base or quality symbol in mate {mate} FASTQ record {number}')
    token=h[1:].split(None,1)[0]
    suffix=token[-2:]
    if suffix in (b'/1',b'/2'):
        if suffix!=f'/{mate}'.encode():
            raise RuntimeError(f'Unexpected mate identifier in mate {mate} record {number}')
        token=token[:-2]
    if not token:
        raise RuntimeError(f'Missing read ID in mate {mate} record {number}')
    return token, quality


def validate_raw_pairs(r1: Path,r2: Path) -> dict:
    print('Checking all raw read pairs, record structure, identifiers, and stored quality values...',flush=True)
    q=[{'distinct':set(),'within_read_variable_reads':0} for _ in range(2)]
    count=0
    with r1.open('rb') as f, r2.open('rb') as g:
        while True:
            a=fq_record(f,1,count+1)
            b=fq_record(g,2,count+1)
            if a is None and b is None:
                break
            if a is None or b is None or a[0]!=b[0]:
                raise RuntimeError(f'Mate counts or read identifiers disagree near pair {count+1}')
            for i,read in enumerate((a,b)):
                scores=q[i]
                quals=read[1]
                scores['distinct'].update(quals)
                if quals and min(quals)!=max(quals):
                    scores['within_read_variable_reads']+=1
            count+=1
            if count%1_000_000==0:
                print(f'  Verified {count:,} raw pairs',flush=True)
    if count==0:
        raise RuntimeError('SRA yielded zero paired reads; stop.')
    quality={}
    for i in range(2):
        v=q[i]
        scores=sorted(k-33 for k in v['distinct'])
        if len(scores)<2 or v['within_read_variable_reads']<=0:
            raise RuntimeError(f'R{i+1} has insufficient stored quality variation; stop and review representation.')
        quality[f'R{i+1}']={'distinct_scores':scores,
                         'within_read_variable_reads':v['within_read_variable_reads']}
        print(f'R{i+1} stored quality values: {scores}',flush=True)
    return {'status':'PASS','archive_validation':'PASS','fastq_structure':'PASS',
            'paired_identifiers':'PASS','read_pairs':count,'quality':quality,
            'ena_sequence_comparison_performed':False}


def ensure_fastq(root: Path, run: str, archive: Path, sra_bin: Path, logs: Path)->tuple[Path,Path]:
    out=root/'data'/'raw'/'normalized_fastq'/run
    out.mkdir(parents=True,exist_ok=True)
    r1=out/f'{run}_1.fastq'; r2=out/f'{run}_2.fastq'
    if r1.is_file() and r2.is_file():
        print(f'Reusing existing extracted FASTQs for {run}; they will be fully revalidated.',flush=True)
        return r1,r2
    if r1.exists() or r2.exists():
        raise RuntimeError(f'Only one FASTQ file exists for {run}; possible partial extraction. Do not overwrite it.')
    reserve=max(75*GIB,archive.stat().st_size*12 + 55*GIB)
    need_space(root,reserve,f'{run} FASTQ extraction')
    temp=root/'data'/'tmp'/'fasterq'/run
    temp.mkdir(parents=True,exist_ok=True)
    run_cmd([str(sra_bin/'fasterq-dump'),'--split-files','--threads','4',
             '--progress','--outdir',str(out),'--temp',str(temp),str(archive)],
            logs/f'{run}_fasterq.log')
    if not r1.is_file() or not r2.is_file():
        raise RuntimeError(f'Missing one or both paired FASTQs for {run}; review extraction output.')
    if list(out.glob(f'{run}_3.fastq')):
        raise RuntimeError('Unexpected third read file. Stop and investigate before filtering.')
    return r1,r2


def fastqc_counts(path: Path) -> int:
    with zipfile.ZipFile(path) as z:
        names=[x for x in z.namelist() if x.endswith('/fastqc_data.txt')]
        if len(names)!=1:
            raise RuntimeError(f'Cannot parse FastQC ZIP: {path}')
        text=z.read(names[0]).decode('utf-8')
    if '##FastQC\t0.12.1' not in text:
        raise RuntimeError(f'Unexpected FastQC version in {path}')
    match=re.search(r'^Total Sequences\t([0-9,]+)$',text,re.M)
    if not match:
        raise RuntimeError(f'FastQC missing Total Sequences: {path}')
    return int(match.group(1).replace(',',''))


def find_normalized_record(root:Path,run:str)->Path|None:
    for folder in sorted((root/'metadata').glob(f'batch_normalized_{run}_*'),reverse=True):
        p=folder/'normalized_record.json'
        if not p.is_file():
            continue
        r=json.loads(p.read_text())
        if r.get('status')!='COMPLETE_REVIEW_REQUIRED' or r.get('run_accession')!=run:
            continue
        rows=r.get('normalized_fastq',[])
        if len(rows)!=2 or {x.get('mate') for x in rows}!={'R1','R2'}:
            raise RuntimeError(f'Previously completed normalized record for {run} is malformed.')
        for row in rows:
            fp=(root/row['path']).resolve()
            if not fp.is_relative_to(root) or not fp.is_file() or fp.stat().st_size!=row.get('bytes'):
                raise RuntimeError(f'Previously verified raw FASTQ changed/missing: {fp}')
            print(f'Checking existing normalized-input SHA-256: {fp.name}',flush=True)
            if sha256(fp)!=row['sha256']:
                raise RuntimeError(f'Previously verified raw FASTQ checksum changed: {fp}')
        fastqc=root/r['qc_directory']/'fastqc'
        pairs=r['batch_validation']['read_pairs']
        for i in (1,2):
            if fastqc_counts(fastqc/f'{run}_{i}_fastqc.zip') != pairs:
                raise RuntimeError('Prior FastQC read count has changed.')
        return folder
    return None


def prepare_normalized_record(root:Path,run:str,archive:Path,r1:Path,r2:Path,
                              qc_env:str,logs:Path)->Path:
    previous=find_normalized_record(root,run)
    if previous:
        print(f'Reusing previously completed full-file validation + FastQC: {previous}',flush=True)
        return previous
    stats=validate_raw_pairs(r1,r2)
    folder=root/'metadata'/f'batch_normalized_{run}_{tag()}'
    qcdir=root/'results'/'qc'/'normalized_batch'/folder.name.replace('batch_normalized_','')
    (qcdir/'fastqc').mkdir(parents=True,exist_ok=False)
    folder.mkdir(parents=True,exist_ok=False)
    run_cmd(conda_cmd(qc_env,'fastqc','--threads','2','--outdir',str(qcdir/'fastqc'),str(r1),str(r2)),
            logs/f'{run}_raw_fastqc.log')
    for i in (1,2):
        count=fastqc_counts(qcdir/'fastqc'/f'{run}_{i}_fastqc.zip')
        if count!=stats['read_pairs']:
            raise RuntimeError(f'{run} FastQC R{i} count does not match validated FASTQ pair count.')
    inputs=[]
    for mate,path in zip(('R1','R2'),(r1,r2)):
        path.chmod(0o444)
        inputs.append({'mate':mate,'path':rel(root,path),'bytes':path.stat().st_size,'sha256':sha256(path)})
    record={'project':PROJECT,'run_accession':run,'sample_alias':RUN_TO_SAMPLE[run][0],
            'status':'COMPLETE_REVIEW_REQUIRED','started_utc':now(),
            'source_format':'NCBI SRA Normalized Format archive',
            'archive_path':rel(root,archive),'archive_sha256':sha256(archive),
            'normalized_fastq':inputs,'batch_validation':stats,
            'qc_directory':rel(root,qcdir),
            'interpretation':'Full-file FASTQ structure/pairing and variable stored qualities checked. '
                             'No ENA-vs-normalized identity comparison performed for this run.'}
    (folder/'normalized_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (folder/'normalized_summary.txt').write_text(
        f'RUN: {run} | {RUN_TO_SAMPLE[run][0]}\n'
        f'Validated raw pairs: {stats["read_pairs"]:,}\n'
        f'Full-file FASTQ structure: PASS\nPaired identifiers: PASS\n'
        f'Archive structural validation: PASS\nStored quality variation: PASS (both mates)\n'
        f'ENA read-sequence equivalence: NOT TESTED for this run\n'
        f'Archive: {rel(root,archive)}\nFastQC: {rel(root,qcdir)}\n')
    return folder


def record_counts(pr:dict, hr:dict)->dict:
    return {'raw_pairs':pr['fragment_accounting']['input_pairs'],
            'qc_pairs':pr['fragment_accounting']['retained_pairs'],
            'human_screened_pairs':hr['fragment_accounting']['retained_pairs_both_unmapped']}


def write_summary(root:Path,folder:Path,rows:list[dict],status:str,error:str=''):
    fields=['run_accession','sample_id','treatment_plant','status','raw_pairs','qc_pairs',
            'human_screened_pairs','preprocessing_report','human_screen_report','message']
    with (folder/'sample_status.tsv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,delimiter='\t',fieldnames=fields)
        writer.writeheader();writer.writerows(rows)
    lines=[f'PROJECT: {PROJECT}',f'STATUS: {status}',f'TIME: {now()}',
           'Normalized SRA -> full-file FASTQ checks -> FastQC -> fastp 1.4.0 -> GRCh38.p14 Bowtie2 2.5.4 + SAMtools 1.22.1',
           'Samples are never pooled. All samples were collected on dry days per original publication.',
           'All filtered single mates are stored separately and are NOT included in the primary paired reads.',
           'Human-screened means neither mate aligned to the specified reference with the stated settings.',
           'No protozoan identification, pathogen validation, infection inference, or prior host-removal verification is claimed.',
           '', 'Per-run status and numeric counts: sample_status.tsv']
    if error:lines.extend(['','ERROR: '+error])
    (folder/'batch_summary.txt').write_text('\n'.join(lines)+'\n')
    with zipfile.ZipFile(folder/'batch_review.zip','w',zipfile.ZIP_DEFLATED) as z:
        for name in ('batch_summary.txt','sample_status.tsv'):
            z.write(folder/name,name)
        for row in rows:
            for key,name in [('preprocessing_report','preprocess_summary.txt'),
                             ('human_screen_report','human_screen_summary.txt')]:
                raw=row.get(key)
                if not raw:continue
                p=(root/raw)/name
                if p.is_file(): z.write(p,f'{row["run_accession"]}/{name}')
    print(f'Batch review report: {rel(root,folder)}/batch_review.zip',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--preflight',action='store_true',help='Inspect setup without downloading/processing reads')
    ap.add_argument('--qc-env',default='kampala-qc')
    ap.add_argument('--align-env',default='kampala-align')
    ap.add_argument('--sra-bin',type=Path,default=Path.home()/'bioinformatics'/'tools'/'sratoolkit.3.4.1-mac-arm64'/'bin')
    args=ap.parse_args()
    root=Path(__file__).resolve().parent.parent
    args.sra_bin=args.sra_bin.expanduser().resolve()
    check_runtime(root,args)
    if args.preflight:
        print('Preflight complete. Run python scripts/10_batch_pipeline.py to execute the workflow.')
        return 0
    logdir=root/'metadata'/f'complete_batch_{tag()}'
    (logdir/'logs').mkdir(parents=True,exist_ok=False)
    rows=[];failure='';status='COMPLETE'
    lockfile=root/'metadata'/'.complete_batch.lock'
    with lockfile.open('a+') as lk:
        try:
            fcntl.flock(lk,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another batch workflow already holds the project lock. Stop; do not overlap runs.')
        try:
            for run in RUNS:
                alias,plant,*_=RUN_TO_SAMPLE[run]
                row={'run_accession':run,'sample_id':alias,'treatment_plant':plant,
                     'status':'STARTED','raw_pairs':'','qc_pairs':'','human_screened_pairs':'',
                     'preprocessing_report':'','human_screen_report':'','message':''}
                rows.append(row)
                print(f'\n{"="*66}\nRUN {run} | {alias} | {plant}\n{"="*66}',flush=True)
                if run==PILOT:
                    priorp=existing_complete(root,f'preprocess_{run}_*','preprocess_record.json',run,'preprocess')
                    priorh=existing_complete(root,f'human_screen_{run}_*','human_screen_record.json',run,'screen')
                    if not priorp or not priorh:
                        raise RuntimeError('RES001 pilot validated records/outputs missing. Stop; do not reprocess automatically.')
                    pdir,pr=priorp;hdir,hr=priorh
                    row.update(record_counts(pr,hr));row.update(
                        status='REUSED_VERIFIED',preprocessing_report=rel(root,pdir),
                        human_screen_report=rel(root,hdir),message='Previously completed pilot; outputs rechecked')
                    write_summary(root,logdir,rows,'RUNNING')
                    continue
                priorp=existing_complete(root,f'preprocess_{run}_*','preprocess_record.json',run,'preprocess')
                priorh=existing_complete(root,f'human_screen_{run}_*','human_screen_record.json',run,'screen')
                if priorh and not priorp:
                    raise RuntimeError(f'Human-screen result exists for {run} but preprocessing record is missing.')
                if priorp and priorh:
                    pdir,pr=priorp;hdir,hr=priorh
                    row.update(record_counts(pr,hr));row.update(
                        status='REUSED_VERIFIED',preprocessing_report=rel(root,pdir),
                        human_screen_report=rel(root,hdir),message='Both completed stages verified')
                    write_summary(root,logdir,rows,'RUNNING')
                    continue
                if priorp:
                    pdir,pr=priorp
                    print(f'Reusing complete fastp output for {run}; human screen is outstanding.',flush=True)
                else:
                    archive=fetch_archive(root,run,args.sra_bin,logdir/'logs')
                    r1,r2=ensure_fastq(root,run,archive,args.sra_bin,logdir/'logs')
                    normalized=prepare_normalized_record(root,run,archive,r1,r2,args.qc_env,logdir/'logs')
                    need_space(root,60*GIB,f'{run} fastp stage')
                    run_cmd(conda_cmd(args.qc_env,'python',str(root/'scripts'/'10_batch_preprocess_reads.py'),
                                      '--normalized',str(normalized),'--run',run),
                            logdir/'logs'/f'{run}_preprocess.log')
                    priorp=existing_complete(root,f'preprocess_{run}_*','preprocess_record.json',run,'preprocess')
                    if not priorp:
                        raise RuntimeError(f'{run} fastp completed but no complete preprocessing record was found.')
                    pdir,pr=priorp
                row.update(raw_pairs=pr['fragment_accounting']['input_pairs'],
                           qc_pairs=pr['fragment_accounting']['retained_pairs'],
                           preprocessing_report=rel(root,pdir),status='FILTERED')
                write_summary(root,logdir,rows,'RUNNING')
                need_space(root,80*GIB,f'{run} human-reference screening')
                run_cmd(conda_cmd(args.align_env,'python',str(root/'scripts'/'09_screen_human_reads.py'),
                                  '--preprocess',str(pdir),'--run',run),
                        logdir/'logs'/f'{run}_human_screen.log')
                priorh=existing_complete(root,f'human_screen_{run}_*','human_screen_record.json',run,'screen')
                if not priorh:
                    raise RuntimeError(f'{run} screening completed but no complete human-screen record was found.')
                hdir,hr=priorh
                retained=hr['fragment_accounting']['retained_pairs_both_unmapped']
                if retained<=0:
                    raise RuntimeError(f'{run} has no retained nonaligning pairs; requires review.')
                row.update(human_screened_pairs=retained,human_screen_report=rel(root,hdir),
                           status='COMPLETE',message='All stages and outputs validated')
                write_summary(root,logdir,rows,'RUNNING')
        except (Exception,KeyboardInterrupt) as exc:
            status='STOPPED'
            failure=f'{type(exc).__name__}: {exc}'
            print(f'\nSTOP: {failure}',flush=True)
            if rows and rows[-1]['status'] in ('STARTED','FILTERED'):
                rows[-1]['status']='STOPPED';rows[-1]['message']=failure
        finally:
            write_summary(root,logdir,rows,status,failure)
            fcntl.flock(lk,fcntl.LOCK_UN)
    if failure:
        print('Do not treat a stopped run as biologically negative. Inspect logs before rerunning.',flush=True)
        return 2
    print('ALL SIX RUNS PRESENT AND HUMAN-SCREENED (RES001 reused; five other runs processed/reused).',flush=True)
    print('Share only the small batch_review.zip, NOT sequence files.',flush=True)
    return 0

if __name__=='__main__':
    try:sys.exit(main())
    except Exception as exc:
        print(f'STOP: {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(2)
