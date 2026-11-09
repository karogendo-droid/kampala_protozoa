#!/usr/bin/env python3
"""Kampala protozoa, step 09: reference-locked human-read screening.

Run from the project scripts folder using Python 3.12, Bowtie2 2.5.4,
and SAMtools 1.22.1. Downloads only the public GRCh38.p14 reference and
its metadata. Local specimen sequences are never uploaded.

Inputs: the completed step-08 preprocess_record.json.
Primary output: intact QC-filtered pairs for which BOTH mates are unmapped
against this human reference under the specified alignment settings.
Unpaired step-08 reads are not processed in this primary analysis.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import fcntl
import gzip
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
import time
import urllib.error
import urllib.request
import zipfile

VERSION = '1.0.0'
PROJECT = 'PRJNA1091428'
ASSEMBLY = 'GCF_000001405.40'
LABEL = 'GRCh38.p14'
STEM = f'{ASSEMBLY}_{LABEL}'
BASE = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/{STEM}/'
FA_GZ = f'{STEM}_genomic.fna.gz'
ASSEMBLY_REPORT = f'{STEM}_assembly_report.txt'
BOWTIE_VERSION = '2.5.4'
SAMTOOLS_VERSION = '1.22.1'
BLOCK = 4 * 1024 * 1024
GIB = 1024**3


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def stamp() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')


def hash_file(p: Path, algorithm: str = 'sha256') -> str:
    h = hashlib.new(algorithm)
    with p.open('rb') as f:
        for b in iter(lambda: f.read(BLOCK), b''):
            h.update(b)
    return h.hexdigest()


def save_json(p: Path, obj: dict) -> None:
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')
    tmp.replace(p)


def inside(root: Path, value: str | Path) -> Path:
    p = Path(value).expanduser()
    p = (root / p).resolve() if not p.is_absolute() else p.resolve()
    if not p.is_relative_to(root.resolve()):
        raise ValueError('A recorded path lies outside this project; review it before proceeding.')
    return p


def tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def load_inputs(root: Path, source: Path, run: str) -> tuple[dict, dict, int]:
    r = json.loads((source / 'preprocess_record.json').read_text())
    if r.get('project') != PROJECT or r.get('run_accession') != run:
        raise ValueError('The preprocessing record is for another project or run.')
    if r.get('status') != 'COMPLETE_REVIEW_REQUIRED' or not r.get('normalized_inputs_unchanged_sha256'):
        raise ValueError('Step 08 is not recorded as successfully completed.')
    n = r.get('fragment_accounting', {}).get('retained_pairs')
    if type(n) is not int or n < 1:
        raise ValueError('Invalid number of QC-filtered input pairs.')
    selected = {}
    for x in r.get('outputs', []):
        g = x.get('group')
        if g not in ('R1_paired', 'R2_paired'):
            continue
        if g in selected:
            raise ValueError('Duplicate paired-output record.')
        p = inside(root, x['path'])
        if not p.is_file() or p.stat().st_size != x.get('bytes'):
            raise ValueError(f'{g}: missing file or unexpected file size.')
        if not re.fullmatch(r'[0-9a-f]{64}', x.get('sha256', '')):
            raise ValueError(f'{g}: invalid input checksum.')
        if r.get('output_counts', {}).get(g, {}).get('reads') != n:
            raise ValueError('Recorded paired-read counts do not reconcile.')
        selected[g] = {**x, 'absolute': p}
    if set(selected) != {'R1_paired', 'R2_paired'}:
        raise ValueError('Exactly two paired inputs are required; do not use the unpaired files.')
    return r, selected, n


def parse_checksums(text: str) -> dict[str, str]:
    out = {}
    for line in text.splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or not re.fullmatch(r'[0-9a-fA-F]{32}', fields[0]):
            continue
        name = fields[1].lstrip('*')
        if name.startswith('./'):
            name = name[2:]
        if name in out and out[name] != fields[0].lower():
            raise ValueError('Conflicting reference MD5 entries.')
        out[name] = fields[0].lower()
    for name in (FA_GZ, ASSEMBLY_REPORT):
        if name not in out:
            raise ValueError(f'Official checksum list does not include {name}.')
    return out


def fetch_small(url: str, path: Path) -> None:
    """Atomic public-metadata retrieval. No specimen data in request."""
    req = urllib.request.Request(url, headers={'User-Agent': 'KampalaProtozoa/1.0', 'Accept-Encoding': 'identity'})
    with urllib.request.urlopen(req, timeout=90) as r:
        if r.geturl().split(':', 1)[0] != 'https':
            raise ValueError('Refusing a non-HTTPS metadata redirect.')
        content = r.read(8 * 1024 * 1024 + 1)
        if len(content) > 8 * 1024 * 1024:
            raise ValueError('Unexpectedly large reference metadata response.')
    temp = path.with_name(path.name + '.part')
    temp.write_bytes(content)
    temp.replace(path)


def download_verified(url: str, path: Path, expected_md5: str, emit) -> None:
    """Resume a public reference download; trust it only after MD5 verification."""
    if path.exists():
        emit(f'Checking cached reference file: {path.name}')
        if hash_file(path, 'md5') != expected_md5:
            raise ValueError(f'Cached reference checksum mismatch: {path.name}. No automatic replacement.')
        return
    part = path.with_name(path.name + '.part')
    for attempt in range(3):
        offset = part.stat().st_size if part.exists() else 0
        headers = {'User-Agent': 'KampalaProtozoa/1.0', 'Accept-Encoding': 'identity'}
        if offset:
            headers['Range'] = f'bytes={offset}-'
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                if not r.geturl().startswith('https://'):
                    raise ValueError('Refusing a non-HTTPS reference redirect.')
                status = r.getcode()
                if status == 206:
                    m = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', r.headers.get('Content-Range', ''))
                    if not m or int(m[1]) != offset:
                        raise ValueError('Invalid partial-download response; retained partial is not trusted.')
                    expected_total, mode = int(m[3]), 'ab'
                elif status == 200:
                    # Server ignored Range: restart safely, rather than append a full response.
                    offset, mode = 0, 'wb'
                    expected_total = int(r.headers['Content-Length']) if r.headers.get('Content-Length') else None
                else:
                    raise ValueError(f'Unexpected reference response: {status}')
                done = offset
                next_notice = done + 128 * 1024**2
                emit(f'Downloading {path.name}; starting at {done / 1e6:.1f} MB.')
                with part.open(mode) as f:
                    while True:
                        b = r.read(BLOCK)
                        if not b:
                            break
                        f.write(b)
                        done += len(b)
                        if done >= next_notice:
                            emit(f'  Received {done / 1e6:.1f} MB of reference data.')
                            next_notice += 128 * 1024**2
                if expected_total is not None and done != expected_total:
                    raise OSError('Reference transfer ended before the advertised length.')
            break
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and part.exists() and hash_file(part, 'md5') == expected_md5:
                break
            if exc.code not in (408, 429, 500, 502, 503, 504) or attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
    if hash_file(part, 'md5') != expected_md5:
        raise ValueError(f'Reference checksum failed: {part.name}. Stop; do not use this partial file.')
    part.replace(path)
    path.chmod(0o444)
    emit(f'  Official MD5 verified: {path.name}')


def assembly_lengths(text: str) -> dict[str, int]:
    if ASSEMBLY not in text or LABEL not in text:
        raise ValueError('The downloaded assembly report does not identify the selected version.')
    if not re.search(r'^#\s*Taxid:\s*9606\s*$', text, re.M):
        raise ValueError('Assembly report does not declare Homo sapiens taxon 9606.')
    lengths = {}
    for line in text.splitlines():
        if not line or line.startswith('#'):
            continue
        a = line.split('\t')
        if len(a) < 9:
            raise ValueError('Unexpected assembly-report column structure.')
        if a[6] in ('na', 'n/a', '-'):
            continue
        if a[6] in lengths:
            raise ValueError('Duplicate RefSeq sequence identifier in assembly report.')
        lengths[a[6]] = int(a[8])
    required = {'NC_000001.11', 'NC_000023.11', 'NC_000024.10', 'NC_012920.1'}
    if not required.issubset(lengths):
        raise ValueError('Expected nuclear and mitochondrial reference sequences are missing.')
    if not 3_000_000_000 < sum(lengths.values()) < 4_000_000_000:
        raise ValueError('Unexpected total reference sequence length.')
    return lengths


def fai_lengths(path: Path) -> dict[str, int]:
    out = {}
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if len(f) < 5 or f[0] in out:
            raise ValueError('Malformed or duplicate FASTA-index record.')
        out[f[0]] = int(f[1])
    return out


def run_command(cmd: list[str], name: str, reports: Path, rec: dict, emit, stdout_path: Path | None = None) -> str:
    emit('Running: ' + shlex.join(map(str, cmd)))
    entry = {'name': name, 'argv': list(map(str, cmd)), 'started_utc': utc()}
    rec['commands'].append(entry)
    err = reports / 'logs' / f'{name}.stderr.txt'
    out = stdout_path or reports / 'logs' / f'{name}.stdout.txt'
    with out.open('wb') as oh, err.open('wb') as eh:
        p = subprocess.Popen(cmd, stdout=oh, stderr=eh)
        try:
            code = p.wait()
        except BaseException:
            p.terminate()
            try:
                p.wait(timeout=10)
            except subprocess.TimeoutExpired:
                p.kill(); p.wait()
            raise
    entry.update(returncode=code, finished_utc=utc())
    if code:
        raise RuntimeError(f'{name} failed (exit {code}); inspect its local log. Do not upload sequence files.')
    return out.read_text(errors='replace') if out.stat().st_size < 4 * 1024**2 else ''


def pipeline(commands: list[list[str]], name: str, reports: Path, rec: dict, emit) -> None:
    emit('Running pipeline: ' + ' | '.join(shlex.join(c) for c in commands))
    entry = {'name': name, 'pipeline': commands, 'started_utc': utc()}
    rec['commands'].append(entry)
    procs = []
    with contextlib.ExitStack() as stack:
        sink = stack.enter_context(open(os.devnull, 'wb'))
        previous = None
        try:
            for i, cmd in enumerate(commands):
                err = stack.enter_context((reports / 'logs' / f'{name}_{i+1}.stderr.txt').open('wb'))
                p = subprocess.Popen(cmd, stdin=previous, stdout=subprocess.PIPE if i < len(commands)-1 else sink, stderr=err)
                if previous is not None:
                    previous.close()
                previous = p.stdout
                procs.append(p)
            codes = [p.wait() for p in procs]
            entry.update(returncodes=codes, finished_utc=utc())
            if any(codes):
                raise RuntimeError(f'{name} failed with exit codes {codes}; inspect local command logs.')
        finally:
            for p in procs:
                if p.poll() is None:
                    p.terminate()
            for p in procs:
                try:
                    p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    p.kill(); p.wait()


def prepare_reference(root: Path, tools: dict, reports: Path, rec: dict, emit) -> tuple[Path, dict]:
    home = root / 'references' / 'human' / STEM
    home.mkdir(parents=True, exist_ok=True)
    with (home / '.reference.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest = home / 'md5checksums.txt'
        if not manifest.exists():
            emit('Retrieving the official reference checksum list.')
            fetch_small(BASE + 'md5checksums.txt', manifest)
        sums = parse_checksums(manifest.read_text())
        for name in (ASSEMBLY_REPORT, FA_GZ):
            download_verified(BASE + name, home / name, sums[name], emit)
        expected = assembly_lengths((home / ASSEMBLY_REPORT).read_text())
        fasta = home / FA_GZ[:-3]
        fa_record = home / 'fasta_record.json'
        if fasta.exists():
            if not fa_record.exists():
                raise ValueError('Existing reference FASTA lacks a completed checksum record; review, do not overwrite.')
            old = json.loads(fa_record.read_text())
            if old.get('source_md5') != sums[FA_GZ] or hash_file(fasta) != old.get('sha256'):
                raise ValueError('Cached reference FASTA does not match its recorded checksum.')
            fasta_sha = old['sha256']
        else:
            emit('Decompressing the verified reference; specimen files are not involved.')
            temp = fasta.with_name(fasta.name + '.partial')
            h = hashlib.sha256()
            with gzip.open(home / FA_GZ, 'rb') as src, temp.open('wb') as dst:
                for block in iter(lambda: src.read(BLOCK), b''):
                    dst.write(block); h.update(block)
            temp.replace(fasta)
            fasta_sha = h.hexdigest()
            save_json(fa_record, {'assembly': ASSEMBLY, 'source_md5': sums[FA_GZ], 'sha256': fasta_sha, 'bytes': fasta.stat().st_size})
            fasta.chmod(0o444)
        run_command([tools['samtools'], 'faidx', str(fasta)], 'reference_faidx', reports, rec, emit)
        actual = fai_lengths(Path(str(fasta)+'.fai'))
        if actual != expected:
            raise ValueError('FASTA sequence names/lengths disagree with the assembly report. No indexing or filtering performed.')
        emit(f'Reference verified: {len(actual)} sequences; {sum(actual.values()):,} bases; nuclear and mitochondrial sequences included.')
        complete = home / f'bowtie2_{BOWTIE_VERSION}_index_record.json'
        suffixes = ['.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l', '.rev.1.bt2l', '.rev.2.bt2l']
        if complete.exists():
            idx = json.loads(complete.read_text())
            if idx.get('fasta_sha256') != fasta_sha or idx.get('version') != BOWTIE_VERSION:
                raise ValueError('Existing reference index is not for the locked FASTA/version.')
            prefix = inside(root, idx['prefix'])
            if len(idx.get('files', [])) != 6:
                raise ValueError('Incomplete reference-index manifest.')
            for s in suffixes:
                p = Path(str(prefix)+s)
                entries = [x for x in idx['files'] if inside(root, x['path']) == p]
                if len(entries) != 1 or not p.is_file() or hash_file(p) != entries[0]['sha256']:
                    raise ValueError('A cached Bowtie2 index file failed verification.')
            emit('Reusing the completed, checksum-verified human reference index.')
        else:
            indexdir = home / f'bowtie2_{BOWTIE_VERSION}_{stamp()}'
            indexdir.mkdir()
            prefix = indexdir / 'human'
            # Force the documented large-index path. The blockwise builder is
            # appropriate for the complete human reference; limit CPU threads.
            cmd = [tools['bowtie2-build'], '--large-index', '--threads', '4', '--packed', '--bmaxdivn', '16', str(fasta), str(prefix)]
            emit('Building the human-reference search index. Progress is written to the reference_build logs.')
            run_command(cmd, 'reference_build', reports, rec, emit)
            files = []
            for s in suffixes:
                p = Path(str(prefix)+s)
                if not p.is_file() or p.stat().st_size == 0:
                    raise ValueError('Bowtie2 did not produce a complete six-file large index.')
                files.append({'path': str(p.relative_to(root)), 'bytes': p.stat().st_size, 'sha256': hash_file(p)})
                p.chmod(0o444)
            idx = {'version': BOWTIE_VERSION, 'fasta_sha256': fasta_sha, 'prefix': str(prefix.relative_to(root)), 'files': files, 'command': cmd, 'created_utc': utc()}
            save_json(complete, idx)
        for f in (manifest, home / ASSEMBLY_REPORT, fa_record, complete):
            shutil.copy2(f, reports / 'provenance' / f.name)
        ref = {'accession': ASSEMBLY, 'name': LABEL, 'taxid': 9606, 'source_url': BASE + FA_GZ,
               'source_md5': sums[FA_GZ], 'fasta_sha256': fasta_sha, 'fasta_bytes': fasta.stat().st_size,
               'reference_sequences': len(actual), 'reference_bases': sum(actual.values()),
               'selection': 'All RefSeq sequences in this assembly report, including nuclear, mitochondrial, alt and patch records where represented.',
               'index_record': str(complete.relative_to(root)), 'index_prefix': str(prefix.relative_to(root))}
        return prefix, ref


def alignment_counts(flagstat: dict, input_pairs: int) -> dict:
    categories = ('QC-passed reads', 'QC-failed reads')
    def count(name: str) -> int:
        vals = [flagstat.get(c, {}).get(name) for c in categories]
        if any(type(v) is not int or v < 0 for v in vals):
            raise ValueError(f'Invalid SAMtools flagstat count: {name}.')
        return sum(vals)
    if count('secondary') or count('supplementary'):
        raise ValueError('Unexpected secondary/supplementary records in single-alignment output.')
    if count('primary') != 2*input_pairs or count('total') != 2*input_pairs or count('paired in sequencing') != 2*input_pairs:
        raise ValueError('Alignment record counts disagree with the processed input pairs.')
    if count('read1') != input_pairs or count('read2') != input_pairs:
        raise ValueError('Alignment mate counts do not match.')
    both_reads = count('with itself and mate mapped')
    if both_reads % 2:
        raise ValueError('Odd count of reads with both mates mapped.')
    both = both_reads // 2
    one = count('singletons')
    neither = input_pairs - both - one
    if neither < 0 or count('primary mapped') != 2*both+one:
        raise ValueError('Alignment categories do not reconcile.')
    return {'input_pairs': input_pairs, 'both_mates_mapped_pairs': both, 'exactly_one_mate_mapped_pairs': one,
            'excluded_pairs_any_mate_mapped': both+one, 'retained_pairs_both_unmapped': neither,
            'retained_percent_of_qc_pairs': 100*neither/input_pairs}


def fq_record(h, number: int) -> tuple[bytes, bytes, bytes] | None:
    a = h.readline()
    if not a:
        return None
    b, c, d = h.readline(), h.readline(), h.readline()
    a, b, c, d = (s.rstrip(b'\r\n') for s in (a,b,c,d))
    if not a.startswith(b'@') or not c.startswith(b'+') or not b or len(b) != len(d):
        raise ValueError(f'Malformed retained FASTQ record at position {number}.')
    token = a[1:].split()[0]
    if token.endswith((b'/1',b'/2')):
        token = token[:-2]
    if not token or not re.fullmatch(b'[ACGTNacgtn]+', b) or any(x<33 or x>126 for x in d):
        raise ValueError(f'Invalid retained FASTQ characters at position {number}.')
    if not 75 <= len(b) <= 151 or b.upper().count(b'N') > 5:
        raise ValueError('A retained read violates the previous sequence-length/N limits.')
    if sum(x<53 for x in d)*5 > len(d):
        raise ValueError('A retained read violates the previous stored-quality limit.')
    return token, b, d


def validate_outputs(r1: Path, r2: Path, expected: int, emit) -> dict:
    n = 0
    stats = [{'reads':0, 'bases':0, 'min_length':None, 'max_length':None} for _ in range(2)]
    with gzip.open(r1,'rb') as a, gzip.open(r2,'rb') as b:
        while True:
            x, y = fq_record(a,n+1), fq_record(b,n+1)
            if x is None and y is None:
                break
            if x is None or y is None or x[0] != y[0]:
                raise ValueError('Retained paired-read identifiers or record counts disagree.')
            n += 1
            for st, rec in zip(stats,(x,y)):
                length = len(rec[1]); st['reads'] += 1; st['bases'] += length
                st['min_length'] = length if st['min_length'] is None else min(length,st['min_length'])
                st['max_length'] = length if st['max_length'] is None else max(length,st['max_length'])
            if n % 1_000_000 == 0:
                emit(f'  Validated {n:,} retained read pairs.')
    if n != expected:
        raise ValueError('Retained FASTQ counts differ from the both-unmapped alignment count.')
    return {'pairs':n, 'R1':stats[0], 'R2':stats[1], 'gzip_fastq_pairing':'PASS'}


def run_main(args, root: Path) -> int:
    source = inside(root,args.preprocess)
    prior, inputs, n = load_inputs(root,source,args.run)
    tools = {x:shutil.which(x) for x in ('bowtie2','bowtie2-build','samtools')}
    if not all(tools.values()):
        raise ValueError('Activate kampala-align with Bowtie2 2.5.4 and SAMtools 1.22.1 installed.')
    if platform.system() == 'Darwin' and platform.machine() != 'arm64':
        raise ValueError('Use native Apple Silicon Python, not an Intel-emulated session.')
    if shutil.disk_usage(root).free < 50*GIB:
        raise ValueError('Less than the conservative 50 GiB free-space reserve; review storage first.')
    tag = f'{args.run}_{stamp()}'
    reports = root/'metadata'/('human_screen_'+tag)
    out = root/'data'/'processed'/'human_screen'/tag
    for p in (reports/'logs',reports/'provenance',out/'working'):
        p.mkdir(parents=True,exist_ok=False)
    (out/'working').chmod(0o700)
    log = reports/'human_screen_log.txt'
    def emit(text):
        print(text,flush=True)
        with log.open('a',encoding='utf-8') as f:
            f.write(text+'\n')
    rec = {'project':PROJECT, 'run_accession':args.run, 'sample_alias':prior.get('sample_alias','NOT_REPORTED'),
           'status':'STARTED','started_utc':utc(),'script_version':VERSION,'script_sha256':hash_file(Path(__file__)),
           'source_folder':str(source.relative_to(root)),'python':sys.version,'architecture':platform.machine(),
           'commands':[],'versions':{},'output_directory':str(out.relative_to(root)),
           'parameters':{'bowtie2_preset':'--very-sensitive-local','phred':33,'threads':4,'seed':42,
                         'reporting':'default single primary alignment per read; mixed/discordant search enabled',
                         'mapq_exclusion_cutoff':None,'keep_flags':13,'exclude_flags':2304,
                         'fragment_rule':'Exclude the complete pair if EITHER mate has a reported human-reference alignment.',
                         'unpaired_inputs_processed':False},
           'human_genetic_analysis_performed':False,'parasite_screening_performed':False}
    success = False
    lock = (out.parent / f'{args.run}.lock').open('a+')
    try:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        emit(f'HUMAN-REFERENCE SCREENING: {args.run} | {rec["sample_alias"]}')
        emit(f'Input: {n:,} QC-filtered pairs. Prior single mates remain separate and unscreened.')
        emit('Only public human-reference files are downloaded. No specimen reads are uploaded.')
        for name,tool in tools.items():
            text = run_command([tool,'--version'],name.replace('-','_')+'_version',reports,rec,emit)
            rec['versions'][name] = text.strip()
        if not re.search(r'\bversion\s+2\.5\.4\b',rec['versions']['bowtie2']):
            raise ValueError('This workflow is locked to Bowtie2 2.5.4.')
        if not re.search(r'\bversion\s+2\.5\.4\b',rec['versions']['bowtie2-build']):
            raise ValueError('Bowtie2 and its index builder must both be version 2.5.4.')
        if not re.search(r'^samtools 1\.22\.1\b',rec['versions']['samtools']):
            raise ValueError('This workflow is locked to SAMtools 1.22.1.')
        for key,r in inputs.items():
            emit(f'Checking input checksum: {key}')
            if hash_file(r['absolute']) != r['sha256']:
                raise ValueError(f'{key} input checksum differs from the completed preprocessing record.')
        rec['inputs'] = [{k:v for k,v in r.items() if k!='absolute'} for r in inputs.values()]
        for name in ('preprocess_record.json','preprocess_summary.txt'):
            if (source/name).is_file():
                shutil.copy2(source/name,reports/'provenance'/name)
        prefix, reference = prepare_reference(root,tools,reports,rec,emit)
        rec['reference'] = reference
        save_json(reports/'human_screen_record.json',rec)
        bam = out/'working'/'all_alignments.bam'
        bowtie = [tools['bowtie2'],'--very-sensitive-local','--phred33','--reorder','--seed','42','-p','4',
                  '-x',str(prefix),'-1',str(inputs['R1_paired']['absolute']),'-2',str(inputs['R2_paired']['absolute'])]
        view = [tools['samtools'],'view','-b','-@','2','-o',str(bam),'-']
        pipeline([bowtie,view],'human_alignment',reports,rec,emit)
        run_command([tools['samtools'],'quickcheck','-v',str(bam)],'alignment_quickcheck',reports,rec,emit)
        flags = reports/'alignment_flagstat.json'
        run_command([tools['samtools'],'flagstat','-@','2','-O','json',str(bam)],'alignment_flagstat',reports,rec,emit,flags)
        counts = alignment_counts(json.loads(flags.read_text()),n)
        rec['fragment_accounting'] = counts
        emit(f'Alignment accounting: {counts["excluded_pairs_any_mate_mapped"]:,} pairs excluded; {counts["retained_pairs_both_unmapped"]:,} pairs eligible to retain.')
        r1 = out/f'{args.run}_R1.human_screened.fastq.gz'
        r2 = out/f'{args.run}_R2.human_screened.fastq.gz'
        anomaly0 = out/'working'/'unexpected_category0.fastq'
        anomaly1 = out/'working'/'unexpected_singletons.fastq'
        # f=13 explicitly requires paired, read-unmapped AND mate-unmapped.
        # Never use --un-conc as a substitute for this both-unmapped rule.
        select = [tools['samtools'],'view','-u','-f','13','-F','2304',str(bam)]
        collate = [tools['samtools'],'collate','-u','-O','-@','2','-',str(out/'working'/'collate')]
        fastq = [tools['samtools'],'fastq','-@','2','-n','-1',str(r1),'-2',str(r2),'-0',str(anomaly0),'-s',str(anomaly1),'-']
        pipeline([select,collate,fastq],'retain_both_unmapped',reports,rec,emit)
        if any(p.exists() and p.stat().st_size for p in (anomaly0,anomaly1)):
            raise ValueError('Unexpected singleton/category-0 records after paired filtering; outputs require investigation.')
        emit('Checking every retained FASTQ pair against the expected count and prior QC limits.')
        rec['output_validation'] = validate_outputs(r1,r2,counts['retained_pairs_both_unmapped'],emit)
        for r in inputs.values():
            if hash_file(r['absolute']) != r['sha256']:
                raise ValueError('A QC-filtered input file changed during screening.')
        rec['inputs_unchanged_sha256'] = True
        rec['outputs'] = []
        for mate,p in (('R1',r1),('R2',r2)):
            rec['outputs'].append({'mate':mate,'path':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':hash_file(p),
                                   'reads':counts['retained_pairs_both_unmapped']})
            p.chmod(0o444)
        tsv(reports/'fragment_accounting.tsv',['metric','value'],[{'metric':k,'value':v} for k,v in counts.items()])
        refshare = {k:v for k,v in reference.items() if k not in ('index_record','index_prefix')}
        save_json(reports/'reference_summary.json',refshare)
        # Remove only the working directory created by THIS invocation after full
        # validation. This is not secure erasure; original inputs are unchanged.
        shutil.rmtree(out/'working')
        rec['temporary_alignment_files_removed'] = True
        rec['status'] = 'COMPLETE_REVIEW_REQUIRED'
        success = True
        summary = '\n'.join([
            'HUMAN-REFERENCE SCREEN COMPLETE - REVIEW REQUIRED',f'BioProject: {PROJECT}',f'Run: {args.run}',
            f'Sample alias: {rec["sample_alias"]}',f'Reference: {LABEL} ({ASSEMBLY}), full downloaded RefSeq assembly',
            f'Reference sequences: {reference["reference_sequences"]:,}; bases: {reference["reference_bases"]:,}',
            f'Bowtie2: {BOWTIE_VERSION}; preset: --very-sensitive-local; SAMtools: {SAMTOOLS_VERSION}',
            f'Input QC-filtered read pairs: {n:,}',f'Pairs with both mates mapped: {counts["both_mates_mapped_pairs"]:,}',
            f'Pairs with exactly one mate mapped: {counts["exactly_one_mate_mapped_pairs"]:,}',
            f'Excluded pairs (either or both mates mapped): {counts["excluded_pairs_any_mate_mapped"]:,}',
            f'Retained pairs (both mates unmapped): {counts["retained_pairs_both_unmapped"]:,}',
            f'Retention relative to QC-filtered pairs: {counts["retained_percent_of_qc_pairs"]:.4f}%',
            'Paired identifiers, FASTQ/gzip integrity and count reconciliation: PASS',
            'Input QC-filtered FASTQ checksums unchanged: YES',
            'New temporary alignment files removed after validation: YES (not a secure-erasure guarantee)',
            '', 'INTERPRETATION:',
            '- Retained = no alignment reported to this reference with these settings; not certified human-free.',
            '- Excluded = human-reference-aligning or ambiguous sequence evidence; not individually proven human molecules.',
            '- Conservative pair exclusion can also remove microbial sequences resembling human DNA.',
            '- This is not a human-genetics, parasite-detection, viability or patient study.',
            '- Deposited host_removed processing remains incompletely documented; original discarded reads cannot be restored.',
            '- Single mates saved by step 08 have NOT undergone this screen and must not be added to these paired outputs.',
            '- Reference and index are reusable for the other eligible runs; no second FastQC pass is required just for this removal step.',
            '- A zero-retained result is a review stop, not a negative parasite finding.',
            '',f'Processed files: {out.relative_to(root)}',f'Reports: {reports.relative_to(root)}',''])
        (reports/'human_screen_summary.txt').write_text(summary)
        emit('\n'+summary)
    except KeyboardInterrupt:
        rec.update(status='INTERRUPTED',error='Interrupted. Do not use partial outputs. Local temporary alignment files may remain.')
        emit('\nSTOP: '+rec['error'])
    except Exception as exc:
        rec.update(status='FAILED',error=f'{type(exc).__name__}: {exc}')
        emit('\nSTOP: '+rec['error'])
    finally:
        rec['finished_utc'] = utc()
        save_json(reports/'human_screen_record.json',rec)
        if not success:
            (reports/'human_screen_summary.txt').write_text('HUMAN-REFERENCE SCREEN NOT COMPLETED\n'+f'Run: {args.run}\nStatus: {rec["status"]}\n'+rec.get('error','')+'\n'+
                f'Reports: {reports.relative_to(root)}\nDo not upload FASTQ/BAM/SAM files. Temporary files may remain locally.\n')
        with zipfile.ZipFile(reports/'human_screen_review.zip','w',zipfile.ZIP_DEFLATED) as z:
            for name in ('human_screen_summary.txt','fragment_accounting.tsv','reference_summary.json','alignment_flagstat.json'):
                if (reports/name).is_file(): z.write(reports/name,name)
        emit(f'Review files: {reports.relative_to(root)}')
        lock.close()
    return 0 if success else 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preprocess',required=True,help='Completed metadata/preprocess_RUN_TIMESTAMP folder.')
    parser.add_argument('--run',default='SRR28434174')
    args = parser.parse_args()
    if not re.fullmatch(r'SRR\d+',args.run):
        raise ValueError('Expected an SRA run accession.')
    script = Path(__file__).resolve()
    root = script.parent.parent
    if script.parent.name!='scripts' or not (root/'metadata').is_dir():
        raise ValueError('Copy the script into the existing project scripts folder before running.')
    return run_main(args,root)


if __name__=='__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'STOP: {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(2)
