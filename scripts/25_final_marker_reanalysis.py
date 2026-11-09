#!/usr/bin/env python3
"""Stage 4: final integrated, contamination-aware, genus-level marker screen.

Uses all six already host-screened FASTQ pairs; PR2 v5.1.1 SSU markers,
peer-reviewed E. bieneusi ITS GenBank accessions, and local competitors.
Never treats a reference-genome Kraken2 hit as confirmation.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PR2 = ROOT / 'references' / 'markers' / 'pr2_v5.1.1'
PR2_FA = PR2 / 'pr2_version_5.1.1_SSU_mothur.fasta'
PR2_TAX = PR2 / 'pr2_version_5.1.1_SSU_mothur.tax'
OLD_GENOMES = ROOT / 'references' / 'kraken2' / 'kampala_curated_v1' / 'staged_reference_library.fna'
OLD_BLAST = ROOT / 'metadata' / 'taxonomy' / 'stage3c_reference_specificity_review.zip'
STAGE1_ZIP = ROOT / 'metadata' / 'taxonomy' / 'classification_v1_review.zip'
OUT = ROOT / 'results' / 'confirmation' / 'stage4_integrated_marker_v1'
REFERENCES = ROOT / 'references' / 'markers' / 'integrated_stage4_v1'
REVIEW = ROOT / 'metadata' / 'taxonomy' / 'stage4_final_marker_review.zip'
RUNS = [
    ('RES001', 'SRR28434174', 'Bugolobi'),
    ('RES002', 'SRR28434173', 'Bugolobi'),
    ('RES003', 'SRR28434172', 'Lubigi'),
    ('RES004', 'SRR28434171', 'Lubigi'),
    ('RES005', 'SRR28434170', 'Bugolobi'),
    ('RES013', 'SRR28434169', 'Lubigi'),
]
TARGETS = ('Giardia', 'Cryptosporidium', 'Entamoeba', 'Cyclospora',
           'Cystoisospora', 'Blastocystis', 'Enterocytozoon', 'Encephalitozoon')
ITS_ACCESSIONS = ('KJ475401', 'KJ475402', 'KJ475403', 'KJ475404')
BACTERIAL_DECOY_ASSEMBLIES = ('GCF_000005845.2_', 'GCF_000006765.1_',
                              'GCF_000009045.1_')
# Fixed analysis thresholds; not adjusted in response to outcome.
MIN_ALIGNED_BASES = 90
MIN_IDENTITY = 0.97
COMPETING_AS_MARGIN = 8
MIN_INDEPENDENT_PAIRS = 3
MIN_DIFFERENT_100BP_START_BINS = 2
MAX_REPORTED_ALIGNMENTS = 10
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def stop(message):
    raise SystemExit('STOP: ' + message)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for data in iter(lambda: handle.read(2 * 1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf8', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def read_fasta(path):
    name, bases = None, []
    with path.open(encoding='utf8', errors='replace') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(bases)
                name, bases = line[1:], []
            else:
                if name is None:
                    stop(f'Malformed FASTA {path}: sequence before first header')
                bases.append(line.upper().replace('U', 'T'))
    if name is not None:
        yield name, ''.join(bases)


def write_fasta_record(handle, ident, seq):
    handle.write('>' + ident + '\n')
    for i in range(0, len(seq), 80):
        handle.write(seq[i:i + 80] + '\n')


def target_genus(lineage):
    levels = [v.strip().split('(')[0] for v in lineage.rstrip(';').split(';')]
    for genus in TARGETS:
        if any(re.fullmatch(re.escape(genus) + r'(?:\d+)?(?:[_ ].*)?', level, re.I)
               for level in levels):
            return genus
    return None


def find_fastqs():
    root = ROOT / 'data' / 'processed' / 'human_screen'
    if not root.is_dir():
        stop(f'Human-screened FASTQ folder missing: {root}')
    data = []
    for sample, run, plant in RUNS:
        mates = []
        for n in (1, 2):
            options = sorted(root.rglob(f'{run}_R{n}.human_screened.fastq.gz'))
            if len(options) != 1:
                stop(f'{sample}: expected exactly ONE screened R{n} file, found {len(options)}')
            if options[0].stat().st_size < 1000:
                stop(f'{sample}: FASTQ appears empty: {options[0]}')
            mates.append(options[0])
        data.append((sample, run, plant, *mates))
    return data


def preflight():
    for tool in ('bowtie2', 'bowtie2-build', 'samtools'):
        if not shutil.which(tool):
            stop(f'{tool} not installed in current environment. Activate kampala-align.')
    for path in (PR2_FA, PR2_TAX, OLD_GENOMES, OLD_BLAST, STAGE1_ZIP):
        if not path.is_file() or path.stat().st_size == 0:
            stop('Required completed-stage reference file not found: ' + str(path))
    with zipfile.ZipFile(OLD_BLAST) as z:
        if 'query_manifest.tsv' not in z.namelist():
            stop('Stage 3C archive lacks flagged-reference manifest')
    with zipfile.ZipFile(STAGE1_ZIP) as z:
        if 'classification_summary.tsv' not in z.namelist():
            stop('Stage 1 review ZIP lacks classification_summary.tsv')
    data = find_fastqs()
    print('PREFLIGHT PASS: 6 paired FASTQ datasets, PR2 marker files,')
    print('  44-genome reference library and 12 flagged-contig manifest found.')
    print('  Required programs: bowtie2 / bowtie2-build / samtools')
    print('  No genomes or FASTQs re-downloaded; no alignment performed.')
    return data



def load_screened_pair_totals():
    with zipfile.ZipFile(STAGE1_ZIP) as z:
        with z.open('classification_summary.tsv') as handle:
            reader = csv.DictReader(io.TextIOWrapper(handle, encoding='utf8'), delimiter='\t')
            out = {r['sample_id']: (r['run_accession'], int(r['input_pairs'])) for r in reader}
    for sample, run, _ in RUNS:
        if sample not in out or out[sample][0] != run or out[sample][1] <= 0:
            stop(f'Missing/mismatched screened read-pair total in Stage 1 ZIP for {sample}')
    return {sample: out[sample][1] for sample, _, _ in RUNS}


def get_its_records():
    """Retrieve only peer-reviewed GenBank ITS sequences (not samples)."""
    cache = REFERENCES / 'published_E_bieneusi_ITS.fasta'
    if cache.exists() and cache.stat().st_size > 0:
        records = list(read_fasta(cache))
        existing = {x.split()[0].split('.')[0]: seq for x, seq in records}
        if all(a in existing and len(existing[a]) >= 200 for a in ITS_ACCESSIONS):
            print('REUSING published ITS marker accessions from', cache)
            return [(a, existing[a]) for a in ITS_ACCESSIONS]
        stop(f'Cached ITS FASTA has missing/short accessions; inspect {cache}')
    print('Downloading four published E. bieneusi ITS references from GenBank...', flush=True)
    records = []
    for accession in ITS_ACCESSIONS:
        url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
               + urllib.parse.urlencode(dict(db='nuccore', id=accession,
                                              rettype='fasta', retmode='text')))
        result = None
        last_error = ''
        for attempt in range(3):
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'Kampala-Protozoa-Study/1.0'})
                with urllib.request.urlopen(req, timeout=50) as response:
                    result = response.read().decode('utf8', errors='replace').strip()
                if not result.startswith('>'):
                    raise ValueError(f'NCBI returned unexpected content for {accession}')
                break
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                last_error = str(exc)
                if attempt < 2:
                    time.sleep(3 * (attempt + 1))
        if not result or not result.startswith('>'):
            stop(f'Could not retrieve published ITS accession {accession}: {last_error}. '
                 'Do not change the predefined accession list to force a result.')
        header = result.splitlines()[0]
        if accession not in header or 'Enterocytozoon bieneusi' not in header:
            stop(f'Taxonomic identity of GenBank {accession} not verified: {header}')
        seq = ''.join(result.splitlines()[1:]).upper().replace('U', 'T')
        if len(seq) < 200 or any(c not in 'ACGTN' for c in seq):
            stop(f'GenBank {accession} has invalid bases or length {len(seq)}')
        records.append((accession, seq))
        print(f'  {accession}: {len(seq)} bases, published E. bieneusi ITS record', flush=True)
        time.sleep(0.4)
    tmp = cache.with_suffix('.fasta.partial')
    with tmp.open('w') as out:
        for acc, seq in records:
            write_fasta_record(out, acc, seq)
    tmp.replace(cache)
    return records


def flagged_contigs():
    with zipfile.ZipFile(OLD_BLAST) as z:
        with z.open('query_manifest.tsv') as handle:
            rows = csv.DictReader(io.TextIOWrapper(handle, encoding='utf8'), delimiter='\t')
            flagged = {r['reference_contig'] for r in rows}
    if len(flagged) != 12:
        stop(f'Expected 12 unique flagged contigs, found {len(flagged)}')
    return flagged


def make_reference_library():
    """All PR2 sequences compete with ITS, known bad contigs and bacterial decoys."""
    REFERENCES.mkdir(parents=True, exist_ok=True)
    its = get_its_records()
    bad = flagged_contigs()
    ref = REFERENCES / 'integrated_marker_competitive.fasta'
    meta = REFERENCES / 'integrated_marker_manifest.tsv'
    tempref = ref.with_suffix('.fasta.partial')
    tempmeta = meta.with_suffix('.tsv.partial')
    tax = {}
    with PR2_TAX.open(encoding='utf8', errors='replace') as f:
        for line in f:
            parts = line.rstrip('\n\r').split('\t')
            if len(parts) >= 2:
                tax[parts[0].split()[0]] = parts[1]
    if len(tax) < 100:
        stop('PR2 taxonomy table unexpectedly small')
    coverage = collections.Counter()
    reference_count = 0
    observed = set()
    qc = collections.Counter()
    with tempref.open('w') as fo, tempmeta.open('w', newline='') as mo:
        fields = ['sequence_id', 'group', 'taxon', 'source_accession', 'length_bases']
        writer = csv.DictWriter(mo, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        for rawid, seq in read_fasta(PR2_FA):
            key = rawid.split()[0]
            lineage = tax.get(key)
            if not lineage:
                stop(f'PR2 FASTA identifier {key} missing taxonomy')
            if len(seq) < 90:
                continue
            reference_count += 1
            ident = f'PR2_{reference_count:07d}'
            genus = target_genus(lineage)
            label = genus if genus else 'OTHER_PR2'
            if genus:
                coverage[genus] += 1
            write_fasta_record(fo, ident, seq)
            writer.writerow(dict(sequence_id=ident, group='SSU', taxon=label,
                                 source_accession=key, length_bases=len(seq)))
        if reference_count < 10000:
            stop(f'PR2 FASTA appears incomplete: {reference_count} records')
        for acc, seq in its:
            ident = 'ITS_' + acc
            write_fasta_record(fo, ident, seq)
            writer.writerow(dict(sequence_id=ident, group='ITS',
                                 taxon='Enterocytozoon', source_accession=acc,
                                 length_bases=len(seq)))
            qc['published_ITS'] += 1
        for rawid, seq in read_fasta(OLD_GENOMES):
            key = rawid.split()[0]
            if key in bad:
                # Flagged material is added as competing DECOY, never a protozoan target.
                marker = 'FLAGGED'
                observed.add(key)
            elif key.startswith(BACTERIAL_DECOY_ASSEMBLIES) or key.startswith('GRCh38_chrM'):
                marker = 'BACKGROUND'
            else:
                continue
            ident = f'DECOY_{qc["decoy_total"]+1:05d}'
            write_fasta_record(fo, ident, seq)
            writer.writerow(dict(sequence_id=ident, group='DECOY',
                                 taxon='DECOY', source_accession=key,
                                 length_bases=len(seq)))
            qc['decoy_total'] += 1
            qc[marker] += 1
    if observed != bad:
        stop('Failed to include these flagged contigs: ' + ', '.join(sorted(bad - observed)))
    if qc['BACKGROUND'] < 3 or any(coverage[g] == 0 for g in TARGETS):
        stop(f'Missing background or target PR2 reference representation: {dict(coverage)}')
    tempref.replace(ref)
    tempmeta.replace(meta)
    print(f'REFERENCES: {reference_count:,} PR2 SSU sequences; '
          f'{qc["published_ITS"]} published ITS; '
          f'{qc["FLAGGED"]} flagged contigs; '
          f'{qc["BACKGROUND"]} bacterial/human background contigs', flush=True)
    return ref, meta, coverage, qc


def index_files(prefix):
    small = [Path(str(prefix) + '.' + x + '.bt2') for x in
             ('1', '2', '3', '4', 'rev.1', 'rev.2')]
    large = [Path(str(prefix) + '.' + x + '.bt2l') for x in
             ('1', '2', '3', '4', 'rev.1', 'rev.2')]
    if all(p.exists() and p.stat().st_size for p in small):
        return small
    if all(p.exists() and p.stat().st_size for p in large):
        return large
    return []


def build_index(ref, threads):
    idx = REFERENCES / 'marker_competitive'
    stamp = REFERENCES / 'index_manifest.json'
    refhash = sha256(ref)
    if index_files(idx):
        if not stamp.is_file():
            stop('Existing index lacks fingerprint manifest; do not guess whether it matches')
        record = json.loads(stamp.read_text())
        if record.get('fasta_sha256') != refhash:
            stop('Reference FASTA differs from indexed FASTA; do not reuse index')
        print('REUSING fingerprint-verified Bowtie2 marker index')
        return idx, refhash
    print('BUILDING independent PR2+ITS+decoy Bowtie2 index...', flush=True)
    if stamp.exists():
        stop('Index fingerprint exists but index files are incomplete. Inspect previous failure.')
    result = subprocess.run(['bowtie2-build', '--threads', str(threads), str(ref), str(idx)],
                            stdout=subprocess.DEVNULL)
    if result.returncode or not index_files(idx):
        stop('Bowtie2 index build failed or index incomplete. No FASTQs were modified.')
    stamp.write_text(json.dumps({'fasta_sha256': refhash,
                                 'date_utc': datetime.now(timezone.utc).isoformat()}, indent=2))
    return idx, refhash


def load_manifest(path):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    m = {r['sequence_id']: r for r in rows}
    if len(m) != len(rows):
        stop('Duplicate internal marker identifiers')
    return m


def aligned_bases_identity(cigar, nm):
    ops = CIGAR.findall(cigar)
    aligned = sum(int(n) for n, op in ops if op in ('M', '=', 'X', 'I'))
    return aligned, 1 - nm / max(1, aligned)


def parse_sam_line(line, manifest):
    a = line.rstrip('\n').split('\t')
    if len(a) < 11:
        stop('Malformed SAM line in marker BAM')
    flag = int(a[1])
    if flag & (4 | 2048):
        return None
    refname = a[2]
    if refname not in manifest:
        stop('Unknown reference sequence in BAM: ' + refname)
    tags = {}
    for token in a[11:]:
        if token.startswith('NM:i:'):
            tags['NM'] = int(token[5:])
        elif token.startswith('AS:i:'):
            tags['AS'] = int(token[5:])
    if not {'NM', 'AS'} <= tags.keys():
        return None
    aligned, identity = aligned_bases_identity(a[5], tags['NM'])
    return dict(qname=a[0], mate=1 if flag & 64 else (2 if flag & 128 else 0),
                ref=refname, pos=int(a[3]), aligned=aligned, identity=identity,
                AS=tags['AS'], taxon=manifest[refname]['taxon'],
                group=manifest[refname]['group'])


def judge_fragment(reads):
    """One fragment classified only when a target genus beats non-target competitors."""
    per_mate = collections.defaultdict(list)
    for r in reads:
        if r['aligned'] >= 75 and r['identity'] >= 0.90:
            per_mate[r['mate']].append(r)
    target_options = collections.defaultdict(list)
    for mate, candidates in per_mate.items():
        if not candidates:
            continue
        competing_top = max(r['AS'] for r in candidates)
        for r in candidates:
            if (r['taxon'] in TARGETS and r['aligned'] >= MIN_ALIGNED_BASES
                    and r['identity'] >= MIN_IDENTITY):
                target_options[r['taxon']].append((mate, r, competing_top))
    if not target_options:
        return 'no_qualifying_target', None, None
    valid = []
    for genus, hits in target_options.items():
        # Must beat ALL comparable alignments to another target genus or any decoy.
        usable = []
        for mate, r, _ in hits:
            competitor_scores = [alt['AS'] for alt in per_mate[mate]
                                 if alt['taxon'] != genus]
            other_best = max(competitor_scores) if competitor_scores else -10**9
            if r['AS'] - other_best >= COMPETING_AS_MARGIN:
                usable.append((mate, r))
        if usable:
            valid.append((genus, usable))
    if len(valid) != 1:
        return 'ambiguous_competing_reference', None, None
    genus, hits = valid[0]
    supported_mates = {mate for mate, _ in hits}
    for mate, others in per_mate.items():
        if mate in supported_mates:
            continue
        # Avoid contradictory mate assignments, even when alignment is imperfect.
        if any(x['identity'] >= MIN_IDENTITY and x['aligned'] >= MIN_ALIGNED_BASES
               and x['taxon'] != genus for x in others):
            return 'conflicting_mates', None, None
    best = max((r for _, r in hits), key=lambda r: (r['AS'], r['aligned']))
    return 'accepted', genus, dict(read=best, mates=len(supported_mates))


def fragments_from_sam(sam_lines, manifest):
    cur = None
    buffer = []
    for line in sam_lines:
        q = line.split('\t', 1)[0]
        if cur is not None and q != cur:
            yield cur, buffer
            buffer = []
        cur = q
        parsed = parse_sam_line(line, manifest)
        if parsed is not None:
            buffer.append(parsed)
    if cur is not None:
        yield cur, buffer


def summarize_bam(bam, manifest, sample, run, plant, auditout):
    found = collections.defaultdict(lambda: dict(pairs=0, two_mate=0, 
        starts=set(), bins=set(), ids=set(), its_pairs=0, ssu_pairs=0))
    status = collections.Counter()
    proc = subprocess.Popen(['samtools', 'view', '-F', '2052', str(bam)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with gzip.open(auditout, 'wt', newline='') as audit:
        w = csv.writer(audit, delimiter='\t')
        w.writerow(['sample_id', 'read_pair_id', 'decision', 'genus',
                    'best_reference', 'best_reference_position_1based',
                    'aligned_bases', 'approx_identity', 'paired_mates_support'])
        for qname, reads in fragments_from_sam(proc.stdout, manifest):
            decision, genus, detail = judge_fragment(reads)
            status[decision] += 1
            if decision != 'accepted':
                if decision != 'no_qualifying_target':
                    w.writerow([sample, qname, decision, '', '', '', '', '', ''])
                continue
            r = detail['read']
            info = found[genus]
            info['pairs'] += 1
            info['two_mate'] += 1 if detail['mates'] == 2 else 0
            info['starts'].add((r['ref'], r['pos']))
            info['bins'].add((r['ref'], (r['pos'] - 1) // 100))
            info['ids'].add(r['ref'])
            info['its_pairs'] += 1 if r['group'] == 'ITS' else 0
            info['ssu_pairs'] += 1 if r['group'] == 'SSU' else 0
            w.writerow([sample, qname, decision, genus, r['ref'], r['pos'],
                        r['aligned'], f'{r["identity"]:.5f}', detail['mates']])
    err = proc.stderr.read()
    if proc.wait():
        stop(f'{sample}: samtools failed reading marker BAM: {err[-600:]}')
    rows = []
    for genus in TARGETS:
        a = found[genus]
        # Conservative and predeclared; only marker-level evidence, never "confirmed".
        enough = (a['pairs'] >= MIN_INDEPENDENT_PAIRS and
                  len(a['bins']) >= MIN_DIFFERENT_100BP_START_BINS and
                  len(a['starts']) >= 2)
        rows.append(dict(sample_id=sample, run_accession=run, plant=plant,
                         genus=genus, accepted_distinct_read_pairs=a['pairs'],
                         two_mate_support_pairs=a['two_mate'],
                         ssu_best_reference_pairs=a['ssu_pairs'],
                         its_best_reference_pairs=a['its_pairs'],
                         distinct_reference_ids=len(a['ids']),
                         distinct_reference_start_sites=len(a['starts']),
                         distinct_100bp_reference_start_bins=len(a['bins']),
                         interpretation='PROVISIONAL_MARKER_SUPPORT' if enough
                                        else 'NOT_CORROBORATED_AT_PRESET_THRESHOLD'))
    return rows, status


def align_sample(row, idx, refhash, threads, manifest):
    sample, run, plant, r1, r2 = row
    here = OUT / sample
    here.mkdir(parents=True, exist_ok=True)
    bam = here / 'competitive_marker.name_sorted.bam'
    audit = here / 'marker_candidate_read_ids.tsv.gz'
    report = here / 'genus_marker_summary.tsv'
    stamp = here / 'completion.json'
    inputstamp = here / 'input_signature.json'
    signature = dict(ref_sha256=refhash, r1_size=r1.stat().st_size,
                     r2_size=r2.stat().st_size, script_version='stage4_v1')
    if stamp.exists():
        prior = json.loads(stamp.read_text())
        if prior.get('input_signature') != signature or not bam.exists() or not report.exists():
            stop(f'{sample}: incomplete or mismatched prior results; do not silently overwrite')
        print(f'REUSING completed marker alignment for {sample}', flush=True)
        with report.open(newline='') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        return rows, prior['read_audit']
    if bam.exists():
        if not inputstamp.exists() or json.loads(inputstamp.read_text()) != signature:
            stop(f'{sample}: incomplete BAM with missing/mismatched input fingerprint; inspect it before rerunning')
        if subprocess.run(['samtools', 'quickcheck', '-v', str(bam)],
                          stdout=subprocess.DEVNULL).returncode:
            stop(f'{sample}: partial BAM failed integrity check; inspect it before rerunning')
        print(f'RESUMING {sample}: validating existing BAM without remapping reads', flush=True)
        rows, counters = summarize_bam(bam, manifest, sample, run, plant, audit)
        write_tsv(report, list(rows[0]), rows)
        stamp.write_text(json.dumps(dict(input_signature=signature, read_audit=dict(counters),
                completed_utc=datetime.now(timezone.utc).isoformat()), indent=2))
        return rows, dict(counters)
    inputstamp.write_text(json.dumps(signature, indent=2) + '\n')
    print(f'ALIGNING {sample} ({plant}) to PR2 + E. bieneusi ITS + contaminants...', flush=True)
    log = here / 'bowtie2_samtools.log'
    temp = here / 'competitive_marker.partial.bam'
    prefix = here / 'sort_tmp'
    cmd = ['bowtie2', '--very-sensitive-local', '-k', str(MAX_REPORTED_ALIGNMENTS),
           '--no-unal', '--threads', str(threads), '-x', str(idx),
           '-1', str(r1), '-2', str(r2)]
    with log.open('w') as err:
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err)
        p2 = subprocess.Popen(['samtools', 'view', '-b', '-F', '4', '-'],
                              stdin=p1.stdout, stdout=subprocess.PIPE, stderr=err)
        p1.stdout.close()
        p3 = subprocess.run(['samtools', 'sort', '-n', '-@', '2', '-m', '256M',
                             '-T', str(prefix), '-o', str(temp), '-'],
                            stdin=p2.stdout, stderr=err)
        p2.stdout.close()
        rc2, rc1 = p2.wait(), p1.wait()
    if rc1 or rc2 or p3.returncode:
        temp.unlink(missing_ok=True)
        stop(f'{sample}: alignment failed (bowtie2={rc1}, view={rc2}, sort={p3.returncode}). '
             f'Inspect {log}')
    if subprocess.run(['samtools', 'quickcheck', '-v', str(temp)],
                      stdout=subprocess.DEVNULL).returncode:
        stop(f'{sample}: incomplete BAM; inspect {log}')
    temp.replace(bam)
    rows, counters = summarize_bam(bam, manifest, sample, run, plant, audit)
    write_tsv(report, list(rows[0]), rows)
    stamp.write_text(json.dumps(dict(input_signature=signature, read_audit=dict(counters),
                                     completed_utc=datetime.now(timezone.utc).isoformat()), indent=2))
    print(f'FINISHED {sample}: '+', '.join(f'{r["genus"]}={r["accepted_distinct_read_pairs"]}'
          for r in rows if int(r['accepted_distinct_read_pairs']) > 0) or 'no qualifying target pairs', flush=True)
    return rows, dict(counters)


def controls(ref, idx, manifest, threads):
    """Very small reference-derived controls (not a diagnostic sensitivity assay)."""
    fa = OUT / 'reference_derived_controls.fastq'
    # Record representative target SSU, published ITS, and reference decoys.
    selected = {}
    desired = set(TARGETS) | {'ITS', 'DECOY', 'OTHER_PR2'}
    for name, seq in read_fasta(ref):
        record = manifest[name]
        role = ('ITS' if record['group'] == 'ITS' else record['taxon'])
        if role not in desired or role in selected:
            continue
        seq = seq[:150].upper()
        if len(seq) >= 100 and set(seq) <= {'A', 'C', 'G', 'T'}:
            selected[role] = (name, seq)
        if len(selected) == len(desired):
            break
    if not any(g in selected for g in TARGETS) or 'DECOY' not in selected:
        stop('No valid reference-derived controls available')
    with fa.open('w') as out:
        for role, (name, seq) in selected.items():
            out.write(f'@CTRL_{role}\n{seq}\n+\n{"I"*len(seq)}\n')
    cmd = ['bowtie2', '--very-sensitive-local', '-k', str(MAX_REPORTED_ALIGNMENTS),
           '--no-unal', '-x', str(idx), '-U', str(fa), '-p', str(min(threads, 2))]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        stop(f'Reference-derived controls failed Bowtie2: {p.stderr[-500:]}')
    results = collections.defaultdict(list)
    for line in p.stdout.splitlines():
        if line.startswith('@'):
            continue
        record = parse_sam_line(line, manifest)
        if record:
            results[record['qname']].append(record)
    report = []
    for role, (name, seq) in selected.items():
        decision, genus, _ = judge_fragment(results['CTRL_' + role])
        report.append(dict(control_class=role, sequence_id=name,
                           decision=decision, assigned_genus=genus or ''))
        if role in ('DECOY', 'OTHER_PR2') and decision == 'accepted':
            stop('Negative in-silico control was assigned to a target genus: ' + role)
    return report


def run(args):
    data = preflight()
    OUT.mkdir(parents=True, exist_ok=True)
    ref, metafile, coverage, qc = make_reference_library()
    index, refhash = build_index(ref, args.threads)
    manifest = load_manifest(metafile)
    ctrl_rows = controls(ref, index, manifest, args.threads)
    all_rows, audit_rows = [], []
    screened = load_screened_pair_totals()
    for row in data:
        rows, audit = align_sample(row, index, refhash, args.threads, manifest)
        for r in rows:
            r['screened_read_pairs'] = screened[row[0]]
            r['accepted_pairs_per_million_screened'] = f'{(int(r["accepted_distinct_read_pairs"])/screened[row[0]]*1e6):.6f}'
        all_rows.extend(rows)
        audit_rows.append(dict(sample_id=row[0], run_accession=row[1], plant=row[2],
                               **audit))
    signals = OUT / 'stage4_genus_marker_evidence.tsv'
    write_tsv(signals, list(all_rows[0]), all_rows)
    auditout = OUT / 'stage4_fragment_decision_counts.tsv'
    keys = ['sample_id', 'run_accession', 'plant', 'accepted',
            'ambiguous_competing_reference', 'conflicting_mates', 'no_qualifying_target']
    write_tsv(auditout, keys, audit_rows)
    controlsout = OUT / 'reference_derived_control_checks.tsv'
    write_tsv(controlsout, list(ctrl_rows[0]), ctrl_rows)
    covpath = OUT / 'stage4_reference_target_coverage.tsv'
    write_tsv(covpath, ['genus','pr2_ssu_reference_sequences'],
              [dict(genus=g,pr2_ssu_reference_sequences=coverage[g]) for g in TARGETS])
    methods = OUT / 'stage4_methods_and_limitations.json'
    methods.write_text(json.dumps(dict(stage='4 final contamination-aware marker reanalysis',
        created_utc=datetime.now(timezone.utc).isoformat(),
        input='Six independently host-screened whole FASTQ pairs; no Kraken2 preselection',
        screened_pair_denominators='Stage 1 classification_v1_review.zip, input_pairs',
        targets=list(TARGETS),
        reference_sources=dict(PR2='5.1.1 SSU full reference set, with non-target PR2 taxa competing',
            E_bieneusi_ITS=list(ITS_ACCESSIONS),
            decoys='12 flagged assembly contigs (Stage 3C) plus 3 bacterial genomes and human chrM'),
        reference_fasta_sha256=refhash,
        bowtie2_args='--very-sensitive-local -k 10 --no-unal paired; primary + secondary retained; name-sorted BAM',
        prespecified_requirements=dict(min_aligned_query_bases=MIN_ALIGNED_BASES,
            minimum_approx_sequence_identity=MIN_IDENTITY,
            competing_alignment_score_margin=COMPETING_AS_MARGIN,
            min_independent_read_pairs=MIN_INDEPENDENT_PAIRS,
            min_distinct_100bp_reference_start_bins=MIN_DIFFERENT_100BP_START_BINS,
            min_distinct_start_positions=2),
        caveats=['This is marker corroboration, not species/subtype confirmation or wet-lab validation.',
                 'PR2 represents SSU; E bieneusi ITS is separately represented by four published records.',
                 'Small shotgun read counts and incomplete background markers limit sensitivity and specificity.',
                 'Bowtie2 -k 10 is heuristic; not every possible alternate mapping is enumerated.',
                 'In-silico controls use reference-derived sequences and cannot measure biological sensitivity.',
                 'No certified extraction or PCR-negative controls are available for secondary analysis.',
                 'No qualifying matches cannot establish organism absence.',
                 'Three samples per treatment plant; collection modes and times differ; descriptive analysis only.']),
        indent=2) + '\n')
    REVIEW.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(REVIEW, 'w', compression=zipfile.ZIP_DEFLATED) as zipout:
        for p in (signals, auditout, controlsout, covpath, methods):
            zipout.write(p, p.name)
    print('\nSUCCESS — STAGE 4 INTEGRATED MARKER REVIEW:', REVIEW, flush=True)
    n = sum(x['interpretation'] == 'PROVISIONAL_MARKER_SUPPORT' for x in all_rows)
    print(f'Marker-supported sample-genus combinations at PRESET thresholds: {n}/48')
    print('WARNING: any marker-supported genus remains provisional, NOT confirmed species or subtype.')
    print('Upload only the review ZIP; do not upload read-level files or BAMs.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    op = p.add_mutually_exclusive_group(required=True)
    op.add_argument('--check-only', action='store_true')
    op.add_argument('--run', action='store_true')
    p.add_argument('--threads', type=int, default=4)
    args = p.parse_args()
    if not 1 <= args.threads <= 8:
        stop('--threads must be between 1 and 8')
    if args.check_only:
        preflight()
    else:
        run(args)


if __name__ == '__main__':
    main()
