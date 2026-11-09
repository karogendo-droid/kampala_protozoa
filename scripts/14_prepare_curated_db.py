#!/usr/bin/env python3
"""Download, check and stage the Kampala curated-v1 Kraken2 references.

Uses only public NCBI reference sequences. Does not read any wastewater FASTQ.
Requires the bundled selected-assemblies TSV and local kraken2 tools.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import tarfile
import time
import urllib.request
import urllib.error
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELECTED = ROOT / 'metadata/taxonomy/curated_v1_selected_assemblies.tsv'
DB = ROOT / 'references/kraken2/kampala_curated_v1'
TAXDUMP_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        print('Reusing:', dest.name, flush=True)
        return
    if not url.startswith('https://ftp.ncbi.nlm.nih.gov/'):
        raise ValueError('Unapproved download host: ' + url)
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + '.part')
    for attempt in range(1, 4):
        try:
            print(f'Downloading {dest.name} (attempt {attempt}/3)', flush=True)
            req = urllib.request.Request(url, headers={'User-Agent': 'KampalaProtozoaReferenceResearch/1.0'})
            with urllib.request.urlopen(req, timeout=90) as response, partial.open('wb') as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            if partial.stat().st_size == 0:
                raise IOError('Server returned an empty file')
            os.replace(partial, dest)
            return
        except (OSError, urllib.error.URLError) as exc:
            print(f'  Download attempt failed: {exc}', flush=True)
            if partial.exists():
                partial.unlink()
            if attempt == 3:
                raise RuntimeError(f'Download failed after three attempts: {url}') from exc
            time.sleep(5)


def load_rows():
    with SELECTED.open(newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    if len(rows) < 40 or len({r['assembly_accession'] for r in rows}) != len(rows):
        raise ValueError('Curated manifest is incomplete or has duplicate assemblies')
    for r in rows:
        if not r['taxid'].isdigit():
            raise ValueError('Non-numeric taxid: ' + str(r))
        if not r['ftp_path'].startswith('https://ftp.ncbi.nlm.nih.gov/genomes/all/'):
            raise ValueError('Unexpected assembly URL: ' + r['ftp_path'])
    return rows


def prepare_taxonomy(tax_ids):
    folder = DB / 'taxonomy'
    folder.mkdir(parents=True, exist_ok=True)
    archive = DB / 'downloads' / 'taxdump.tar.gz'
    download(TAXDUMP_URL, archive)
    for name in ('names.dmp', 'nodes.dmp', 'merged.dmp'):
        path = folder / name
        if not path.exists() or path.stat().st_size == 0:
            print('Extracting taxonomy:', name, flush=True)
            with tarfile.open(archive, 'r:gz') as bundle:
                member = bundle.getmember(name)
                with bundle.extractfile(member) as src, path.with_suffix('.part').open('wb') as out:
                    shutil.copyfileobj(src, out)
            os.replace(path.with_suffix('.part'), path)
    valid = set()
    with (folder / 'nodes.dmp').open(encoding='utf-8') as fh:
        for line in fh:
            valid.add(line.split('|', 1)[0].strip())
    missing = sorted(tax_ids - valid)
    if missing:
        raise RuntimeError('Taxonomy IDs absent from downloaded NCBI taxdump: ' + ','.join(missing))
    print(f'All {len(tax_ids)} selected taxids present in NCBI taxonomy.', flush=True)
    return archive


def fasta_records(path):
    """Yield (NCBI original header, sequence lines) while validating gzip structure."""
    with gzip.open(path, 'rt', encoding='ascii') as fh:
        header = None
        seq = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    if not seq:
                        raise ValueError(f'Empty FASTA sequence in {path.name}')
                    yield header, seq
                header = line[1:]
                seq = []
            else:
                if header is None:
                    raise ValueError(f'Malformed FASTA (sequence before header): {path.name}')
                seq.append(line)
        if header is None or not seq:
            raise ValueError(f'No complete FASTA sequences: {path.name}')
        yield header, seq


def prepare():
    rows = load_rows()
    manifest_hash = sha256(SELECTED)
    DB.mkdir(parents=True, exist_ok=True)
    downloads = DB / 'downloads'
    downloads.mkdir(exist_ok=True)
    archives = []
    for i, r in enumerate(rows, 1):
        accession = r['assembly_accession']
        base = r['ftp_path'].rstrip('/').split('/')[-1]
        archive = downloads / f'{base}_genomic.fna.gz'
        url = r['ftp_path'].rstrip('/') + '/' + archive.name
        print(f'[{i}/{len(rows)}] {accession} — {r["organism_name"]}', flush=True)
        download(url, archive)
        # Full gzip CRC check and sequence-length verification happens during staging.
        archives.append(archive)

    tax_archive = prepare_taxonomy({r['taxid'] for r in rows} | {'9606'})
    staged = DB / 'staged_reference_library.fna'
    marker = DB / 'staged_library.complete.json'
    if staged.exists() and marker.exists():
        prior = json.loads(marker.read_text())
        if prior.get('selected_manifest_sha256') != manifest_hash:
            raise RuntimeError('Existing library was built from a different selected-manifest: review before rebuilding')
        if prior.get('staged_sha256') == sha256(staged):
            print('Verified complete staged library; will not rewrite it.', flush=True)
            return
        raise RuntimeError('Existing staged library checksum differs from complete marker; do not overwrite automatically')
    if staged.exists():
        raise RuntimeError('Unmarked existing staged FASTA: inspect before proceeding: ' + str(staged))

    tmp = staged.with_suffix('.fna.part')
    measures = []
    try:
        with tmp.open('w', encoding='ascii', newline='\n') as out:
            for r, archive in zip(rows, archives):
                acc = r['assembly_accession']
                total = 0
                contigs = 0
                for contigs, (header, seq_lines) in enumerate(fasta_records(archive), 1):
                    original_id = header.split()[0].replace('|', '_')
                    out.write(f'>{acc}_{contigs}|kraken:taxid|{r["taxid"]} original={original_id}\n')
                    for seq in seq_lines:
                        total += len(seq)
                        out.write(seq.upper() + '\n')
                expected = r['genome_size']
                if expected.isdigit() and (total < int(expected) * 0.70 or total > int(expected) * 1.30):
                    raise RuntimeError(f'Unexpected genome FASTA size for {acc}: expected {expected}, found {total}')
                print(f'  FASTA validated: {acc} {total:,} bases, {contigs} contigs', flush=True)
                measures.append(dict(accession=acc, organism=r['organism_name'], genome_bp=total,
                                     contig_count=contigs, original_gzip_sha256=sha256(archive)))
            # Residual human mitochondrial matches have a decoy; nuclear host reads were already screened.
            human = ROOT / 'references/human/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna'
            if not human.exists():
                raise FileNotFoundError('Expected existing GRCh38.p14 FASTA missing: ' + str(human))
            found = False
            collect = False
            with human.open('rt', encoding='ascii') as inf:
                for line in inf:
                    if line.startswith('>'):
                        collect = line[1:].split()[0] == 'NC_012920.1'
                        if collect:
                            found = True
                            out.write('>GRCh38_chrM|kraken:taxid|9606 public_human_mitochondrial_reference\n')
                    elif collect:
                        out.write(line.upper())
            if not found:
                raise RuntimeError('Could not locate human mitochondrial sequence NC_012920.1 in existing reference')
        os.replace(tmp, staged)
    finally:
        if tmp.exists():
            tmp.unlink()

    report = ROOT / 'metadata/taxonomy/curated_v1_download_validation.tsv'
    with report.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(measures[0]), delimiter='\t')
        w.writeheader(); w.writerows(measures)
    data = dict(selected_manifest_sha256=manifest_hash, staged_sha256=sha256(staged),
                taxonomy_archive_sha256=sha256(tax_archive), assembly_count=len(rows),
                source='NCBI genome FTP catalogs and verified NCBI bacterial reference pages',
                human_mt='GRCh38.p14 NC_012920.1 from previously downloaded RefSeq reference')
    marker.write_text(json.dumps(data, indent=2) + '\n')
    print(f'READY: {len(rows)} checked assemblies plus human chrM in {staged}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true', help='Validate selected assembly manifest but do not download')
    args = parser.parse_args()
    if args.check_only:
        r = load_rows()
        print(f'MANIFEST PASS: {len(r)} unique assemblies, URLs and NCBI taxonomy IDs well formed')
    else:
        prepare()
