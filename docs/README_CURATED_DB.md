# Kampala enteric-protozoa: curated Kraken2 discovery database, v1

## Goal and scope

This is a **targeted, exploratory protozoan-screening database**, not a comprehensive wastewater microbial profiling database. The package was selected from the nine NCBI assembly metadata files provided in `kampala_reference_review.zip`, supplemented with three individually verified common bacterial RefSeq reference assemblies.

The manifest has **44 selected assembly accessions**: 23 target assemblies, 15 non-target close-relative assemblies, three fungi and three bacterial background genomes. It also adds the **human mitochondrial genome** from the previously downloaded public GRCh38.p14 reference. The earlier human reference screen is not repeated. All selected sequences are assigned explicit NCBI taxonomy IDs and downloaded by accession/version from the NCBI genome FTP service. No GenBank/RefSeq identical assembly duplicates are intentionally included.

**Limits:** Cystoisospora belli lacks a suitable reference genome in the four catalogs checked and remains unassessed by this assembly-based database. Blastocystis subtype 5 is not included. The database represents only a narrow fraction of wastewater bacterial, fungal and environmental-eukaryote diversity, so an assignment to a parasite is preliminary and a negative result does not mean absence. Stage 3 must independently validate candidate matches using competitive alignments, breadth of coverage and unique genomic regions. The included reference-derived smoke test checks operation only: it is *not* independent sensitivity or specificity validation.

## Files added to the project

- `metadata/taxonomy/curated_v1_selected_assemblies.tsv`: 44 selected assemblies, species/subtype, NCBI taxid, accession, download directory, and reason.
- `scripts/14_prepare_curated_db.py`: downloads, checks and combines the reference genomes; obtains NCBI taxonomy; adds explicit `kraken:taxid` identifiers; checks output and hashes.
- `scripts/14_build_curated_db.sh`: checks tools and disk space; calls the preparation script; builds the masked Kraken2 database with four threads and an 8-GB hash-table cap; runs an operational smoke test.
- `scripts/15_smoke_test_curated_db.py`: classifies a tiny set of reads cut from the same included references and reports their Kraken2 assignments. Not an independent test.

## Instructions on the Mac

Start in the **existing** `kampala-taxanomy` environment (note its spelling):

```bash
cd ~/bioinformatics/kampala_protozoa
conda activate kampala-taxanomy
unzip -n ~/Downloads/kampala_curated_db_v1.zip -d .
```

First check that low-complexity masking is available:

```bash
command -v dustmasker
```

If nothing is printed, install NCBI BLAST+ once:

```bash
conda install -n kampala-taxanomy -c conda-forge -c bioconda blast -y
```

Then **one command** downloads, validates, builds and tests the database:

```bash
bash scripts/14_build_curated_db.sh
```

This command does **not** classify wastewater FASTQs. We will write and run a separate stage-1 taxonomic classification script for all six samples after reviewing the database build summary, followed by distinct abundance, independent-validation and between-plant comparison scripts.

## How the main commands work

1. `python scripts/14_prepare_curated_db.py` reads the selected-assembly TSV, downloads each `*_genomic.fna.gz` genome from NCBI, validates the compressed FASTAs, gets NCBI `taxdump.tar.gz`, and checks all taxonomic identifiers against `nodes.dmp`. It creates one FASTA with headers such as `>GCF_000002435.2_1|kraken:taxid|5741` and saves input hashes and metrics. It uses only public references; it does not upload or read wastewater specimens.
2. `kraken2-build --add-to-library ... --db ...` adds the assembled library and uses `dustmasker` to mask low-complexity DNA. Masking reduces some ambiguous short-sequence matches.
3. `kraken2-build --build --threads 4 --max-db-size 8000000000 --db ...` constructs a compact minimizer lookup table used to assign sequencing reads to taxonomy. The `--max-db-size` limit caps the estimated *hash table* size, not the total RAM, total disk usage or number of downloaded genomes. If downsampling is required, sensitivity can decrease; inspect the build log.
4. `python scripts/15_smoke_test_curated_db.py` checks the three built database files and runs Kraken2 on reference-derived 150-bp test sequences. These are in-database examples, NOT held-out independent validation.

## Output locations

- Database: `references/kraken2/kampala_curated_v1/`
- Downloaded public reference FASTAs: `references/kraken2/kampala_curated_v1/downloads/`
- Build log: `metadata/taxonomy/curated_v1_kraken_build.log`
- Verified genome checksums: `metadata/taxonomy/curated_v1_download_validation.tsv`
- Self-check results: `metadata/taxonomy/curated_v1_reference_selfcheck.tsv`

## Safety and restart

The scripts do not read, modify or remove the original `data/raw`, `data/processed`, `results/qc`, or previous screening outputs. Existing downloaded reference genomes are reused; completed staged FASTA files have checksum markers. If the network drops, rerunning the same Bash command reuses completed public downloads. Do not execute two instances at once. After a crash during `kraken2-build --add-to-library` or `--build`, inspect the log; automatic checkpointing is not guaranteed.

**Scientific interpretation:** Reference accession availability alone does not establish organism detection. We still need independently verified wastewater hits before claiming any protozoan presence, abundance, or plant differences. Number of samples is six, all from dry-season sampling, with collection-point differences between Bugolobi and Lubigi. No seasonal analysis is possible.
