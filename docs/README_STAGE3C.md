# Kampala enteric protozoa — Stage 3C: NCBI reference-hotspot specificity screen

## Goal
Your Stage 3B locus audit found a repeated hotspot for *Enterocytozoon bieneusi*, better-distributed alignment support for *Blastocystis* ST6, and heterogeneous hotspots for *Cyclospora cayetanensis*. The next practical check is to learn whether the most commonly matched **public reference genome loci** also closely match other organisms in a much broader reference database. This is **not a patient-read confirmation test** and **not definitive proof of species identity**.

## Privacy and data
**Only reference-genome segments, already public at NCBI, are submitted to NCBI BLAST. No wastewater FASTQ, Kraken2 read ID, sequencing read, or BAM data are transmitted.** The script uses local Stage 3B metrics and the local curated FASTA. It uses the remote NCBI `core_nt` database and does not download it to your Mac.

The `core_nt` search can fail owing to internet or NCBI service problems. The script does not change your existing Kraken2/Bowtie2 outputs. If it fails, the reference FASTA and metadata remain on your computer for a manual/alternative BLAST search.

## Where to install
Unzip from your Downloads into your current project root:

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_stage3c_reference_screen.zip -d .
conda activate kampala-taxanomy
```

The script is installed in `scripts/20_screen_reference_hotspots.py`.

## Preflight
```bash
python scripts/20_screen_reference_hotspots.py --check-only
```

If it says `PREFLIGHT PASS`, you may continue. This step does not submit anything.

## Run
Ensure `blastn` exists:

```bash
command -v blastn
```

If it is not found, install BLAST+ in the active environment:

```bash
conda install -n kampala-taxanomy -c conda-forge -c bioconda blast -y
```

Run the targeted independent-reference screen:

```bash
python scripts/20_screen_reference_hotspots.py --run
```

This performs one remote BLAST search using up to four separate 1,000-bp genomic loci per organism. The results will appear under `results/confirmation/stage3c_reference_specificity_v1/` and a summary ZIP will appear at `metadata/taxonomy/stage3c_reference_specificity_review.zip`.

Upload that ZIP for scientific interpretation. It contains summary metadata and BLAST hit data but **no wastewater reads or raw sequences**.

## What the results mean
- `query_manifest.tsv`: the public reference contigs, coordinates, base composition, and number of sample-specific bin touches.
- `reference_hotspot_blast_summary.tsv`: best database hits and whether a competing non-genus hit occurred close to the top alignment score.
- `top20_reference_blast_hits.tsv`: top 20 BLAST hits by bit score for each query. These top hits are a *screen* only; they may omit other relevant matches.
- `stage3c_methods.json`: method, caveats, and provenance.
- The original FASTA query (local only) and complete BLAST result TSV (local only) are retained for follow-up.

**Important limitations:** This tests sequence specificity of selected **reference** loci, not whether the wastewater reads exactly match the same taxon. Detecting alternative high-scoring hits can weaken confidence in an annotation; finding none does not prove it correct. Specificity of the actual read alignments and marker coverage must still be considered before any organism is reported as confirmed. *Blastocystis* subtype claims should ultimately be cross-checked against SSU-rRNA subtype barcode resources, e.g., PubMLST. The result supports transparent sensitivity analysis before the final descriptive Bugolobi–Lubigi comparison.
