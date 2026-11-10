# Kampala wastewater — Stage 3A: Competitive alignment verification

**This is a new, separate analysis stage.** It does not redo Kraken2 classification (Stage 1), RPM calculations (Stage 2), or comparison between Bugolobi and Lubigi (Stage 4). It uses original human-screened paired FASTQs and the Stage 1 read-level taxonomic assignments already on your Mac. Scripts remain in the project's `scripts/` folder.

## What the stage can—and cannot—prove

Kraken2 classifies DNA sequences using matching short sequence patterns. Stage 3A uses a *different algorithm*, **Bowtie2 alignment**, to ask whether candidate paired-end fragments align convincingly against parasite genomes **while also allowing alignment to the non-target reference sequences we already included**.

We check: (1) how many fragments align to a genome, (2) how many align with high mapping quality, (3) whether both mates support the same reference, (4) breadth of coverage across the reference, (5) how many distinct 1,000-base genomic intervals contain support.

**This is conditional computational validation, not final independent confirmation.** Candidate reads were chosen using Kraken2 and are aligned against the *same* limited reference genomes, just with an orthogonal algorithm. High mapping quality means unique **within this panel**, not uniquely attributable in the complete environment. We will examine the Stage 3A results and, where evidence is substantial, perform Stage 3B using independent markers or a more comprehensive external comparison (such as BLAST) before claiming confirmed species. No automatic positive/negative species call is made in Stage 3A.

## Installation

Download `kampala_stage3a_alignment.zip` into your Mac's Downloads folder, then:

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_stage3a_alignment.zip -d .
conda activate kampala-align
```

`unzip -n` preserves existing files; it installs only `scripts/18_validate_candidate_alignments.py` plus this README.

## Step 1 — verify prerequisites, without processing

```bash
python scripts/18_validate_candidate_alignments.py --check-only
```

This checks the existence and relationships of the six original paired FASTQs, the Kraken2 assignments and Stage 2 tables, the curated reference genome FASTA, taxonomy records and the `bowtie2` and `samtools` programs. A pass means the file setup is ready; it does **not** classify or align reads.

## Step 2 — prepare the competitive reference index (separate operation)

```bash
python scripts/18_validate_candidate_alignments.py --build-index-only
```

`bowtie2-build --threads 4` builds a searchable index of our **44 curated reference genomes and the human mitochondrial decoy**. It does *not* download any data, reanalyze FASTQs, or repeat the Kraken2 classification. It creates a checkpoint with the reference FASTA SHA-256 checksum and refuses to overwrite an existing, unverified index.

The source FASTA is `references/kraken2/kampala_curated_v1/staged_reference_library.fna`. The index will be under `references/validation_bowtie2/curated_v1/`.

## Step 3 — independently align all candidate read pairs, in one batch

```bash
python scripts/18_validate_candidate_alignments.py --run
```

Each run is processed separately. Script 18 automatically:

1. Reads `results/taxonomy/kraken2_curated_v1_baseline/<run>/assignments.tsv.gz`. It selects read IDs whose Kraken2 taxonomy assignment belongs to **Giardia, Cryptosporidium, Entamoeba, Cyclospora, Blastocystis, Enterocytozoon, or Encephalitozoon**. The genus-level selection counts must match Stage 2's previous reports.
2. Scans the corresponding original `*.human_screened.fastq.gz` files. It writes *only selected candidate read pairs* to new files (`candidate_R1.fastq.gz` and `candidate_R2.fastq.gz`), preserving pairs and validating original read counts. These are **new copies**, not modifications to originals.
3. Runs Bowtie2 against the **combined** reference library (`--very-sensitive-local`, up to three alignments per fragment, four threads) and pipes the alignments into SAMtools to produce a sorted and indexed BAM file. The comparator references are retained; mapping is *competitive*, not one genome at a time.
4. Examines **primary mapped alignments**. A high-quality alignment must have mapping quality **MAPQ ≥20**, an aligned segment of at least **75 bases**, and estimated alignment identity **≥95%** based on `NM` edit distance. These are *screening criteria*, not an established clinical detection threshold.
5. Counts distinct read-pair identifiers with qualifying alignments, properly paired support, covered reference bases, number of 1-kb genomic bins with coverage and reference breadth percentage.
6. Writes a sample-by-assembly summary and a small shareable ZIP. It never uploads or removes private sequence files.

### Key command options explained

- `--very-sensitive-local`: Bowtie2 looks for sensitive **local** alignments, allowing clipped ends (e.g., after adapter/quality filtering).
- `-k 3`: reports up to three alignments, so some ambiguous matches are exposed instead of forcing one unique genome.
- `-1`, `-2`: FASTQ files containing the two mates from the *same* sequencing fragments.
- `samtools view -b -F 4`: converts SAM records to BAM while removing unmapped records.
- `samtools sort`: orders BAM alignments by reference coordinate to support coverage assessment and indexing.
- `samtools index`: makes a BAM index to permit later investigation of mapped genomic regions.
- `samtools view -F 2308`: excludes unmapped (4), secondary (256) and supplementary (2048) alignments for the primary-evidence summary.
- `MAPQ`: Bowtie2's mapping confidence within the *reference panel*; it is **not** a probability the organism was truly present.
- **Genome breadth**: number of distinct reference bases touched by qualifying alignments ÷ genome length. Deep coverage of one repeated segment is not the same as broad genomic evidence.

### Results and sharing

All detailed results remain on your Mac under:

`results/confirmation/stage3a_curated_v1/`

The ZIP created for upload is:

`metadata/taxonomy/stage3a_alignment_review.zip`

**Upload only this small ZIP.** Do not upload `candidate_R1.fastq.gz`, `candidate_R2.fastq.gz`, `*.bam`, `*.bai` or Stage 1 `assignments.tsv.gz`; those files contain actual sequence data or read identifiers.

The process checks for previously finished runs and refuses to overwrite incomplete outputs. If it stops, share the Terminal error before starting another copy. If it says a completed run is being reused, it verifies the previous completion checkpoint.

## Interpretation before Stage 4

- A genome with numerous high-quality pair alignments **spread across multiple regions** deserves further species-specific validation.
- A genome with a large Kraken2 count but almost no convincing alignments, or alignments confined to a tiny region, is a warning for misleading classification.
- Closely related protozoa may share highly similar loci, resulting in low MAPQ; this is inconclusive, not proof of absence.
- Missing *Cystoisospora belli* genome coverage means this stage **cannot assess that organism**.
- The outputs do **not** estimate parasite abundance, environmental viability, infection rates, seasonal differences or plant treatment efficiency.

**Validation note:** Script syntax and local parsing/extraction logic were tested with synthetic records. Bowtie2 and SAMtools are not installed in the file-building environment, so actual genome indexing and mapping must be validated by the script's own checks on your Mac. Do not interpret any species as confirmed based only on Stage 3A.
