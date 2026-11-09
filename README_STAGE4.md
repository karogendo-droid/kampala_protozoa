# Kampala wastewater project — Stage 4 final integrated marker reanalysis

**Purpose:** perform one final, evidence-first evaluation of whether shotgun reads from PRJNA1091428 support *genus-level* enteric protozoan marker sequences. This is **not** designed to manufacture positive detections. Zero qualifying hits is a valid outcome. Do not infer pathogen absence, infectivity, human infection, source prevalence, or a difference between treatment plants.

## What changes relative to Stages 1–3F

- Reuses **all six human-screened paired FASTQs**; does not restrict attention to Kraken2 candidate reads and does not repeat downloads, read trimming, host screening, or genome-based Kraken2 classification.
- Uses **PR2 v5.1.1 SSU-rRNA** as a broad marker library (including *all other PR2 taxa* as competing references), rather than contaminated candidate whole-genome reference contigs.
- Adds four **peer-reviewed GenBank *Enterocytozoon bieneusi* ITS** records, KJ475401–KJ475404, to complement PR2. The accession species names and basic lengths are checked before use. ITS is not provided by an SSU-only collection.
- Competes with the **12 flagged assembly contigs** identified by Stage 3C and representative complete bacterial and human mitochondrial sequences already present in the curated genome library. Flagged contigs are treated as **negative competitor sequences**, not enteric protozoan references.
- Builds a new Bowtie2 marker index and screens all six samples. Reports up to 10 alternate alignments per read/mate instead of assuming the first hit is unique.
- Only counts a candidate genus when its marker alignment is at least 90 aligned query bases with ≥97% approximate identity and beats other target genera, non-target PR2 taxa, and the included decoys by ≥8 Bowtie2 alignment-score points. If a mate strongly supports a conflicting genus/decoy, the pair is rejected.
- A *provisional marker-support* label requires ≥3 accepted independent paired-fragment IDs, ≥2 reference start positions, and ≥2 distinct 100-bp reference start bins **per sample and genus**. These are conservative prespecified operational criteria, **not a calibrated diagnostic assay**.
- Includes basic reference-derived positive/negative pipeline checks. These **are not** extraction controls, sequencing controls, estimates of analytic sensitivity, or proof of reference specificity in all wastewater organisms.

**Limitations:** Bowtie2 `-k 10` search is heuristic, and this library does not contain all possible environmental competitors; marker hits remain *provisional*. 18S generally does not establish *Blastocystis* ST6 or distinguish all near-neighbor *Cryptosporidium* taxa. Four ITS sequences do not cover every *E. bieneusi* genotype. Some marker loci may have zero coverage in shotgun metagenomes. The prior PR2 screen found only 14 low-quality target alignments; no positive outcome is presumed. There are just three non-identical sampling observations per facility, so the study is descriptive only.

## Installation

1. Download `kampala_stage4_integrated_marker.zip` into the Mac's `~/Downloads` folder.
2. Open Terminal and run:

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_stage4_integrated_marker.zip -d .
conda activate kampala-align
```

This installs **one script** under `scripts/25_final_marker_reanalysis.py`. Previous scripts and data are not overwritten.

## Check then run

```bash
python scripts/25_final_marker_reanalysis.py --check-only
```

The check verifies exactly one R1 and R2 screened FASTQ for each run, the retained Stage 3E PR2 reference files, the original 44-genome reference FASTA, the Stage 3C flagged-contig review ZIP, the Stage 1 read-pair denominator report, and Bowtie2/SAMtools. No reads are reprocessed during preflight.

If it passes:

```bash
python scripts/25_final_marker_reanalysis.py --run --threads 4
```

On its first execution, the script downloads **only four small public *E. bieneusi* ITS marker records**. It then builds the separate marker index and processes the six sample pairs sequentially. The RNA/ITS reference downloads require connectivity to NCBI E-utilities; failure is reported explicitly. No sequencing reads are uploaded to NCBI.

Intermediate results and read-pair identifiers stay locally in:

`results/confirmation/stage4_integrated_marker_v1/`

The reusable marker reference and index are in:

`references/markers/integrated_stage4_v1/`

The small, shareable output is:

`metadata/taxonomy/stage4_final_marker_review.zip`

Upload **only this ZIP** for analysis. It contains:

- `stage4_genus_marker_evidence.tsv`: 48 sample-by-target-genus records, accepted fragment counts, SSU vs ITS support, number of distinct reference locations, accepted pairs per million originally screened pairs, **provisional** evidence labels.
- `stage4_fragment_decision_counts.tsv`: accepted, ambiguous, conflicting, and non-target-mapped read-pair decisions. This counts *read pairs with mapped marker/decoy alignments*, **not every screened pair**.
- `reference_derived_control_checks.tsv`: reference-derived simulated read classification; **not** an empirical sensitivity assessment.
- `stage4_reference_target_coverage.tsv`: target genus PR2 coverage.
- `stage4_methods_and_limitations.json`: reference hash, scoring/filtering thresholds, dataset limits and interpretation constraints.

The script preserves all older classifications and stage reports. If interrupted after an integrity-checked sample BAM is written, it can reuse that BAM if its stored input fingerprint still matches. For an unexplained mismatch or partial index, do not delete the folder; share the exact error.

## Interpretation of output

`PROVISIONAL_MARKER_SUPPORT` = the **prespecified marker-level mapping rule** was satisfied. It is still **not confirmed species/subtype identification**. `NOT_CORROBORATED_AT_PRESET_THRESHOLD` = the method found insufficient discriminatory reads; **not absence**. Species/subtype, viability, transmission, public-health prevalence, and differential treatment-plant risk remain unproven.

## Sources for marker choice

- PR2 database: https://pr2-database.org/ (version 5.1.1)
- Published *E. bieneusi* ITS accessions KJ475401–KJ475404: https://journals.asm.org/doi/10.1128/AEM.00807-14
- *Blastocystis* barcode typing and allele database: https://pubmlst.org/organisms/blastocystis-spp/more-information
