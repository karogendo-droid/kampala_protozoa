# Kampala Protozoa — Stage 2: exploratory candidate-signal normalization

This is a *separate* Stage 2 script. It does not run taxonomic classification, download references, remap reads, estimate genome equivalents, or compare treatment plants.

## Install

Unzip the package into `~/bioinformatics/kampala_protozoa`. It installs only `scripts/17_assess_candidate_signal_abundance.py` and this README. All scripts stay in `scripts/`.

## Run

```bash
cd ~/bioinformatics/kampala_protozoa
conda activate kampala-taxanomy
python scripts/17_assess_candidate_signal_abundance.py \
  --review-zip metadata/taxonomy/classification_v1_review.zip
```

If the review ZIP was saved elsewhere, supply its actual path after `--review-zip`. No extra Python packages are required.

## What it computes

1. Reads `classification_summary.tsv`, `target_taxid_signals.tsv`, all six reports and `classification_methods.json` **from the classification ZIP**.
2. Verifies six distinct runs, input counts, reported classified and unclassified counts, and taxon counts in per-run Kraken2 reports.
3. Calculates candidate signal in **assigned read pairs per million input human-screened read pairs (RPM)**. For instance, 100 assignments among 10 million input pairs = 10 RPM.
4. Keeps **Kraken2 genus *clade* counts** and exact-NCBI-taxon **direct assignment counts** in **separate** TSV files. Clade and direct counts must not be added to each other.
5. Saves output in `results/abundance/stage2_v1/` and a small shareable `metadata/taxonomy/abundance_v1_review.zip`.

## Command explanation

- `python`: invokes Python in the current shell; this script needs only Python standard libraries.
- `scripts/17_assess_candidate_signal_abundance.py`: path to our Stage 2 analysis script.
- `--review-zip`: tells it where to find the Stage 1 aggregate reports.
- `--out-dir`: optional alternative output folder if you intentionally repeat the analysis; the default existing nonempty folder is protected against overwrites.

## Scientific limitations

RPM is **not biological abundance**. It does not correct for genome length, DNA extraction efficiency, recovery or database completeness. Our targeted database is far from a comprehensive wastewater background; cross-species assignments are possible. The output reports preliminary sequence-assignment signals and **not confirmed organisms**, clinical infections or parasite viability. Formal treatment-plant comparisons belong to Stage 4, *after* candidate confirmation in Stage 3.
