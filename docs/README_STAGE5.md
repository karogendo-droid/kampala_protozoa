# Stage 5 — Verified results synthesis and manuscript draft

Project: Enteric Protozoa Kampala · BioProject PRJNA1091428  
Reporting version: 8 November 2026

## What has been completed

Nine preserved review ZIPs from the actual completed analyses were reconciled. Stage 5 does not rerun read retrieval, preprocessing, human screening, Kraken2, BLAST, or Bowtie2. It turns archived aggregate reports into publication-oriented summary tables and figures, and a working Word manuscript draft.

This is a reporting step, not a new protozoan-detection analysis.

## Main conclusion

Six samples were analyzed (three Bugolobi, three Lubigi; 68,288,156 screened pairs). Preliminary species/genus assignments were affected by reference-specificity concerns. The final integrated marker screen found three isolated accepted marker pairs across 48 sample–genus combinations, but **none** met the predetermined corroboration rule. This does not demonstrate that protozoa are absent.

## Project organization

- `scripts/26_synthesize_stage5.py`: standalone, Python-standard-library-only reporting script.
- `metadata/stage5_inputs/`: copies of the nine archived review ZIPs already uploaded in conversation. No FASTQ, human reads or BAM files are included.
- `results/stage5/tables/`: TSV results tables, ready for R, Excel or Word.
- `results/stage5/figures/`: editable vector SVG figures plus PNG versions.
- `results/stage5/stage5_results_summary.md`: cross-checked narrative and reporting boundaries.
- `results/stage5/stage5_manifest.json`: SHA-256 source hashes and integrity record.
- `metadata/taxonomy/stage5_results_review.zip`: small output archive of the reporting tables, vectors and summary.
- `manuscript/Kampala_Protozoa_Stage5_Manuscript_Draft.docx`: editable draft with numbered citations, tables, figures and author-action placeholders.

## Install on your Mac

Download `kampala_stage5_reporting.zip` and save it in your Downloads folder. In Terminal:

```bash
cd ~/bioinformatics/kampala_protozoa
conda activate kampala-taxanomy
unzip -n ~/Downloads/kampala_stage5_reporting.zip -d .
```

The `-n` option avoids overwriting existing files. All scripts remain in `scripts/`.

## Reproduce Stage 5 on your Mac (optional)

```bash
python scripts/26_synthesize_stage5.py
```

Expected terminal ending:

```
STAGE 5 COMPLETE: all 9 source reviews reconciled
Six samples: 68,288,156 screened pairs; 2,512,116 narrow-panel Kraken2 classified
Flagged leading-candidate pairs: 6,251/9,855
Isolated marker pairs: 3; corroborated target sample–genus combinations: 0/48
```

This command reads only the embedded review ZIP files; it does **not** access the internet or reprocess sequencing reads. It refreshes its own generated tables, SVGs, summary, manifest and output ZIP in `results/stage5/` and `metadata/taxonomy/`, not other stages.

## Files to review for the manuscript

1. Read the Word draft and correct `[AUTHOR ACTION: ...]` placeholders before submitting to any journal.
2. Review Table 1 (six sample processing metrics), Table 2 (flagged alignment evidence), Table 3 (marker corroboration by genus), Figures 1–3 and supplementary TSV files.
3. Preserve the distinction between initial Kraken2 candidates and the final uncorroborated taxonomic evidence.
4. The original broader concept proposal included bacterial-indicator concordance. That analysis has not been performed and is **not** claimed as completed.
5. These data do not support an inference about community prevalence, seasonality, treatment efficacy or a difference in parasite occurrence between the plants.
6. In-silico control checks are not wet-lab controls; especially note ambiguous Cyclospora/Cystoisospora outcomes.

## Scientific notes

- The Kraken2 database was narrowly targeted and is not suitable for determining total community composition.
- Stage 3A and Kraken2 used the same underlying assemblies; Stage 3A was not an independent reference validation.
- The Stage 3D flagged-contig exclusion was conservative and limited to 12 specific contigs; it does not assess total genome database contamination.
- PR2 v5.1.1 is an SSU marker resource. The integrated screen also used four E. bieneusi ITS sequences. Marker non-corroboration in shotgun data is **not** a diagnostic negative.
- Manuscript wording and citations are a working draft requiring author approval and journal-specific review.
