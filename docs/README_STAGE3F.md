# Stage 3F — Diagnose the zero PR2 target hits

This short diagnostic **reuses six existing Stage 3E BAM files**; it does **not** re-download references, rerun Bowtie2 or change sequencing files.

## Install
From the existing project root on Mac:

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_stage3f_pr2_diagnostic.zip -d .
conda activate kampala-align
```

## Run

```bash
python scripts/24_diagnose_pr2_zero_targets.py --check-only
python scripts/24_diagnose_pr2_zero_targets.py --run
```

## Outputs

Upload `metadata/taxonomy/stage3f_pr2_zero_audit_review.zip`. This holds two TSVs and a methods JSON file, and no sequence data.

`pr2_threshold_diagnostic_by_sample_genus.tsv` provides counts of **primary mapped read alignments** to each target genus under progressively stronger filtering: any MAPQ, MAPQ >=20, alignment >=75bp, >=75bp plus identity >=90%, >=75bp plus identity >=97%, and all three criteria (MAPQ >=20, >=75bp, identity >=97%). These are nested checks, not independent counts. Two reads of one pair can contribute two alignments.

`pr2_alignment_totals_by_sample.tsv` lists total mapped primary alignments (any MAPQ), total with MAPQ >=20, and high-MAPQ non-target alignments. These help explain why Stage 3E reported zero qualifying target markers.

**Interpretation:** A low-MAPQ target alignment is only a possible sequence match, not confirmation. If there are no target alignments, or only ambiguous low-quality alignments, no organism detection claim should be made. The dataset's molecular sampling and reference database may limit sensitivity.
