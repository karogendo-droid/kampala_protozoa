# Stage 1 — taxonomic classification of all six Kampala wastewater samples

This is **only Stage 1**. No abundance estimation, independent confirmation, statistical comparison, or inference of infections is performed.

## Install

1. Save the attached `kampala_taxonomy_stage1.zip` in the Mac's `Downloads` folder.
2. Open Terminal:

```bash
cd ~/bioinformatics/kampala_protozoa
conda activate kampala-taxanomy
unzip -n ~/Downloads/kampala_taxonomy_stage1.zip -d .
```

3. Verify the input files and the previously built database without processing reads:

```bash
python scripts/16_classify_all_samples.py --check-only
```

4. If you see `PREFLIGHT PASS`, classify all six samples automatically:

```bash
python scripts/16_classify_all_samples.py
```

The pipeline classifies six samples **sequentially**, never mixes reads between plants or samples, retains the original screened FASTQs untouched, and produces one Kraken2 report and one gzip-compressed read-pair assignment file per sample. It checks that all Kraken2 pair counts match the completed human-screen records. It stores SHA-256 fingerprints of the reference database files and requires a `BUILD_COMPLETE` marker. A stopped partial run is not silently overwritten; ask for help if you see an `INCOMPLETE` or `STOP` message. Successfully completed samples will be reused on rerun after verifying their recorded checksums.

## Explanation of the Kraken2 command

The core command, run by Python for each sample, is functionally:

```bash
kraken2 \
  --db references/kraken2/kampala_curated_v1 \
  --threads 4 \
  --paired \
  --gzip-compressed \
  --confidence 0.0 \
  --minimum-hit-groups 2 \
  --report reports/EXAMPLE.report.txt \
  R1.human_screened.fastq.gz R2.human_screened.fastq.gz
```

The Python script captures standard output and compresses it to `assignments.tsv.gz` as it is generated, rather than writing a huge uncompressed TSV. `--db` points to the curated reference collection. `--threads 4` limits parallel threads. `--paired` instructs Kraken2 to classify pairs as one fragment (both mates considered). `--gzip-compressed` reads `.fastq.gz` inputs without uncompressing them to disk. `--confidence 0.0` is a **sensitive discovery pass** at Kraken2's default confidence, *not* a stringent positive-detection threshold. `--minimum-hit-groups 2` retains the default requirement for two hit groups. `--report` writes per-taxon assignment counts. The read-level output is kept locally for Stage 3 verification.

No `--classified-out` or `--unclassified-out` FASTQs are created. Saved single mates from fastp are **not** included. Existing normalized SRA archives and human-screened FASTQs remain untouched.

## Interpreting results

- **Classified pairs** means Kraken2 assigned a taxonomic identifier in this **limited 44-assembly + chrM database**.
- **Unclassified pairs** means Kraken2 found no adequate match in this specific database; it does *not* mean that no microbe was present.
- `target_taxid_signals.tsv` contains **preliminary clade-level and directly assigned pair counts** for taxids in the selected-assembly manifest. Clade counts include descendants, can overlap among parent and strain IDs, and **must not be added across taxids**. These are not confirmed detections or abundance estimates.
- The database covers a tiny subset of wastewater taxa and omits *Cystoisospora belli* nuclear references. False-positive and missed assignments are possible. We must confirm biologically interesting candidate hits by independent competitive alignment and coverage analysis in Stage 3.
- All six published sequenced samples are dry-season samples. No statistical claims or treatment-efficiency conclusions follow from this stage.

## Output locations

`results/taxonomy/kraken2_curated_v1_baseline/` contains:

- `classification_summary.tsv`: classified / unclassified pair counts per sample.
- `target_taxid_signals.tsv`: exploratory target taxid clade and direct count signals.
- `classification_methods.json`: database fingerprints and run settings.
- `SRR.../kraken2.report.txt`: per-sample Kraken2 taxonomic report.
- `SRR.../assignments.tsv.gz`: compressed Kraken2 per-read-pair assignment output **kept locally; do not upload**.
- `SRR.../completion.json`: per-sample completion record and checksums.

At successful completion the script creates a **small review ZIP** containing only taxonomy reports and aggregate tables:

`metadata/taxonomy/classification_v1_review.zip`

Upload that ZIP to ChatGPT for interpretation and to plan Stage 2. Do **not** upload FASTQs, the read-assignment files, or original sequencing data.

## Known limitations

The database's prior reference-derived smoke test was operational only; it was not held-out sensitivity/specificity validation. This stage is **discovery**, not diagnostic confirmation. Confirmatory analysis and deeper database coverage review remain essential before making parasite-presence claims.
