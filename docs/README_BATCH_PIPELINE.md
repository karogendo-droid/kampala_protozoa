# Kampala enteric protozoa — validated preprocessing batch

This bundle orchestrates the same analysis settings already piloted on
SRR28434174 (RES001) for the other five runs of PRJNA1091428. It does **not**
perform protozoan classification yet.

## Files

- `scripts/10_batch_pipeline.py` — downloads/validates missing normalized SRA
  archives, extracts paired FASTQ, checks **every** FASTQ pair and variable
  stored quality scores, runs prefilter FastQC, performs baseline preprocessing,
  human screening, and writes a consolidated report. The five nonpilot runs
  are handled sequentially to control RAM/disk use. Can be restarted.
- `scripts/10_batch_preprocess_reads.py` — generalized copy of the previously
  reviewed `08_preprocess_reads.py`, maintaining exactly the same fastp settings
  and numeric/gzip/paired-read validations. Added a separate independently
  verified input path for runs that have not undergone ENA-vs-normalized
  sequence comparison. No comparison is falsely claimed.
- `scripts/run_complete_batch.sh` — one-command launcher.

The existing `scripts/09_screen_human_reads.py` is **reused unchanged**. Its
SHA-256 is checked so a different version cannot run silently.

## Install (without overwriting originals)

1. Save `kampala_complete_batch_pipeline.zip` in your Mac's **Downloads** folder.
2. Open Terminal and run:

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_complete_batch_pipeline.zip -d .
```

3. Check dependencies and manifest, without processing any reads:

```bash
bash scripts/run_complete_batch.sh --preflight
```

4. Start the **complete** workflow (one command):

```bash
bash scripts/run_complete_batch.sh
```

Leave this Terminal open while the workflow runs. This is a batch job on **your
own Mac**, not a background task being executed by ChatGPT. If the job stops,
read the diagnostic output; after correcting the issue, run the same command
again. Never assume stopped runs were negative for parasites.

## What gets reused

- **RES001 (SRR28434174):** validated pilot preprocessing and human-screening
  outputs are checksummed and reused; not downloaded or reprocessed.
- **RES002 (SRR28434173):** existing 1.6-GiB normalized SRA archive and
  extracted uncompressed paired FASTQs are reused after revalidation.
- **RES003, RES004, RES005, RES013:** download missing normalized SRA archives,
  extract, quality check, fastp filter, human screen.
- If any other run has a **fully documented** completed preprocessing or human
  screening record, the batch validates its paired outputs before skipping it.

## Controls and provenance

- Sequence data stay inside this project's local folders. No user FASTQs,
  alignments, or human reads are uploaded to external servers.
- Uses SRA Toolkit 3.4.1 at the previously identified **manual** install path.
- Requires existing Conda environments `kampala-qc` (fastp 1.4.0,
  FastQC 0.12.1, MultiQC) and `kampala-align` (Bowtie2 2.5.4,
  SAMtools 1.22.1). These versions are checked before any processing.
- fastp: paired overlap/Nextera adapter; four-base sliding windows Q20;
  minimum length 75 bp; no more than 20% bases below Q20; at most 5 Ns;
  poly-G/poly-X, merging, correction and deduplication OFF. Passing single
  mates are saved separately and **not** included in primary paired output.
- Human screen: reuse full GRCh38.p14 RefSeq Bowtie2 index, `--very-sensitive-local`;
  exclude the entire pair when either mate aligns. The script validates all
  retained paired FASTQ records and reconciles SAMtools counts. It cannot
  prove that a read is human-derived or that retained reads are entirely
  human-free.
- Each sample is processed separately and kept separate. No wet-season
  samples were sequenced in the published data. Previously deposited files
  were already labelled host-removed and the initial host-removal method
  remains insufficiently documented.
- Preflight enforces **50+ GiB** at the existing human-screen stage; the
  wrapper imposes larger pre-stage reserves. It stops if storage drops below
  these thresholds, without automatically deleting your files.
- Failed or interrupted steps can leave partial outputs and temporary files;
  their cleanup is **never automatic**. Do not manually delete them without
  identifying the failed stage.

## Reports

A new directory is created with each run:

```text
metadata/complete_batch_YYYYMMDDTHHMMSSffffffZ/
    batch_summary.txt
    sample_status.tsv
    batch_review.zip
    logs/
```

After the run completes, find the latest review file using:

```bash
find metadata -name batch_review.zip -print
```

Upload only `batch_review.zip` to ChatGPT, **never FASTQ, BAM, SAM, or raw
FastQC archives**. Per-sample QC, fastp, and human-screening records remain
in the established `metadata/`, `results/qc/`, and `data/processed/` folders.

**Caution:** Successful preprocessing does not establish enteric protozoa
presence. Detection, independent confirmation, interpretation and study
statistics are separate next-stage analyses.
