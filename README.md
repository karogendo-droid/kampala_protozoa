# Enteric protozoan signals in Kampala wastewater

**Reproducible shotgun-metagenomics workflow | Kampala, Uganda**

[NCBI BioProject PRJNA1091428](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1091428) · [Analysis scripts](scripts/) · [Sample manifest](metadata/analysis_sample_manifest.csv) · [Additional quality-assurance documentation](docs/ADDITIONAL_QA.md)

## About the study

This repository contains the bioinformatics scripts, selected public metadata, and stage-specific instructions used for a secondary analysis of published shotgun metagenomes from untreated wastewater and fecal-sludge inputs associated with the **Bugolobi** and **Lubigi** treatment facilities in Kampala, Uganda.

**Study title:** *Enteric Protozoan Detection in Wastewater from Two Treatment Plants in Kampala, Uganda: A Secondary Analysis of Published Shotgun Metagenomes.*

**Public-health question:** Can already-published wastewater shotgun-metagenomic data provide preliminary evidence of enteric protozoan-associated genetic signals relevant to community-level surveillance?

The **primary analysis** comprises Kraken2 taxonomic classification, descriptive comparison of protozoan-associated assignments, and independent reference-based ribosomal RNA/ITS marker analysis of the same sequencing reads. The analyses are exploratory and do not estimate infection prevalence or demonstrate confirmed pathogen circulation.

## Public sequencing data

Six runs from [NCBI BioProject PRJNA1091428](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1091428) were included. The [sample manifest](metadata/analysis_sample_manifest.csv) links the original study identifiers to their accession numbers and recorded sampling contexts.

| Sample | SRA run | Facility | Collection source | Sampling condition |
|---|---|---|---|---|
| RES001 | SRR28434174 | Bugolobi | Emptier truck | Dry morning |
| RES002 | SRR28434173 | Bugolobi | Environmental | Dry morning |
| RES003 | SRR28434172 | Lubigi | Environmental | Dry morning |
| RES004 | SRR28434171 | Lubigi | Emptier truck | Dry morning |
| RES005 | SRR28434170 | Bugolobi | Environmental | Dry evening |
| RES013 | SRR28434169 | Lubigi | Environmental | Dry evening |

There are **three samples per facility**. Collection sources differ, so these data support **descriptive, sample-level comparisons**, not unbiased comparisons of plant catchments, organism prevalence, pathogen concentration, or disease incidence. The original sequencing reads are *not* stored in this repository; download them from NCBI/ENA using the public accessions.

Additional run accessions and summaries are available under [`metadata/PRJNA1091428_20261101T094644565419Z/`](metadata/PRJNA1091428_20261101T094644565419Z/).

## Primary analysis workflow

| Step | Analysis | Code and instructions |
|---|---|---|
| **1** | **Retrieve sequencing data and review metadata.** Audit the six accessions and download/validate the sequencing runs. | `scripts/01_fetch_run_metadata.py`–`07_retrieve_normalized_reads.py` |
| **2** | **Quality control and preprocessing.** Assess read quality and perform adapter/quality processing as documented. | `scripts/04_run_raw_qc.py`, `08_preprocess_reads.py`, `10_batch_preprocess_reads.py`, `10_batch_pipeline.py`, `run_complete_batch.sh` · [Batch guide](docs/README_BATCH_PIPELINE.md) |
| **3** | **Human-reference screening.** Exclude paired reads mapping to the selected human reference before downstream taxonomic classification. | `scripts/09_screen_human_reads.py` · [Batch guide](docs/README_BATCH_PIPELINE.md) |
| **4** | **Kraken2 taxonomic classification.** Select references, build and check the curated database, then classify screened samples. | `scripts/11_find_reference_candidates.py`–`16_classify_all_samples.py` · [Database guide](docs/README_CURATED_DB.md) · [Classification guide](docs/README_STAGE1.md) |
| **5** | **Describe protozoan-associated signals.** Summarize taxonomic assignments by sample and compare Bugolobi and Lubigi descriptively, accounting for the number of screened read pairs. | `scripts/17_assess_candidate_signal_abundance.py` · [Signal summary guide](docs/README_STAGE2.md) |
| **6** | **Independent ribosomal RNA and ITS marker analysis.** Examine the same metagenomic reads against independently curated SSU/18S rRNA references (PR2) and *Enterocytozoon bieneusi* ITS reference sequences, using documented support thresholds. | `scripts/23_screen_pr2_18s.py`, `24_diagnose_pr2_zero_targets.py`, `25_final_marker_reanalysis.py` · [PR2 guide](docs/README_STAGE3E.md) · [Threshold diagnostics](docs/README_STAGE3F.md) · [Integrated marker guide](docs/README_STAGE4.md) |
| **7** | **Synthesize findings and public-health implications.** Integrate initial assignments and marker evidence for descriptive reporting. | `scripts/26_synthesize_stage5.py` · [Reporting guide](docs/README_STAGE5.md) |

**Script numbering reflects the historical development sequence**, not a seven-command installation or execution interface. Scripts 18–22 represent additional checks conducted during development and remain available for transparency; see [additional quality assurance](docs/ADDITIONAL_QA.md). Stage-specific guides contain the detailed commands, required inputs and outputs. A fresh-machine, end-to-end rerun has **not** yet been independently verified.

## Set up the software

The original analyses were performed on an **Apple M1 Mac with 32 GiB of RAM**. Running the workflow requires substantial additional disk space for compressed and decompressed reads, reference files, indices, and temporary outputs.

Main command-line tools included **NCBI SRA Toolkit, FastQC, MultiQC, fastp, Bowtie2, SAMtools, Kraken2, and Python**, with different software environments for different stages.

```bash
git clone https://github.com/karogendo-droid/kampala_protozoa.git
cd kampala_protozoa

# Create environments as needed; consult the corresponding YAML first.
conda env create -f envs/kampala-qc.yml
conda env create -f envs/kampala-align.yml
```

The [`envs/`](envs/) directory also includes a fastp-equipped QC environment. The original interactive environment was named `kampala-taxanomy`, which is not necessarily identical to either shared environment YAML. The commands shown above are setup examples, **not a tested universal installation recipe**. Follow the stage-specific README and the actual program arguments for each script.

## Reference resources and analysis conventions

The first Kraken2 implementation used a **purpose-curated collection of 44 selected reference assemblies plus a human mitochondrial sequence**, not a comprehensive wastewater metagenome database. This database composition limits the biological specificity of its initial taxonomic labels.

The reference selection process is documented in:

- [`metadata/taxonomy/protozoa_target_manifest.tsv`](metadata/taxonomy/protozoa_target_manifest.tsv)
- [`metadata/taxonomy/curated_v1_selected_assemblies.tsv`](metadata/taxonomy/curated_v1_selected_assemblies.tsv)
- [`metadata/taxonomy/reference_shortlist_for_review.tsv`](metadata/taxonomy/reference_shortlist_for_review.tsv)

For descriptive comparisons, **reads per million (RPM)** refers to the number of Kraken2-assigned *read-pair fragments* divided by the total human-screened read-pair fragments for the relevant sample or pooled facility, multiplied by one million. These values are **classification-signal metrics**, not measures of organism concentration, viable pathogen burden, or infection prevalence. Keep the same read-pair counting conventions when reproducing the analysis.

For independent marker analysis, document the exact PR2 version, *E. bieneusi* ITS references, competing/decoy references, alignment filters, and the fragment-support criteria used in the original scripts. Marker analysis uses the **same sequencing dataset** with a different reference system; it is not independent wet-laboratory confirmation.

## Interpretation of the study signals

The initial Kraken2 analysis generated candidate taxonomic assignments to enteric protozoan genera, including *Enterocytozoon*, *Blastocystis*, *Cyclospora*, *Entamoeba*, *Cryptosporidium*, and *Giardia*. These labels are **preliminary sequence signals**, not organism confirmations.

The independent marker analyses recovered a small number of isolated *Blastocystis*- and *Entamoeba*-associated fragments, but **no target genus met the predefined marker corroboration thresholds** across the samples. Accordingly, the findings may motivate follow-up surveillance, but they do not establish community infection, pathogen viability, or confirmed circulation. A failure to reach a marker threshold does not demonstrate absence of the organism.

Additional sequence/reference quality assessments undertaken during development materially qualify some of the initial classifications. These are retained in the code and [supplementary QA documentation](docs/ADDITIONAL_QA.md), rather than framed as separate primary study objectives.

## Repository organization

```text
kampala_protozoa/
├── README.md                   # Entry point and primary analysis workflow
├── README_BATCH_PIPELINE.md    # Data preparation and batch details
├── README_CURATED_DB.md        # Reference selection and Kraken2 database
├── README_STAGE*.md            # Original stage-specific instructions
├── scripts/                    # Original numbered Python/Bash scripts
├── envs/                       # Conda environment specifications
├── metadata/                   # Selected public run and reference manifests
└── docs/
    └── ADDITIONAL_QA.md        # Additional assessment/provenance information
```

Original data, generated FASTQ/BAM files, human-derived read outputs, local paths, logs, transient analysis archives, and private-location records should **not** be committed. Examine staged changes using `git diff --cached` before every public push. A `.gitignore` does not remove files that Git already tracks.

## Reproducibility status

This repository contains the **analysis source code, selected public manifests, environment definitions, and stage documentation**. It has **not yet passed an independent, clean-environment, end-to-end reproducibility test**. To fully reproduce a published result, researchers should additionally capture version-pinned references, hashes, commands, environment specifications, and the input/output provenance of each step. Exact result tables should be compared with documented expected outputs before declaring a successful replication.

## Data citation, pipeline citation and license

**Source dataset:** [NCBI BioProject PRJNA1091428](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1091428). Cite the original data-generating investigators and associated publication, where verified.

**Pipeline:** Cite this repository and a stable release/commit identifier. Add the manuscript citation and DOI after publication.

A software license has not been specified by this README. Confirm code ownership and third-party licensing before adding a `LICENSE` file or redistributing third-party reference data.
