# Final analysis and reproducibility

## Study

Investigation of Enteric Protozoan Genetic Signals in Kampala Wastewater:
Implications for Diarrhoeal Disease Surveillance

## Research objectives

1. Identify and quantify candidate enteric protozoan genetic signals.
2. Describe their distribution across wastewater specimens.
3. Evaluate molecular evidence supporting their taxonomic identification.

## Data source

NCBI BioProject: PRJNA1091428

Six paired-end sequencing runs:
SRR28434169–SRR28434174.

Three specimens originated from Bugolobi and three from Lubigi,
Kampala, Uganda.

## Analysis stages

### 1. Read preparation

The original project scripts document sequence retrieval,
quality preprocessing and human-read screening.

Final human-screened read pairs: 68,288,156.

### 2. Genomic reference preparation

The final genomic screening used 55 distinct microbial and
eukaryotic assemblies and an additional human competitor reference.

Relevant scripts:
- 10_download_enteric_refseq.py
- 11_download_enteric_expansion.py
- 12_prepare_kraken_references.py
- 13_extract_legacy_unique.py
- 14_prepare_human_competitor.py

The final assembly manifest and reference checksums must be
documented alongside this workflow.

### 3. Preliminary taxonomic screening

Relevant scripts:
- 15_build_and_screen_pilot.sh
- 16_screen_remaining.sh

Publication output:
publication_results/genus_counts_all_six.csv

Total candidate read pairs: 89,444.

### 4. Candidate genomic assessment

Relevant scripts:
- 17_extract_candidate_pairs.py
- 18_competitive_alignment.sh
- 19_summarize_alignment_support.py
- 20_corrected_alignment.sh
- 21_genomic_distribution.py
- 22_export_blast_candidates.py

Publication outputs:
- corrected_genus_support.csv
- genomic_distribution_summary.csv

Qualifying genomic-alignment pairs: 10,149.

These alignments represent within-reference-panel support,
not confirmed organism detection.

### 5. Whole-metagenome marker assessment

This analysis examined the complete human-screened read sets,
independently of Kraken2 candidate selection.

Original analysis scripts and documentation from the integrated
marker investigation remain preserved in the repository.

Publication outputs:
publication_results/whole_metagenome_marker/

Three isolated target-associated read pairs were retained:
two Blastocystis-associated and one Entamoeba-associated.

No specimen-genus combination met the predefined
multi-fragment corroboration criterion.

### 6. Candidate-specific marker assessment

Relevant scripts:
- 23_new_marker_validation.sh
- 24_finalize_new_marker_evidence.py
- 25_final_marker_taxonomy.py

Publication outputs:
- final_corrected_marker_evidence.csv
- new_PR2_marker_evidence.csv

This analysis evaluated Kraken2-selected candidate reads
and used different acceptance criteria from the
whole-metagenome marker assessment.

No pairs met the stringent genus-concordant PR2 criterion.

## Interpretation

The investigation identified preliminary genetic signals associated
with enteric protozoa and microsporidia, but the combined evidence
did not establish corroborated organism-level detections.

These results do not demonstrate pathogen absence, organism
viability, community infection prevalence or infections missed
by existing surveillance programmes.

The public-health objective is to evaluate whether enteric protozoa
merit further investigation as additional wastewater surveillance
targets for diarrhoeal diseases in Uganda.

## Reproducibility considerations

The publication_results directory contains selected summary outputs.

Full reproduction requires the original sequencing reads, reference
sequences, final assembly manifest, software dependencies and
associated analytical parameters.

The stage descriptions identify the relevant scripts but are not
a verified end-to-end executable workflow.

Large FASTQ files, reference databases and intermediate alignments
are intentionally excluded from this repository.

## Final genomic reference inventory

The 55 non-human assemblies used in the final genomic screening
are documented in:

metadata/taxonomy/final_screen_reference_inventory.tsv

The inventory combines:
- 28 assemblies from enteric_combined_v1.
- 27 non-overlapping assemblies recovered from legacy_unique_v1.

The 28-assembly collection was downloaded and subjected to
source-checksum verification during reference preparation.

The 27 recovered assemblies were identified and inventoried
from an earlier local reference archive. Their original-source
checksums were not independently revalidated.

The human competitor reference was included separately and
is not counted among the 55 assemblies.

The screening database was constructed under:

references/kraken2/kampala_screen_v3

Reference sequence files and the compiled Kraken2 database
are excluded from GitHub because of their size.
