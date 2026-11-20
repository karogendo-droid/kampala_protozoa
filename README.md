# Investigation of Enteric Protozoan Genetic Signals in Kampala Wastewater: Implications for Diarrhoeal Disease Surveillance

## Study overview

Wastewater and environmental surveillance (WES) provides a complementary approach to infectious disease monitoring. However, many surveillance applications have primarily targeted bacterial and viral pathogens, while enteric protozoa remain comparatively underrepresented.

Enteric protozoa, including Cryptosporidium spp., Giardia duodenalis, Entamoeba histolytica and Cyclospora cayetanensis, contribute to diarrhoeal disease and may shed genetic material into sewage through human faecal inputs.

This study investigated enteric protozoan-associated genetic signals in wastewater from Kampala, Uganda, and examined their potential relevance to broader surveillance of diarrhoeal diseases.

## Study objectives

1. Identify and quantify candidate enteric protozoan-associated genetic sequences in wastewater shotgun metagenomes.
2. Describe the distribution of candidate genetic signals across the investigated wastewater specimens.
3. Evaluate the molecular evidence supporting the taxonomic identification of the candidate signals.

## Study data

The investigation used six publicly available paired-end shotgun metagenomes from NCBI BioProject PRJNA1091428.

The specimens originated from wastewater treatment facilities at Bugolobi and Lubigi in Kampala, Uganda.

Each facility contributed:
- Two environmental wastewater specimens.
- One vacuum-truck faecal-sludge specimen.

All six specimens represented dry-condition sampling.

The original sequencing data were generated and deposited by Kanyerezi and colleagues.

BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1091428/

Run accessions: SRR28434169–SRR28434174.

## Overview of the analysis

Sequencing reads underwent quality preprocessing and human-sequence screening.

Candidate enteric protozoan and microsporidian sequences were investigated using taxonomic classification, competitive genomic alignment, genomic-distribution assessment, nucleotide similarity searches and complementary marker-based analyses.

The genomic screening reference collection contained 55 distinct microbial and eukaryotic assemblies, with an additional human competitor reference.

The investigation included Giardia, Cryptosporidium, Entamoeba, Cyclospora, Blastocystis, Enterocytozoon and Encephalitozoon.

Enterocytozoon and Encephalitozoon were considered separately as microsporidia rather than conventional protozoa.

## Principal findings

A total of 91,612,337 original paired-end read pairs were available, of which 68,288,156 remained following preprocessing and human-sequence screening.

Initial taxonomic classification assigned 89,444 read pairs to the seven investigated genera.

Of these, 10,149 pairs (11.3%) satisfied the predefined genus-concordant genomic-alignment criteria.

The genomic distribution and nucleotide similarity assessments identified uncertainty concerning the specificity of some preliminary assignments.

Whole-metagenome marker screening retained three isolated target-associated read pairs, comprising two Blastocystis-associated pairs and one Entamoeba-associated pair.

None of the investigated specimen–genus combinations met the predefined marker corroboration criteria.

The complementary marker assessment of Kraken2-selected candidate reads also identified no pairs satisfying its stringent genus-concordance requirements.

These results represent preliminary genetic signals rather than confirmed organism-level detections.

## Public-health relevance

Enteric protozoa contribute to gastrointestinal disease but remain comparatively underrepresented in wastewater surveillance.

The study addresses whether these organisms warrant consideration as additional targets within wastewater-based surveillance of diarrhoeal diseases in Uganda.

The findings do not demonstrate confirmed protozoan occurrence, organism viability, community infection prevalence or infections missed by existing surveillance programmes.

Prospective investigations using validated organism-specific molecular assays, repeated wastewater sampling and complementary clinical surveillance are needed to determine whether protozoan monitoring could provide additional information for public-health decision-making.

## Repository organization

- `scripts/`: Analysis scripts, including preprocessing, taxonomic screening and molecular assessment.
- `metadata/`: Specimen information and selected reference documentation.
- `publication_results/`: Compact summary tables supporting the finalized manuscript.
- `docs/`: Analysis documentation and methodological notes.
- `envs/`: Software environment information.

The large sequencing files, reference databases and intermediate analysis outputs are not distributed through this repository.

## Publication results

The `publication_results/` directory contains the selected summary tables from the finalized analysis.

These tables should be interpreted as sequence-level evidence and not as estimates of parasite abundance, organism viability or community infection prevalence.

## Data availability

The original sequencing data are available through NCBI BioProject PRJNA1091428.

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1091428/

## Acknowledgement

We acknowledge Kanyerezi S., Nabisubi P., Kebirungi G. and their co-authors for generating and publicly depositing the original wastewater shotgun metagenomic data.

Source publication:

Kanyerezi S, et al. Metagenomics insights into the microbial resistome and virulome composition of Kampala's wastewater. Open Research Africa. 2024;7:8.

https://doi.org/10.12688/openresafrica.15040.1

## Supporting documentation

- [Final analysis and reproducibility](docs/FINAL_ANALYSIS.md)
- [Final reference inventory](metadata/taxonomy/final_screen_reference_inventory.tsv)
- [Reference input manifests](metadata/taxonomy/reference_inputs/README.md)
- [Specimen and sequencing-run manifest](metadata/analysis_sample_manifest.csv)
- [Software environments](envs/)
- [Publication-supporting results](publication_results/)

## Reproducibility

The repository includes analysis scripts and selected summary outputs. Full computational reproduction additionally requires the archived sequencing reads, corresponding reference sequences, software environments and reference manifests.

Analysis-specific documentation should be consulted before executing individual scripts.
