# Additional quality-assurance and provenance analyses

This document describes **supplementary analytical checks performed during development** of the Kampala wastewater enteric-protozoa workflow. These checks are not independent public-health study objectives, but are retained because they affect the credibility and interpretation of the exploratory taxonomic assignments.

## Why these analyses are retained

The initial Kraken2 database was purpose-curated and represented only a limited number of competing organisms. Therefore, assigning reads to a protozoan reference does not by itself establish that the sequence originated from the intended organism. Additional comparisons were used to assess the specificity and genomic distribution of some candidate assignments. They must remain available to anyone seeking to reproduce the full analytical history and interpretation.

## Historical supplementary scripts

| Script | Purpose | Original documentation |
|---|---|---|
| `scripts/18_validate_candidate_alignments.py` | Competitive alignment of candidate reads against the curated reference collection | [Stage 3A](../README_STAGE3A.md) |
| `scripts/19_audit_candidate_loci.py` | Inspect genomic breadth and concentration of qualifying alignments | [Stage 3B](../README_STAGE3B.md) |
| `scripts/20_screen_reference_hotspots.py` | Check selected reference sequence hotspots using remote sequence comparison | [Stage 3C](../README_STAGE3C.md) |
| `scripts/21_screen_reference_hotspots_https.py` | HTTPS implementation of remote hotspot checking | [Stage 3C HTTPS](../README_STAGE3C_HTTPS.md) |
| `scripts/22_audit_flagged_reference_contigs.py` | Quantify assignments involving reference contigs flagged by specificity checks | [Stage 3D](../README_STAGE3D.md) |

## Important observations

The first competitive-alignment stage produced qualifying read pairs for selected assemblies, including *Enterocytozoon bieneusi* (8,149), *Blastocystis* ST6 (970), and *Cyclospora cayetanensis* (736). These figures were generated using the original curated reference collection and are **not independently confirmed organism counts**.

Further examination found substantial hotspot concentration and reference-specificity concerns. Of the qualifying pairs above, 5,222/8,149 (64.1%) for *E. bieneusi*, 556/970 (57.3%) for *Blastocystis* ST6, and 473/736 (64.3%) for *Cyclospora* occurred on contigs flagged by the reference review. The fact that some pairs remained after excluding these contigs does **not** establish their specificity or validate organism detection.

Separate marker-based analyses using PR2 SSU/18S sequences and *E. bieneusi* ITS references recovered isolated candidate fragments but did not meet the predefined multi-fragment/multi-locus support thresholds for any target genus. These results limit interpretation of the initial classification signals, while not demonstrating that the biological organisms were absent from wastewater.

## Reporting and reproducibility

The seven-step [primary workflow](../README.md#primary-analysis-workflow) emphasizes the public-health question and its planned taxonomic and marker analyses. The supplementary checks above form part of the **actual analysis provenance**, and their material implications should be disclosed when discussing detection confidence. Keep the historical scripts and result logs accessible where they can be shared safely, without redistributing raw reads, human-derived sequences, credentials, or private information.
