# Kampala wastewater — Stage 3D: contamination sensitivity audit

## Why this step is required

Stage 3C found that the 12 selected hotspots in reference assemblies labelled *Enterocytozoon bieneusi* H348, *Blastocystis* ST6 and *Cyclospora cayetanensis* have strong BLAST matches to non-target bacteria or human DNA. In the Stage 3A competitive alignment, those suspect reference sequences were labelled as protozoa. Consequently, previously large Kraken2 and Bowtie2 read counts could have reflected bacterial or human reads that matched those sequences.

**Important:** these findings do not prove that all the underlying assemblies or all other protozoa calls are contaminated or false. Stage 3C is a limited region-level screen, not a whole-genome contamination assessment.

This short, **read-only** diagnostic counts qualifying Stage 3A read pairs (MAPQ≥20, aligned length≥75bp, approximate identity≥95%) with alignments to the 12 suspicious reference contigs, then reports how many remain on other contigs of the three selected target assemblies. It does not download references, reprocess FASTQ reads, rebuild databases, or perform treatment plant comparisons. Whole-contig exclusion is precautionary: the foreign match may span only part of a contig. Original outputs are never overwritten.

## Installation

Download `kampala_stage3d_contamination_audit.zip` to your Downloads folder. From your existing project root:

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_stage3d_contamination_audit.zip -d .
conda activate kampala-align
```

All Python code installs to `scripts/22_audit_flagged_reference_contigs.py`.

## First: verify files

The uploaded Stage 3C ZIP should already exist locally at `metadata/taxonomy/stage3c_reference_specificity_review.zip`. The six Stage 3A BAMs and metrics must also be present locally; the script checks this automatically.

```bash
python scripts/22_audit_flagged_reference_contigs.py --check-only
```

## Then run the audit

```bash
python scripts/22_audit_flagged_reference_contigs.py --run
```

The command uses `samtools view` to examine **existing** BAMs only. It reconciles the counts with the original Stage 3A metrics before writing any output. No FASTQ sequences or read identifiers are included in the shareable ZIP.

## Upload this review ZIP

`metadata/taxonomy/stage3d_contamination_audit_review.zip`

The outputs include `affected_species_by_sample.tsv`, `flagged_reference_contigs.tsv`, `stage3d_methods.json` and `stage3d_summary.txt`.

### Scientific interpretation

- **Flagged pairs** = read pairs with ≥1 qualifying alignment anywhere on one of the suspect contigs. These were previously counted as support for their labelled protozoan assembly, but their specificity is not established.
- **Remaining pairs** = pairs supported by other contigs of the same target assembly, excluding any pair also matching a flagged contig. These remain **unconfirmed**, not confirmed protozoan reads.
- Counts in the two columns sum to the previous Stage 3A counts for each species/sample.
- Other protozoan reference assemblies have not yet undergone this contamination screen. This diagnostic is not a replacement for broader reference re-curation, independent markers, or a cleaned-database reclassification.

**Do not compare the two plants biologically using uncorrected RPM values**. We need reference remediation before trusting that those signals represent protozoa rather than abundant environmental bacteria.
