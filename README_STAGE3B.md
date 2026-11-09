# Kampala enteric protozoa — Stage 3B locus distribution audit

This is a focused audit of the Stage 3A alignment evidence. It is **not** taxonomic confirmation from an independent sequence source. The existing 44-assembly panel remains incomplete, and reads were selected based on Kraken2 results. A later independent marker/broader-decoy check is still required.

## Install and run

Run these commands inside your project. The ZIP contains the script under the `scripts/` folder.

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_stage3b_locus_audit.zip -d .
conda activate kampala-align
python scripts/19_audit_candidate_loci.py --check-only
python scripts/19_audit_candidate_loci.py --run
```

The last command reads existing BAM files in `results/confirmation/stage3a_curated_v1/`, checks all six samples, reproduces Stage 3A high-quality assignment counts, and calculates whether read-pair alignments are concentrated in one or a few 1,000-base reference regions. It does not rerun Kraken2, trim reads, or download new reference genomes. It does not change your existing data.

## What is being measured and why

* `high_quality_distinct_pairs` — distinct fragments that pass the same MAPQ >=20, aligned length >=75 bp, and approximate identity >=95% rules as Stage 3A.
* `total_1kb_bins` — distinct 1,000-base regions touched by qualifying alignments.
* `top_1kb_bin_pair_fraction` — fraction of qualifying fragments touching the single most frequently supported region. High values can indicate localized/conserved/repetitive sequence effects.
* `top_3_bins_union_pair_fraction` — fraction of fragments touching any of the three most-supported bins, counting each fragment once.
* `top_contig_pair_fraction` — fraction of fragments that align to the most-supported reference contig.
* `distinct_reference_start_sites` and `largest_exact_start_count` — whether many fragments start at identical coordinates, which may warrant duplication or artifact review.

These are exploratory diagnostic measurements, **not pass/fail confirmation thresholds**. Low genome breadth can reflect true low concentration, degraded DNA, conserved loci, or artifacts. This analysis alone cannot distinguish all these causes.

## Outputs

```
results/confirmation/stage3b_locus_audit_v1/locus_concentration_metrics.tsv
results/confirmation/stage3b_locus_audit_v1/top_reference_loci.tsv
results/confirmation/stage3b_locus_audit_v1/stage3b_locus_methods.json
metadata/taxonomy/stage3b_locus_review.zip
```

Upload only `metadata/taxonomy/stage3b_locus_review.zip` for aggregate review. Do not upload raw FASTQ or BAM files.

## Troubleshooting

If preflight reports missing `candidate_mapped.sorted.bam`, do not rerun previous stages or delete any data. Share the message so we can locate the existing Stage 3A results. The script refuses to overwrite prior Stage 3B output.
