# Stage 3E: Independent PR2 small-subunit rRNA marker screen

This step provides an independent **reference-source** cross-check for candidate protists after the contamination audit found suspect genome contigs. It screens all six **host-screened paired FASTQs**, rather than recycling reads selected by the first Kraken2 database. It **does not** imply organism or subtype confirmation.

## Install into the existing Kampala project

The ZIP contains `scripts/23_screen_pr2_18s.py`. Extract it at `~/bioinformatics/kampala_protozoa`, preserving the `scripts/` subfolder. Nothing previously produced is overwritten.

```bash
cd ~/bioinformatics/kampala_protozoa
unzip -n ~/Downloads/kampala_stage3e_pr2_marker_screen.zip -d .
conda activate kampala-align
python scripts/23_screen_pr2_18s.py --check-only
python scripts/23_screen_pr2_18s.py --run --threads 4
```

If `--check-only` reports missing `curl`, `bowtie2`, `bowtie2-build`, or `samtools`, install the missing executable in the alignment Conda environment before continuing. If the script reports multiple candidate FASTQs for the same run/mate, resolve that explicitly rather than guessing which files are original.

The script downloads the published PR2 v5.1.1 SSU Mothur FASTA and taxonomy tables, builds a separate Bowtie2 index in `references/markers/pr2_v5.1.1/`, and maps all six samples. It retains BAM/logs in `results/confirmation/stage3e_pr2_18s_v1/` and creates a small report ZIP at `metadata/taxonomy/stage3e_pr2_marker_review.zip`.

Unlike the previous genome database, PR2 is an independently curated **SSU (mainly 18S rRNA) marker** database. We require MAPQ >=20, aligned query length >=75bp, and approximate identity >=97%; those are conservative alignment filters, not definitive taxonomic validation. Genus-level read-pair counts are exploratory. The species/subtype cannot be confirmed from this analysis alone, and no matches do not establish absence because marker coverage and taxonomic representation can be incomplete. *Enterocytozoon bieneusi* ITS, *Blastocystis* 18S barcode subtype, and *Cyclospora* diagnostic targets may need specific follow-up.

Reference: https://github.com/pr2database/pr2database/releases/tag/v5.1.1 and https://pr2-database.org/

The downloader uses GitHub release assets; if NCBI or GitHub connectivity is interrupted, re-running is safe. Existing downloaded files, indexes, and completed BAMs are reused. The program does not upload wastewater reads to an external server.

## Share results

Upload only: `metadata/taxonomy/stage3e_pr2_marker_review.zip`.

IMPORTANT: We are **not** claiming that any organism is present or absent based on previous Kraken2 reads, remaining contigs, or this marker screening alone. Final plant comparisons will use only appropriately verified evidence and retain limitations for small n=3/plant and different sampling conditions.
