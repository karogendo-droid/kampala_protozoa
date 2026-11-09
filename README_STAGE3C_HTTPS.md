# Kampala wastewater protozoa — Stage 3C HTTPS alternative

## Why

The original `scripts/20_screen_reference_hotspots.py --run` tried BLAST+'s remote dispatcher and failed. Your tests showed `https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi` is reachable. This helper uses the documented NCBI BLAST Common HTTPS URL API instead.

**No wastewater reads are transmitted.** Only the twelve already-prepared 801–1000 bp public genome reference windows are submitted. The script reuses script 20's selected locus manifest, interpretation and review ZIP packaging.

## Installation

From the project directory:

```bash
cd ~/bioinformatics/kampala_protozoa
conda activate kampala-taxanomy
unzip -n ~/Downloads/kampala_stage3c_https_alternative.zip -d .
```

All executable Python scripts remain in `scripts/`.

## One quick check and one run

```bash
python scripts/21_screen_reference_hotspots_https.py --check-only
python scripts/21_screen_reference_hotspots_https.py --run
```

Optionally add `--email your_email@example.org` to identify a contact for NCBI API usage. This script sends no private sequencing material.

The job is submitted once and assigned an RID. Subsequent runs **resume** the RID instead of starting another job. Status requests are at least 65 seconds apart as required by NCBI; an interrupted run can be resumed using the identical command. `--restart` submits a *new* job and should be reserved for invalid or expired RIDs.

Output review ZIP:

`metadata/taxonomy/stage3c_reference_specificity_review.zip`

## Files written

- `results/confirmation/stage3c_reference_specificity_v1/ncbi_blast_https_job.json`: saved RID and inputs for resumability.
- `.../reference_hotspots_core_nt_https.json`: original NCBI machine-readable response.
- `.../reference_hotspots_core_nt_blast.tsv`: converted taxon-annotated BLAST matches.
- `.../reference_hotspot_blast_summary.tsv`: one summary per reference window.
- `.../top20_reference_blast_hits.tsv`: leading matches for each window.
- `metadata/taxonomy/stage3c_reference_specificity_review.zip`: share this file in ChatGPT.

### Important limitation

A genomic reference window finding matches in NCBI core_nt is **a reference specificity screen**, not a separate comparison of wastewater read sequences. We still cannot confirm protozoan species from this step alone. Taxonomic breadth of NCBI core_nt is incomplete, and multi-taxon descriptions are not independent observations.

Official NCBI documentation: https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html and https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html
