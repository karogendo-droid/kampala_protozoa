# Genomic reference input manifests

This directory documents the reference assembly collections
used to prepare the final Kampala wastewater genomic screening
database.

## Reference collections

| File | Assemblies | Description |
|---|---:|---|
| refseq_16.tsv | 16 | Original RefSeq download manifest |
| expansion_12.tsv | 12 | Additional GenBank assembly download manifest |
| combined_28.tsv | 28 | Combined newly prepared reference collection |
| legacy_27.tsv | 27 | Additional assemblies recovered from an earlier local reference archive |

The final non-human reference collection contained 55 distinct
assembly accessions: 28 newly prepared and 27 recovered.

The complete combined inventory is available at:

metadata/taxonomy/final_screen_reference_inventory.tsv

The human GRCh38.p14 competitor was included separately and
is not counted among the 55 assemblies.

## Reference provenance

The 28 newly prepared assemblies were subjected to source
checksum verification during reference preparation.

The 27 recovered assemblies were checked against the earlier
local reference inventory, including accession identifiers,
taxonomy assignments and sequence counts.

Their original-source checksums were not independently
revalidated.

## Reproduction requirements

The reference preparation scripts expect their input manifests
in the corresponding directories under:

references/assemblies/

For example, combined_28.tsv corresponds to:

references/assemblies/enteric_combined_v1/manifest.tsv

The download manifests similarly correspond to the
enteric_refseq_v1 and enteric_expansion_v1 directories.

The recovered legacy reference collection depends on an
earlier locally constructed Kraken2 reference archive.
Its sequences cannot be reconstructed from the published
manifest alone.

The manifests contain reference metadata, not genomic
sequence files.

The final Kraken2 screening database additionally included
a human competitor reference, whereas the Bowtie2 genomic
validation index used the 55 non-human assemblies.

The reference panel was targeted rather than a comprehensive
representation of wastewater organisms.

## Scientific interpretation

Inclusion of an assembly in the reference database does not
establish that the corresponding organism was present in
the sampled wastewater.

Taxonomic assignments require additional assessment of
alignment specificity and molecular corroboration.
