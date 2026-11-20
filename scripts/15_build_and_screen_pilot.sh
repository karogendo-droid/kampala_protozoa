#!/bin/bash
set -euo pipefail

DB="references/kraken2/kampala_screen_v3"
OUT="results/kraken2_screen_v3"

mkdir -p "$DB/taxonomy" "$OUT"

# Reuse the verified NCBI taxonomy.
TAX="references/kraken2/kampala_expanded_https_v2/taxonomy"
cp "$TAX/names.dmp" "$TAX/nodes.dmp" "$DB/taxonomy/"

# Use the human reference we have already prepared.
HUMAN="references/assemblies/competing_v1/human_GRCh38_kraken.fna"

if [ ! -s "$HUMAN" ]; then
    gzip -dc "${HUMAN}.gz" > "$HUMAN"
fi

# Add the 28 newly downloaded assemblies.
echo "Adding 28 reference assemblies..."
kraken2-build --add-to-library \
    references/assemblies/enteric_combined_v1/enteric_targets_28.fna \
    --db "$DB"

# Add the 27 recovered legacy assemblies.
echo "Adding 27 additional assemblies..."
kraken2-build --add-to-library \
    references/assemblies/legacy_unique_v1/legacy_unique_27.fna \
    --db "$DB"

# Add human as a competing reference.
echo "Adding human reference..."
kraken2-build --add-to-library "$HUMAN" \
    --no-masking --db "$DB"

# Build the database with a controlled hash-table size.
echo "Building Kraken2 database..."
kraken2-build --build \
    --threads 4 \
    --max-db-size 8000000000 \
    --db "$DB"

echo "PASS: Kraken2 database constructed."

# Locate the previously processed pilot sample.
R1=$(find data/processed/human_screen -type f \
    -name 'SRR28434174_R1.human_screened.fastq.gz' \
    -print -quit)

R2=$(find data/processed/human_screen -type f \
    -name 'SRR28434174_R2.human_screened.fastq.gz' \
    -print -quit)

if [ -z "$R1" ] || [ -z "$R2" ]; then
    echo "ERROR: Pilot FASTQ files not found."
    exit 1
fi

echo "Classifying SRR28434174..."

kraken2 \
    --db "$DB" \
    --threads 4 \
    --paired \
    --gzip-compressed \
    --report-minimizer-data \
    --report "$OUT/SRR28434174.report" \
    --output "$OUT/SRR28434174.kraken" \
    "$R1" "$R2"

echo "PASS: Pilot classification completed."

echo "Candidate protozoan classifications:"
grep -Ei \
    'Giardia|Cryptosporidium|Entamoeba|Cyclospora|Blastocystis|Enterocytozoon|Encephalitozoon' \
    "$OUT/SRR28434174.report" || true
