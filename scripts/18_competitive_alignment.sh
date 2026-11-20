#!/bin/bash
set -euo pipefail

REF1="references/assemblies/enteric_combined_v1/enteric_targets_28.fna"
REF2="references/assemblies/legacy_unique_v1/legacy_unique_27.fna"

INDEX_DIR="references/validation_bowtie2/expanded_v3"
INDEX="$INDEX_DIR/competitive_55"

INPUT="results/protozoa_validation_v3"
OUTPUT="$INPUT/competitive_alignment"

mkdir -p "$INDEX_DIR" "$OUTPUT"

for TOOL in bowtie2-build bowtie2 samtools; do
    if ! command -v "$TOOL" >/dev/null 2>&1; then
        echo "ERROR: $TOOL is not installed in this environment."
        exit 1
    fi
done

# Build the alignment index once.
if [[ ! -s "${INDEX}.1.bt2" && ! -s "${INDEX}.1.bt2l" ]]; then
    echo "Building competitive Bowtie2 index..."

    bowtie2-build \
        --threads 4 \
        --packed \
        "$REF1,$REF2" \
        "$INDEX" \
        > "$OUTPUT/index_build.log" 2>&1

    echo "PASS: Alignment index constructed."
else
    echo "Existing alignment index found."
fi

# Align candidate reads from all six samples.
for NUMBER in 69 70 71 72 73 74
do
    RUN="SRR284341${NUMBER}"

    echo ""
    echo "========== Aligning $RUN =========="

    R1="$INPUT/${RUN}_candidate_R1.fastq.gz"
    R2="$INPUT/${RUN}_candidate_R2.fastq.gz"

    BAM="$OUTPUT/${RUN}.bam"

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "ERROR: Candidate FASTQ missing for $RUN"
        exit 1
    fi

    if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
        echo "Previously completed alignment found. Skipping."
        continue
    fi

    bowtie2 \
        --very-sensitive \
        -k 5 \
        -p 4 \
        -x "$INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        2> "$OUTPUT/${RUN}.bowtie2.log" \
        | samtools sort \
            -@ 2 \
            -m 512M \
            -o "$OUTPUT/${RUN}.tmp.bam" -

    mv "$OUTPUT/${RUN}.tmp.bam" "$BAM"

    samtools index "$BAM"

    samtools flagstat "$BAM" \
        > "$OUTPUT/${RUN}.flagstat.txt"

    echo "PASS: $RUN alignment completed."
done

echo ""
echo "SUCCESS: Competitive alignment completed for all six samples."
