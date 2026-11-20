#!/bin/bash
set -euo pipefail

REF="references/markers/integrated_stage4_v1/integrated_marker_competitive.fasta"
IDX="references/markers/marker_validation_v3/marker_index"
IN="results/protozoa_validation_v3"
OUT="results/protozoa_validation_v3/new_marker_validation"

mkdir -p "$(dirname "$IDX")" "$OUT"

if [[ ! -s "$REF" ]]; then
    echo "ERROR: Marker reference not found."
    exit 1
fi

# Build the marker index only if it does not exist.
if [[ ! -s "${IDX}.1.bt2" && ! -s "${IDX}.1.bt2l" ]]; then
    echo "Preparing marker alignment index..."

    bowtie2-build --threads 4 \
        "$REF" "$IDX" \
        > "$OUT/index.log" 2>&1

    echo "Marker index ready."
else
    echo "Reusing marker index."
fi

# Align the candidate read pairs from all six samples.
for N in 69 70 71 72 73 74
do
    RUN="SRR284341${N}"
    BAM="$OUT/${RUN}.bam"

    if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
        echo "$RUN already completed."
        continue
    fi

    echo "Marker validation: $RUN"

    bowtie2 \
        --very-sensitive \
        --no-unal \
        -p 4 \
        -x "$IDX" \
        -1 "$IN/${RUN}_candidate_R1.fastq.gz" \
        -2 "$IN/${RUN}_candidate_R2.fastq.gz" \
        2> "$OUT/${RUN}.log" \
        | samtools sort -@ 2 -m 512M \
            -o "$OUT/${RUN}.tmp.bam" -

    samtools quickcheck "$OUT/${RUN}.tmp.bam"

    mv "$OUT/${RUN}.tmp.bam" "$BAM"
    samtools index "$BAM"

    echo "PASS: $RUN"
done

echo "PASS: All six new marker alignments completed."
