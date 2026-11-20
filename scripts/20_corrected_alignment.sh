#!/bin/bash
set -euo pipefail

INDEX="references/validation_bowtie2/expanded_v3/competitive_55"
INPUT="results/protozoa_validation_v3"
OUT="$INPUT/competitive_alignment/default_mode"

mkdir -p "$OUT"

for N in 69 70 71 72 73 74; do
    RUN="SRR284341${N}"

    if [[ -s "$OUT/${RUN}.bam" && -s "$OUT/${RUN}.bam.bai" ]]; then
        echo "$RUN already completed; skipping."
        continue
    fi

    echo "Aligning $RUN using default Bowtie2 reporting..."

    bowtie2 \
        --very-sensitive \
        -p 4 \
        -x "$INDEX" \
        -1 "$INPUT/${RUN}_candidate_R1.fastq.gz" \
        -2 "$INPUT/${RUN}_candidate_R2.fastq.gz" \
        2> "$OUT/${RUN}.log" \
        | samtools sort -@ 2 -m 512M \
            -o "$OUT/${RUN}.tmp.bam" -

    mv "$OUT/${RUN}.tmp.bam" "$OUT/${RUN}.bam"

    samtools quickcheck "$OUT/${RUN}.bam"
    samtools index "$OUT/${RUN}.bam"

    echo "PASS: $RUN"
done

echo "All six corrected alignments completed."
