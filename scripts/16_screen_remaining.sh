#!/bin/bash
set -euo pipefail

DB="references/kraken2/kampala_screen_v3"
OUT="results/kraken2_screen_v3"

mkdir -p "$OUT"

for RUN in SRR28434169 SRR28434170 SRR28434171 SRR28434172 SRR28434173
do
    echo ""
    echo "========== Processing $RUN =========="

    R1=$(find data/processed/human_screen -type f \
        -name "${RUN}_R1.human_screened.fastq.gz" -print -quit)

    R2=$(find data/processed/human_screen -type f \
        -name "${RUN}_R2.human_screened.fastq.gz" -print -quit)

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "ERROR: FASTQ files missing for $RUN"
        exit 1
    fi

    if [[ -s "$OUT/${RUN}.report" && -s "$OUT/${RUN}.kraken" ]]; then
        echo "Existing complete output found. Skipping."
        continue
    fi

    kraken2 \
        --db "$DB" \
        --threads 4 \
        --paired \
        --gzip-compressed \
        --report-minimizer-data \
        --report "$OUT/${RUN}.report.tmp" \
        --output "$OUT/${RUN}.kraken.tmp" \
        "$R1" "$R2"

    mv "$OUT/${RUN}.report.tmp" "$OUT/${RUN}.report"
    mv "$OUT/${RUN}.kraken.tmp" "$OUT/${RUN}.kraken"

    echo "PASS: $RUN completed."
done

echo ""
echo "ALL FIVE ADDITIONAL SAMPLES COMPLETED."
