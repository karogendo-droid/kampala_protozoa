#!/usr/bin/env bash
# One-command first-pass, curated protozoa-inclusive Kraken2 reference database build.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
DB="$ROOT/references/kraken2/kampala_curated_v1"
LOG="$ROOT/metadata/taxonomy/curated_v1_kraken_build.log"
mkdir -p "$DB" "$(dirname "$LOG")"
for cmd in python kraken2 kraken2-build dustmasker; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "STOP: missing executable '$cmd' in the current Conda environment." >&2
    if [[ "$cmd" == "dustmasker" ]]; then
      echo "Install NCBI BLAST+: conda install -n kampala-taxanomy -c conda-forge -c bioconda blast -y" >&2
    fi
    exit 2
  fi
done
python scripts/14_prepare_curated_db.py --check-only
FREE_KIB="$(df -Pk "$ROOT" | awk 'NR==2{print $4}')"
if [[ "$FREE_KIB" -lt 36700160 ]]; then
  echo "STOP: at least 35 GiB free space required for staged files and build." >&2
  exit 2
fi
if [[ -f "$DB/BUILD_COMPLETE" ]]; then
  echo "Existing completed database found: $DB"
  python scripts/15_smoke_test_curated_db.py --db "$DB"
  exit 0
fi
python scripts/14_prepare_curated_db.py 2>&1 | tee -a "$LOG"
# The marker protects resume after successful addition. A partial prior add is not silently reused.
if [[ ! -f "$DB/LIBRARY_ADDED" ]]; then
  echo 'Adding reference FASTA to Kraken2 (dustmasker will handle low-complexity regions)...'
  kraken2-build --add-to-library "$DB/staged_reference_library.fna" --db "$DB" 2>&1 | tee -a "$LOG"
  touch "$DB/LIBRARY_ADDED"
fi
if [[ ! -f "$DB/BUILD_COMPLETE" ]]; then
  echo 'Building Kraken2 hash table with an 8 GB hash-table cap...'
  kraken2-build --build --threads 4 --max-db-size 8000000000 --db "$DB" 2>&1 | tee -a "$LOG"
  for file in hash.k2d opts.k2d taxo.k2d; do
    if [[ ! -s "$DB/$file" ]]; then
      echo "STOP: Kraken2 database missing required file $file" >&2
      exit 2
    fi
  done
  echo 'Database files verified.'
  touch "$DB/BUILD_COMPLETE"
fi
python scripts/15_smoke_test_curated_db.py --db "$DB" 2>&1 | tee -a "$LOG"
echo "Database ready for separate taxonomic classification stage: $DB"
