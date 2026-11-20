#!/usr/bin/env python3

"""Reconstruct and verify two publication summary tables."""

import csv
from collections import Counter
from functools import lru_cache
from pathlib import Path

RUNS = [f"SRR284341{i}" for i in range(69, 75)]

TARGETS = {
    27971: "Enterocytozoon",
    12967: "Blastocystis",
    5758: "Entamoeba",
    5806: "Cryptosporidium",
    44417: "Cyclospora",
    5740: "Giardia",
    6033: "Encephalitozoon",
}

GENERA = list(TARGETS.values())

KRAKEN = Path("results/kraken2_screen_v3")
DISTRIBUTION = Path(
    "results/protozoa_validation_v3/"
    "competitive_alignment/default_mode/"
    "genomic_distribution_summary.csv"
)

NODES = Path(
    "references/kraken2/kampala_screen_v3/"
    "taxonomy/nodes.dmp"
)

PUBLICATION = Path("publication_results")
OUTPUT = Path("results/publication_table_rebuild")

# Load NCBI taxonomy relationships.
parents = {}

with NODES.open() as f:
    for line in f:
        fields = line.split("|")
        parents[int(fields[0].strip())] = int(
            fields[1].strip()
        )


@lru_cache(maxsize=None)
def target_genus(taxid):
    seen = set()

    while taxid not in seen:
        if taxid in TARGETS:
            return TARGETS[taxid]

        seen.add(taxid)
        parent = parents.get(taxid)

        if parent is None or parent == taxid:
            break

        taxid = parent

    return None


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["Run"] + GENERA,
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def verify(filename):
    rebuilt = read_csv(OUTPUT / filename)
    published = read_csv(PUBLICATION / filename)

    if rebuilt != published:
        raise SystemExit(
            f"MISMATCH: {filename} differs "
            "from the published table."
        )

    print(f"PASS: {filename} matches publication.")


# Stage 1: Reconstruct genus-clade counts from Kraken2.
print("\n1. Checking preliminary taxonomic assignments")

preliminary = []

for run in RUNS:
    counts = Counter()

    source = KRAKEN / f"{run}.kraken"

    with source.open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t", 3)

            if fields[0] != "C":
                continue

            genus = target_genus(int(fields[2]))

            if genus:
                counts[genus] += 1

    row = {"Run": run}

    for genus in GENERA:
        row[genus] = counts[genus]

    preliminary.append(row)

    print(run, ":", sum(counts.values()))


# Stage 2: Summarize final qualifying genomic pairs.
print("\n2. Checking qualifying genomic alignments")

distribution_rows = read_csv(DISTRIBUTION)

support = {}

for row in distribution_rows:
    key = (row["Run"], row["Genus"])

    if key in support:
        raise SystemExit(f"Duplicate result: {key}")

    support[key] = int(row["Supporting_pairs"])

expected_keys = {
    (run, genus)
    for run in RUNS
    for genus in GENERA
}

if set(support) != expected_keys:
    raise SystemExit(
        "Unexpected or missing run-genus combinations."
    )

corrected = []

for run in RUNS:
    row = {"Run": run}

    for genus in GENERA:
        row[genus] = support[(run, genus)]

    corrected.append(row)


# Stage 3: Save rebuilt outputs outside tracked files.
OUTPUT.mkdir(parents=True, exist_ok=True)

write_csv(
    OUTPUT / "genus_counts_all_six.csv",
    preliminary
)

write_csv(
    OUTPUT / "corrected_genus_support.csv",
    corrected
)


# Stage 4: Compare with publication copies.
print("\n3. Verifying publication tables")

verify("genus_counts_all_six.csv")
verify("corrected_genus_support.csv")

preliminary_total = sum(
    row[g] for row in preliminary for g in GENERA
)

corrected_total = sum(
    row[g] for row in corrected for g in GENERA
)

print("\nFINAL SUMMARY")
print("Preliminary candidate pairs:", preliminary_total)
print("Qualifying genomic pairs:", corrected_total)

if preliminary_total != 89444 or corrected_total != 10149:
    raise SystemExit("Totals differ from the manuscript.")

print("\nPASS: Publication summary verification completed.")
print("Rebuilt tables:", OUTPUT)
