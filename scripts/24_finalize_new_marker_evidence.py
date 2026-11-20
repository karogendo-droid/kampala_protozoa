import csv
import gzip
import re
import subprocess
from pathlib import Path
from collections import defaultdict, Counter
from functools import lru_cache

BASE = Path("results/protozoa_validation_v3/new_marker_validation")
KRAKEN = Path("results/kraken2_screen_v3")
PR2 = Path("references/markers/pr2_v5.1.1")
NODES = Path("references/kraken2/kampala_screen_v3/taxonomy/nodes.dmp")

TARGETS = {
    27971: "Enterocytozoon",
    12967: "Blastocystis",
    5758: "Entamoeba",
    5806: "Cryptosporidium",
    44417: "Cyclospora",
    5740: "Giardia",
    6033: "Encephalitozoon"
}

# Locate the existing PR2 taxonomy.
tax_files = sorted(PR2.rglob("pr2_version_5.1.1_SSU_mothur.tax*"))

if not tax_files:
    print("PR2 taxonomy file not found.")
    print("Available taxonomy-like files:")
    for p in sorted(PR2.rglob("*tax*"))[:20]:
        print(p)
    raise SystemExit("Cannot assign marker genera without PR2 taxonomy.")

tax_file = tax_files[0]
opener = gzip.open if tax_file.suffix == ".gz" else open

print("Using PR2 taxonomy:", tax_file, flush=True)

pr2_genus = {}
with opener(tax_file, "rt") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t", 1)
        if len(parts) != 2:
            continue

        identifier, annotation = parts
        matches = [
            g for g in TARGETS.values()
            if re.search(
                r"(?<![A-Za-z])" + re.escape(g) +
                r"(?![A-Za-z])",
                annotation,
                re.IGNORECASE
            )
        ]

        pr2_genus[identifier] = (
            matches[0] if len(matches) == 1 else None
        )

parents = {}
with NODES.open() as f:
    for line in f:
        p = line.split("|")
        parents[int(p[0].strip())] = int(p[1].strip())

@lru_cache(None)
def kraken_genus(taxid):
    visited = set()
    while taxid not in visited:
        if taxid in TARGETS:
            return TARGETS[taxid]
        visited.add(taxid)
        parent = parents.get(taxid)
        if parent is None or parent == taxid:
            break
        taxid = parent
    return None

def clean_id(value):
    value = value.split()[0]
    if value.endswith(("/1", "/2")):
        return value[:-2]
    return value

rows = []
genera = list(TARGETS.values())

for n in range(69, 75):
    run = f"SRR284341{n}"
    print("Processing", run, flush=True)

    bam = BASE / f"{run}.bam"

    output = subprocess.run(
        ["samtools", "view", "-f", "2", "-F", "2308", str(bam)],
        capture_output=True,
        text=True,
        check=True
    )

    pairs = defaultdict(dict)
    other_records = 0

    for line in output.stdout.splitlines():
        p = line.split("\t")
        rid = clean_id(p[0])
        flag = int(p[1])
        ref = p[2]
        pos = int(p[3])
        mapq = int(p[4])

        if ref not in pr2_genus:
            other_records += 1
            continue

        mate = 1 if flag & 64 else 2 if flag & 128 else 0

        if mate:
            pairs[rid][mate] = (
                pr2_genus[ref], ref, pos, mapq
            )

    # Require both mates to match the same target genus
    # using PR2 markers with numerical MAPQ >= 10.
    eligible = {}

    for rid, mates in pairs.items():
        if set(mates) != {1, 2}:
            continue

        a, b = mates[1], mates[2]

        if a[0] is None or a[0] != b[0]:
            continue

        if not (10 <= a[3] < 255 and 10 <= b[3] < 255):
            continue

        eligible[rid] = (a[0], a, b)

    # Confirm agreement with the new Kraken2 assignment.
    accepted = {}
    with (KRAKEN / f"{run}.kraken").open() as f:
        for line in f:
            p = line.rstrip("\n").split("\t", 3)
            rid = clean_id(p[1])

            if rid not in eligible or p[0] != "C":
                continue

            genus, a, b = eligible[rid]

            if kraken_genus(int(p[2])) == genus:
                accepted[rid] = (genus, a, b)

    for g in genera:
        selected = [
            value for value in accepted.values()
            if value[0] == g
        ]

        references = set()
        starts = set()
        bins = set()

        for _, a, b in selected:
            for x in (a, b):
                references.add(x[1])
                starts.add((x[1], x[2]))
                bins.add((x[1], (x[2]-1)//100))

        rows.append({
            "Run": run,
            "Genus": g,
            "Concordant_marker_pairs": len(selected),
            "Distinct_PR2_references": len(references),
            "Distinct_start_positions": len(starts),
            "Distinct_100bp_bins": len(bins)
        })

    print(
        "  PR2 genus-concordant pairs:",
        len(accepted),
        "| Non-PR2 alignment records:",
        other_records
    )

output_file = BASE / "new_PR2_marker_evidence.csv"

with output_file.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print("\n=== FINAL NEW PR2 MARKER EVIDENCE ===")
print("Genus | Supporting pairs | References | 100-bp bins")

for genus in genera:
    subset = [r for r in rows if r["Genus"] == genus]
    print(
        genus,
        "|", sum(r["Concordant_marker_pairs"] for r in subset),
        "|", sum(r["Distinct_PR2_references"] for r in subset),
        "|", sum(r["Distinct_100bp_bins"] for r in subset)
    )

print("\nPASS: New PR2 marker evidence summarized.")
print("Saved:", output_file)
