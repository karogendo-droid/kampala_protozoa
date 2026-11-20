import csv
import re
import subprocess
from pathlib import Path
from collections import Counter, defaultdict
from functools import lru_cache

base = Path("results/protozoa_validation_v3/competitive_alignment/default_mode")
kraken = Path("results/kraken2_screen_v3")
taxfile = Path("references/kraken2/kampala_screen_v3/taxonomy/nodes.dmp")

targets = {
    27971: "Enterocytozoon",
    12967: "Blastocystis",
    5758: "Entamoeba",
    5806: "Cryptosporidium",
    44417: "Cyclospora",
    5740: "Giardia",
    6033: "Encephalitozoon"
}

parents = {}
with taxfile.open() as f:
    for line in f:
        p = line.split("|")
        parents[int(p[0].strip())] = int(p[1].strip())

@lru_cache(None)
def genus(taxid):
    seen = set()
    while taxid not in seen:
        if taxid in targets:
            return targets[taxid]
        seen.add(taxid)
        parent = parents.get(taxid)
        if parent is None or parent == taxid:
            break
        taxid = parent
    return None

def clean_id(x):
    x = x.split()[0]
    return x[:-2] if x.endswith(("/1", "/2")) else x

rows = []
all_windows = defaultdict(Counter)
all_contigs = defaultdict(set)

for n in range(69, 75):
    run = f"SRR284341{n}"
    print("Checking", run, flush=True)

    assigned = {}
    with (kraken / f"{run}.kraken").open() as f:
        for line in f:
            p = line.rstrip("\n").split("\t", 3)
            if p[0] == "C":
                g = genus(int(p[2]))
                if g:
                    assigned[clean_id(p[1])] = g

    pairs = defaultdict(dict)

    process = subprocess.Popen(
        ["samtools", "view", "-f", "2", "-F", "2308",
         str(base / f"{run}.bam")],
        stdout=subprocess.PIPE,
        text=True
    )

    for line in process.stdout:
        p = line.split("\t", 6)
        rid = clean_id(p[0])
        if rid not in assigned:
            continue

        mapq = int(p[4])
        if not 10 <= mapq < 255:
            continue

        match = re.search(r"kraken:taxid\|(\d+)", p[2])
        if not match:
            raise RuntimeError("Missing reference taxid")

        if genus(int(match.group(1))) != assigned[rid]:
            continue

        flag = int(p[1])
        mate = 1 if flag & 64 else 2 if flag & 128 else 0

        if mate:
            if mate in pairs[rid]:
                raise RuntimeError("Duplicate primary mate: " + rid)
            pairs[rid][mate] = (p[2], int(p[3]))

    if process.wait() != 0:
        raise RuntimeError("samtools failed: " + run)

    by_genus = defaultdict(Counter)

    for rid, mates in pairs.items():
        if 1 not in mates or 2 not in mates:
            continue

        ref1, pos1 = mates[1]
        ref2, pos2 = mates[2]

        if ref1 != ref2:
            continue

        g = assigned[rid]

        # A 1-kb genomic interval, uniquely identified
        # by reference contig and interval number.
        interval = (ref1, (min(pos1, pos2) - 1) // 1000)

        by_genus[g][interval] += 1
        all_windows[g][interval] += 1
        all_contigs[g].add(ref1)

    for g in targets.values():
        counter = by_genus[g]
        count = sum(counter.values())

        rows.append({
            "Run": run,
            "Genus": g,
            "Supporting_pairs": count,
            "Distinct_1kb_intervals": len(counter),
            "Distinct_reference_contigs": len(
                {key[0] for key in counter}
            ),
            "Top_interval_percent": round(
                100 * max(counter.values()) / count, 2
            ) if count else 0
        })

out = base / "genomic_distribution_summary.csv"

with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print("\n=== GENOMIC DISTRIBUTION SUMMARY ===")
print("Genus | Supporting pairs | Contigs | 1-kb intervals | Top 10 intervals (%)")

for g in targets.values():
    counter = all_windows[g]
    total = sum(counter.values())

    if not total:
        print(g, "| 0 | 0 | 0 | 0")
        continue

    top10 = sum(sorted(counter.values(), reverse=True)[:10])

    print(
        g, "|", total,
        "|", len(all_contigs[g]),
        "|", len(counter),
        "|", round(100 * top10 / total, 2)
    )

print("\nSaved:", out)
