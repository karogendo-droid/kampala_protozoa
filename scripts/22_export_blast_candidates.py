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
    12967: "Blastocystis",
    27971: "Enterocytozoon",
    44417: "Cyclospora"
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

counts = defaultdict(Counter)
examples = {}

for n in range(69, 75):
    run = f"SRR284341{n}"
    assigned = {}

    with (kraken / f"{run}.kraken").open() as f:
        for line in f:
            p = line.split("\t", 3)
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
        p = line.rstrip("\n").split("\t")
        rid = clean_id(p[0])
        if rid not in assigned:
            continue

        flag = int(p[1])
        mapq = int(p[4])
        if not 10 <= mapq < 255:
            continue

        match = re.search(r"kraken:taxid\|(\d+)", p[2])
        if not match:
            continue

        if genus(int(match.group(1))) != assigned[rid]:
            continue

        mate = 1 if flag & 64 else 2 if flag & 128 else 0
        if mate:
            pairs[rid][mate] = (
                p[2], int(p[3]), p[9], mapq
            )

    if process.wait() != 0:
        raise RuntimeError(f"samtools failed: {run}")

    for rid, mates in pairs.items():
        if set(mates) != {1, 2}:
            continue

        a, b = mates[1], mates[2]
        if a[0] != b[0]:
            continue

        g = assigned[rid]
        window = (a[0], (min(a[1], b[1]) - 1) // 1000)
        key = (g, window)

        counts[g][window] += 1
        score = a[3] + b[3]

        if key not in examples or score > examples[key][0]:
            examples[key] = (score, run, rid, a[2], b[2])

output = base / "blast_candidate_hotspots"
output.mkdir(exist_ok=True)

with (output / "hotspots.tsv").open("w", newline="") as tab:
    writer = csv.writer(tab, delimiter="\t")
    writer.writerow([
        "Genus", "Rank", "Reference",
        "Interval_start", "Supporting_pairs", "Example_run"
    ])

    for g in targets.values():
        fasta = output / f"{g}_representative_reads.fasta"

        with fasta.open("w") as out:
            for rank, (window, number) in enumerate(
                counts[g].most_common(3), 1
            ):
                ref, index = window
                score, run, rid, seq1, seq2 = examples[(g, window)]

                writer.writerow([
                    g, rank, ref,
                    index * 1000 + 1, number, run
                ])

                out.write(
                    f">{g}_hotspot{rank}_R1 {run}\n{seq1}\n"
                    f">{g}_hotspot{rank}_R2 {run}\n{seq2}\n"
                )

        print(f"{g}: {fasta}")

print("\nPASS: Representative BLAST sequences prepared.")
print("Hotspot details:", output / "hotspots.tsv")
