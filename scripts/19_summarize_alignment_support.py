import csv
import re
import subprocess
from pathlib import Path
from functools import lru_cache

ROOT = Path("results/kraken2_screen_v3")
ALIGN = Path("results/protozoa_validation_v3/competitive_alignment")
TAX = Path("references/kraken2/kampala_screen_v3/taxonomy/nodes.dmp")

TARGETS = {
    27971: "Enterocytozoon",
    12967: "Blastocystis",
    5758: "Entamoeba",
    5806: "Cryptosporidium",
    44417: "Cyclospora",
    5740: "Giardia",
    6033: "Encephalitozoon"
}

genera = list(TARGETS.values())
parents = {}

with TAX.open() as f:
    for line in f:
        fields = line.split("|")
        parents[int(fields[0].strip())] = int(fields[1].strip())

@lru_cache(maxsize=None)
def genus_of(taxid):
    seen = set()

    while taxid not in seen:
        if taxid in TARGETS:
            return TARGETS[taxid]

        seen.add(taxid)
        parent = parents.get(taxid)

        if parent is None or parent == taxid:
            break

        taxid = parent

    return "Other"

def clean_id(value):
    value = value.split()[0]
    if value.endswith(("/1", "/2")):
        value = value[:-2]
    return value

rows = []

for number in range(69, 75):
    run = f"SRR284341{number}"
    print(f"Analyzing {run}...", flush=True)

    kraken_file = ROOT / f"{run}.kraken"
    bam = ALIGN / f"{run}.bam"

    assigned = {}

    with kraken_file.open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t", 3)

            if fields[0] != "C":
                continue

            genus = genus_of(int(fields[2]))

            if genus != "Other":
                assigned[clean_id(fields[1])] = genus

    alignments = {}

    process = subprocess.Popen(
        ["samtools", "view", "-F", "4", str(bam)],
        stdout=subprocess.PIPE,
        text=True
    )

    for line in process.stdout:
        fields = line.split("\t", 6)

        read_id = clean_id(fields[0])

        if read_id not in assigned:
            continue

        flag = int(fields[1])
        reference = fields[2]
        mapq = int(fields[4])

        match = re.search(r"kraken:taxid\|(\d+)", reference)

        if not match:
            raise RuntimeError(
                f"Reference has no taxonomy ID: {reference}"
            )

        mapped_genus = genus_of(int(match.group(1)))

        state = alignments.setdefault(
            read_id,
            {"primary": {}, "genera": set()}
        )

        state["genera"].add(mapped_genus)

        # Ignore secondary and supplementary alignments
        # when counting primary mapped read pairs.
        if flag & 0x100 or flag & 0x800:
            continue

        mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0

        if mate:
            if mate in state["primary"]:
                raise RuntimeError(
                    f"Duplicate primary alignment: {read_id}"
                )

            state["primary"][mate] = (
                mapped_genus,
                mapq,
                bool(flag & 0x2)
            )

    if process.wait() != 0:
        raise RuntimeError(f"samtools failed for {run}")

    for genus in genera:
        ids = [rid for rid, g in assigned.items() if g == genus]

        any_mapped = 0
        both_mapped = 0
        genus_agree = 0
        proper_agree = 0
        proper_agree_mapq10 = 0
        competing_genus = 0

        for rid in ids:
            state = alignments.get(rid)

            if not state:
                continue

            primary = state["primary"]

            if primary:
                any_mapped += 1

            if state["genera"] - {genus}:
                competing_genus += 1

            if 1 not in primary or 2 not in primary:
                continue

            both_mapped += 1

            a = primary[1]
            b = primary[2]

            if a[0] != genus or b[0] != genus:
                continue

            genus_agree += 1

            if a[2] and b[2]:
                proper_agree += 1

                if a[1] >= 10 and b[1] >= 10:
                    proper_agree_mapq10 += 1

        rows.append({
            "Run": run,
            "Genus": genus,
            "Kraken_candidate_pairs": len(ids),
            "At_least_one_primary_mapped": any_mapped,
            "Both_mates_primary_mapped": both_mapped,
            "Both_mates_match_Kraken_genus": genus_agree,
            "Proper_pairs_matching_genus": proper_agree,
            "Proper_matching_pairs_MAPQ10": proper_agree_mapq10,
            "Pairs_with_other_genus_alignment": competing_genus
        })

output = ALIGN / "genus_alignment_support.csv"

with output.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print("\n=== ALIGNMENT SUPPORT SUMMARY ===")
print("Genus | Candidates | Proper genus-matching pairs | MAPQ >= 10")

for genus in genera:
    subset = [r for r in rows if r["Genus"] == genus]

    print(
        genus,
        "|", sum(r["Kraken_candidate_pairs"] for r in subset),
        "|", sum(r["Proper_pairs_matching_genus"] for r in subset),
        "|", sum(r["Proper_matching_pairs_MAPQ10"] for r in subset)
    )

print("\nPASS: All six BAM files summarized.")
print("Detailed results:", output)
