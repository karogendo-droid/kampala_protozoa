import csv
import hashlib
import re
import subprocess
from pathlib import Path
from collections import Counter, defaultdict
from functools import lru_cache

ROOT = Path("references/markers")
PR2 = ROOT / "pr2_v5.1.1"
COMBINED = ROOT / "integrated_stage4_v1/integrated_marker_competitive.fasta"
BASE = Path("results/protozoa_validation_v3/new_marker_validation")
KRAKEN = Path("results/kraken2_screen_v3")
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
GENERA = list(TARGETS.values())

def fasta_records(path):
    name, parts = None, []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            yield name, "".join(parts)

def digest(sequence):
    sequence = sequence.upper().replace("-", "").replace(".", "")
    return hashlib.sha1(sequence.encode()).hexdigest()

def tax_genus(annotation):
    for field in annotation.strip().split(";"):
        field = field.strip().replace(" ", "_").lower()
        for genus in GENERA:
            g = genus.lower()
            if field == g or field.startswith(g + "_"):
                return genus
    return "Non-target"

def clean_id(value):
    value = value.split()[0]
    return value[:-2] if value.endswith(("/1", "/2")) else value

print("1. Reading existing marker alignments...", flush=True)

pairs_by_run = {}
needed_refs = set()
audit = {}

for n in range(69, 75):
    run = f"SRR284341{n}"
    bam = BASE / f"{run}.bam"
    pairs = defaultdict(dict)
    counts = Counter()

    process = subprocess.Popen(
        ["samtools", "view", "-f", "2", "-F", "2308", str(bam)],
        stdout=subprocess.PIPE,
        text=True
    )

    for line in process.stdout:
        p = line.split("\t", 6)
        ref = p[2]
        counts["all"] += 1

        if not ref.startswith("PR2_"):
            counts["other"] += 1
            continue

        counts["pr2"] += 1
        needed_refs.add(ref)

        flag = int(p[1])
        mate = 1 if flag & 64 else 2 if flag & 128 else 0
        mapq = int(p[4])

        if mate and 10 <= mapq < 255:
            rid = clean_id(p[0])
            if mate in pairs[rid]:
                raise RuntimeError("Duplicate primary mate: " + rid)
            pairs[rid][mate] = (ref, int(p[3]))

    if process.wait() != 0:
        raise RuntimeError("samtools failed for " + run)

    pairs_by_run[run] = pairs
    audit[run] = counts

print("PR2 reference IDs to identify:", len(needed_refs), flush=True)

print("2. Matching combined reference sequences...", flush=True)

hash_to_refs = defaultdict(set)

for ref, sequence in fasta_records(COMBINED):
    if ref in needed_refs:
        hash_to_refs[digest(sequence)].add(ref)

found_refs = set().union(*hash_to_refs.values()) if hash_to_refs else set()

if found_refs != needed_refs:
    raise SystemExit(
        f"ERROR: {len(needed_refs - found_refs)} mapped PR2 "
        "references were not found in the combined FASTA."
    )

print("3. Recovering original PR2 taxonomy...", flush=True)

taxonomy = {}

tax_path = PR2 / "pr2_version_5.1.1_SSU_mothur.tax"
fasta_path = PR2 / "pr2_version_5.1.1_SSU_mothur.fasta"

with tax_path.open() as f:
    for line in f:
        fields = line.rstrip("\n").split("\t", 1)
        if len(fields) == 2:
            taxonomy[fields[0]] = fields[1]

ref_labels = defaultdict(set)

for original_id, sequence in fasta_records(fasta_path):
    h = digest(sequence)
    if h not in hash_to_refs:
        continue

    label = tax_genus(taxonomy.get(original_id, ""))

    for ref in hash_to_refs[h]:
        ref_labels[ref].add(label)

matched_refs = set(ref_labels)
print("Mapped PR2 IDs:", len(needed_refs))
print("Matched to original taxonomy:", len(matched_refs))

if not needed_refs or len(matched_refs) < 0.9 * len(needed_refs):
    raise SystemExit(
        "STOP: Insufficient exact sequence-to-taxonomy matches. "
        "Do not interpret zero counts as biological negatives."
    )

# If identical sequences have different taxonomic labels,
# do not assign them to any one target genus.
ref_genus = {
    ref: next(iter(labels))
    for ref, labels in ref_labels.items()
    if len(labels) == 1 and next(iter(labels)) in GENERA
}

print("4. Checking marker and Kraken2 genus agreement...", flush=True)

parents = {}
with NODES.open() as f:
    for line in f:
        p = line.split("|")
        parents[int(p[0].strip())] = int(p[1].strip())

@lru_cache(None)
def kraken_genus(taxid):
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

rows = []

for run, pairs in pairs_by_run.items():
    eligible = {}

    for rid, mates in pairs.items():
        if set(mates) != {1, 2}:
            continue

        ref1, pos1 = mates[1]
        ref2, pos2 = mates[2]

        if ref1 != ref2:
            continue

        genus = ref_genus.get(ref1)

        if genus:
            eligible[rid] = (
                genus, ref1, pos1, pos2
            )

    concordant = {}

    if eligible:
        with (KRAKEN / f"{run}.kraken").open() as f:
            for line in f:
                p = line.rstrip("\n").split("\t", 3)

                if p[0] != "C":
                    continue

                rid = clean_id(p[1])

                if rid in eligible:
                    genus = eligible[rid][0]

                    if kraken_genus(int(p[2])) == genus:
                        concordant[rid] = eligible[rid]

    for genus in GENERA:
        marker_matches = {
            rid: x for rid, x in eligible.items()
            if x[0] == genus
        }

        accepted = [
            x for x in concordant.values()
            if x[0] == genus
        ]

        bins = set()

        for _, ref, p1, p2 in accepted:
            bins.add((ref, (p1 - 1) // 100))
            bins.add((ref, (p2 - 1) // 100))

        rows.append({
            "Run": run,
            "Genus": genus,
            "Marker_supported_pairs": len(marker_matches),
            "Kraken_marker_concordant_pairs": len(accepted),
            "Distinct_100bp_bins": len(bins)
        })

output = BASE / "final_corrected_marker_evidence.csv"

with output.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print("\n=== FINAL CORRECTED MARKER EVIDENCE ===")
print("Genus | Marker pairs | Kraken-concordant pairs | 100-bp bins")

for genus in GENERA:
    subset = [r for r in rows if r["Genus"] == genus]
    print(
        genus,
        "|", sum(r["Marker_supported_pairs"] for r in subset),
        "|", sum(r["Kraken_marker_concordant_pairs"] for r in subset),
        "|", sum(r["Distinct_100bp_bins"] for r in subset)
    )

print("\n=== PR2 ALIGNMENT AUDIT ===")
for run, c in audit.items():
    print(run, "| PR2:", c["pr2"], "| Other:", c["other"])

print("\nSaved:", output)
print("Analysis complete. No alignments repeated.")
