import gzip
from pathlib import Path
from collections import Counter
from functools import lru_cache

root = Path("data/processed/human_screen")
results = Path("results/kraken2_screen_v3")
output = Path("results/protozoa_validation_v3")
output.mkdir(parents=True, exist_ok=True)

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

taxonomy = Path(
    "references/kraken2/kampala_screen_v3/"
    "taxonomy/nodes.dmp"
)

with taxonomy.open() as f:
    for line in f:
        parts = line.split("|")
        parents[int(parts[0].strip())] = int(parts[1].strip())

@lru_cache(maxsize=None)
def target_genus(taxid):
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

def read_id(header):
    identifier = header.split()[0].lstrip("@")
    if identifier.endswith(("/1", "/2")):
        identifier = identifier[:-2]
    return identifier

for number in range(69, 75):
    run = f"SRR284341{number}"
    print(f"\nProcessing {run}", flush=True)

    classifications = results / f"{run}.kraken"

    r1_files = list(
        root.rglob(f"{run}_R1.human_screened.fastq.gz")
    )
    r2_files = list(
        root.rglob(f"{run}_R2.human_screened.fastq.gz")
    )

    if len(r1_files) != 1 or len(r2_files) != 1:
        raise RuntimeError(f"FASTQ files not uniquely found: {run}")

    if not classifications.is_file():
        raise RuntimeError(f"Kraken output missing: {run}")

    candidates = {}

    with classifications.open() as f:
        for line in f:
            fields = line.split("\t", 3)

            if fields[0] != "C":
                continue

            genus = target_genus(int(fields[2]))

            if genus:
                candidates[read_id(fields[1])] = genus

    print(f"Candidate read pairs: {len(candidates):,}", flush=True)

    out1 = output / f"{run}_candidate_R1.fastq.gz"
    out2 = output / f"{run}_candidate_R2.fastq.gz"

    tmp1 = output / f"{run}_candidate_R1.fastq.gz.partial"
    tmp2 = output / f"{run}_candidate_R2.fastq.gz.partial"

    counts = Counter()
    total_pairs = 0

    with gzip.open(r1_files[0], "rt") as f1, \
         gzip.open(r2_files[0], "rt") as f2, \
         gzip.open(tmp1, "wt", compresslevel=3) as o1, \
         gzip.open(tmp2, "wt", compresslevel=3) as o2:

        while True:
            a = [f1.readline() for _ in range(4)]
            b = [f2.readline() for _ in range(4)]

            if not a[0] and not b[0]:
                break

            if not all(a) or not all(b):
                raise RuntimeError(f"Incomplete FASTQ record: {run}")

            id1 = read_id(a[0])
            id2 = read_id(b[0])

            if id1 != id2:
                raise RuntimeError(
                    f"Read-pair ID mismatch: {id1} / {id2}"
                )

            total_pairs += 1

            genus = candidates.get(id1)

            if genus:
                o1.writelines(a)
                o2.writelines(b)
                counts[genus] += 1

    if sum(counts.values()) != len(candidates):
        raise RuntimeError(
            f"Extraction mismatch for {run}: "
            f"{sum(counts.values())} extracted, "
            f"{len(candidates)} expected"
        )

    tmp1.replace(out1)
    tmp2.replace(out2)

    print(f"Original read pairs examined: {total_pairs:,}")
    print(f"Successfully extracted: {sum(counts.values()):,}")

    for genus, count in sorted(counts.items()):
        print(f"  {genus}: {count:,}")

print("\nPASS: Candidate extraction completed for all six runs.")
