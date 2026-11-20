import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

root = Path("references")
old = root / "kraken2/kampala_curated_v1/library/added"
combined = root / "assemblies/enteric_combined_v1/manifest.tsv"
output_dir = root / "assemblies/legacy_unique_v1"
output_dir.mkdir(parents=True, exist_ok=True)

with combined.open() as f:
    new_accessions = {
        row["accession"]
        for row in csv.DictReader(f, delimiter="\t")
    }

expected = Counter()
taxids = defaultdict(set)

with (old / "prelim_map.txt").open() as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 3 or fields[0] != "TAXID":
            continue

        match = re.match(r"(GC[AF]_\d+\.\d+)", fields[1])

        if match:
            accession = match.group(1)

            if accession not in new_accessions:
                expected[accession] += 1
                taxids[accession].add(fields[2])

if len(expected) != 27:
    raise SystemExit(
        f"Expected 27 old-only assemblies, found {len(expected)}"
    )

for accession, assigned_taxids in taxids.items():
    if len(assigned_taxids) != 1:
        raise SystemExit(
            f"Inconsistent taxonomy for {accession}"
        )

names = {}

taxonomy_file = (
    root / "kraken2/kampala_expanded_https_v2/"
    "taxonomy/names.dmp"
)

with taxonomy_file.open() as f:
    for line in f:
        parts = [part.strip() for part in line.split("|")]

        if len(parts) >= 4 and parts[3] == "scientific name":
            names[parts[0]] = parts[1]

output = output_dir / "legacy_unique_27.fna"
temporary = output_dir / "legacy_unique_27.fna.partial"

observed = Counter()
current = None
keep = False

try:
    with (old / "gPFqW1vMOo.fna").open() as source, \
         temporary.open("w") as destination:

        for line in source:

            if line.startswith(">"):
                identifier = line[1:].split()[0]
                match = re.match(
                    r"(GC[AF]_\d+\.\d+)", identifier
                )

                current = match.group(1) if match else None
                keep = current in expected

                if keep:
                    tag = re.search(
                        r"(?:^|\|)kraken:taxid\|(\d+)",
                        identifier
                    )

                    if not tag or tag.group(1) not in taxids[current]:
                        raise ValueError(
                            f"Taxonomy mismatch: {identifier}"
                        )

                    observed[current] += 1

            if keep:
                destination.write(line)

    if observed != expected:
        for accession in sorted(expected):
            if observed[accession] != expected[accession]:
                print(
                    accession,
                    "expected:", expected[accession],
                    "observed:", observed[accession]
                )
        raise RuntimeError("Sequence-count verification failed")

    temporary.replace(output)

except Exception:
    temporary.unlink(missing_ok=True)
    raise

manifest = output_dir / "manifest.tsv"

with manifest.open("w") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "accession", "taxid", "organism", "sequence_count"
    ])

    for accession in sorted(expected):
        taxid = next(iter(taxids[accession]))
        writer.writerow([
            accession,
            taxid,
            names.get(taxid, "Unknown"),
            expected[accession]
        ])

print("\n========== EXTRACTION SUMMARY ==========")
print("Additional assemblies:", len(expected))
print("Extracted sequences:", f"{sum(observed.values()):,}")
print("Taxonomy mismatches: 0")
print("Sequence-count mismatches: 0")
print("\nPASS: Legacy reference extraction completed.")
print("FASTA:", output)
print("Manifest:", manifest)
