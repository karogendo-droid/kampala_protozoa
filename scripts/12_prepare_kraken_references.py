import csv
import gzip
from pathlib import Path
from collections import Counter

root = Path("references/assemblies")

manifest = root / "enteric_combined_v1/manifest.tsv"

output = root / "enteric_combined_v1/enteric_targets_28.fna"
temporary = output.with_suffix(".fna.partial")

genus_counts = Counter()
sequence_count = 0
total_bases = 0

with manifest.open() as f:
    records = list(csv.DictReader(f, delimiter="\t"))

if len(records) != 28:
    raise SystemExit("Expected 28 reference assemblies.")

try:
    with temporary.open("w") as out:

        for number, record in enumerate(records, 1):

            accession = record["accession"]
            taxid = record["taxid"]
            organism = record["organism"]

            source = (
                root / record["collection"] / record["filename"]
            )

            print(
                f"[{number}/28] Processing {organism}",
                flush=True
            )

            contig_count = 0
            has_sequence = False

            with gzip.open(source, "rt") as genome:

                for line in genome:

                    if line.startswith(">"):

                        if contig_count > 0 and not has_sequence:
                            raise ValueError(
                                f"Empty sequence in {accession}"
                            )

                        original_id = line[1:].split()[0]

                        new_id = (
                            f"kraken:taxid|{taxid}|"
                            f"{accession}_{original_id}"
                        )

                        out.write(
                            f">{new_id} {organism}\n"
                        )

                        contig_count += 1
                        sequence_count += 1
                        has_sequence = False

                    else:

                        sequence = line.strip()

                        if sequence:

                            if contig_count == 0:
                                raise ValueError(
                                    f"Sequence before FASTA header: "
                                    f"{accession}"
                                )

                            out.write(sequence + "\n")

                            total_bases += len(sequence)
                            has_sequence = True

            if contig_count == 0 or not has_sequence:
                raise ValueError(
                    f"Missing or empty FASTA sequence: {accession}"
                )

            genus = organism.split()[0]
            genus_counts[genus] += 1

    temporary.replace(output)

except Exception:
    temporary.unlink(missing_ok=True)
    raise

print("\n========== PREPARATION SUMMARY ==========")
print("Reference assemblies:", len(records))
print("Total FASTA sequences:", sequence_count)
print("Total nucleotide bases:", f"{total_bases:,}")

print("\nAssemblies by genus:")

for genus, count in sorted(genus_counts.items()):
    print(f"{genus}: {count}")

print("\nPASS: All 28 references prepared successfully.")
print("Output:", output)
