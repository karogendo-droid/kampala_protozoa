import gzip
from pathlib import Path

source = Path(
    "references/human/GCF_000001405.40_GRCh38.p14/"
    "GCF_000001405.40_GRCh38.p14_genomic.fna"
)

folder = Path("references/assemblies/competing_v1")
folder.mkdir(parents=True, exist_ok=True)

output = folder / "human_GRCh38_kraken.fna.gz"
temporary = folder / "human_GRCh38_kraken.fna.gz.partial"

if not source.is_file():
    raise SystemExit("ERROR: Human reference FASTA not found.")

if output.exists():
    raise SystemExit(
        "Output already exists. Existing file was not modified."
    )

sequences = 0
bases = 0
identifiers = set()

try:
    with source.open("rt") as inp, gzip.open(
        temporary, "wt", compresslevel=3
    ) as out:

        for line in inp:

            if line.startswith(">"):
                original = line[1:].strip()
                identifier = original.split()[0]

                if identifier in identifiers:
                    raise ValueError(
                        f"Duplicate sequence ID: {identifier}"
                    )

                identifiers.add(identifier)

                out.write(
                    f">{identifier}|kraken:taxid|9606 "
                    f"{original}\n"
                )

                sequences += 1

            else:
                sequence = line.strip()

                if sequence:
                    if sequences == 0:
                        raise ValueError(
                            "Sequence encountered before FASTA header."
                        )

                    bases += len(sequence)
                    out.write(sequence + "\n")

    if sequences == 0 or bases == 0:
        raise ValueError("No valid reference sequences found.")

    temporary.replace(output)

except Exception:
    temporary.unlink(missing_ok=True)
    raise

print("\n========== HUMAN REFERENCE SUMMARY ==========")
print("Organism: Homo sapiens")
print("NCBI taxonomy ID: 9606")
print("FASTA sequences:", f"{sequences:,}")
print("Total bases:", f"{bases:,}")
print("Output:", output)
print("\nPASS: Human reference prepared successfully.")
