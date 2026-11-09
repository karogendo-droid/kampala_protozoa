import csv
from collections import Counter
from pathlib import Path

ROOT = Path.home() / "bioinformatics/kampala_protozoa"
META = ROOT / "metadata/taxonomy"

TARGET_FILE = META / "protozoa_target_manifest.tsv"

CATALOGS = {
    "GenBank": META / "genbank_protozoa_assembly_summary.txt",
    "RefSeq": META / "refseq_protozoa_assembly_summary.txt",
}

OUTPUT = META / "protozoa_assembly_candidates.tsv"
SUMMARY = META / "protozoa_coverage_summary.tsv"


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def read_assembly_catalog(path, source):
    records = []

    with path.open(encoding="utf-8") as file:

        header = None

        for line in file:
            if line.startswith("#assembly_accession\t"):
                header = line[1:].rstrip("\r\n").split("\t")
                break

        if header is None:
            raise ValueError(f"Assembly header not found: {path}")

        for row in csv.DictReader(
            file, fieldnames=header, delimiter="\t"
        ):
            if not row.get("assembly_accession"):
                continue

            row["source"] = source
            records.append(row)

    return records


def matches_target(organism_name, target_name):
    organism = organism_name.strip().casefold()
    target = target_name.strip().casefold()

    if target == "blastocystis sp.":
        return organism.startswith("blastocystis ")

    return organism == target or organism.startswith(target + " ")


def main():

    targets = read_tsv(TARGET_FILE)

    assemblies = []

    for source, path in CATALOGS.items():
        records = read_assembly_catalog(path, source)
        print(f"{source}: {len(records)} assembly records loaded")
        assemblies.extend(records)

    fields = [
        "target_group",
        "target_organism",
        "priority",
        "source",
        "assembly_accession",
        "organism_name",
        "taxid",
        "species_taxid",
        "assembly_level",
        "version_status",
        "genome_size",
        "seq_rel_date",
        "gbrs_paired_asm",
        "paired_asm_comp",
        "ftp_path",
    ]

    summary_fields = [
        "target_organism",
        "priority",
        "genbank_count",
        "refseq_count",
        "status",
    ]

    with OUTPUT.open("w", newline="") as out, \
         SUMMARY.open("w", newline="") as summary:

        writer = csv.DictWriter(
            out, fieldnames=fields, delimiter="\t"
        )
        summary_writer = csv.DictWriter(
            summary, fieldnames=summary_fields, delimiter="\t"
        )

        writer.writeheader()
        summary_writer.writeheader()

        for target in targets:

            name = target["organism"]

            matches = [
                a for a in assemblies
                if matches_target(a.get("organism_name", ""), name)
                and a.get("version_status", "").lower() == "latest"
            ]

            counts = Counter(a["source"] for a in matches)

            for assembly in matches:
                record = {
                    "target_group": target["target_group"],
                    "target_organism": name,
                    "priority": target["priority"],
                    "source": assembly["source"],
                }

                for field in fields[4:]:
                    record[field] = assembly.get(field, "")

                writer.writerow(record)

            summary_writer.writerow({
                "target_organism": name,
                "priority": target["priority"],
                "genbank_count": counts["GenBank"],
                "refseq_count": counts["RefSeq"],
                "status": (
                    "CANDIDATES_FOUND"
                    if matches else "NOT_FOUND_IN_THESE_CATALOGS"
                ),
            })

            print(
                f"{name}: "
                f"GenBank={counts['GenBank']}, "
                f"RefSeq={counts['RefSeq']}"
            )

    print("\nCandidate table:", OUTPUT)
    print("Coverage summary:", SUMMARY)


if __name__ == "__main__":
    main()
