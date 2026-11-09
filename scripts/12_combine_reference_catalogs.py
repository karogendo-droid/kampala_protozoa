import csv
from collections import Counter
from pathlib import Path

ROOT = Path.home() / "bioinformatics/kampala_protozoa"
META = ROOT / "metadata/taxonomy"

CATALOGS = [
    ("GenBank", "protozoa"),
    ("RefSeq", "protozoa"),
    ("GenBank", "fungi"),
    ("RefSeq", "fungi"),
]

ALIASES = {
    "Cystoisospora belli": [
        "Cystoisospora belli",
        "Isospora belli",
    ],
}

DETAIL = META / "combined_reference_candidates.tsv"
SUMMARY = META / "combined_reference_coverage.tsv"


def load_catalog(database, group):
    path = META / f"{database.lower()}_{group}_assembly_summary.txt"

    with path.open(encoding="utf-8", newline="") as file:
        for line in file:
            if line.startswith("#assembly_accession\t"):
                header = line[1:].rstrip("\r\n").split("\t")
                break
        else:
            raise ValueError(f"No assembly header in {path}")

        records = []

        for row in csv.DictReader(
            file, fieldnames=header, delimiter="\t"
        ):
            if row.get("version_status") != "latest":
                continue

            row["database"] = database
            row["genome_group"] = group
            records.append(row)

    return records


def matches(organism, target):
    organism = organism.casefold().strip()

    if target == "Blastocystis sp.":
        return organism.startswith("blastocystis ")

    names = ALIASES.get(target, [target])

    return any(
        organism == name.casefold()
        or organism.startswith(name.casefold() + " ")
        for name in names
    )


def main():
    with (META / "protozoa_target_manifest.tsv").open(
        encoding="utf-8", newline=""
    ) as file:
        targets = list(csv.DictReader(file, delimiter="\t"))

    assemblies = []

    for database, group in CATALOGS:
        rows = load_catalog(database, group)
        assemblies.extend(rows)
        print(f"{database} {group}: {len(rows)} current records")

    fields = [
        "target", "priority", "database", "genome_group",
        "assembly_accession", "organism_name", "taxid",
        "species_taxid", "assembly_level", "genome_size",
        "genome_rep", "version_status", "seq_rel_date",
        "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
    ]

    coverage_fields = [
        "target", "protozoa_genbank", "protozoa_refseq",
        "fungi_genbank", "fungi_refseq", "total_records", "status",
    ]

    with DETAIL.open("w", newline="") as detail_file, \
         SUMMARY.open("w", newline="") as summary_file:

        detail_writer = csv.DictWriter(
            detail_file, fieldnames=fields, delimiter="\t"
        )
        summary_writer = csv.DictWriter(
            summary_file, fieldnames=coverage_fields, delimiter="\t"
        )

        detail_writer.writeheader()
        summary_writer.writeheader()

        for target in targets:
            name = target["organism"]

            found = [
                record for record in assemblies
                if matches(record.get("organism_name", ""), name)
            ]

            counts = Counter(
                (r["genome_group"], r["database"]) for r in found
            )

            for record in found:
                output = {
                    "target": name,
                    "priority": target["priority"],
                }

                for field in fields[2:]:
                    output[field] = record.get(field, "")

                detail_writer.writerow(output)

            summary_writer.writerow({
                "target": name,
                "protozoa_genbank": counts["protozoa", "GenBank"],
                "protozoa_refseq": counts["protozoa", "RefSeq"],
                "fungi_genbank": counts["fungi", "GenBank"],
                "fungi_refseq": counts["fungi", "RefSeq"],
                "total_records": len(found),
                "status": (
                    "CANDIDATES_FOUND" if found
                    else "NOT_FOUND_IN_FOUR_CATALOGS"
                ),
            })

            print(f"{name}: {len(found)} candidate records")

    print("\nSaved:", DETAIL)
    print("Saved:", SUMMARY)


if __name__ == "__main__":
    main()
