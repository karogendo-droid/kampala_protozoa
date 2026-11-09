import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata" / "taxonomy"

def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def read_catalog(path):
    with path.open(newline="", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#assembly_accession\t"):
                header = line[1:].rstrip("\r\n").split("\t")
                break
        else:
            raise ValueError(f"Missing assembly header: {path}")

        return list(csv.DictReader(
            f, fieldnames=header, delimiter="\t"
        ))

def number(value):
    try:
        return int(value)
    except (ValueError, TypeError):
        return 0

catalogs = {}

for database in ("GenBank", "RefSeq"):
    for group in ("protozoa", "fungi"):
        path = META / (
            f"{database.lower()}_{group}_assembly_summary.txt"
        )

        for record in read_catalog(path):
            accession = record.get("assembly_accession", "")
            catalogs[(database, group, accession)] = record

targets = read_tsv(META / "protozoa_target_manifest.tsv")
candidates = read_tsv(
    META / "combined_reference_candidates.tsv"
)

fields = [
    "target", "priority", "database", "genome_group",
    "assembly_accession", "organism_name", "taxid",
    "infraspecific_name", "isolate", "asm_name",
    "refseq_category", "assembly_level", "genome_rep",
    "genome_size", "contig_count", "scaffold_count",
    "seq_rel_date", "gbrs_paired_asm", "paired_asm_comp",
    "ftp_path", "review_rank", "review_status"
]

by_target = defaultdict(list)

for candidate in candidates:
    key = (
        candidate["database"],
        candidate["genome_group"],
        candidate["assembly_accession"]
    )

    if key not in catalogs:
        raise ValueError(f"Missing catalog record: {key}")

    original = catalogs[key]

    result = dict(candidate)

    for field in fields:
        if field not in ("review_rank", "review_status"):
            if field in original:
                result[field] = original[field]

    result["review_status"] = "NOT_YET_SELECTED"
    by_target[result["target"]].append(result)

def ranking(record):
    levels = {
        "Complete Genome": 0,
        "Chromosome": 1,
        "Scaffold": 2,
        "Contig": 3,
    }

    contigs = number(record.get("contig_count"))
    if contigs == 0:
        contigs = 10**12

    return (
        0 if record.get("genome_rep") == "Full" else 1,
        levels.get(record.get("assembly_level"), 9),
        0 if record.get("database") == "RefSeq" else 1,
        contigs,
        record["assembly_accession"]
    )

def counterpart_key(record):
    accession = record["assembly_accession"]
    paired = record.get("gbrs_paired_asm", "")

    if (
        record.get("paired_asm_comp") == "identical"
        and paired.startswith(("GCA_", "GCF_"))
    ):
        return tuple(sorted((accession, paired)))

    return (accession,)

all_records = []
shortlist = []

for target in targets:
    name = target["organism"]
    records = sorted(by_target[name], key=ranking)
    seen_pairs = set()

    for rank, record in enumerate(records, start=1):
        record["review_rank"] = rank
        all_records.append(record)

        pair = counterpart_key(record)
        if pair not in seen_pairs and len(
            [r for r in shortlist if r["target"] == name]
        ) < 5:
            seen_pairs.add(pair)
            shortlist.append(record)

    print(f"\n{name}: {len(records)} candidate records")

    for record in [r for r in shortlist if r["target"] == name]:
        print(
            "  ", record["assembly_accession"],
            "|", record["assembly_level"],
            "| contigs:", record.get("contig_count", "NA"),
            "| genome bp:", record.get("genome_size", "NA")
        )

for filename, records in (
    ("reference_quality_review.tsv", all_records),
    ("reference_shortlist_for_review.tsv", shortlist),
):
    with (META / filename).open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f, fieldnames=fields, delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(records)

print("\nReview tables saved in metadata/taxonomy/")
print("No genomes downloaded or selected.")
