import csv
import hashlib
import subprocess
from pathlib import Path

folder = Path("references/assemblies/enteric_refseq_v1")
manifest = folder / "manifest.tsv"

def md5sum(path):
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def download(url, path):
    subprocess.run([
        "curl", "-fL",
        "--retry", "5",
        "--retry-all-errors",
        "--retry-delay", "10",
        "--connect-timeout", "30",
        "--max-time", "1800",
        "-o", str(path),
        url
    ], check=True)

with manifest.open() as f:
    records = list(csv.DictReader(f, delimiter="\t"))

success = []
failed = []

for number, record in enumerate(records, 1):
    accession = record["accession"]
    organism = record["organism"]

    filename = record["fasta_url"].split("/")[-1]
    destination = folder / filename
    checksum_file = folder / f"{accession}_md5checksums.txt"

    print(
        f"\n[{number}/{len(records)}] {organism}",
        flush=True
    )

    try:
        download(record["md5_url"], checksum_file)

        expected = None

        for line in checksum_file.read_text().splitlines():
            parts = line.split()

            if len(parts) >= 2 and \
               Path(parts[-1].lstrip("*")).name == filename:
                expected = parts[0].lower()
                break

        if expected is None:
            raise RuntimeError("NCBI checksum not found")

        if destination.exists():
            if md5sum(destination) == expected:
                print("PASS: Already downloaded and verified.")
                success.append(accession)
                continue

        verified = False

        for attempt in range(1, 4):
            print(f"Download attempt {attempt}/3", flush=True)

            temporary = folder / (filename + ".part")

            try:
                download(record["fasta_url"], temporary)

                if md5sum(temporary) == expected:
                    temporary.replace(destination)
                    print("PASS: MD5 verified.", flush=True)
                    verified = True
                    break

                print("WARNING: MD5 mismatch.", flush=True)

            except subprocess.CalledProcessError:
                print("WARNING: Connection interrupted.", flush=True)

            finally:
                temporary.unlink(missing_ok=True)

        if verified:
            success.append(accession)
        else:
            failed.append(accession)

    except (subprocess.CalledProcessError, RuntimeError) as error:
        print(f"FAILED: {error}", flush=True)
        failed.append(accession)

print("\n========== SUMMARY ==========")
print(f"Total assemblies: {len(records)}")
print(f"Verified genomes: {len(success)}")
print(f"Failed genomes:   {len(failed)}")

if failed:
    print("Incomplete accessions:")
    for accession in failed:
        print(" -", accession)
    print("\nRun this script again to retry missing files.")
    raise SystemExit(1)

print("\nSUCCESS: All reference genomes verified.")
