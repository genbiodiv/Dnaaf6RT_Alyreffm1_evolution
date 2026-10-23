#!/usr/bin/env python3
"""Download fixed, accession-versioned Lgsn loci and record their provenance."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "data" / "synteny" / "species.tsv"
RAW = ROOT / "data" / "synteny" / "loci"
META = ROOT / "data" / "synteny" / "metadata"
USER_AGENT = "Dnaaf6rt-synteny-reanalysis/1.0 (contact: rneme@uninorte.edu.co)"


def request_bytes(url: str, attempts: int = 5) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=90) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    META.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(CONFIG.open(encoding="utf-8"), delimiter="\t"))
    gene_ids = ",".join(row["lgsn_gene_id"] for row in rows)
    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "gene", "retmode": "json", "id": gene_ids})
    )
    summary_path = META / "ncbi_gene_esummary.json"
    summary_path.write_bytes(request_bytes(summary_url))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    provenance = []
    for row in rows:
        gene = summary["result"][row["lgsn_gene_id"]]
        genomic = gene["genomicinfo"][0]
        observed_accession = genomic["chraccver"]
        if observed_accession != row["sequence_accession"]:
            raise RuntimeError(
                f"Pinned accession mismatch for {row['scientific_name']}: "
                f"{row['sequence_accession']} != {observed_accession}"
            )
        observed_start = min(int(genomic["chrstart"]), int(genomic["chrstop"])) + 1
        observed_end = max(int(genomic["chrstart"]), int(genomic["chrstop"])) + 1
        observed_strand = "+" if int(genomic["chrstart"]) < int(genomic["chrstop"]) else "-"
        pinned = (int(row["lgsn_start_1based"]), int(row["lgsn_end_1based"]), row["lgsn_strand"])
        observed = (observed_start, observed_end, observed_strand)
        if pinned != observed:
            raise RuntimeError(
                f"Pinned Lgsn coordinates mismatch for {row['scientific_name']}: "
                f"{pinned} != {observed}"
            )
        current_assembly = gene.get("locationhist", [{}])[0].get("assemblyaccver")
        if current_assembly != row["assembly_accession"]:
            raise RuntimeError(
                f"Pinned assembly mismatch for {row['scientific_name']}: "
                f"{row['assembly_accession']} != {current_assembly}"
            )
        start = int(row["region_start_1based"])
        end = int(row["region_end_1based"])
        params = {
            "db": "nuccore",
            "id": row["sequence_accession"],
            "seq_start": str(start),
            "seq_stop": str(end),
            "rettype": "fasta",
            "strand": "1",
        }
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
        output = RAW / f"{row['species_slug']}.{row['sequence_accession']}.{start}-{end}.fasta"
        output.write_bytes(request_bytes(url))
        if not output.read_text(encoding="utf-8").startswith(">"):
            raise RuntimeError(f"NCBI did not return FASTA for {row['scientific_name']}")
        provenance.append(
            {
                "species_slug": row["species_slug"],
                "scientific_name": row["scientific_name"],
                "assembly_accession": row["assembly_accession"],
                "sequence_accession": row["sequence_accession"],
                "region_start_1based": start,
                "region_end_1based": end,
                "retrieved_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
                "sha256": sha256(output),
                "source_url": url,
                "local_file": str(output.relative_to(ROOT)),
            }
        )
        time.sleep(0.4)

    columns = list(provenance[0])
    with (META / "locus_provenance.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(provenance)
    manifest = {
        "retrieval_date_utc": dt.datetime.now(dt.timezone.utc).date().isoformat(),
        "species_count": len(rows),
        "anchor_gene": "Lgsn",
        "window": "60 kb on each side of the pinned NCBI Lgsn gene interval",
        "gene_summary_url": summary_url,
        "gene_summary_sha256": sha256(summary_path),
    }
    (META / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Downloaded {len(rows)} accession-versioned loci to {RAW}")


if __name__ == "__main__":
    main()
