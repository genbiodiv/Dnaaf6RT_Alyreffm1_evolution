#!/usr/bin/env python3
"""Fail-fast consistency checks for the archived Alyreffm1-202 results."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "expression"
PLOTS = ROOT / "plots" / "expression"


def read_tsv(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    summary = json.loads((RESULTS / "validation_summary.json").read_text(encoding="utf-8"))
    assert summary["current_transcripts"] == 10
    assert summary["transcripts_using_inserted_first_exon"] == 6
    assert summary["adult_biological_replicates_with_exact_202_junction_reads"] == 2
    assert summary["exact_202_junction_reads_by_replicate"] == {"ENCFF001JPX": 16, "ENCFF001JPZ": 20}
    assert summary["pacbio_flnc_reads_screened"] == 129_471
    assert summary["pacbio_reads_with_shared_promoter_junction"] == 3
    assert summary["pacbio_reads_with_unique_202_junction"] == 0
    assert summary["candidate_orf_length_aa"] == 127

    junctions = read_tsv("junction_kmer_summary.tsv")
    specific = [row for row in junctions if row["junction"] == "alyreffm1_202_exon2_to_exon3"]
    assert [int(row["exact_24nt_hits"]) for row in specific] == [16, 20]
    assert all(int(row["unique_read_names"]) == int(row["exact_24nt_hits"]) for row in specific)

    lengths = read_tsv("junction_kmer_length_validation.tsv")
    specific_40 = [
        row
        for row in lengths
        if row["junction"] == "alyreffm1_202_exon2_to_exon3" and row["junction_kmer_length"] == "40"
    ]
    assert [int(row["exact_hits"]) for row in specific_40] == [12, 16]
    assert all(int(row["contiguous_occurrences_in_140kb_mm10_locus"]) == 0 for row in lengths)

    peptide = "".join(
        line.strip()
        for line in (RESULTS / "sequences" / "Alyreffm1-202.longest_ORF.faa").read_text().splitlines()
        if not line.startswith(">")
    )
    assert len(peptide) == 127

    required_plots = [
        PLOTS / "expression_evidence_summary.png",
        PLOTS / "alyreffm1_202_evidence_map.png",
    ]
    assert all(path.stat().st_size > 10_000 for path in required_plots)

    print("PASS\tannotation\t10 transcripts; 6 use the insertion-derived first exon")
    print("PASS\tRNA-seq\t16 and 20 exact Alyreffm1-202 junction reads")
    print("PASS\tPacBio\t3 promoter-junction reads; 0 full -202 junction chains")
    print("PASS\tORF\t127-aa candidate reconstructed")
    print("PASS\tplots\tdiagnostic PNG outputs present")


if __name__ == "__main__":
    main()

