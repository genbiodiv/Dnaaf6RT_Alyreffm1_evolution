#!/usr/bin/env python3
"""Scan FASTQ reads for exact sequence spanning diagnostic splice junctions.

The input may be an uncompressed FASTQ path or ``-`` for standard input.  Exact
24-mers use 12 nt from each side of a junction; both orientations are tested.
This deliberately conservative scan yields directly inspectable reads and does
not infer abundance from ambiguous transcript assignment.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


MOTIFS = {
    "shared_promoter_exon1_to_exon2": "TCAGCCACAGAGTTCATACTTCTT",
    "alyreffm1_202_exon2_to_exon3": "CGCCAACACCAGGGTTCGGCATCA",
}


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq", default="-", help="Uncompressed FASTQ or '-' for stdin")
    parser.add_argument("--source", required=True, help="Stable source/accession label")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    handle = sys.stdin if args.fastq == "-" else open(args.fastq, encoding="utf-8")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    total = 0
    try:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not (sequence and plus and quality):
                raise ValueError("Truncated FASTQ record")
            total += 1
            for junction, motif in MOTIFS.items():
                rc = reverse_complement(motif)
                orientation = "+" if motif in sequence else "-" if rc in sequence else None
                if orientation:
                    rows.append(
                        {
                            "source": args.source,
                            "read_name": header.strip().removeprefix("@"),
                            "junction": junction,
                            "orientation": orientation,
                            "motif_24nt": motif if orientation == "+" else rc,
                            "read_length": len(sequence),
                            "read_sequence": sequence,
                        }
                    )
    finally:
        if handle is not sys.stdin:
            handle.close()

    columns = [
        "source",
        "read_name",
        "junction",
        "orientation",
        "motif_24nt",
        "read_length",
        "read_sequence",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"FASTQ_READS\t{total}")
    print(f"MATCHES\t{len(rows)}")
    for junction in MOTIFS:
        print(f"{junction}\t{sum(row['junction'] == junction for row in rows)}")


if __name__ == "__main__":
    main()
