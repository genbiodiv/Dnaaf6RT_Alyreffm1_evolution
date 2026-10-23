#!/usr/bin/env python3
"""Regenerate diagnostic expression plots from the archived result tables.

This script does not infer new observations and does not create manuscript
figures.  Every plotted value is read from a TSV produced by the full public-
data workflow.  The purpose of the plots is to make the evidence easy to audit.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/alyreffm1_expression_plots")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "expression"
PLOTS = ROOT / "plots" / "expression"


def read_tsv(name: str) -> list[dict[str, str]]:
    """Return a result table as a list of dictionaries."""
    with (RESULTS / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def plot_summary() -> None:
    """Plot four independent summaries used to inspect expression evidence."""
    timecourse = read_tsv("fantom5_testis_timecourse.tsv")
    junctions = read_tsv("junction_kmer_summary.tsv")
    rsem = read_tsv("encode_rsem_alyreffm1.tsv")
    ribo = read_tsv("riboseq_region_counts.tsv")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4), constrained_layout=True)

    ax = axes[0, 0]
    stages = [row["stage"] for row in timecourse]
    x = list(range(len(stages)))
    ax.plot(
        x,
        [float(row["alyreffm1_promoter_tpm"]) for row in timecourse],
        "o-",
        color="#245a9a",
        lw=2,
        label="Alyreffm1 direction",
    )
    ax.plot(
        x,
        [float(row["dnaaf6rt_promoter_tpm"]) for row in timecourse],
        "s-",
        color="#b24a3a",
        lw=2,
        label="Dnaaf6rt direction",
    )
    ax.set_yscale("symlog", linthresh=0.2)
    ax.set_xticks(x, stages, rotation=45, ha="right")
    ax.set_ylabel("Combined CAGE TPM (symlog)")
    ax.set_title("FANTOM5 testis time course", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    sources = sorted({row["source"] for row in junctions})
    shared = [
        int(next(row["exact_24nt_hits"] for row in junctions if row["source"] == source and row["junction"].startswith("shared")))
        for source in sources
    ]
    unique = [
        int(next(row["exact_24nt_hits"] for row in junctions if row["source"] == source and row["junction"].startswith("alyreffm1_202")))
        for source in sources
    ]
    width = 0.36
    labels = ["Replicate 1", "Replicate 2"]
    ax.bar([i - width / 2 for i in range(len(labels))], shared, width, color="#e5a44a", label="shared exon 1→2")
    ax.bar([i + width / 2 for i in range(len(labels))], unique, width, color="#6a4c93", label="-202-specific exon 2→3")
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylabel("Reads with an exact 24-nt junction k-mer")
    ax.set_title("Exact RNA-seq junction scan", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    rows_202 = [row for row in rsem if row["transcript_id"].split(".")[0] == "ENSMUST00000135245"]
    rows_202.sort(key=lambda row: (row["sample_group"], row["file_accession"]))
    values = [float(row["TPM"]) for row in rows_202]
    labels = [f"{row['file_accession']}\n{row['replicate']}" for row in rows_202]
    ax.bar(range(len(values)), values, color=["#bdbdbd" if value == 0 else "#2a9d8f" for value in values])
    ax.set_xticks(range(len(values)), labels, rotation=35, ha="right", fontsize=7)
    ax.set_ylabel("Alyreffm1-202 RSEM TPM")
    ax.set_title("ENCODE transcript quantification", loc="left", fontweight="bold")

    ax = axes[1, 1]
    region_order = ["5UTR_exon1", "ORF_exon1", "ORF_exon2", "ORF_exon3_with_stop", "3UTR_exon3"]
    region_labels = ["5′ UTR", "ORF exon 1", "ORF exon 2", "ORF exon 3/stop", "3′ UTR"]
    samples = sorted({row["sample"] for row in ribo})
    width = 0.36
    for sample_index, sample in enumerate(samples):
        values = [
            float(next(row["relative_density_within_transcript"] for row in ribo if row["sample"] == sample and row["region"] == region))
            for region in region_order
        ]
        offsets = [index + (sample_index - 0.5) * width for index in range(len(region_order))]
        ax.bar(offsets, values, width, label=f"Replicate {sample_index + 1}", color=("#5b8e7d", "#bc8da0")[sample_index])
    ax.axhline(1, color="#555555", lw=0.8, ls="--")
    ax.set_xticks(range(len(region_order)), region_labels, rotation=30, ha="right")
    ax.set_ylabel("Relative RPF density within transcript")
    ax.set_title("Ribosome-footprint distribution", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTS / "expression_evidence_summary.png", dpi=220)
    plt.close(fig)


def plot_transcript_evidence() -> None:
    """Plot the annotated RNA, candidate ORF, and direct junction support."""
    orfs = read_tsv("alyreffm1_202_orfs.tsv")
    longest = max(orfs, key=lambda row: int(row["protein_length_aa"]))
    junctions = read_tsv("junction_kmer_summary.tsv")
    long_reads = read_tsv("gse137490_long_read_junction_support.tsv")

    shared = [int(row["exact_24nt_hits"]) for row in junctions if row["junction"].startswith("shared")]
    unique = [int(row["exact_24nt_hits"]) for row in junctions if row["junction"].startswith("alyreffm1_202")]
    promoter_reads = sum(row["supports_shared_promoter_junction"] == "yes" for row in long_reads)
    unique_long_reads = sum(row["supports_unique_202_junction"] == "yes" for row in long_reads)

    fig, ax = plt.subplots(figsize=(12, 4.6), constrained_layout=True)
    ax.set_xlim(0, 910)
    ax.set_ylim(-1.25, 2.3)
    ax.set_yticks([])
    ax.set_xlabel("Position in spliced Alyreffm1-202 transcript (nt)")
    ax.set_title("Diagnostic map of Alyreffm1-202 evidence", loc="left", fontweight="bold")

    exon_specs = [(1, 209, "Insertion-derived exon 1", "#e9a03b"), (210, 431, "Exon 2", "#4c78a8"), (432, 889, "-202-specific exon 3", "#4c78a8")]
    for start, end, label, color in exon_specs:
        ax.add_patch(Rectangle((start, 1.25), end - start + 1, 0.34, color=color, alpha=0.9))
        ax.text((start + end) / 2, 1.67, label, ha="center", va="bottom")
    ax.text(1, 1.05, "Annotated RNA: 889 nt", ha="left", color="#333333")

    orf_start = int(longest["start_nt_1based"])
    stop_start = int(longest["stop_codon_start_nt_1based"])
    ax.add_patch(Rectangle((orf_start, 0.48), stop_start + 3 - orf_start, 0.30, color="#6a4c93", alpha=0.9))
    ax.text(orf_start, 0.32, "ATG", ha="left", color="#4f366e")
    ax.text(stop_start + 2, 0.32, longest["stop_codon"], ha="right", color="#4f366e")
    ax.text((orf_start + stop_start) / 2, 0.88, f"Candidate ORF: {longest['protein_length_aa']} aa", ha="center", color="#4f366e", fontweight="bold")

    ax.add_patch(FancyArrowPatch((195, -0.05), (225, -0.05), connectionstyle="arc3,rad=-0.55", arrowstyle="-", lw=2.2, color="#e59f3a"))
    ax.add_patch(FancyArrowPatch((415, -0.05), (448, -0.05), connectionstyle="arc3,rad=-0.55", arrowstyle="-", lw=2.2, color="#6a4c93"))
    ax.text(210, -0.52, f"Shared junction: {' / '.join(map(str, shared))} reads", ha="center", color="#9a641f")
    ax.text(432, -0.52, f"-202-specific junction: {' / '.join(map(str, unique))} reads", ha="center", color="#4f366e")
    ax.text(1, -1.02, f"PacBio testis: {promoter_reads} reads support the first-exon junction; {unique_long_reads} support the complete -202-specific junction chain.", ha="left")

    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTS / "alyreffm1_202_evidence_map.png", dpi=220)
    plt.close(fig)


def main() -> None:
    plot_summary()
    plot_transcript_evidence()
    print(f"Diagnostic plots written to {PLOTS}")


if __name__ == "__main__":
    main()
