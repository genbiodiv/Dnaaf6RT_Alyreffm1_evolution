#!/usr/bin/env python3
"""Reanalyse the focal Dnaaf6/Dnaaf6rt/Alyreffm1 retroposition event.

Coordinates are 1-based and inclusive unless a BED file is being parsed.
The analysis intentionally separates three questions:
1) Was there a processed retrotransposition?
2) Did it copy a bidirectional regulatory segment?
3) How is that segment used by current Alyreffm1 transcript models?

Only accessioned sequence and annotation files under ``data/event`` are read.
Tables are written to ``results/event`` and diagnostic plots to ``plots/event``.
No manuscript figure is assembled by this script.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/retrotransposition_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle
from Bio import Align, SeqIO


MM_CHR1_REGION_START = 31_155_362
MM_CHRX_REGION_START = 138_955_000

LEFT_TSD = (31_261_459, 31_261_472)
INSERTION = (31_261_473, 31_263_369)
RIGHT_TSD = (31_263_370, 31_263_383)
DNAAF6RT = (31_261_920, 31_263_368)
COPIED_SOURCE = (138_963_238, 138_963_709)
COPIED_TARGET = (31_261_483, 31_261_938)
NUP62CL_END = 138_963_456
DNAAF6_START = 138_963_628

ALYREF_TSS = [(31_261_685, 31_261_701), (31_261_708, 31_261_724), (31_261_725, 31_261_758)]
DNAAF6RT_TSS = [(31_261_863, 31_261_873), (31_261_889, 31_261_902)]


@dataclass
class AlignmentSummary:
    score: float
    matches: int
    ungapped_bases: int
    alignment_columns: int
    identity_ungapped_pct: float
    query_start_0: int
    query_end_0: int
    target_start_0: int
    target_end_0: int


def fasta_records(path: Path) -> dict[str, str]:
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}


def record_with_prefix(records: dict[str, str], prefix: str) -> str:
    for key, value in records.items():
        if key.startswith(prefix):
            return value
    raise KeyError(f"No FASTA record begins with {prefix!r}")


def genomic_slice(sequence: str, region_start: int, start: int, end: int) -> str:
    return sequence[start - region_start : end - region_start + 1]


def aligner(mode: str = "local") -> Align.PairwiseAligner:
    obj = Align.PairwiseAligner()
    obj.mode = mode
    obj.match_score = 2
    obj.mismatch_score = -2
    obj.open_gap_score = -5
    obj.extend_gap_score = -1
    return obj


def summarize_alignment(query: str, target: str, mode: str = "local") -> AlignmentSummary:
    aln = aligner(mode).align(query, target)[0]
    qcoords, tcoords = aln.coordinates
    matches = 0
    ungapped = 0
    columns = 0
    for idx in range(len(qcoords) - 1):
        dq = int(qcoords[idx + 1] - qcoords[idx])
        dt = int(tcoords[idx + 1] - tcoords[idx])
        columns += max(dq, dt)
        if dq == dt:
            qblock = query[int(qcoords[idx]) : int(qcoords[idx + 1])]
            tblock = target[int(tcoords[idx]) : int(tcoords[idx + 1])]
            matches += sum(a == b for a, b in zip(qblock, tblock))
            ungapped += dq
    return AlignmentSummary(
        score=float(aln.score),
        matches=matches,
        ungapped_bases=ungapped,
        alignment_columns=columns,
        identity_ungapped_pct=100.0 * matches / ungapped,
        query_start_0=int(qcoords[0]),
        query_end_0=int(qcoords[-1]),
        target_start_0=int(tcoords[0]),
        target_end_0=int(tcoords[-1]),
    )


def direct_identity(left: str, right: str) -> tuple[int, int, float]:
    if len(left) != len(right):
        raise ValueError("Direct identity requires equal-length sequences")
    matches = sum(a == b for a, b in zip(left, right))
    return matches, len(left) - matches, 100.0 * matches / len(left)


def find_tsd(
    sequence: str,
    region_start: int,
    left_window: tuple[int, int],
    right_window: tuple[int, int],
    min_length: int = 8,
    max_length: int = 30,
) -> tuple[int, int, str]:
    best: tuple[int, int, str] | None = None
    for left_coord in range(left_window[0], left_window[1] + 1):
        left_idx = left_coord - region_start
        for right_coord in range(right_window[0], right_window[1] + 1):
            right_idx = right_coord - region_start
            length = 0
            while (
                length < max_length
                and left_idx + length < len(sequence)
                and right_idx + length < len(sequence)
                and sequence[left_idx + length] == sequence[right_idx + length]
            ):
                length += 1
            if length >= min_length:
                candidate = (left_coord, right_coord, sequence[left_idx : left_idx + length])
                if best is None or len(candidate[2]) > len(best[2]):
                    best = candidate
    if best is None:
        raise RuntimeError("No target-site duplication found")
    return best


def overlaps(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return max(left[0], right[0]) <= min(left[1], right[1])


def write_tsv(path: Path, rows: Iterable[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_bed(path: Path) -> list[dict[str, str]]:
    """Read the six-column BED records used for the archived liftOver audit."""
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, score, strand = line.rstrip().split("\t")
            rows.append(
                {
                    "chrom": chrom,
                    "start_0based": start,
                    "end_0based": end,
                    "feature": name,
                    "score": score,
                    "strand": strand,
                }
            )
    return rows


def summarize_liftover(liftover_dir: Path) -> list[dict[str, object]]:
    """Convert the archived mm39-to-rn7 liftOver files into one audit table."""
    source = read_bed(liftover_dir / "mm39_target_features.bed")
    mapped = {row["feature"]: row for row in read_bed(liftover_dir / "rn7_target_features.bed")}

    reasons: dict[str, str] = {}
    pending_reason = "unmapped"
    with (liftover_dir / "rn7_target_features.unmapped").open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                pending_reason = line.removeprefix("#").strip()
                continue
            fields = line.split("\t")
            reasons[fields[3]] = pending_reason

    rows: list[dict[str, object]] = []
    for item in source:
        target = mapped.get(item["feature"])
        rows.append(
            {
                "feature": item["feature"],
                "mm39_interval_0based_halfopen": (
                    f"{item['chrom']}:{item['start_0based']}-{item['end_0based']}"
                ),
                "status": "mapped" if target else "unmapped",
                "rn7_interval_0based_halfopen": (
                    f"{target['chrom']}:{target['start_0based']}-{target['end_0based']}"
                    if target
                    else ""
                ),
                "rn7_strand": target["strand"] if target else "",
                "liftOver_note": reasons.get(item["feature"], ""),
            }
        )
    return rows


def arrow_gene(ax, start: int, end: int, y: float, direction: int, color: str, label: str) -> None:
    width = end - start
    if direction > 0:
        arrow = FancyArrow(start, y, width, 0, width=0.13, head_width=0.26, head_length=max(width * 0.08, 18),
                           length_includes_head=True, color=color)
    else:
        arrow = FancyArrow(end, y, -width, 0, width=0.13, head_width=0.26, head_length=max(width * 0.08, 18),
                           length_includes_head=True, color=color)
    ax.add_patch(arrow)
    ax.text((start + end) / 2, y + 0.22, label, ha="center", va="bottom", fontsize=9)


def plot_event_reconstruction(output_png: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 8.2), constrained_layout=True)

    ax = axes[0]
    ax.set_title("A. Locus fuente en chrX (GRCm39)", loc="left", fontweight="bold")
    ax.axvspan(COPIED_SOURCE[0], COPIED_SOURCE[1], color="#f4a261", alpha=0.28, label="segmento copiado")
    arrow_gene(ax, 138_963_238, NUP62CL_END, 0.55, -1, "#457b9d", "extremo 5′ de Nup62cl")
    arrow_gene(ax, DNAAF6_START, 138_963_709, 0.55, 1, "#6a4c93", "exón 1 de Dnaaf6")
    ax.add_patch(Rectangle((NUP62CL_END + 1, 0.28), DNAAF6_START - NUP62CL_END - 2, 0.12,
                           color="#e9c46a", alpha=0.9))
    ax.text((NUP62CL_END + DNAAF6_START) / 2, 0.18, "región promotora bidireccional (171 pb)",
            ha="center", va="top", fontsize=8)
    ax.set_xlim(138_963_180, 138_963_760)
    ax.set_ylim(0, 1.15)
    ax.set_yticks([])
    ax.set_xlabel("coordenada chrX")
    ax.ticklabel_format(style="plain", axis="x", useOffset=False)
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    ax.set_title("B. Inserto en chr1: retrogen + promotor divergente", loc="left", fontweight="bold")
    ax.axvspan(INSERTION[0], INSERTION[1], color="#f4a261", alpha=0.20)
    ax.add_patch(Rectangle((LEFT_TSD[0], 0.18), LEFT_TSD[1] - LEFT_TSD[0] + 1, 0.14, color="#2a9d8f"))
    ax.add_patch(Rectangle((RIGHT_TSD[0], 0.18), RIGHT_TSD[1] - RIGHT_TSD[0] + 1, 0.14, color="#2a9d8f"))
    ax.text(LEFT_TSD[0], 0.07, "TSD 14 pb", ha="left", fontsize=8, color="#146c60")
    ax.text(RIGHT_TSD[1], 0.07, "TSD 14 pb", ha="right", fontsize=8, color="#146c60")
    arrow_gene(ax, 31_261_557, 31_261_765, 0.72, -1, "#457b9d", "exón 1 de Alyreffm1")
    arrow_gene(ax, DNAAF6RT[0], DNAAF6RT[1], 0.72, 1, "#6a4c93", "Dnaaf6rt (1 exón)")
    for start, end in ALYREF_TSS:
        ax.axvspan(start, end, ymin=0.49, ymax=0.62, color="#1d3557", alpha=0.8)
    for start, end in DNAAF6RT_TSS:
        ax.axvspan(start, end, ymin=0.49, ymax=0.62, color="#9b2226", alpha=0.8)
    ax.text(31_261_715, 0.48, "TSS Alyreffm1 ←", ha="center", fontsize=8, color="#1d3557")
    ax.text(31_261_885, 0.48, "→ TSS Dnaaf6rt", ha="center", fontsize=8, color="#9b2226")
    ax.set_xlim(31_261_400, 31_263_450)
    ax.set_ylim(0, 1.15)
    ax.set_yticks([])
    ax.set_xlabel("coordenada chr1")
    ax.ticklabel_format(style="plain", axis="x", useOffset=False)

    ax = axes[2]
    ax.set_title("C. Prueba de sintenia mm39 → rn7", loc="left", fontweight="bold")
    ax.axis("off")
    evidence = [
        ("Flancos del inserto", "chr1 ratón → chr9 rata; quedan adyacentes", "conservados"),
        ("Dnaaf6rt + inserto", "sin mapeo ortólogo en chr9 de rata", "no resueltos por liftOver"),
        ("Alyreffm1-202 exones 2–3", "chr1 ratón → chr9 rata", "locus receptor"),
        ("Alyreffm1-202 exón 1", "chr1 ratón → chrX rata", "locus fuente"),
    ]
    xcols = [0.02, 0.31, 0.84]
    headers = ["Elemento", "Resultado liftOver", "Lectura"]
    for x, header in zip(xcols, headers):
        ax.text(x, 0.90, header, transform=ax.transAxes, fontweight="bold", fontsize=9)
    for idx, row in enumerate(evidence):
        y = 0.72 - idx * 0.18
        for x, value in zip(xcols, row):
            ax.text(x, y, value, transform=ax.transAxes, fontsize=9, va="top")
        ax.plot([0.02, 0.98], [y - 0.055, y - 0.055], transform=ax.transAxes, color="#dddddd", lw=0.7)

    fig.suptitle("Reconstrucción actualizada del evento Dnaaf6 → Dnaaf6rt / Alyreffm1", fontsize=14, fontweight="bold")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def plot_transcripts(gene: dict, output_png: Path) -> None:
    transcripts = sorted(gene["Transcript"], key=lambda item: item["display_name"])
    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    insertion_color = "#f4a261"
    ax.axvspan(INSERTION[0], INSERTION[1], color=insertion_color, alpha=0.22,
               label="inserto retrotranspuesto")
    for row, transcript in enumerate(transcripts):
        y = len(transcripts) - row
        exons = sorted((int(e["start"]), int(e["end"])) for e in transcript["Exon"])
        ax.plot([exons[0][0], exons[-1][1]], [y, y], color="#999999", lw=0.8, zorder=1)
        for start, end in exons:
            first_exon = max(exons, key=lambda value: value[1])
            color = insertion_color if (start, end) == first_exon and overlaps((start, end), COPIED_TARGET) else "#457b9d"
            ax.add_patch(Rectangle((start, y - 0.14), max(end - start + 1, 600), 0.28,
                                   color=color, zorder=2))
        translation = transcript.get("Translation")
        if translation:
            ax.plot([translation["start"], translation["end"]], [y, y], color="#111111", lw=4,
                    solid_capstyle="butt", zorder=3)
        label = transcript["display_name"].replace("Alyreffm1-", "-")
        suffix = "  proteína 241 aa" if translation else ""
        ax.text(31_030_500, y, label + suffix, ha="right", va="center", fontsize=8)
    ax.set_xlim(31_030_000, 31_265_000)
    ax.set_ylim(0.3, len(transcripts) + 1.1)
    ax.set_yticks([])
    ax.set_xlabel("chr1, GRCm39")
    ax.ticklabel_format(style="plain", axis="x", useOffset=False)
    ax.set_title("Isoformas actuales de Alyreffm1: el inserto aporta un exón 5′ alternativo",
                 loc="left", fontweight="bold")
    ax.text(31_037_600, 0.55, "barra negra: CDS traducido", fontsize=8, ha="center")
    ax.legend(frameon=False, loc="lower right")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.project_root.resolve()
    data = root / "data" / "event"
    results = root / "results" / "event"
    plots = root / "plots" / "event"
    results.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    reference = data / "reference"
    # Reuse the accessioned mouse locus from the multispecies module instead of
    # storing a second, overlapping copy of the same GRCm39 sequence.
    chr1 = record_with_prefix(
        fasta_records(
            root
            / "data"
            / "synteny"
            / "loci"
            / "mus_musculus.NC_000067.7.31155362-31303806.fasta"
        ),
        "NC_000067.7",
    )
    chrx = record_with_prefix(fasta_records(reference / "mouse_GRCm39_chrX_138955000_139015000.fasta"), "NC_000086.8")

    parent_rna = record_with_prefix(fasta_records(reference / "dnaaf6_parent_rna.fna"), "NM_177921.3")
    parent_cds = record_with_prefix(fasta_records(reference / "dnaaf6_parent_cds.fna"), "NM_177921.3")
    parent_protein = record_with_prefix(fasta_records(reference / "dnaaf6_parent_protein.faa"), "NP_808589.1")
    retro_cds = record_with_prefix(fasta_records(reference / "dnaaf6rt_cds.fna"), "NM_029062.2")
    retro_protein = record_with_prefix(fasta_records(reference / "dnaaf6rt_protein.faa"), "NP_083338.1")

    target_window_start, target_window_end = 31_261_000, 31_263_600
    parent_to_target = summarize_alignment(
        parent_rna,
        genomic_slice(chr1, MM_CHR1_REGION_START, target_window_start, target_window_end),
    )
    source_window_start, source_window_end = 138_962_500, 138_964_000
    promoter_target_start, promoter_target_end = 31_260_500, 31_262_500
    promoter_alignment = summarize_alignment(
        genomic_slice(chrx, MM_CHRX_REGION_START, source_window_start, source_window_end),
        genomic_slice(chr1, MM_CHR1_REGION_START, promoter_target_start, promoter_target_end),
    )
    cds_matches, cds_changes, cds_identity = direct_identity(parent_cds, retro_cds)
    protein_matches, protein_changes, protein_identity = direct_identity(parent_protein, retro_protein)
    tsd_left, tsd_right, tsd_seq = find_tsd(
        chr1,
        MM_CHR1_REGION_START,
        (31_261_400, 31_261_520),
        (31_263_320, 31_263_450),
    )

    metrics = [
        {"analysis": "parent_rna_vs_target_genome", "metric": "aligned_query_nt", "value": parent_to_target.ungapped_bases},
        {"analysis": "parent_rna_vs_target_genome", "metric": "matches", "value": parent_to_target.matches},
        {"analysis": "parent_rna_vs_target_genome", "metric": "identity_ungapped_pct", "value": f"{parent_to_target.identity_ungapped_pct:.2f}"},
        {"analysis": "parent_rna_vs_target_genome", "metric": "local_alignment_score", "value": f"{parent_to_target.score:.1f}"},
        {"analysis": "parent_cds_vs_retrogene_cds", "metric": "length_nt", "value": len(parent_cds)},
        {"analysis": "parent_cds_vs_retrogene_cds", "metric": "matches", "value": cds_matches},
        {"analysis": "parent_cds_vs_retrogene_cds", "metric": "substitutions", "value": cds_changes},
        {"analysis": "parent_cds_vs_retrogene_cds", "metric": "identity_pct", "value": f"{cds_identity:.2f}"},
        {"analysis": "parent_protein_vs_retrogene_protein", "metric": "length_aa", "value": len(parent_protein)},
        {"analysis": "parent_protein_vs_retrogene_protein", "metric": "matches", "value": protein_matches},
        {"analysis": "parent_protein_vs_retrogene_protein", "metric": "differences", "value": protein_changes},
        {"analysis": "parent_protein_vs_retrogene_protein", "metric": "identity_pct", "value": f"{protein_identity:.2f}"},
        {"analysis": "source_promoter_vs_target", "metric": "source_interval", "value": f"chrX:{source_window_start + promoter_alignment.query_start_0}-{source_window_start + promoter_alignment.query_end_0 - 1}"},
        {"analysis": "source_promoter_vs_target", "metric": "target_interval", "value": f"chr1:{promoter_target_start + promoter_alignment.target_start_0}-{promoter_target_start + promoter_alignment.target_end_0 - 1}"},
        {"analysis": "source_promoter_vs_target", "metric": "matches", "value": promoter_alignment.matches},
        {"analysis": "source_promoter_vs_target", "metric": "ungapped_bases", "value": promoter_alignment.ungapped_bases},
        {"analysis": "source_promoter_vs_target", "metric": "identity_ungapped_pct", "value": f"{promoter_alignment.identity_ungapped_pct:.2f}"},
        {"analysis": "target_site_duplication", "metric": "sequence", "value": tsd_seq},
        {"analysis": "target_site_duplication", "metric": "length_bp", "value": len(tsd_seq)},
        {"analysis": "target_site_duplication", "metric": "left_interval", "value": f"chr1:{tsd_left}-{tsd_left + len(tsd_seq) - 1}"},
        {"analysis": "target_site_duplication", "metric": "right_interval", "value": f"chr1:{tsd_right}-{tsd_right + len(tsd_seq) - 1}"},
    ]
    write_tsv(results / "summary_metrics.tsv", metrics, ["analysis", "metric", "value"])

    gene = json.loads((data / "annotation" / "ENSMUSG00000086727.lookup.json").read_text(encoding="utf-8"))
    transcript_rows = []
    for transcript in sorted(gene["Transcript"], key=lambda item: item["display_name"]):
        exons = sorted((int(e["start"]), int(e["end"])) for e in transcript["Exon"])
        first_exon = max(exons, key=lambda value: value[1])
        translation = transcript.get("Translation")
        transcript_rows.append({
            "transcript_id": transcript["id"],
            "display_name": transcript["display_name"],
            "biotype": transcript["biotype"],
            "exon_count": len(exons),
            "five_prime_exon": f"chr1:{first_exon[0]}-{first_exon[1]}(-)",
            "uses_insertion_derived_exon": "yes" if overlaps(first_exon, COPIED_TARGET) else "no",
            "has_translation": "yes" if translation else "no",
            "translation_length_aa": translation.get("length", "") if translation else "",
        })
    write_tsv(
        results / "alyreffm1_transcripts.tsv",
        transcript_rows,
        ["transcript_id", "display_name", "biotype", "exon_count", "five_prime_exon",
         "uses_insertion_derived_exon", "has_translation", "translation_length_aa"],
    )

    liftover_rows = summarize_liftover(data / "liftover")
    write_tsv(
        results / "liftover_summary.tsv",
        liftover_rows,
        [
            "feature",
            "mm39_interval_0based_halfopen",
            "status",
            "rn7_interval_0based_halfopen",
            "rn7_strand",
            "liftOver_note",
        ],
    )

    validation = {
        "target_site_duplication": tsd_seq,
        "target_site_duplication_length_bp": len(tsd_seq),
        "parent_vs_retrogene_cds_identity_pct": round(cds_identity, 2),
        "source_block_identity_ungapped_pct": round(promoter_alignment.identity_ungapped_pct, 2),
        "alyreffm1_transcript_models": len(transcript_rows),
        "models_using_insertion_derived_first_exon": sum(
            row["uses_insertion_derived_exon"] == "yes" for row in transcript_rows
        ),
        "liftOver_mapped_features": sum(row["status"] == "mapped" for row in liftover_rows),
        "liftOver_unmapped_features": sum(row["status"] == "unmapped" for row in liftover_rows),
    }
    (results / "validation_summary.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )

    plot_event_reconstruction(plots / "event_reconstruction_diagnostic.png")
    plot_transcripts(gene, plots / "alyreffm1_transcript_model_audit.png")
    print(f"Analysis complete. Results: {results}")


if __name__ == "__main__":
    main()
