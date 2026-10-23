#!/usr/bin/env python3
"""Auditable multispecies synteny analysis of the Dnaaf6rt insertion.

All downloaded loci are oriented in the transcriptional direction of the
independently annotated Lgsn ortholog. Coordinates written as genomic are
1-based and inclusive; internal coordinates are 0-based and half-open.

Inputs are under ``data/synteny``. Tables and extracted sequences are written
to ``results/synteny``; two diagnostic plots are written to ``plots/synteny``.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/dnaaf6rt_synteny_matplotlib")

import edlib
import mappy as mp
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "dnaaf6rt-synteny-v1"
import matplotlib.pyplot as plt
from Bio import Align, SeqIO
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "data" / "synteny"
RAW = CONFIG / "loci"
RESULTS = ROOT / "results" / "synteny"
DERIVED = RESULTS / "derived"
PLOTS = ROOT / "plots" / "synteny"
ALIGNMENTS = RESULTS / "alignments"


@dataclass
class MapHit:
    feature: str
    query_length: int
    q_start: int
    q_end: int
    r_start: int
    r_end: int
    strand: int
    matching_bases: int
    alignment_block: int
    mapq: int
    cs: str

    @property
    def coverage(self) -> float:
        return (self.q_end - self.q_start) / self.query_length

    @property
    def identity(self) -> float:
        return self.matching_bases / self.alignment_block if self.alignment_block else 0.0


@dataclass
class PairwiseSummary:
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    query_span: int
    ungapped_bases: int
    matches: int
    identity: float
    score: float
    aligned_query: str
    aligned_target: str

    @property
    def query_coverage(self) -> float:
        return self.query_span


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    fields = columns or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_path(slug: str) -> Path:
    hits = list(RAW.glob(f"{slug}.*.fasta"))
    if len(hits) != 1:
        raise RuntimeError(f"Expected one FASTA for {slug}, found {len(hits)}")
    return hits[0]


def fasta_sequence(path: Path) -> str:
    return str(next(SeqIO.parse(path, "fasta")).seq).upper()


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(textwrap.wrap(sequence, width))


def local_aligner() -> Align.PairwiseAligner:
    obj = Align.PairwiseAligner()
    obj.mode = "local"
    obj.match_score = 2
    obj.mismatch_score = -2
    obj.open_gap_score = -5
    obj.extend_gap_score = -1
    return obj


def global_aligner() -> Align.PairwiseAligner:
    obj = local_aligner()
    obj.mode = "global"
    return obj


def summarize_pairwise(alignment: Align.Alignment, query: str, target: str) -> PairwiseSummary:
    qcoords, tcoords = alignment.coordinates
    matches = 0
    ungapped = 0
    for idx in range(len(qcoords) - 1):
        dq = int(qcoords[idx + 1] - qcoords[idx])
        dt = int(tcoords[idx + 1] - tcoords[idx])
        if dq == dt:
            qblock = query[int(qcoords[idx]) : int(qcoords[idx + 1])]
            tblock = target[int(tcoords[idx]) : int(tcoords[idx + 1])]
            matches += sum(left == right for left, right in zip(qblock, tblock))
            ungapped += dq
    return PairwiseSummary(
        query_start=int(qcoords[0]),
        query_end=int(qcoords[-1]),
        target_start=int(tcoords[0]),
        target_end=int(tcoords[-1]),
        query_span=int(qcoords[-1] - qcoords[0]),
        ungapped_bases=ungapped,
        matches=matches,
        identity=matches / ungapped if ungapped else 0.0,
        score=float(alignment.score),
        aligned_query=alignment[0],
        aligned_target=alignment[1],
    )


def map_query(index: mp.Aligner, feature: str, query: str, expected: tuple[int, int]) -> MapHit | None:
    candidates = []
    for hit in index.map(query, cs=True):
        if hit.strand != 1:
            continue
        midpoint = (hit.r_st + hit.r_en) / 2
        if expected[0] <= midpoint <= expected[1]:
            candidates.append(hit)
    if not candidates:
        return None
    hit = max(candidates, key=lambda item: (item.mlen, item.q_en - item.q_st, item.mapq))
    return MapHit(
        feature=feature,
        query_length=len(query),
        q_start=int(hit.q_st),
        q_end=int(hit.q_en),
        r_start=int(hit.r_st),
        r_end=int(hit.r_en),
        strand=int(hit.strand),
        matching_bases=int(hit.mlen),
        alignment_block=int(hit.blen),
        mapq=int(hit.mapq),
        cs=str(hit.cs or ""),
    )


def oriented_to_genomic(row: dict[str, str], raw_length: int, start: int, end: int) -> tuple[int, int]:
    region_start = int(row["region_start_1based"])
    if row["lgsn_strand"] == "+":
        raw_start, raw_end = start, end
    else:
        raw_start, raw_end = raw_length - end, raw_length - start
    return region_start + raw_start, region_start + raw_end - 1


def project_reference_interval(
    aligned_reference: str, aligned_target: str, start: int, end: int
) -> str:
    projected: list[str] = []
    position = 0
    for ref_base, target_base in zip(aligned_reference, aligned_target):
        if ref_base != "-":
            if start <= position < end:
                projected.append(target_base)
            position += 1
    return "".join(projected).replace("-", "")


def best_target_motif(reference: str, target: str) -> tuple[str, int, float]:
    result = edlib.align(reference, target, mode="HW", task="locations")
    if result["editDistance"] < 0 or not result["locations"]:
        return "", -1, 0.0
    start, end = result["locations"][0]
    motif = target[start : end + 1]
    identity = 1.0 - result["editDistance"] / max(len(reference), len(motif))
    return motif, int(result["editDistance"]), identity


def map_row(
    species: dict[str, str], hit: MapHit | None, raw_length: int
) -> dict[str, object]:
    if hit is None:
        return {
            "species_slug": species["species_slug"],
            "scientific_name": species["scientific_name"],
            "feature": "",
            "query_length": "",
            "query_start_0based": "",
            "query_end_0based": "",
            "oriented_target_start_0based": "",
            "oriented_target_end_0based": "",
            "genomic_start_1based": "",
            "genomic_end_1based": "",
            "coverage": "",
            "identity": "",
            "mapq": "",
            "cs": "",
        }
    genomic_start, genomic_end = oriented_to_genomic(species, raw_length, hit.r_start, hit.r_end)
    return {
        "species_slug": species["species_slug"],
        "scientific_name": species["scientific_name"],
        "feature": hit.feature,
        "query_length": hit.query_length,
        "query_start_0based": hit.q_start,
        "query_end_0based": hit.q_end,
        "oriented_target_start_0based": hit.r_start,
        "oriented_target_end_0based": hit.r_end,
        "genomic_start_1based": min(genomic_start, genomic_end),
        "genomic_end_1based": max(genomic_start, genomic_end),
        "coverage": f"{hit.coverage:.6f}",
        "identity": f"{hit.identity:.6f}",
        "mapq": hit.mapq,
        "cs": hit.cs,
    }


def plot_synteny(rows: list[dict[str, object]], png: Path) -> None:
    ordered = list(reversed(rows))
    fig, ax = plt.subplots(figsize=(11.2, 7.2), constrained_layout=True)
    colors = {"present": "#D97706", "absent": "#64748B", "ambiguous": "#B91C1C"}
    for idx, row in enumerate(ordered):
        y = idx
        gene_start = int(row["lgsn_oriented_start_0based"])
        gene_end = int(row["lgsn_oriented_end_0based"])
        site_start = int(row["site_oriented_start_0based"])
        site_end = int(row["site_oriented_end_0based"])
        gene_len = gene_end - gene_start
        rel_start = site_start - gene_end
        rel_end = site_end - gene_end
        ax.plot([-gene_len, rel_end], [y, y], color="#CBD5E1", lw=1.2, zorder=1)
        ax.add_patch(plt.Rectangle((-gene_len, y - 0.16), gene_len, 0.32, color="#334155", zorder=2))
        state = str(row["state"])
        if state == "present":
            ax.add_patch(
                plt.Rectangle(
                    (rel_start, y - 0.22), max(rel_end - rel_start, 20), 0.44,
                    color=colors[state], zorder=3
                )
            )
        else:
            ax.scatter([(rel_start + rel_end) / 2], [y], marker="|", s=260,
                       linewidths=2.2, color=colors[state], zorder=3)
        state_label = {"present": "presente", "absent": "ausente", "ambiguous": "ambiguo"}[state]
        ax.text(rel_end + 1100, y, f"{row['site_gap_bp']} pb; {state_label}", va="center", fontsize=8.5)
    ax.axvline(0, color="#94A3B8", lw=0.8, ls="--")
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels([str(row["scientific_name"]) for row in ordered], fontstyle="italic")
    ax.set_xlabel("Posición relativa al extremo 3′ de Lgsn en la orientación normalizada (pb)")
    ax.set_title("Sintenia del sitio receptor de Dnaaf6rt", loc="left", fontweight="bold", pad=16)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E2E8F0", lw=0.6)
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=240)
    plt.close(fig)


def plot_phylogeny(png: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw={"width_ratios": [1.3, 1]}, constrained_layout=True)
    ax = axes[0]
    leaves = [
        (0, "Rattini  ·  0/2 presentes", "absent"),
        (1, "Arvicanthini  ·  0/2 presentes", "absent"),
        (2, "Apodemini  ·  1/1 presente", "present"),
        (3, "Praomyini  ·  1/1 presente", "present"),
        (4, "Murini  ·  3/3 presentes", "present"),
    ]
    xpos = {"root": 0.4, "arv": 1.35, "apod": 2.25, "pm": 3.15, "leaf": 4.3}
    # Ladderized topology: (Rattini,(Arvicanthini,(Apodemini,(Praomyini,Murini))))
    nodes = [(xpos["root"], 2.0), (xpos["arv"], 2.5), (xpos["apod"], 3.0), (xpos["pm"], 3.5)]
    ax.plot([xpos["root"], xpos["leaf"]], [0, 0], color="#475569", lw=1.6)
    ax.plot([xpos["root"], xpos["root"]], [0, 2.5], color="#475569", lw=1.6)
    ax.plot([xpos["root"], xpos["arv"]], [2.5, 2.5], color="#475569", lw=1.6)
    ax.plot([xpos["arv"], xpos["leaf"]], [1, 1], color="#475569", lw=1.6)
    ax.plot([xpos["arv"], xpos["arv"]], [1, 3.0], color="#475569", lw=1.6)
    ax.plot([xpos["arv"], xpos["apod"]], [3.0, 3.0], color="#D97706", lw=5.0, solid_capstyle="butt")
    ax.plot([xpos["apod"], xpos["leaf"]], [2, 2], color="#475569", lw=1.6)
    ax.plot([xpos["apod"], xpos["apod"]], [2, 3.5], color="#475569", lw=1.6)
    ax.plot([xpos["apod"], xpos["pm"]], [3.5, 3.5], color="#475569", lw=1.6)
    ax.plot([xpos["pm"], xpos["leaf"]], [3, 3], color="#475569", lw=1.6)
    ax.plot([xpos["pm"], xpos["leaf"]], [4, 4], color="#475569", lw=1.6)
    ax.plot([xpos["pm"], xpos["pm"]], [3, 4], color="#475569", lw=1.6)
    for y, label, state in leaves:
        color = "#D97706" if state == "present" else "#64748B"
        ax.scatter([xpos["leaf"]], [y], s=55, color=color, zorder=4)
        ax.text(xpos["leaf"] + 0.12, y, label, va="center", fontsize=9)
    ax.text((xpos["arv"] + xpos["apod"]) / 2, 3.17, "ganancia inferida",
            color="#B45309", ha="center", fontsize=9, fontweight="bold")
    ax.set_xlim(0.1, 6.5)
    ax.set_ylim(-0.55, 4.55)
    ax.axis("off")
    ax.set_title("A. Reconstrucción de máxima parsimonia", loc="left", fontweight="bold")

    ax = axes[1]
    ax.axvspan(9.62, 12.30, color="#FED7AA", alpha=0.45, label="envolvente de los nodos")
    ax.plot([10.84, 11.70], [0.55, 0.55], color="#D97706", lw=9, solid_capstyle="butt")
    ax.scatter([10.84, 11.70], [0.55, 0.55], color="#7C2D12", s=45, zorder=3)
    ax.text(10.84, 0.72, "10,84 Ma\nApodemini", ha="center", fontsize=9)
    ax.text(11.70, 0.72, "11,70 Ma\nArvicanthini", ha="center", fontsize=9)
    ax.text(11.27, 0.37, "ventana puntual de la inserción", ha="center", color="#92400E", fontsize=9)
    ax.set_xlim(13.0, 9.0)
    ax.set_ylim(0.05, 1.15)
    ax.set_yticks([])
    ax.set_xlabel("Millones de años antes del presente")
    ax.set_title("B. Acotación temporal", loc="left", fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E2E8F0", lw=0.6)
    fig.suptitle("Origen filogenético de la inserción Dnaaf6rt", fontweight="bold", fontsize=14)
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=240)
    plt.close(fig)


def main() -> None:
    for directory in (DERIVED, RESULTS, PLOTS, ALIGNMENTS, RESULTS / "sequences"):
        directory.mkdir(parents=True, exist_ok=True)
    species_rows = read_tsv(CONFIG / "species.tsv")
    feature_rows = {row["feature"]: row for row in read_tsv(CONFIG / "mouse_features.tsv")}
    mouse_row = next(row for row in species_rows if row["species_slug"] == "mus_musculus")
    mouse_raw = fasta_sequence(fasta_path("mus_musculus"))
    mouse_region_start = int(mouse_row["region_start_1based"])

    def mouse_feature(name: str) -> str:
        row = feature_rows[name]
        start = int(row["start_1based"]) - mouse_region_start
        end = int(row["end_1based"]) - mouse_region_start + 1
        return mouse_raw[start:end]

    features = {
        "left_flank_5kb": mouse_feature("left_flank_5kb"),
        "insertion": mouse_feature("insertion"),
        "copied_regulatory_block": mouse_feature("copied_regulatory_block"),
        "dnaaf6rt": mouse_feature("dnaaf6rt"),
        "right_flank_5kb": mouse_feature("right_flank_5kb"),
    }
    left_tsd = mouse_feature("left_tsd")
    right_tsd = mouse_feature("right_tsd")
    if left_tsd != right_tsd:
        raise RuntimeError("Mouse TSD copies are not identical")
    core_site = left_tsd + features["insertion"] + right_tsd
    query_fasta = []
    for name, sequence in features.items():
        query_fasta.append(f">{name}\n{wrap_fasta(sequence)}")
    (DERIVED / "mouse_queries.fasta").write_text("\n".join(query_fasta) + "\n", encoding="utf-8")

    presence_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    tsd_rows: list[dict[str, object]] = []
    site_fastas: list[str] = []
    insertion_fastas: list[str] = []
    empty_fastas: list[str] = []

    for species in species_rows:
        raw_path = fasta_path(species["species_slug"])
        raw = fasta_sequence(raw_path)
        raw_length = len(raw)
        oriented = reverse_complement(raw) if species["lgsn_strand"] == "-" else raw
        raw_gene_start = int(species["lgsn_start_1based"]) - int(species["region_start_1based"])
        raw_gene_end = int(species["lgsn_end_1based"]) - int(species["region_start_1based"]) + 1
        if species["lgsn_strand"] == "+":
            gene_start, gene_end = raw_gene_start, raw_gene_end
        else:
            gene_start, gene_end = raw_length - raw_gene_end, raw_length - raw_gene_start
        expected = (gene_end + 5_000, min(len(oriented), gene_end + 50_000))
        index = mp.Aligner(seq=oriented, k=9, w=5, best_n=20, min_chain_score=20)
        if not index:
            raise RuntimeError(f"Failed to build minimap2 index for {species['scientific_name']}")

        left_candidates = [
            ("left_anchor_1kb", features["left_flank_5kb"][-1000:]),
            ("left_anchor_2kb", features["left_flank_5kb"][-2000:]),
            ("left_flank_5kb", features["left_flank_5kb"]),
        ]
        left_hit = None
        for name, query in left_candidates:
            candidate = map_query(index, name, query, expected)
            if candidate is not None and candidate.coverage >= 0.70 and candidate.q_end >= 0.85 * len(query):
                left_hit = candidate
                break
            if left_hit is None or (candidate is not None and candidate.matching_bases > left_hit.matching_bases):
                left_hit = candidate
        right_candidates = [
            ("right_anchor_1kb", features["right_flank_5kb"][:1000]),
            ("right_anchor_2kb", features["right_flank_5kb"][:2000]),
            ("right_flank_5kb", features["right_flank_5kb"]),
        ]
        right_hit = None
        for name, query in right_candidates:
            candidate = map_query(index, name, query, expected)
            if candidate is not None and candidate.coverage >= 0.70 and candidate.q_start <= 0.15 * len(query):
                right_hit = candidate
                break
            if right_hit is None or (candidate is not None and candidate.matching_bases > right_hit.matching_bases):
                right_hit = candidate
        if left_hit is None or right_hit is None:
            raise RuntimeError(f"Could not recover both syntenic anchors for {species['scientific_name']}")
        collinear = left_hit.r_end < right_hit.r_start
        site_start, site_end = left_hit.r_end, right_hit.r_start
        site = oriented[site_start:site_end] if collinear else ""
        site_n = site.count("N")
        site_gap = len(site)

        insertion_hit = map_query(index, "insertion", features["insertion"], expected)
        retro_hit = map_query(index, "dnaaf6rt", features["dnaaf6rt"], expected)
        for hit in (left_hit, insertion_hit, retro_hit, right_hit):
            if hit is not None:
                mapping_rows.append(map_row(species, hit, raw_length))

        # Refine the complete insertion and its 5-prime block near the syntenic site.
        if retro_hit is not None:
            window_start = max(0, retro_hit.r_start - 900)
            window_end = min(len(oriented), retro_hit.r_end + 200)
        else:
            window_start = max(0, site_start - 600)
            window_end = min(len(oriented), site_end + 600)
        window = oriented[window_start:window_end]
        insertion_alignment = local_aligner().align(features["insertion"], window)[0]
        insertion_summary = summarize_pairwise(insertion_alignment, features["insertion"], window)
        insertion_target_start = window_start + insertion_summary.target_start
        insertion_target_end = window_start + insertion_summary.target_end

        if retro_hit is not None:
            reg_window_start = max(0, retro_hit.r_start - 900)
            reg_window_end = min(len(oriented), retro_hit.r_start + 200)
        else:
            reg_window_start = max(0, site_start - 600)
            reg_window_end = min(len(oriented), site_end + 600)
        reg_window = oriented[reg_window_start:reg_window_end]
        regulatory_alignment = local_aligner().align(features["copied_regulatory_block"], reg_window)[0]
        regulatory_summary = summarize_pairwise(
            regulatory_alignment, features["copied_regulatory_block"], reg_window
        )
        regulatory_coverage = regulatory_summary.query_span / len(features["copied_regulatory_block"])

        insertion_coverage = insertion_summary.query_span / len(features["insertion"])
        retro_coverage = retro_hit.coverage if retro_hit else 0.0
        retro_identity = retro_hit.identity if retro_hit else 0.0
        present = (
            collinear
            and insertion_coverage >= 0.85
            and insertion_summary.identity >= 0.75
            and retro_coverage >= 0.85
            and retro_identity >= 0.75
            and site_n == 0
        )
        absent = (
            collinear
            and site_gap <= 250
            and retro_hit is None
            and (insertion_hit is None or insertion_hit.coverage < 0.25)
            and site_n == 0
        )
        state = "present" if present else "absent" if absent else "ambiguous"
        left_resolved = left_hit.coverage >= 0.70 and left_hit.q_end >= 0.85 * left_hit.query_length
        right_resolved = right_hit.coverage >= 0.70 and right_hit.q_start <= 0.15 * right_hit.query_length
        boundary_status = "resolved" if left_resolved and right_resolved else "partial"
        site_gstart, site_gend = oriented_to_genomic(species, raw_length, site_start, site_end)
        insert_gstart, insert_gend = oriented_to_genomic(
            species, raw_length, insertion_target_start, insertion_target_end
        )

        row = {
            "species_slug": species["species_slug"],
            "scientific_name": species["scientific_name"],
            "tribe": species["tribe"],
            "assembly_accession": species["assembly_accession"],
            "sequence_accession": species["sequence_accession"],
            "state": state,
            "boundary_status": boundary_status,
            "site_gap_bp": site_gap,
            "site_n_bases": site_n,
            "site_oriented_start_0based": site_start,
            "site_oriented_end_0based": site_end,
            "site_genomic_start_1based": min(site_gstart, site_gend),
            "site_genomic_end_1based": max(site_gstart, site_gend),
            "insertion_coverage": f"{insertion_coverage:.6f}" if insertion_hit else "",
            "insertion_identity": f"{insertion_summary.identity:.6f}" if insertion_hit else "",
            "insertion_genomic_start_1based": min(insert_gstart, insert_gend) if present else "",
            "insertion_genomic_end_1based": max(insert_gstart, insert_gend) if present else "",
            "dnaaf6rt_coverage": f"{retro_coverage:.6f}" if retro_hit else "",
            "dnaaf6rt_identity": f"{retro_identity:.6f}" if retro_hit else "",
            "regulatory_coverage": f"{regulatory_coverage:.6f}",
            "regulatory_identity": f"{regulatory_summary.identity:.6f}",
            "lgsn_oriented_start_0based": gene_start,
            "lgsn_oriented_end_0based": gene_end,
            "distance_lgsn_to_site_bp": site_start - gene_end,
        }
        presence_rows.append(row)

        qc_rows.append(
            {
                "species_slug": species["species_slug"],
                "scientific_name": species["scientific_name"],
                "fasta_file": str(raw_path.relative_to(ROOT)),
                "fasta_sha256": sha256(raw_path),
                "locus_length_bp": raw_length,
                "locus_n_bases": raw.count("N"),
                "locus_n_fraction": f"{raw.count('N') / raw_length:.8f}",
                "lgsn_strand_in_assembly": species["lgsn_strand"],
                "left_anchor_feature": left_hit.feature,
                "left_anchor_coverage": f"{left_hit.coverage:.6f}",
                "left_anchor_identity": f"{left_hit.identity:.6f}",
                "right_anchor_feature": right_hit.feature,
                "right_anchor_coverage": f"{right_hit.coverage:.6f}",
                "right_anchor_identity": f"{right_hit.identity:.6f}",
                "anchors_collinear": str(collinear).lower(),
                "boundary_status": boundary_status,
                "site_n_bases": site_n,
            }
        )

        site_fastas.append(
            f">{species['species_slug']}|{species['sequence_accession']}|site_gap={site_gap}|state={state}\n"
            f"{wrap_fasta(site)}"
        )
        if present:
            insertion_homolog = oriented[insertion_target_start:insertion_target_end]
            insertion_fastas.append(
                f">{species['species_slug']}|{species['sequence_accession']}|"
                f"oriented={insertion_target_start}-{insertion_target_end}|"
                f"genomic={min(insert_gstart, insert_gend)}-{max(insert_gstart, insert_gend)}\n"
                f"{wrap_fasta(insertion_homolog)}"
            )
            alignment_path = ALIGNMENTS / f"{species['species_slug']}.insertion_vs_mouse.fasta"
            alignment_path.write_text(
                f">mus_musculus_insertion|aligned_query={insertion_summary.query_start}-{insertion_summary.query_end}\n"
                f"{wrap_fasta(insertion_summary.aligned_query)}\n"
                f">{species['species_slug']}|aligned_target={insertion_target_start}-{insertion_target_end}\n"
                f"{wrap_fasta(insertion_summary.aligned_target)}\n",
                encoding="utf-8",
            )
        else:
            empty_fastas.append(
                f">{species['species_slug']}|{species['sequence_accession']}|"
                f"genomic={min(site_gstart, site_gend)}-{max(site_gstart, site_gend)}\n"
                f"{wrap_fasta(site)}"
            )

        if present:
            core_alignment = global_aligner().align(core_site, site)[0]
            projected_left = project_reference_interval(core_alignment[0], core_alignment[1], 0, len(left_tsd))
            projected_right = project_reference_interval(
                core_alignment[0], core_alignment[1], len(core_site) - len(right_tsd), len(core_site)
            )
            if projected_left and projected_right:
                edit = edlib.align(projected_left, projected_right, mode="NW")["editDistance"]
                tsd_identity = 1.0 - edit / max(len(projected_left), len(projected_right))
                display = f"{100 * tsd_identity:.1f}%"
            else:
                tsd_identity = ""
                display = "límite parcial"
            tsd_rows.append(
                {
                    "species_slug": species["species_slug"],
                    "scientific_name": species["scientific_name"],
                    "state": state,
                    "left_tsd_or_empty_site": projected_left or "no resuelta",
                    "right_tsd_or_best_motif": projected_right or "no resuelta",
                    "tsd_identity": f"{tsd_identity:.6f}" if tsd_identity != "" else "",
                    "tsd_identity_display": display,
                    "empty_site_length": "",
                    "edit_distance_to_mouse_tsd": "",
                }
            )
        else:
            motif, edit, identity = best_target_motif(left_tsd, site)
            tsd_rows.append(
                {
                    "species_slug": species["species_slug"],
                    "scientific_name": species["scientific_name"],
                    "state": state,
                    "left_tsd_or_empty_site": site,
                    "right_tsd_or_best_motif": motif,
                    "tsd_identity": f"{identity:.6f}",
                    "tsd_identity_display": f"{100 * identity:.1f}% vs. TSD de ratón",
                    "empty_site_length": len(site),
                    "edit_distance_to_mouse_tsd": edit,
                }
            )

    write_tsv(RESULTS / "presence_absence.tsv", presence_rows)
    write_tsv(RESULTS / "locus_qc.tsv", qc_rows)
    write_tsv(RESULTS / "feature_mappings.tsv", mapping_rows)
    write_tsv(RESULTS / "tsd_and_empty_sites.tsv", tsd_rows)
    (RESULTS / "sequences" / "orthologous_site_intervals.fasta").write_text(
        "\n".join(site_fastas) + "\n", encoding="utf-8"
    )
    (RESULTS / "sequences" / "insertion_homologs.fasta").write_text(
        "\n".join(insertion_fastas) + "\n", encoding="utf-8"
    )
    (RESULTS / "sequences" / "ancestral_empty_sites.fasta").write_text(
        "\n".join(empty_fastas) + "\n", encoding="utf-8"
    )
    dating = [
        {
            "inference": "maximum_parsimony_gain",
            "older_boundary_node": "Arvicanthini_vs_Apodemini_Praomyini_Murini",
            "older_boundary_ma": "11.70",
            "younger_boundary_node": "Apodemini_vs_Praomyini_Murini",
            "younger_boundary_ma": "10.84",
            "point_estimate_wording": "approximately 11 Ma",
            "conservative_outer_envelope_ma": "9.62-12.30",
            "citation": "Thomas et al. 2025; doi:10.1093/gbe/evaf017",
        }
    ]
    write_tsv(RESULTS / "dating_summary.tsv", dating)
    plot_synteny(presence_rows, PLOTS / "orthologous_site_length_and_state.png")
    plot_phylogeny(PLOTS / "presence_absence_phylogeny.png")

    validation = {
        "species": len(presence_rows),
        "present": sum(row["state"] == "present" for row in presence_rows),
        "absent": sum(row["state"] == "absent" for row in presence_rows),
        "ambiguous": sum(row["state"] == "ambiguous" for row in presence_rows),
        "all_sites_without_N": all(int(row["site_n_bases"]) == 0 for row in presence_rows),
        "all_anchors_collinear": all(row["anchors_collinear"] == "true" for row in qc_rows),
        "mouse_tsd": left_tsd,
        "thresholds": {
            "presence_min_insertion_coverage": 0.85,
            "presence_min_identity": 0.75,
            "absence_max_site_gap_bp": 250,
        },
    }
    (RESULTS / "validation_summary.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    manifest = {
        "analysis_date": "2026-08-29",
        "python": platform.python_version(),
        "packages": {
            "biopython": importlib.metadata.version("biopython"),
            "edlib": importlib.metadata.version("edlib"),
            "mappy": importlib.metadata.version("mappy"),
            "matplotlib": importlib.metadata.version("matplotlib"),
        },
        "script_sha256": sha256(Path(__file__)),
        "config_sha256": {
            "species.tsv": sha256(CONFIG / "species.tsv"),
            "mouse_features.tsv": sha256(CONFIG / "mouse_features.tsv"),
        },
        "coordinate_conventions": {
            "internal": "0-based half-open",
            "reported_genomic": "1-based inclusive",
        },
        "validation": validation,
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if validation["ambiguous"] != 0 or not validation["all_sites_without_N"] or not validation["all_anchors_collinear"]:
        raise RuntimeError(f"Validation failed: {validation}")
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
