#!/usr/bin/env python3
"""Fail-fast scientific and packaging checks for the compact repository.

The checks are intentionally explicit.  They verify that the recomputed tables
contain the reported observations, that all included Python source parses, and
that no manuscript-figure assets or random/simulated-data code entered the
minimal archive.
"""

from __future__ import annotations

import ast
import csv
import json
import struct
import zlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metric(rows: list[dict[str, str]], analysis: str, name: str) -> str:
    return next(row["value"] for row in rows if row["analysis"] == analysis and row["metric"] == name)


def record(checks: list[dict[str, str]], name: str, condition: bool, detail: str) -> None:
    if not condition:
        raise AssertionError(f"{name}: {detail}")
    checks.append({"check": name, "status": "PASS", "detail": detail})


def inspect_python_sources() -> tuple[int, list[str]]:
    """Parse source and flag imports/calls that could generate fake observations."""
    files = sorted((ROOT / "scripts").glob("*.py"))
    findings: list[str] = []
    for path in files:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == "random" or alias.name.startswith("numpy.random"):
                        findings.append(f"{path.name}: import {alias.name}")
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if module == "random" or module.startswith("numpy.random"):
                    findings.append(f"{path.name}: from {module} import ...")
            elif isinstance(node, ast.Call):
                name = ""
                if isinstance(node.func, ast.Name):
                    name = node.func.id
                elif isinstance(node.func, ast.Attribute):
                    name = node.func.attr
                if name.lower() in {"simulate", "simulation", "make_mock_data", "make_dummy_data"}:
                    findings.append(f"{path.name}: call {name}()")
    return len(files), findings


def valid_png(path: Path) -> bool:
    """Validate PNG chunk boundaries and CRC values without extra packages."""
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        return False
    position = 8
    while position < len(data):
        if position + 12 > len(data):
            return False
        length = struct.unpack(">I", data[position : position + 4])[0]
        chunk_type = data[position + 4 : position + 8]
        end = position + 12 + length
        if end > len(data):
            return False
        payload = data[position + 8 : position + 8 + length]
        observed_crc = struct.unpack(">I", data[position + 8 + length : end])[0]
        expected_crc = zlib.crc32(chunk_type + payload) & 0xFFFFFFFF
        if observed_crc != expected_crc:
            return False
        position = end
        if chunk_type == b"IEND":
            return position == len(data)
    return False


def main() -> None:
    checks: list[dict[str, str]] = []

    # Structural reconstruction.
    event_metrics = read_tsv(RESULTS / "event" / "summary_metrics.tsv")
    record(
        checks,
        "exact_target_site_duplication",
        metric(event_metrics, "target_site_duplication", "sequence") == "AAAAATGATGCTTA"
        and metric(event_metrics, "target_site_duplication", "length_bp") == "14",
        "14-bp sequence AAAAATGATGCTTA recovered from GRCm39",
    )
    record(
        checks,
        "processed_copy_identity",
        metric(event_metrics, "parent_cds_vs_retrogene_cds", "length_nt") == "657"
        and metric(event_metrics, "parent_cds_vs_retrogene_cds", "identity_pct") == "93.15",
        "657-nt CDS comparison; 93.15% identity",
    )
    record(
        checks,
        "source_homologous_block",
        metric(event_metrics, "source_promoter_vs_target", "ungapped_bases") == "452"
        and metric(event_metrics, "source_promoter_vs_target", "identity_ungapped_pct") == "80.31",
        "452 aligned ungapped bases; 80.31% identity",
    )
    transcripts = read_tsv(RESULTS / "event" / "alyreffm1_transcripts.tsv")
    record(
        checks,
        "transcript_models",
        len(transcripts) == 10 and sum(row["uses_insertion_derived_exon"] == "yes" for row in transcripts) == 6,
        "10 current Ensembl models; 6 use the insertion-derived first exon",
    )
    lift = {row["feature"]: row for row in read_tsv(RESULTS / "event" / "liftover_summary.tsv")}
    record(
        checks,
        "mouse_rat_liftover",
        lift["candidate_insertion"]["status"] == "unmapped"
        and lift["Dnaaf6rt"]["status"] == "unmapped"
        and lift["Alyreffm1_202_exon2"]["rn7_interval_0based_halfopen"].startswith("chr9:")
        and lift["Alyreffm1_202_exon3"]["rn7_interval_0based_halfopen"].startswith("chr9:")
        and lift["Alyreffm1_202_exon1"]["rn7_interval_0based_halfopen"].startswith("chrX:"),
        "recipient exons map to rat chr9; the insertion lacks an orthologous chr9 liftOver mapping; first exon maps to rat chrX",
    )

    # Multispecies bracket.
    synteny = json.loads((RESULTS / "synteny" / "validation_summary.json").read_text(encoding="utf-8"))
    record(
        checks,
        "nine_species_state_calls",
        (synteny["species"], synteny["present"], synteny["absent"], synteny["ambiguous"]) == (9, 5, 4, 0),
        "9 loci: 5 present, 4 absent, 0 ambiguous",
    )
    dating = read_tsv(RESULTS / "synteny" / "dating_summary.tsv")[0]
    record(
        checks,
        "phylogenetic_age_bracket",
        dating["older_boundary_ma"] == "11.70" and dating["younger_boundary_ma"] == "10.84",
        "physical insertion bracketed between 11.70 and 10.84 Ma",
    )

    # Expression evidence.  This validates archived public-data derivatives;
    # the optional full workflow regenerates them from downloaded source files.
    expression = json.loads((RESULTS / "expression" / "validation_summary.json").read_text(encoding="utf-8"))
    record(
        checks,
        "isoform_specific_short_reads",
        expression["exact_202_junction_reads_by_replicate"] == {"ENCFF001JPX": 16, "ENCFF001JPZ": 20},
        "16 and 20 exact -202-specific junction reads in two biological replicates",
    )
    record(
        checks,
        "long_read_limit",
        expression["pacbio_flnc_reads_screened"] == 129471
        and expression["pacbio_reads_with_shared_promoter_junction"] == 3
        and expression["pacbio_reads_with_unique_202_junction"] == 0,
        "129,471 FLNC reads screened; 3 shared-junction reads; no full -202 junction chain",
    )
    record(
        checks,
        "candidate_orf_boundary",
        expression["candidate_orf_length_aa"] == 127,
        "127-aa ORF reconstructed; translation remains unproven",
    )

    # Repository scope and code checks.
    forbidden_suffixes = {".pptx", ".svg", ".ai", ".eps"}
    forbidden_parts = {"05_editable_figures", "drafts", "powerpoint", "source_data"}
    forbidden_files = [
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
        and (path.suffix.lower() in forbidden_suffixes or any(part in forbidden_parts for part in path.parts))
    ]
    record(checks, "no_manuscript_figure_assets", not forbidden_files, "no PPTX/SVG/AI/EPS or figure-source directories")

    python_count, synthetic_findings = inspect_python_sources()
    record(checks, "python_source_parses", python_count >= 6, f"{python_count} Python scripts parsed with ast")
    record(checks, "no_simulated_data_code", not synthetic_findings, "no random-data imports or simulation/mock-data calls")

    absolute_paths = []
    for path in (ROOT / "scripts").iterdir():
        if path.is_file() and path.suffix in {".py", ".sh"} and path.name != Path(__file__).name:
            text = path.read_text(encoding="utf-8")
            if "/workspace/scratch/" in text or "/home/rafik/" in text:
                absolute_paths.append(path.name)
    record(checks, "portable_paths", not absolute_paths, "no workspace- or user-specific absolute paths")

    required_plots = [
        PLOTS / "event" / "event_reconstruction_diagnostic.png",
        PLOTS / "event" / "alyreffm1_transcript_model_audit.png",
        PLOTS / "synteny" / "orthologous_site_length_and_state.png",
        PLOTS / "synteny" / "presence_absence_phylogeny.png",
        PLOTS / "expression" / "expression_evidence_summary.png",
        PLOTS / "expression" / "alyreffm1_202_evidence_map.png",
    ]
    record(
        checks,
        "diagnostic_plots",
        all(path.is_file() and path.stat().st_size > 10_000 and valid_png(path) for path in required_plots),
        "six evidence/QC plots regenerated as complete, CRC-valid PNG files",
    )

    output = RESULTS / "validation_report.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)
    print(f"PASS: {len(checks)} checks. Report: {output}")


if __name__ == "__main__":
    main()
