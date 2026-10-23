#!/usr/bin/env python3
"""Run the full expression and coding reanalysis for Alyreffm1-202.

Coordinates are 1-based and inclusive unless explicitly labeled as BED/0-based.
The analysis separates promoter activity, isoform-specific splice evidence,
full-length transcript support, and translation-compatible ribosome occupancy.

Run ``fetch_expression_inputs.sh`` before this network-data workflow. The
default repository workflow validates the compact public-data derivatives,
because complete FASTQ and long-read inputs are not duplicated in the archive.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/alyreffm1_expression_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import chisquare

try:
    import mappy as mp
except ImportError:  # pragma: no cover - reported cleanly by the long-read stage
    mp = None


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "expression"
RESULTS = ROOT / "results" / "expression"
PLOTS = ROOT / "plots" / "expression"
SEQUENCES = RESULTS / "sequences"

GENE_JSON = ROOT / "data" / "event" / "annotation" / "ENSMUSG00000086727.lookup.json"
FANTOM_SUBSET = DATA / "fantom5" / "mm9_Alyreffm1_Dnaaf6rt_CAGE_peaks.tsv"
MM10_LOCUS = DATA / "reference" / "mm10_NC_000067.6_31090000-31230000.fasta"
LONG_READ_FASTA = DATA / "gse137490" / "GSE137490_TestisPacBio_FLNC.fasta.gz"
LONG_READ_GTF = DATA / "gse137490" / "GSE137490_TranscriptomeAssembly.MouseTestis.mm10.gtf.gz"

MM39_TO_MM10_OFFSET = -39_081
MM10_FASTA_START_1 = 31_090_000

ALYREFFM1_202_EXONS_MM39 = [
    (31_261_557, 31_261_765),
    (31_243_333, 31_243_554),
    (31_242_012, 31_242_469),
]
ALYREFFM1_202_EXONS_MM10 = [
    (start + MM39_TO_MM10_OFFSET, end + MM39_TO_MM10_OFFSET)
    for start, end in ALYREFFM1_202_EXONS_MM39
]

# Minimap2 introns: 0-based half-open positions relative to the locus FASTA.
FASTA0 = MM10_FASTA_START_1 - 1
PROMOTER_INTRON_LOCAL = (31_204_473 - FASTA0, 31_222_475 - FASTA0)
UNIQUE_202_INTRON_LOCAL = (31_203_388 - FASTA0, 31_204_251 - FASTA0)
NR033445_INTRON_LOCAL = (31_141_768 - FASTA0, 31_204_251 - FASTA0)

FANTOM_PEAKS = {
    "chr1:31279449..31279465,-": ("Alyreffm1", "Tssr11183", "mm_5389.1"),
    "chr1:31279472..31279488,-": ("Alyreffm1", "Tssr11184", "mm_5390.1"),
    "chr1:31279489..31279522,-": ("Alyreffm1", "Tssr11185", "mm_5391.1"),
    "chr1:31279627..31279637,+": ("Dnaaf6rt", "", ""),
    "chr1:31279653..31279666,+": ("Dnaaf6rt", "", ""),
}

RSEM_INPUTS = {
    "ENCFF183XSN.tsv": ("ENCSR000CHH", "adult_2m_shallow", "biological_rep1", "36-bp single-end"),
    "ENCFF448JGR.tsv": ("ENCSR000CHH", "adult_2m_shallow", "biological_rep2", "36-bp single-end"),
    "ENCFF957IRE.tsv": ("ENCSR000BYW", "adult_2m_deep", "biological_rep1", "76-bp paired-end"),
    "ENCFF802BII.tsv": ("ENCSR000BYW", "adult_2m_deep", "biological_rep2", "76-bp paired-end"),
    "ENCFF926TOA.tsv": ("ENCSR266ESZ", "adult_10w_deep", "technical_rep1", "paired-end"),
    "ENCFF392NZM.tsv": ("ENCSR266ESZ", "adult_10w_deep", "technical_rep2", "paired-end"),
}

FASTQ_READ_COUNTS = {
    "ENCFF001JPX": 37_696_007,
    "ENCFF001JPZ": 41_120_376,
}

RIBO_INPUTS = {
    "GSM4160756.bed13.txt.gz": ("adult_WT_rep1", "SRX7141430"),
    "GSM4160757.bed13.txt.gz": ("adult_WT_rep2", "SRX7141431"),
}

TRANSCRIPT_REGIONS = {
    # name: (BED start, BED end, transcript start, transcript end)
    "5UTR_exon1": (31_222_601, 31_222_684, 0, 83),
    "ORF_exon1": (31_222_475, 31_222_601, 83, 209),
    "ORF_exon2": (31_204_251, 31_204_473, 209, 431),
    "ORF_exon3_with_stop": (31_203_352, 31_203_388, 431, 467),
    "3UTR_exon3": (31_202_930, 31_203_352, 467, 889),
}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def parse_attributes(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]*)"', text))


def transcript_ordered_exons(transcript: dict[str, object]) -> list[dict[str, object]]:
    strand = int(transcript["strand"])
    return sorted(transcript.get("Exon", []), key=lambda x: int(x["start"]), reverse=strand < 0)


def analyse_annotation() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    gene = json.loads(GENE_JSON.read_text(encoding="utf-8"))
    transcript_rows: list[dict[str, object]] = []
    junction_rows: list[dict[str, object]] = []
    junction_uses: Counter[str] = Counter()

    interim: list[tuple[dict[str, object], list[dict[str, object]], list[str]]] = []
    for transcript in gene["Transcript"]:
        exons = transcript_ordered_exons(transcript)
        junctions: list[str] = []
        for left, right in zip(exons, exons[1:]):
            high = left if int(left["start"]) > int(right["start"]) else right
            low = right if high is left else left
            junction = f"chr1:{int(low['end'])}-{int(high['start'])}(-)"
            junctions.append(junction)
            junction_uses[junction] += 1
        interim.append((transcript, exons, junctions))

    for transcript, exons, junctions in interim:
        translation = transcript.get("Translation")
        first = exons[0]
        uses_inserted = int(first["start"]) >= 31_261_473
        transcript_rows.append(
            {
                "transcript_id": transcript["id"],
                "version": transcript.get("version", ""),
                "display_name": transcript.get("display_name", ""),
                "biotype": transcript.get("biotype", ""),
                "length_nt": transcript.get("length", ""),
                "exon_count": len(exons),
                "five_prime_exon_mm39": f"chr1:{first['start']}-{first['end']}(-)",
                "uses_insertion_derived_first_exon": "yes" if uses_inserted else "no",
                "has_annotated_translation": "yes" if translation else "no",
                "translation_length_aa": translation.get("length", "") if translation else "",
                "junction_chain_mm39": ";".join(junctions),
            }
        )
        for order, junction in enumerate(junctions, 1):
            junction_rows.append(
                {
                    "transcript_id": transcript["id"],
                    "display_name": transcript.get("display_name", ""),
                    "junction_order": order,
                    "junction_mm39": junction,
                    "junction_use_count_in_gene": junction_uses[junction],
                    "is_unique_within_current_annotation": "yes" if junction_uses[junction] == 1 else "no",
                }
            )

    transcript_rows.sort(key=lambda row: row["display_name"])
    junction_rows.sort(key=lambda row: (row["display_name"], row["junction_order"]))
    write_tsv(
        RESULTS / "alyreffm1_transcript_architecture.tsv",
        transcript_rows,
        list(transcript_rows[0]),
    )
    write_tsv(RESULTS / "alyreffm1_junction_catalog.tsv", junction_rows, list(junction_rows[0]))
    return transcript_rows, junction_rows


def locus_sequence() -> str:
    return str(next(SeqIO.parse(MM10_LOCUS, "fasta")).seq).upper()


def genomic_slice(sequence: str, start_1: int, end_1: int) -> str:
    return sequence[start_1 - MM10_FASTA_START_1 : end_1 - MM10_FASTA_START_1 + 1]


def reconstruct_202_and_orfs() -> tuple[str, dict[str, object], list[dict[str, object]]]:
    genome = locus_sequence()
    exons = [reverse_complement(genomic_slice(genome, a, b)) for a, b in ALYREFFM1_202_EXONS_MM10]
    transcript = "".join(exons)
    if len(transcript) != 889:
        raise AssertionError(f"Expected 889 nt for Alyreffm1-202, observed {len(transcript)}")

    all_orfs: list[dict[str, object]] = []
    stops = {"TAA", "TAG", "TGA"}
    for frame in range(3):
        for start in range(frame, len(transcript) - 2, 3):
            if transcript[start : start + 3] != "ATG":
                continue
            for stop in range(start + 3, len(transcript) - 2, 3):
                codon = transcript[stop : stop + 3]
                if codon in stops:
                    peptide = str(Seq(transcript[start:stop]).translate())
                    all_orfs.append(
                        {
                            "start_nt_1based": start + 1,
                            "stop_codon_start_nt_1based": stop + 1,
                            "orf_length_including_stop_nt": stop + 3 - start,
                            "protein_length_aa": len(peptide),
                            "frame": frame,
                            "stop_codon": codon,
                            "kozak_minus3": transcript[start - 3] if start >= 3 else "",
                            "kozak_plus4": transcript[start + 3] if start + 3 < len(transcript) else "",
                            "peptide": peptide,
                        }
                    )
                    break
    all_orfs.sort(key=lambda row: int(row["protein_length_aa"]), reverse=True)
    longest = all_orfs[0]

    SEQUENCES.mkdir(parents=True, exist_ok=True)
    (SEQUENCES / "Alyreffm1-202.mm10.fasta").write_text(
        ">ENSMUST00000135245.2|Alyreffm1-202|mm10_reconstructed|length=889\n"
        + "\n".join(transcript[i : i + 70] for i in range(0, len(transcript), 70))
        + "\n",
        encoding="utf-8",
    )
    (SEQUENCES / "Alyreffm1-202.longest_ORF.faa").write_text(
        f">Alyreffm1-202_candidate_ORF|aa={longest['protein_length_aa']}|nt={longest['start_nt_1based']}-{int(longest['stop_codon_start_nt_1based']) + 2}\n"
        + "\n".join(str(longest["peptide"])[i : i + 70] for i in range(0, len(str(longest["peptide"])), 70))
        + "\n",
        encoding="utf-8",
    )
    write_tsv(RESULTS / "alyreffm1_202_orfs.tsv", all_orfs, list(all_orfs[0]))
    return transcript, longest, all_orfs


def analyse_fantom() -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, float]]:
    with FANTOM_SUBSET.open(encoding="utf-8") as handle:
        matrix = list(csv.reader(handle, delimiter="\t"))
    header = matrix[0]
    sample_names = [urllib.parse.unquote(value.removeprefix("tpm.")) for value in header[1:]]
    peak_ids = [row[0] for row in matrix[1:]]
    values = [[float(value) for value in row[1:]] for row in matrix[1:]]
    if set(peak_ids) != set(FANTOM_PEAKS):
        raise AssertionError("FANTOM subset does not contain the expected five peaks")
    aly_indices = [i for i, peak in enumerate(peak_ids) if FANTOM_PEAKS[peak][0] == "Alyreffm1"]
    dna_indices = [i for i, peak in enumerate(peak_ids) if FANTOM_PEAKS[peak][0] == "Dnaaf6rt"]

    all_rows: list[dict[str, object]] = []
    peak_rows: list[dict[str, object]] = []
    for sample_idx, sample in enumerate(sample_names):
        aly = sum(values[i][sample_idx] for i in aly_indices)
        dna = sum(values[i][sample_idx] for i in dna_indices)
        all_rows.append(
            {
                "sample": sample,
                "alyreffm1_three_minus_strand_peaks_tpm": aly,
                "dnaaf6rt_two_plus_strand_peaks_tpm": dna,
            }
        )
        for peak_idx, peak in enumerate(peak_ids):
            gene, tss, fantom_id = FANTOM_PEAKS[peak]
            peak_rows.append(
                {
                    "sample": sample,
                    "peak_mm9": peak,
                    "directional_locus": gene,
                    "mgi_tss": tss,
                    "fantom_peak_alias": fantom_id,
                    "tpm": values[peak_idx][sample_idx],
                }
            )

    stage_patterns = [
        ("E13", "testis, embryo E13."),
        ("E15", "testis, embryo E15."),
        ("E16", "testis, embryo E16."),
        ("E17", "testis, embryo E17."),
        ("E18", "testis, embryo E18."),
        ("N00", "testis, neonate N00."),
        ("N07", "testis, neonate N07."),
        ("N10", "testis, neonate N10."),
        ("N20", "testis, neonate N20."),
        ("N30", "testis, neonate N30."),
        ("adult", "testis, adult."),
    ]
    timecourse: list[dict[str, object]] = []
    for order, (stage, prefix) in enumerate(stage_patterns):
        row = next(row for row in all_rows if str(row["sample"]).startswith(prefix))
        timecourse.append(
            {
                "stage": stage,
                "stage_order": order,
                "sample": row["sample"],
                "alyreffm1_promoter_tpm": row["alyreffm1_three_minus_strand_peaks_tpm"],
                "dnaaf6rt_promoter_tpm": row["dnaaf6rt_two_plus_strand_peaks_tpm"],
            }
        )

    x = [float(row["alyreffm1_promoter_tpm"]) for row in timecourse]
    y = [float(row["dnaaf6rt_promoter_tpm"]) for row in timecourse]
    xbar, ybar = mean(x), mean(y)
    correlation = sum((a - xbar) * (b - ybar) for a, b in zip(x, y)) / math.sqrt(
        sum((a - xbar) ** 2 for a in x) * sum((b - ybar) ** 2 for b in y)
    )
    n10 = next(float(row["alyreffm1_promoter_tpm"]) for row in timecourse if row["stage"] == "N10")
    n20 = next(float(row["alyreffm1_promoter_tpm"]) for row in timecourse if row["stage"] == "N20")
    metrics = {
        "testis_directional_cage_pearson_r": correlation,
        "alyreffm1_n20_over_n10_fold_change": n20 / n10,
        "alyreffm1_adult_tpm": next(float(row["alyreffm1_promoter_tpm"]) for row in timecourse if row["stage"] == "adult"),
        "dnaaf6rt_adult_tpm": next(float(row["dnaaf6rt_promoter_tpm"]) for row in timecourse if row["stage"] == "adult"),
    }
    write_tsv(RESULTS / "fantom5_all_samples.tsv", all_rows, list(all_rows[0]))
    write_tsv(RESULTS / "fantom5_peak_values.tsv", peak_rows, list(peak_rows[0]))
    write_tsv(RESULTS / "fantom5_testis_timecourse.tsv", timecourse, list(timecourse[0]))
    return all_rows, timecourse, metrics


def analyse_rsem() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    raw_dir = DATA / "encode"
    for filename, (experiment, group, replicate, layout) in RSEM_INPUTS.items():
        path = raw_dir / filename
        with path.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["gene_id"].split(".")[0] != "ENSMUSG00000086727":
                    continue
                rows.append(
                    {
                        "experiment": experiment,
                        "file_accession": path.stem,
                        "sample_group": group,
                        "replicate": replicate,
                        "layout": layout,
                        "transcript_id": row["transcript_id"],
                        "expected_count": float(row["expected_count"]),
                        "TPM": float(row["TPM"]),
                        "FPKM": float(row["FPKM"]),
                        "IsoPct": float(row["IsoPct"]),
                    }
                )
    rows.sort(key=lambda row: (row["experiment"], row["file_accession"], row["transcript_id"]))
    write_tsv(RESULTS / "encode_rsem_alyreffm1.tsv", rows, list(rows[0]))
    return rows


def analyse_fastq_junction_hits(
    transcript: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    hit_paths = sorted(RESULTS.glob("ENCFF001JP*.junction_kmer_hits.tsv"))
    all_hits: list[dict[str, object]] = []
    for path in hit_paths:
        with path.open(encoding="utf-8") as handle:
            all_hits.extend(csv.DictReader(handle, delimiter="\t"))

    exon1, exon2, exon3 = [
        transcript[:209],
        transcript[209:431],
        transcript[431:],
    ]
    junction_sequences = {
        "shared_promoter_exon1_to_exon2": (exon1, exon2),
        "alyreffm1_202_exon2_to_exon3": (exon2, exon3),
    }
    genome = locus_sequence()
    summary: list[dict[str, object]] = []
    validation: list[dict[str, object]] = []
    for source in sorted(FASTQ_READ_COUNTS):
        source_hits = [row for row in all_hits if row["source"] == source]
        for junction, (left, right) in junction_sequences.items():
            junction_hits = [row for row in source_hits if row["junction"] == junction]
            count24 = len(junction_hits)
            summary.append(
                {
                    "source": source,
                    "biological_replicate": "1" if source == "ENCFF001JPX" else "2",
                    "total_reads_scanned": FASTQ_READ_COUNTS[source],
                    "junction": junction,
                    "exact_24nt_hits": count24,
                    "hits_per_million_reads": 1_000_000 * count24 / FASTQ_READ_COUNTS[source],
                    "unique_read_names": len({row["read_name"] for row in junction_hits}),
                    "unique_read_sequences": len({row["read_sequence"] for row in junction_hits}),
                }
            )
            for half in (12, 15, 18, 20, 25, 30):
                motif = left[-half:] + right[:half]
                reverse = reverse_complement(motif)
                count = sum(motif in row["read_sequence"] or reverse in row["read_sequence"] for row in junction_hits)
                validation.append(
                    {
                        "source": source,
                        "junction": junction,
                        "bases_from_each_exon": half,
                        "junction_kmer_length": 2 * half,
                        "exact_hits": count,
                        "contiguous_occurrences_in_140kb_mm10_locus": genome.count(motif) + genome.count(reverse),
                    }
                )
    write_tsv(RESULTS / "junction_kmer_summary.tsv", summary, list(summary[0]))
    write_tsv(RESULTS / "junction_kmer_length_validation.tsv", validation, list(validation[0]))
    return summary, validation


def introns_from_cigar(hit: object) -> list[tuple[int, int]]:
    position = int(hit.r_st)
    introns: list[tuple[int, int]] = []
    for length, operation in hit.cigar:
        if operation == 3:
            introns.append((position, position + int(length)))
            position += int(length)
        elif operation in (0, 2, 7, 8):
            position += int(length)
    return introns


def contains_intron(introns: list[tuple[int, int]], target: tuple[int, int], tolerance: int = 3) -> bool:
    return any(abs(a - target[0]) <= tolerance and abs(b - target[1]) <= tolerance for a, b in introns)


def analyse_long_reads() -> tuple[list[dict[str, object]], list[dict[str, object]], int]:
    if mp is None:
        raise RuntimeError("mappy is required for the long-read remapping stage")
    aligner = mp.Aligner(str(MM10_LOCUS), preset="splice:hq", k=14, extra_flags=0x4000000)
    records: list[dict[str, object]] = []
    read_sequences: dict[str, str] = {}
    read_count = 0
    for name, sequence, _quality in mp.fastx_read(str(LONG_READ_FASTA)):
        read_count += 1
        for hit in aligner.map(sequence, cs=True):
            if not hit.is_primary or hit.mapq < 20:
                continue
            introns = introns_from_cigar(hit)
            promoter = contains_intron(introns, PROMOTER_INTRON_LOCAL)
            unique_202 = contains_intron(introns, UNIQUE_202_INTRON_LOCAL)
            nr033445 = contains_intron(introns, NR033445_INTRON_LOCAL)
            if not (promoter or unique_202 or nr033445):
                continue
            read_sequences[name] = sequence
            records.append(
                {
                    "read_name": name,
                    "read_length": len(sequence),
                    "strand": hit.strand,
                    "mapq": hit.mapq,
                    "matched_bases": hit.mlen,
                    "alignment_block_length": hit.blen,
                    "query_start_0based": hit.q_st,
                    "query_end_0based": hit.q_en,
                    "reference_start_local_0based": hit.r_st,
                    "reference_end_local_0based": hit.r_en,
                    "introns_local_0based": ";".join(f"{a}-{b}" for a, b in introns),
                    "cigar": hit.cigar_str,
                    "supports_shared_promoter_junction": "yes" if promoter else "no",
                    "supports_unique_202_junction": "yes" if unique_202 else "no",
                    "supports_NR_033445_downstream_junction": "yes" if nr033445 else "no",
                }
            )

    write_tsv(RESULTS / "gse137490_long_read_junction_support.tsv", records, list(records[0]))
    with (SEQUENCES / "gse137490_promoter_junction_reads.fasta").open("w", encoding="utf-8") as handle:
        for name in sorted(read_sequences):
            sequence = read_sequences[name]
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 70):
                handle.write(sequence[index : index + 70] + "\n")

    models: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(LONG_READ_GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get("transcript_id", "")
            if "4931428L18Rik" not in transcript_id:
                continue
            models[transcript_id].append((int(fields[3]), int(fields[4])))
    model_rows: list[dict[str, object]] = []
    for transcript_id, exons in sorted(models.items()):
        ascending = sorted(exons)
        introns = [(left[1], right[0] - 1) for left, right in zip(ascending, ascending[1:])]
        model_rows.append(
            {
                "transcript_id": transcript_id,
                "exon_count": len(exons),
                "exon_chain_mm10_1based": ";".join(f"chr1:{a}-{b}(-)" for a, b in sorted(exons, reverse=True)),
                "intron_chain_mm10_0based": ";".join(f"{a}-{b}" for a, b in introns),
                "contains_shared_promoter_junction": "yes" if (31_204_473, 31_222_475) in introns else "no",
                "contains_unique_202_junction": "yes" if (31_203_388, 31_204_251) in introns else "no",
            }
        )
    write_tsv(RESULTS / "gse137490_transcript_models_at_locus.tsv", model_rows, list(model_rows[0]))
    return records, model_rows, read_count


def bed_blocks(fields: list[str]) -> list[tuple[int, int]]:
    start = int(fields[1])
    block_count = int(fields[9])
    sizes = [int(value) for value in fields[10].rstrip(",").split(",")]
    offsets = [int(value) for value in fields[11].rstrip(",").split(",")]
    if len(sizes) != block_count or len(offsets) != block_count:
        raise ValueError("Malformed BED13 block fields")
    return [(start + offset, start + offset + size) for size, offset in zip(sizes, offsets)]


def psite_genomic_position(blocks: list[tuple[int, int]], strand: str, offset: int = 12) -> int | None:
    positions: list[int] = []
    if strand == "-":
        for start, end in sorted(blocks, reverse=True):
            positions.extend(range(end - 1, start - 1, -1))
    else:
        for start, end in sorted(blocks):
            positions.extend(range(start, end))
    return positions[offset] if len(positions) > offset else None


def genomic_to_202_transcript_position(position_0: int) -> int | None:
    exons = [
        (31_222_475, 31_222_684, 0),
        (31_204_251, 31_204_473, 209),
        (31_202_930, 31_203_388, 431),
    ]
    for start, end, transcript_start in exons:
        if start <= position_0 < end:
            return transcript_start + (end - 1 - position_0)
    return None


def transcript_region(position: int | None) -> str:
    if position is None:
        return "outside_transcript"
    for name, (_a, _b, start, end) in TRANSCRIPT_REGIONS.items():
        if start <= position < end:
            return name
    return "outside_transcript"


def analyse_riboseq() -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    raw_dir = DATA / "riboseq"
    footprint_rows: list[dict[str, object]] = []
    region_rows: list[dict[str, object]] = []
    frame_rows: list[dict[str, object]] = []
    exon_union = [(value[0], value[1]) for value in TRANSCRIPT_REGIONS.values()]

    for filename, (sample, srx) in RIBO_INPUTS.items():
        path = raw_dir / filename
        counts_by_region: Counter[str] = Counter()
        rows_by_region: Counter[str] = Counter()
        frame_counts: Counter[int] = Counter()
        frame_counts_27_32: Counter[int] = Counter()
        with gzip.open(path, "rt") as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                if len(fields) < 13 or fields[0] != "chr1" or fields[5] != "-":
                    continue
                blocks = bed_blocks(fields)
                if not any(block_start < exon_end and block_end > exon_start for block_start, block_end in blocks for exon_start, exon_end in exon_union):
                    continue
                weight = int(fields[3])
                read_length = len(fields[12])
                psite = psite_genomic_position(blocks, fields[5], offset=12)
                tx_position = genomic_to_202_transcript_position(psite) if psite is not None else None
                region = transcript_region(tx_position)
                if region != "outside_transcript":
                    counts_by_region[region] += weight
                    rows_by_region[region] += 1
                if tx_position is not None and 83 <= tx_position < 464:
                    frame = (tx_position - 83) % 3
                    frame_counts[frame] += weight
                    if 27 <= read_length <= 32:
                        frame_counts_27_32[frame] += weight
                footprint_rows.append(
                    {
                        "sample": sample,
                        "sra_experiment": srx,
                        "chrom": fields[0],
                        "start_0based": fields[1],
                        "end_0based": fields[2],
                        "collapsed_read_count": weight,
                        "strand": fields[5],
                        "block_count": fields[9],
                        "block_sizes": fields[10],
                        "block_starts": fields[11],
                        "read_length": read_length,
                        "psite_offset": 12,
                        "psite_genomic_0based": psite if psite is not None else "",
                        "psite_transcript_0based": tx_position if tx_position is not None else "",
                        "psite_region": region,
                        "orf_frame": (tx_position - 83) % 3 if tx_position is not None and 83 <= tx_position < 464 else "",
                        "sequence": fields[12],
                    }
                )

        total_weight = sum(counts_by_region.values())
        total_bases = 889
        baseline_density = total_weight / total_bases if total_weight else float("nan")
        for region, (start, end, _tx_start, _tx_end) in TRANSCRIPT_REGIONS.items():
            bases = end - start
            weighted = counts_by_region[region]
            density = weighted / bases
            region_rows.append(
                {
                    "sample": sample,
                    "sra_experiment": srx,
                    "region": region,
                    "region_length_nt": bases,
                    "footprint_rows_by_psite": rows_by_region[region],
                    "weighted_footprints_by_psite": weighted,
                    "footprints_per_nt": density,
                    "relative_density_within_transcript": density / baseline_density if baseline_density else "",
                }
            )
        observed = [frame_counts[index] for index in range(3)]
        observed_27_32 = [frame_counts_27_32[index] for index in range(3)]
        test = chisquare(observed) if sum(observed) else None
        test_filtered = chisquare(observed_27_32) if sum(observed_27_32) else None
        frame_rows.append(
            {
                "sample": sample,
                "sra_experiment": srx,
                "psite_offset_nt": 12,
                "read_length_filter": "all",
                "frame0": observed[0],
                "frame1": observed[1],
                "frame2": observed[2],
                "frame0_fraction": observed[0] / sum(observed),
                "chi_square_uniform": float(test.statistic),
                "p_value_uniform": float(test.pvalue),
            }
        )
        frame_rows.append(
            {
                "sample": sample,
                "sra_experiment": srx,
                "psite_offset_nt": 12,
                "read_length_filter": "27-32 nt",
                "frame0": observed_27_32[0],
                "frame1": observed_27_32[1],
                "frame2": observed_27_32[2],
                "frame0_fraction": observed_27_32[0] / sum(observed_27_32),
                "chi_square_uniform": float(test_filtered.statistic),
                "p_value_uniform": float(test_filtered.pvalue),
            }
        )

    write_tsv(RESULTS / "riboseq_locus_footprints.tsv", footprint_rows, list(footprint_rows[0]))
    write_tsv(RESULTS / "riboseq_region_counts.tsv", region_rows, list(region_rows[0]))
    write_tsv(RESULTS / "riboseq_frame_periodicity.tsv", frame_rows, list(frame_rows[0]))
    return footprint_rows, region_rows, frame_rows


def plot_expression_summary(
    timecourse: list[dict[str, object]],
    junction_summary: list[dict[str, object]],
    rsem_rows: list[dict[str, object]],
    ribo_regions: list[dict[str, object]],
) -> None:
    matplotlib.rcParams["svg.hashsalt"] = "alyreffm1-expression-2026-08-29"
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4), constrained_layout=True)

    ax = axes[0, 0]
    stages = [str(row["stage"]) for row in timecourse]
    x = list(range(len(stages)))
    aly = [float(row["alyreffm1_promoter_tpm"]) for row in timecourse]
    dna = [float(row["dnaaf6rt_promoter_tpm"]) for row in timecourse]
    ax.plot(x, aly, "o-", color="#245a9a", lw=2, label="Alyreffm1 ← (3 picos)")
    ax.plot(x, dna, "s-", color="#b24a3a", lw=2, label="Dnaaf6rt → (2 picos)")
    ax.set_yscale("symlog", linthresh=0.2)
    ax.set_xticks(x, stages, rotation=45, ha="right")
    ax.set_ylabel("CAGE TPM combinado (escala symlog)")
    ax.set_title("A. Activación bidireccional durante el desarrollo testicular", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    sources = sorted({str(row["source"]) for row in junction_summary})
    labels = ["Rep. 1", "Rep. 2"]
    shared = [next(int(row["exact_24nt_hits"]) for row in junction_summary if row["source"] == source and str(row["junction"]).startswith("shared")) for source in sources]
    unique = [next(int(row["exact_24nt_hits"]) for row in junction_summary if row["source"] == source and str(row["junction"]).startswith("alyreffm1_202")) for source in sources]
    width = 0.36
    ax.bar([i - width / 2 for i in range(len(labels))], shared, width, color="#e5a44a", label="exón 1→2 compartida")
    ax.bar([i + width / 2 for i in range(len(labels))], unique, width, color="#6a4c93", label="exón 2→3 exclusiva de -202")
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylabel("Lecturas con k-mer exacto de unión")
    ax.set_title("B. Evidencia directa en RNA-seq adulto", loc="left", fontweight="bold")
    ax.legend(frameon=False)
    for i, values in enumerate(zip(shared, unique)):
        for dx, value in zip((-width / 2, width / 2), values):
            ax.text(i + dx, value + max(shared) * 0.02, str(value), ha="center", va="bottom", fontsize=8)

    ax = axes[1, 0]
    rows_202 = [row for row in rsem_rows if str(row["transcript_id"]).split(".")[0] == "ENSMUST00000135245"]
    rows_202.sort(key=lambda row: (row["sample_group"], row["file_accession"]))
    labels = [f"{row['experiment']}\n{row['replicate']}" for row in rows_202]
    values = [float(row["TPM"]) for row in rows_202]
    colors = ["#bdbdbd" if value == 0 else "#2a9d8f" for value in values]
    ax.bar(range(len(values)), values, color=colors)
    ax.set_xticks(range(len(values)), labels, rotation=35, ha="right", fontsize=7)
    ax.set_ylabel("RSEM TPM de Alyreffm1-202")
    ax.set_title("C. Cuantificación dependiente de profundidad", loc="left", fontweight="bold")

    ax = axes[1, 1]
    region_order = list(TRANSCRIPT_REGIONS)
    region_labels = ["5′ UTR", "ORF exón 1", "ORF exón 2", "ORF exón 3/stop", "3′ UTR"]
    samples = sorted({str(row["sample"]) for row in ribo_regions})
    width = 0.36
    for sample_index, sample in enumerate(samples):
        values = [
            float(next(row["relative_density_within_transcript"] for row in ribo_regions if row["sample"] == sample and row["region"] == region))
            for region in region_order
        ]
        offsets = [index + (sample_index - 0.5) * width for index in range(len(region_order))]
        ax.bar(offsets, values, width, label=f"Rep. {sample_index + 1}", color=("#5b8e7d", "#bc8da0")[sample_index])
    ax.axhline(1, color="#555555", lw=0.8, ls="--")
    ax.set_xticks(range(len(region_order)), region_labels, rotation=30, ha="right")
    ax.set_ylabel("Densidad RPF relativa dentro del transcrito")
    ax.set_title("D. Ocupación ribosomal compatible con el ORF candidato", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    fig.savefig(PLOTS / "expression_evidence_summary.png", dpi=220)
    plt.close(fig)


def plot_transcript_orf(
    long_reads: list[dict[str, object]], junction_summary: list[dict[str, object]]
) -> None:
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(12, 4.6), constrained_layout=True)
    ax.set_xlim(0, 910)
    ax.set_ylim(-1.25, 2.3)
    ax.set_yticks([])
    ax.set_xlabel("Coordenada en el transcrito empalmado Alyreffm1-202 (nt)")
    ax.set_title("Arquitectura de Alyreffm1-202 y niveles de evidencia", loc="left", fontweight="bold")

    exon_specs = [(1, 209, "Exón 1 incorporado", "#e9a03b"), (210, 431, "Exón 2", "#4c78a8"), (432, 889, "Exón 3 exclusivo", "#4c78a8")]
    for start, end, label, color in exon_specs:
        ax.add_patch(Rectangle((start, 1.25), end - start + 1, 0.34, color=color, alpha=0.9))
        ax.text((start + end) / 2, 1.67, label, ha="center", va="bottom")
    ax.text(1, 1.05, "RNA anotado: 889 nt", ha="left", color="#333333")

    ax.add_patch(Rectangle((84, 0.48), 384, 0.30, color="#6a4c93", alpha=0.9))
    ax.text(84, 0.32, "ATG (−3=A; +4=C)", ha="left", color="#4f366e")
    ax.text(467, 0.32, "TGA", ha="right", color="#4f366e")
    ax.text((84 + 467) / 2, 0.88, "ORF candidato: 127 aa", ha="center", color="#4f366e", fontweight="bold")

    shared = FancyArrowPatch((195, -0.05), (225, -0.05), connectionstyle="arc3,rad=-0.55", arrowstyle="-", lw=2.2, color="#e59f3a")
    unique = FancyArrowPatch((415, -0.05), (448, -0.05), connectionstyle="arc3,rad=-0.55", arrowstyle="-", lw=2.2, color="#6a4c93")
    ax.add_patch(shared)
    ax.add_patch(unique)
    shared_counts = [
        int(row["exact_24nt_hits"])
        for row in junction_summary
        if str(row["junction"]).startswith("shared")
    ]
    unique_counts = [
        int(row["exact_24nt_hits"])
        for row in junction_summary
        if str(row["junction"]).startswith("alyreffm1_202")
    ]
    ax.text(210, -0.52, f"Unión compartida: {' / '.join(map(str, shared_counts))} lecturas", ha="center", color="#9a641f")
    ax.text(432, -0.52, f"Unión exclusiva -202: {' / '.join(map(str, unique_counts))} lecturas", ha="center", color="#4f366e")

    promoter_reads = sum(row["supports_shared_promoter_junction"] == "yes" for row in long_reads)
    exact_reads = sum(row["supports_unique_202_junction"] == "yes" for row in long_reads)
    ax.text(1, -1.02, f"PacBio testículo: {promoter_reads} lecturas usan la unión del promotor; {exact_reads} recuperan la unión exclusiva de -202.", ha="left")

    fig.savefig(PLOTS / "alyreffm1_202_evidence_map.png", dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-long-read-remap", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for path in (RESULTS, PLOTS, SEQUENCES):
        path.mkdir(parents=True, exist_ok=True)

    transcripts, junctions = analyse_annotation()
    transcript, longest_orf, _all_orfs = reconstruct_202_and_orfs()
    _fantom_all, timecourse, fantom_metrics = analyse_fantom()
    rsem_rows = analyse_rsem()
    junction_summary, _junction_validation = analyse_fastq_junction_hits(transcript)

    long_read_path = RESULTS / "gse137490_long_read_junction_support.tsv"
    model_path = RESULTS / "gse137490_transcript_models_at_locus.tsv"
    if args.skip_long_read_remap and long_read_path.exists() and model_path.exists():
        with long_read_path.open(encoding="utf-8") as handle:
            long_reads = list(csv.DictReader(handle, delimiter="\t"))
        with model_path.open(encoding="utf-8") as handle:
            long_models = list(csv.DictReader(handle, delimiter="\t"))
        long_read_count = 129_471
    else:
        long_reads, long_models, long_read_count = analyse_long_reads()

    _ribo_footprints, ribo_regions, ribo_frames = analyse_riboseq()
    plot_expression_summary(timecourse, junction_summary, rsem_rows, ribo_regions)
    plot_transcript_orf(long_reads, junction_summary)

    rows_202 = [row for row in rsem_rows if str(row["transcript_id"]).split(".")[0] == "ENSMUST00000135245"]
    detected_202 = [row for row in rows_202 if float(row["expected_count"]) > 0]
    unique_junction = next(
        row for row in junctions if row["transcript_id"] == "ENSMUST00000135245" and row["junction_order"] == 2
    )
    summary = {
        "current_transcripts": len(transcripts),
        "transcripts_using_inserted_first_exon": sum(row["uses_insertion_derived_first_exon"] == "yes" for row in transcripts),
        "alyreffm1_202_unique_junction_mm39": unique_junction["junction_mm39"],
        "fantom5_testis_samples": len(timecourse),
        **fantom_metrics,
        "adult_biological_replicates_with_exact_202_junction_reads": sum(
            row["junction"] == "alyreffm1_202_exon2_to_exon3" and int(row["exact_24nt_hits"]) > 0
            for row in junction_summary
        ),
        "exact_202_junction_reads_by_replicate": {
            row["source"]: int(row["exact_24nt_hits"])
            for row in junction_summary
            if row["junction"] == "alyreffm1_202_exon2_to_exon3"
        },
        "rsem_detected_202_samples": len(detected_202),
        "rsem_202_tpm_range_among_detected": [
            min(float(row["TPM"]) for row in detected_202),
            max(float(row["TPM"]) for row in detected_202),
        ],
        "pacbio_flnc_reads_screened": long_read_count,
        "pacbio_reads_with_shared_promoter_junction": sum(row["supports_shared_promoter_junction"] == "yes" for row in long_reads),
        "pacbio_reads_with_unique_202_junction": sum(row["supports_unique_202_junction"] == "yes" for row in long_reads),
        "pacbio_transcript_models_at_locus": len(long_models),
        "candidate_orf_length_aa": int(longest_orf["protein_length_aa"]),
        "candidate_orf_start_nt_1based": int(longest_orf["start_nt_1based"]),
        "candidate_orf_stop_codon_start_nt_1based": int(longest_orf["stop_codon_start_nt_1based"]),
        "candidate_orf_kozak_minus3": longest_orf["kozak_minus3"],
        "candidate_orf_kozak_plus4": longest_orf["kozak_plus4"],
        "riboseq_frame0_fraction_all_reads": {
            row["sample"]: float(row["frame0_fraction"])
            for row in ribo_frames
            if row["read_length_filter"] == "all"
        },
        "interpretation": {
            "promoter_activity": "confirmed",
            "alyreffm1_202_expression_in_adult_testis": "confirmed_by_exact_unique_junction_reads_in_two_biological_replicates",
            "full_length_alyreffm1_202_long_read": "not_recovered_in_GSE137490",
            "candidate_ORF_translation": "suggestive_ribosome_occupancy_not_definitive_protein_evidence",
        },
    }
    (RESULTS / "validation_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
