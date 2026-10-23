# Dnaaf6rt–Alyreffm1: compact reproducibility package

Version 1.2.0 (23 October 2026).

This repository contains only the material needed to reproduce or audit the
computational evidence for the focal murine retroposition event. Manuscript
figure sources, PowerPoint files, captions, layout code and duplicated inputs
have been removed.

All biological observations come from accessioned public sequence,
annotation or transcriptomic resources. The workflow does not simulate data,
generate random observations or replace missing inputs with synthetic values.

## Directory map

```text
data/       public inputs retained in the compact archive
scripts/    analyses, download helpers, plotting and validation
results/    tables, sequences and machine-readable summaries
plots/      diagnostic/evidence plots (not manuscript figures)
docs/       methods, provenance, result map and interpretation limits
```

Root files provide the environment (`requirements.txt`, `environment.yml`),
citation metadata and checksums.

## Quick start

Python 3.12 is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
bash scripts/run_all.sh
```

To select another Python interpreter:

```bash
RETROEVENT_PYTHON=/absolute/path/to/python bash scripts/run_all.sh
```

The default run:

1. recomputes the sequence-level event reconstruction;
2. recomputes the nine-species synteny and presence/absence analysis;
3. regenerates the diagnostic plots;
4. audits the compact expression outputs derived from public datasets; and
5. runs explicit scientific, code and packaging checks.

A successful run ends with:

```text
All analyses and validation checks completed successfully.
```

## Optional full expression rerun

The complete FASTQ, PacBio FLNC and Ribo-seq files are too large to duplicate
in this repository. Their accessions, source URLs and observed hashes are in
[`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md). To download them and
recompute the expression tables:

```bash
bash scripts/run_all.sh --full-expression
```

This mode requires network access and substantially more time and disk space.
The download script fails if a required public file cannot be obtained; it
does not fabricate a replacement.

## Evidence boundary

The package directly supports:

- an intronless `Dnaaf6`-derived copy at the recipient locus;
- an exact 14-bp target-site duplication (`AAAAATGATGCTTA`);
- a 5′ insertion block homologous to the divergent `Nup62cl–Dnaaf6` source
  region;
- a phylogenetic insertion bracket of 10.84–11.70 Ma from nine murine loci;
- present-day transcription initiation within the insertion-derived block;
- exact short-read support for the `Alyreffm1-202`-specific splice junction.

It does **not** establish the exact ancestral RNA intermediate, prove that an
intact ancestral promoter was carried in that RNA, date the onset of
opposite-strand transcription, recover a complete `Alyreffm1-202` molecule in
the tested long-read dataset, or demonstrate an isoform-specific protein.

See [`docs/RESULTS_GUIDE.md`](docs/RESULTS_GUIDE.md) for the file supporting
each conclusion and [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md) for the claim
boundaries.

## Coordinate conventions

Narrative genomic coordinates are 1-based and inclusive. BED files are
0-based and half-open. Assemblies and accessions are recorded beside every
input in the provenance documentation.

## Before public release

Select and add a software/data licence appropriate for the project. No licence
has been imposed in this archive because that choice belongs to the authors.
