# Validation scope

`scripts/05_validate_package.py` performs fail-fast checks of the central
reported values, code syntax and package contents. It verifies:

- the exact 14-bp TSD;
- parent/retrogene and source-block sequence measurements;
- the current transcript-model counts;
- the retained mouse→rat liftOver pattern;
- nine species with 5 present, 4 absent and 0 ambiguous calls;
- the 11.70–10.84-Ma phylogenetic bracket;
- the exact short-read and long-read support counts;
- the 127-aa ORF boundary;
- the presence of all six diagnostic PNG plots;
- parseable Python source and portable paths; and
- absence of random-data imports, simulation/mock-data calls, PowerPoint,
  SVG/AI/EPS artwork and manuscript-figure source directories.

This audit has two reproducibility levels:

| Module | Default run |
|---|---|
| focal sequence reconstruction | fully recomputed from included FASTA/JSON/BED inputs |
| multispecies synteny | fully recomputed from the nine included loci |
| expression | compact public-data derivatives and exact positive-read tables are validated; plots are regenerated |
| full expression | optional network workflow redownloads and recomputes the large public inputs |

Passing validation means that the repository reproduces the stated numerical
observations within this scope. It is not a substitute for independent
biological replication and does not remove the interpretation limits listed in
`LIMITATIONS.md`.

## Release verification record

The compact workflow was executed end to end on 23 October 2026 with Python
3.12.14, Biopython 1.88, edlib 1.3.9.post1, mappy 2.30, Matplotlib 3.10.5 and
SciPy 1.17.0. All 15 checks passed. The principal structural, synteny and
expression tables were compared with the previously audited v1.1.0 outputs;
their scientific fields were unchanged (the only expected table difference was
the updated relative path of the reorganized FASTA files). All six PNG files
also passed chunk-boundary and CRC validation.

The optional large-file expression download/reanalysis was not repeated during
this packaging step. Its existing public-data derivatives were audited, and the
full retrieval/reanalysis path remains executable and checksum guarded.
