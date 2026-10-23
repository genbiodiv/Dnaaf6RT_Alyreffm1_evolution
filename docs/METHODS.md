# Computational methods

## 1. Focal event reconstruction

The GRCm39 recipient interval (`NC_000067.7`) and source interval
(`NC_000086.8`) are read from the included FASTA files. `Dnaaf6` RNA/CDS and
protein records are compared with the recipient interval and the annotated
`Dnaaf6rt` products. Local pairwise alignment uses Biopython with match `+2`,
mismatch `−2`, gap open `−5` and gap extension `−1`. CDS and protein identity
are also counted directly for equal-length records.

The target-site duplication is found by an exact sequence search across fixed
left and right windows at the GRCm39 insertion boundary, considering matches
of 8–30 bp and reporting the longest match. The source-homologous 5′ block is
measured by local alignment of the fixed chromosome-X and chromosome-1
windows.

The current Ensembl lookup for `ENSMUSG00000086727` is parsed to classify the
first exon and translation status of each transcript. A first exon is called
insertion-derived when it overlaps the defined source-homologous target block.

The mm39→rn7 results were produced with UCSC `liftOver` and the
`mm39ToRn7.over.chain.gz` chain. The compact workflow parses the archived BED
output and unmapped records into one explicit table. Because the chain file
and binary are external UCSC resources, this step is an audit of their retained
output rather than a bundled reimplementation of liftOver.

## 2. Multispecies synteny

Nine accession-versioned murine intervals are anchored on independently
annotated `Lgsn` orthologs. Each interval spans 60 kb on both sides of the
pinned `Lgsn` coordinates and is oriented in the `Lgsn` transcriptional
direction. The species, assemblies, sequence accessions and coordinates are in
`data/synteny/species.tsv`.

Mouse left/right flanks, the complete insertion, the 5′ source-homologous
block and `Dnaaf6rt` are mapped to each locus with minimap2 through `mappy`.
The complete insertion and 5′ block are then refined by local pairwise
alignment.

An insertion is called present when:

- both boundary anchors are collinear;
- insertion coverage is at least 0.85;
- insertion ungapped identity is at least 0.75;
- `Dnaaf6rt` coverage and identity are each at least 0.85 and 0.75,
  respectively; and
- the orthologous site contains no ambiguous bases.

Absence requires collinear resolved boundaries, an interval of at most 250 bp,
no qualifying `Dnaaf6rt` hit, no insertion hit covering 25% or more of the
query, and no ambiguous bases.

The mouse interval containing the left TSD, insertion and right TSD is globally
aligned to each orthologous site. Presence/absence states are placed on the
accepted murine species tree under a single-gain maximum-parsimony model. The
11.70- and 10.84-Ma node ages are taken from Thomas et al. (2025,
<https://doi.org/10.1093/gbe/evaf017>); the resulting interval is a
phylogenetic bracket, not a locus-specific molecular-clock estimate.

## 3. Expression evidence

The expression analysis separates five questions that should not be collapsed
into one claim:

1. **Promoter-region activity.** Five focal FANTOM5 CAGE peaks are extracted
   from the phase 1+2 mouse matrix and summed by direction across the testis
   developmental series.
2. **Isoform-specific RNA evidence.** Two adult-testis ENCODE FASTQ files are
   scanned for exact 24-nt splice-junction k-mers (12 nt from each exon), in
   both orientations. Longer 30–60-nt junction motifs are checked within the
   positive reads, and the motifs are checked against the 140-kb mm10 locus to
   exclude contiguous genomic matches.
3. **Probabilistic transcript quantification.** ENCODE RSEM transcript tables
   are filtered for current `Alyreffm1` transcript identifiers; this evidence
   is kept separate from the exact junction scan.
4. **Full-length support.** PacBio FLNC reads from GSE137490 are remapped to the
   mm10 locus with splice-aware minimap2. Reads are classified by exact intron
   boundaries for the shared first-exon junction and the `-202`-specific
   downstream junction.
5. **Coding compatibility.** Open reading frames are enumerated from the
   reconstructed annotated RNA. Ribo-seq P-site distributions are summarized
   by transcript region and frame. This is treated as compatibility with
   translation, not proof of a stable protein.

The compact default run validates the retained result tables and exact-hit
records. `scripts/run_all.sh --full-expression` downloads the accessioned
large inputs and executes the complete workflow.

## 4. Diagnostic plots

Plots are generated directly from analysis results and serve only as visual
quality-control summaries. No plotting script constructs manuscript panels or
contains values absent from the tabular outputs.
