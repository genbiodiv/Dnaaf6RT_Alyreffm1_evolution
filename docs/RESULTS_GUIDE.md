# Result and evidence map

| Question | Primary result | Script | Diagnostic plot | Supported reading |
|---|---|---|---|---|
| Is the copy RNA-mediated? | `results/event/summary_metrics.tsv` | `01_reconstruct_event.py` | `plots/event/event_reconstruction_diagnostic.png` | processed `Dnaaf6`-derived copy with high CDS identity |
| Is there a target-site duplication? | `results/event/summary_metrics.tsv` | `01_reconstruct_event.py` | same plot | exact 14-bp duplication `AAAAATGATGCTTA` |
| Is the 5′ block source-derived? | `results/event/summary_metrics.tsv` | `01_reconstruct_event.py` | same plot | 452 ungapped aligned bases at 80.31% identity to the chromosome-X source region |
| What maps to rat? | `results/event/liftover_summary.tsv` | `01_reconstruct_event.py` | same plot | recipient exons map to rat chr9; the inserted block lacks an orthologous chr9 liftOver mapping; the focal first exon maps to rat chrX |
| Which current isoforms use the block? | `results/event/alyreffm1_transcripts.tsv` | `01_reconstruct_event.py` | `plots/event/alyreffm1_transcript_model_audit.png` | 6 of 10 current models use an insertion-derived first exon |
| Which lineages carry the insertion? | `results/synteny/presence_absence.tsv` | `02_analyse_synteny.py` | both synteny plots | present in five loci and absent in four, with no ambiguous calls |
| When did the DNA insertion occur? | `results/synteny/dating_summary.tsv` | `02_analyse_synteny.py` | `plots/synteny/presence_absence_phylogeny.png` | physical insertion bracketed between 11.70 and 10.84 Ma |
| Is the focal block transcriptionally active? | `results/expression/fantom5_testis_timecourse.tsv` | full expression workflow | `plots/expression/expression_evidence_summary.png` | bidirectional CAGE activity rises after postnatal day 10 |
| Is `Alyreffm1-202` detected specifically? | `results/expression/junction_kmer_summary.tsv` and two exact-hit tables | `scan_junction_kmers.py` | same plot | 16 and 20 exact isoform-specific junction reads in two biological replicates |
| Was a complete focal isoform recovered by long reads? | `results/expression/gse137490_long_read_junction_support.tsv` | full expression workflow | `plots/expression/alyreffm1_202_evidence_map.png` | no; three reads support only the shared first-exon junction |
| Is a protein demonstrated? | ORF and Ribo-seq result tables | full expression workflow | both expression plots | no; a 127-aa ORF and compatible occupancy are present, but protein evidence is not definitive |

The final executable audit is `results/validation_report.tsv`. It is created by
`scripts/05_validate_package.py` after all other steps succeed.
