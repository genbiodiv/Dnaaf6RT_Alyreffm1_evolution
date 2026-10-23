# Scripts

Run the compact workflow from the repository root:

```bash
bash scripts/run_all.sh
```

Execution order:

1. `01_reconstruct_event.py` — focal sequence comparisons, TSD, transcript
   classification and liftOver audit.
2. `02_analyse_synteny.py` — nine-species locus mapping, insertion-state calls
   and age bracket.
3. `04_plot_expression_evidence.py` — diagnostic plots from archived public-
   data derivatives.
4. `03_validate_expression.py` — fixed consistency checks for expression
   evidence.
5. `05_validate_package.py` — cross-module scientific and packaging checks.

`fetch_synteny_loci.py` redownloads the nine pinned NCBI loci.
`fetch_expression_inputs.sh` and `03_analyse_expression_full.py` implement the
optional full expression rerun. `scan_junction_kmers.py` performs the exact,
strand-aware junction scan used on the ENCODE FASTQ files.

All paths are resolved relative to the repository; no user-specific path is
required.
