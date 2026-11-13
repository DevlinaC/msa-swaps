# Figure coverage and regeneration

Numbering follows the November 12 SI. Commands below run from the repository root. Plot regeneration uses saved data; it does not rerun model inference or reproduce unavailable analyses. S1–S4 retain the supplied plotting layout. S4's historical short labels are explained explicitly in the methods key.

| Display | Included material | Current reproducibility |
| --- | --- | --- |
| Figure 1 | Preview recovered from the current manuscript, with document/media hashes | Final editable/high-resolution artwork awaited; old story image is historical |
| Figure 2 | Model-labelled swapped TM-score table at ≥0.6, historical contact examples and revised separation analysis in S2 | Final exact panel source/selection needs confirmation |
| Figure 3 | Saved structural-score tables and single-sequence records | Query-stratified permutation KS implementation and final plotting source missing |
| Figures 4–5 | Local candidate case-study ensembles inventoried in earlier review | Final run selections, labels, contact matrices and plotting sources required before packaging |
| S1 | 34 target renders, matrices, coordinates, PyMOL sessions, assembly script | Rebuild A, B and combined layout from saved renders/matrices |
| S2 | 100/100/150 GREMLIN results and rank-tolerant MSA Transformer results | Replot complete 362-swap distributions at ≥6/12/24 |
| S3 | Continuous-score Pearson comparisons over top-L unions | Replot both methods; undefined correlations excluded as recorded |
| S4 | Top-L precision/Jaccard, selection and mapping audits, source scoring code | Replot A, B and combined; rank-1 predicted references, not experimental structures |
| S5 | Figure/legend only in local SI, not redistributed here | Full-ensemble calculation table/code and control aggregation missing |
| S6 | Corrected 512-row/unweighted Pairformer package, all 396 maps | Replay analysis and four-metric plots; newer than the three-panel figure in local SI |
| S7 | GA95/GA95-L45Y figure visible in SI | Final inference outputs, exact L45Y inputs and fold-count sources not located |
| S8 | BCCIP/RfaH/NusG figure visible in SI; candidate runs locally | Final run mapping and matched contact matrices/analysis required |
| S9 A–C | Chimera prediction analysis, representative structures, editable renders and panel code | Rebuild prediction panels; filenames migrated from historical S7 to S9 |
| S9 D–E | NMR artwork visible in SI | Raw/processed experimental sources and exact construct mapping awaited |
| Table S1 | CSV/XLSX, sequences, selected chains, CATH and ECOD provenance | Table provided; original ECOD integration script still needs recovery/documentation |

## Commands

```bash
# Verify imported data and key scientific invariants.
python scripts/validate_publication.py

# Reconcile saved main structural scores; no new TM-align.
python scripts/summarize_swapped_structures.py

# S1 A (fold cartoons), B (TM-score heatmaps), and combined figure.
python scripts/plot_s1.py

# S2 and S3, separately or together.
python scripts/plot_s2_s3.py --figure 2
python scripts/plot_s2_s3.py --figure 3

# S4 A (precision), B (ordinary Jaccard), and combined figure.
python scripts/plot_s4.py

# Current Pairformer four-panel figure, individual panels and two rows.
python workflows/pairformer/scripts/plot_summary.py --root workflows/pairformer

# S9 prediction panels A, B and C, then combined A–C.
python workflows/chimera/scripts/export_panels.py --panel A
python workflows/chimera/scripts/export_panels.py --panel B
python workflows/chimera/scripts/export_panels.py --panel C
python workflows/chimera/scripts/plot_chimera_s7.py

# Run all of the above; the original chimera script name is retained for traceability.
python scripts/reproduce_figures.py --figure all
```

S1–S4 write PNG/PDF/SVG into `figures/S1` through `figures/S4`. Pairformer writes into `workflows/pairformer/results/analysis/`, including `panels/`. Chimera writes `Figure_S9.*` and `Panel_*/Figure_S9[A-C].*` under `workflows/chimera/`; imported `Figure_S7.*` files are historical copies, not GA95 figures.

## Reanalysis versus plotting

Pairformer calculations can be independently replayed from the bundled matrices on a CPU. To preserve the imported validation records, copy `workflows/pairformer/` to a fresh working folder, run `analyze_full.py --root <copy>`, then `analyze_rank_tolerant.py --root <copy>` and `plot_summary.py --root <copy>`. Preparation and GPU inference have additional external dependencies documented in that workflow.

The original S2–S4 calculation scripts are in `provenance/source_scripts/contact_analysis/`. They retain the original project paths and scientific checks. They are supplied as auditable source, **not claimed to run against this compact package without the original raw-data layout**. The portable `scripts/plot_*.py` adapters change input/output paths only; their input tables are imported byte-for-byte. Full raw GREMLIN/MSA Transformer matrix and structure archives still need a durable public archive.

## What code should be shared?

Share the code and source data for every quantitative plot, its preprocessing/statistical calculation and the exact commands. A final high-resolution image alone cannot establish how observations were filtered or counted. For manually arranged structure panels, share the selected structures, residue selections, editable PyMOL scenes, renders and layout source. A schematic such as Figure 1 can be supplied as editable vector artwork without a Python generator. Final publication typography need not be perfectly automated if the scientific content is reproducible and manual steps are documented.
