# Figure coverage and regeneration

Numbering follows the November 12 SI. Commands below run from the repository root. Plot regeneration uses saved data; it does not rerun model inference or reproduce unavailable analyses. S1–S4 retain the supplied plotting layout. S4's historical short labels are explained explicitly in the methods key.

| Display | Included material | Current reproducibility |
| --- | --- | --- |
| Figure 1 | Editable PowerPoint, user-supplied PDF export and verified PNG, with source hashes | README image rendered from the PDF; original illustrative calculations not rerun |
| Figure 2 | Packaged supplied structural/contact tables, original plotting scripts and reviewed exports | Replot supplied query/donor panels; historical fixed-100 GREMLIN selection explicitly retained pending reconciliation |
| Figure 3 | Packaged swapped/single-sequence tables, Figure 3A plotting command and recovered statistical source | Replot pooled distributions; AF2 statistical grouping correction and D/p reconciliation remain pending |
| Figure 4 | [PDF/PNG, editable layout and PyMOL assets](../figures/Figure_4/README.md), eight AF2 structures, eight contact maps, four scripts and verified CSV/PNG/SVG results | All four numerical comparisons pass; full layout supplied separately; wording/naming reconciliation remains |
| Figure 5 | [Four reviewed AF2 mutant ensembles](../figures/Figure_5/README.md) (100 models), four alignments, per-residue confidence and configurations | Validate input/model provenance; fold labels, denominators and feature-selected ensembles remain pending |
| S1 | 34 target renders, matrices, coordinates, PyMOL sessions, assembly script | Rebuild A, B and combined layout from saved renders/matrices |
| S2 | 100/100/150 GREMLIN results and rank-tolerant MSA Transformer results | Replot complete 362-swap distributions at ≥6/12/24 |
| S3 | Continuous-score Pearson comparisons over top-L unions | Replot both methods; undefined correlations excluded as recorded |
| S4 | Top-L precision/Jaccard, selection and mapping audits, source scoring code | Replot A, B and combined; rank-1 predicted references, not experimental structures |
| S5 | [Supplied SI PDF/PNG and caption](../figures/S5/README.md) | Full-ensemble calculation table/code and control aggregation missing |
| S6 | Corrected 512-row/unweighted Pairformer package, all 396 maps | Replay analysis and four-metric plots; newer than the three-panel figure in local SI |
| S7 | [GA95/GA95-L45Y SI PNG and caption](../figures/S7/README.md) | Final inference outputs, exact L45Y inputs and fold-count sources not located |
| S8 | Two RfaH/NusG AF3 subsamples per target with matched maps and replayable metrics; other candidates locally | Final panel/run mapping and remaining AF3 case-study counterparts required |
| S9 A–C | Chimera prediction analysis, representative structures, editable renders and panel code | Rebuild prediction panels; filenames migrated from historical S7 to S9 |
| S9 D–E | [Full SI artwork including NMR](../figures/S9/README.md) | Raw/processed experimental sources and exact construct mapping awaited |
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

## Newly integrated collaborator calculations

These commands are separate from `reproduce_figures.py --figure all`, which retains the existing S1–S4/S6/S9 plotting scope. They validate/replay supplied case-study calculations rather than assemble the full main Figure 4 layout:

```bash
python -m pip install -r requirements-structure.txt
python scripts/validate_collaborator_import.py
python scripts/reproduce_case_studies.py
```

See the [case-study guide](../workflows/case_studies/README.md) for output paths and repeat runs, the [mutant records](../data/publication/Figure_5/README.md) for unassigned classifications, and the [remaining-source guide](MISSING_DATA.md) for scientific decisions still pending.

## Main Figures 2 and 3A supplied plots

```bash
python scripts/plot_main_figures.py --figure all
```

This separate command reproduces the supplied Figure 2 query/donor panels and Figure 3A pooled distributions. It uses the [six packaged CSVs](../data/publication/Figures_2_3/README.md), which are all eligible for Git, and writes to `outputs/main_figures/`. Use a fresh `--out-dir` for a repeat run. Figure 2 retains the historical contact selection; no Figure 3 statistical test is executed. The `reproduce_figures.py --figure all` command above retains its S1–S4/S6/S9 scope.

## Supplied SI reference artwork

[Supporting-data review](SUPPORTING_DATA_REVIEW_20261113.md) records the S1–S9 comparison and unresolved caption/data discrepancies. Figures S5–S9 now include source SI artwork under `figures/S*/SI_reference/`, and all nine figure folders contain the supplied caption transcription. These reference files are separate from generated outputs. In particular, [S6](../figures/S6/README.md) preserves the older three-panel image, [S8](../figures/S8/README.md) still needs exact AF3 model/run mapping, and [S9](../figures/S9/README.md) includes NMR artwork without its underlying experimental data.
