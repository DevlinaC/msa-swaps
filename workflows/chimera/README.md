# Figure S9 NusG–RfaH chimera

Use `Figure_S9.pdf` for the combined prediction panels A–C; the complete S9 including NMR is not generated here. PNG (350 dpi) and editable SVG exports are also available. The manuscript-format [legend](Figure_S9_legend.md) describes panels A–C. The current SI places BCCIP/RfaH/NusG in S8. These exports cover only prediction panels A–C of S9; [Full SI artwork including NMR](../../figures/S9/README.md) is now included separately; underlying experimental data and plotting sources are still required.

The figure uses actual predicted coordinates, burgundy #AC0035 for the β-sheet CTD, gold #DAA520 for the α-helical-state key, and grey for the remaining residues. The representative structures were selected by the original model ranking, independently of CTD fold: AF2 by whole-model pLDDT; AF3 by ranking_score. No additional predictions were generated.

## Individual display items (12 November 2026)

- [Panel A: construct schematic](Panel_A_construct/README.md), with mapping data.
- [Panel B: structures](Panel_B_structures/README.md), with individual renders, representative-model data and editable PyMOL session.
- [Panel C: ensemble results](Panel_C_results/README.md), with separate state/confidence plots and all per-model data.
- [Scripts and reproduction instructions](scripts/README.md), including exact original scripts and portable copies.
- [Canonical analysis data](data/README.md).

Each complete panel is exported as PDF, editable SVG and 350-dpi PNG. The two C plots also have separate SVG/PNG exports. The original combined figure and analysis data are preserved. Structure coordinates remain in `structures/`.

## Panel-generation commands

See [commands for each panel and the combined figure](scripts/README.md#commands-for-each-panel-and-the-combined-figure) for copy-and-paste commands to export A, B, and C individually, then generate the combined Figure S9.

## Results

| Measurement | AF2 / ColabFold | AF3 |
|---|---:|---:|
| Predictions examined | 25 | 25 |
| β-sheet CTD | 25/25 | 25/25 |
| α-helical hairpin CTD | 0/25 | 0/25 |
| CTD mean pLDDT ≥70 | 18/25 | 25/25 |
| Median CTD mean pLDDT | 72.38 | 86.95 |
| CTD mean pLDDT range | 61.80–84.30 | 85.72–89.68 |
| Minimum CTD TM-score to method's rank-one model | 0.974 | 0.982 |

All 50 coordinate sequences match the supplied 183-residue FASTA. The analyzed CTD is chimera residues 131–183, exactly identical to wild-type RfaH residues 110–162 used in the manuscript. The longest exact C-terminal match extends from chimera 124–183 to RfaH 103–162. This establishes the CTD mapping; it does not identify an exact splice boundary where parental sequence identities overlap. Panel A marks the analyzed CTD boundary, not an independently verified engineering junction.

PyMOL DSS assigns 34 strand residues and 0–3 helical residues within the 53-residue AF2 CTDs, and 34–35 strand residues and zero helical residues in the AF3 CTDs. The rank-one cartoons and very high CTD-only TM-align similarity across all samples support the common β-sheet fold assignment. A short three-residue helix in some AF2 CTDs does not constitute the alternative α-helical hairpin. State counts include lower-confidence predictions. Confidence is the unweighted mean CA pLDDT over the stated residue interval. Samples are repeated predictions of one sequence, not independent biological replicates.

The result supports the narrower claim that this chimera does not reproduce the α-helical RfaH CTD under the tested conditions. It is consistent with dependence on sequence context. Because fresh MSAs were generated for the chimera, this experiment alone does not isolate the effect of the N-terminal sequence from the accompanying MSA changes, or establish the physical fold of the chimera.

## Files and reproduction

- `data/models.json`: all 50 model identities, source paths/hashes, ranks, CTD secondary structure, confidence and within-method CTD TM-scores.
- `data/summary.json` and `data/domain_mapping.json`: ensemble summaries and exact residue mapping.
- `structures/`: unchanged rank-one coordinate files; raw ensembles remain under `new_af3_runs/chimera_20261106/results/` in the project.
- `editable/Figure_S9_structures.pse`: PyMOL session with full-model and CTD scenes. Full cartoons share a CTD alignment; differing domain arrangements are present in the predictions.
- `provenance/`: input/output hashes, analysis choices, software versions and validation.

Reproduce using `Manuscript/review_20261106/analyze_chimera_s7.py`, then `render_chimera_s7.py`, then `plot_chimera_s7.py` from the complete project. The local `/opt/anaconda3/bin/python3` environment supplies Biopython, NumPy, tmtools, PyMOL and Matplotlib. Rendering and analysis are bounded to one thread. Exported images are cropped only to remove white margins; coordinates are not deformed.

AF2/ColabFold used five models × five seeds, three recycles and max MSA 512:5120. AF3 used five seeds × five samples and ten recycles. Both runs omitted templates and used separately generated MSAs. See the original run directory for job/configuration provenance; no new Biowulf computation was required for this figure.
