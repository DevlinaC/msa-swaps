# MSA-Pairformer manuscript run data guide

Use this project for the corrected swap analysis: **up to 512 MSA sequences,
query weighting disabled**. The earlier `msa_pairformer_full` project used
256 sequences with query bias enabled and is retained as a historical comparison.

Biowulf project root:

```text
/data/chakravartyd2/for_lydia_project/msa_pairformer_512_unweighted/
```

## Matrices and input alignments

The 362 swapped predictions are in `results/swap_<query>_<donor>/`.
The 34 proper reference predictions are in `results/proper_<protein>/`.
Each folder contains:

- `contacts.npy`: full square NumPy contact-score matrix, preferred for analysis.
- `full_contacts.csv`: identical numeric matrix with no headers or row labels.
- `COMPLETE.json`: checksum, input identity, settings, and inference provenance.

For example, `swap_1gb1_2mm4` uses query 1gb1 and donor MSA 2mm4. Reversing these
IDs denotes a different experiment. Scores are model-predicted Cβ contact scores
in [0,1], not distances, binary contacts, or experimental observations.
Matrix rows/columns follow the ungapped first sequence in `inputs/<id>.fas`.
Index zero denotes its first residue; indices are not PDB residue numbers.

The exact selected alignments are in `inputs/`. Each `.fas` has a matching
`.json` containing source paths and hashes, selected row indices, query sequence,
and donor mapping. `inputs/manifest.json` covers all 396 inputs.

For comparisons, use the recorded homolog mapping: swap/query indices are
`start:stop` (zero-based, stop exclusive), and proper-donor indices are those
same indices plus `offset`. Both comparisons use the same overlap. For
`swap_1gb1_2mm4`, swap indices 0–55 map to donor indices 1–56. Do not blindly
compare the full matrices elementwise or replace the registration by cropping.

## Analysis and figures

| Path | Contents |
| --- | --- |
| `results/analysis/comparisons_with_rank_tolerance.csv` | All four metrics, 2,172 swap/reference/separation records. |
| `results/analysis/summary_with_rank_tolerance.csv` | Means and medians for 24 metric/reference/separation groups. |
| `results/analysis/comparisons.csv` | Ordinary top-L Jaccard and both Pearson correlations. |
| `results/rank_tolerant/comparisons.csv` | Rank-tolerant Jaccard with selection/credit counts. |
| `results/analysis/mapping_audit.csv` | Homolog registration and overlap boundaries. |
| `results/analysis/paired_direction.csv` | Donor/query preference counts for the original three metrics. |
| `results/rank_tolerant/paired_direction.csv` | Donor/query preference counts for rank-tolerant Jaccard. |
| `results/analysis/contact_comparison.*` | Main four-panel figure: ordinary Jaccard, rank-tolerant Jaccard, all-pair Pearson, top-L-union Pearson. |
| `results/analysis/figure_legend.txt` | Written figure legend, including the 512-sequence cap and disabled query weighting. |
| `results/analysis/panels/upper_row_jaccard.*` | Separate upper row. |
| `results/analysis/panels/lower_row_pearson.*` | Separate lower row. |
| `results/analysis/panels/panel_A_*` through `panel_D_*` | Individual panels for figure assembly. |
| `results/rank_tolerant/jaccard_comparison.*` | Ordinary versus rank-tolerant sensitivity figure. |
| `results/configuration_comparison/` | Paired differences from the preserved 256-sequence, biased run. |

Figures are provided as PNG, vector PDF, and SVG with editable text. Red denotes
query references and gold donor references. All 362 swaps appear in each group;
black bars mark medians. Jitter is horizontal only. Configuration details appear
in the written figure legend (`results/analysis/figure_legend.txt`) and README;
plot titles identify the metrics.

## Scripts and reproducibility

All code is in `scripts/`. `prepare_full.py` selects rows and registers donors;
`run_full.py` generates predictions; `analyze_full.py` computes ordinary Jaccard
and Pearson; `analyze_rank_tolerant.py` computes the historical rank-tolerant
metric; `plot_summary.py` produces combined tables and all main panel exports.
`compare_configurations.py` compares old and new outputs, and `build_inventory.py`
creates the file indexes. Batch scripts and shared helpers also live here.

`README.md` summarizes the result; `PIPELINE_DETAILS.md` explains execution and
restart. The inference environment and model-weight cache remain in the sibling
`msa_pairformer_pilot` folder on Biowulf. They are not required to read matrices.

Check `results/RUN.json`, `results/WEIGHTING_AUDIT.json`, and the analysis
`VALIDATED.json` files for settings, uniform-weight checks, and numerical
validation. `provenance/protocol_reconciliation.json` records the manuscript
passage and why the original run was replaced for manuscript use. Source MSAs
remain under `/data/AF_interpretation/msa_swaps/`.

## Inventories and sharing

- [Swapped matrix index](inventory/swapped_matrices.tsv): all 362 swaps, query/donor
  IDs, dimensions, depths, paths, and checksums.
- [All matrix index](inventory/all_matrices.tsv): all 396 predictions.
- [CSV file list](inventory/swapped_csv_files.txt) and
  [NumPy file list](inventory/swapped_npy_files.txt): one relative path per swap.
- [Complete file inventory](inventory/files.tsv), [directory listing](inventory/filesystem.txt),
  and [checksums](inventory/SHA256SUMS).

All inventory paths are relative to this project root. Share the whole project
folder for analysis reproducibility, or the desired matrix format plus the
matrix index and this guide for matrix-only sharing. Caches and generated
inventory files are excluded from the main inventory snapshot; generated index
hashes are listed separately in `inventory/VALIDATED.json`.
