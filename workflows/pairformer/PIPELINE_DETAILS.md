# Manuscript matched MSA-Pairformer pipeline

This run reconciles the swap analysis with the manuscript's maximum of 512
sequences and disabled query weighting. It supersedes the separate
`../msa_pairformer_full/` run (256 sequences, query bias enabled) for manuscript
use. The earlier run remains preserved for comparison.

## Inputs and settings

There are 34 proper and 362 swapped inputs, including eight corrected donor-2mm4
alignments. Original files and hashes come from `source_inventory.json` under
`/data/AF_interpretation/msa_swaps/` on Biowulf. Selected inputs are written to
`inputs/` with the query retained as row zero. Deterministic maximum-mean-Hamming
diversity selection retains at most 512 rows; shallower alignments keep every
row. No HHfilter, deduplication, or residue cropping is applied during inference.
The pre-existing homolog-based donor mappings and four comparison definitions
are retained. RfaH/NusG's five random subsamples in the manuscript refer to that
separate case study, not these 362 swaps.

The pinned source, checkpoint hashes, environment, Cβ contact head, BF16
precision, seed 42, and absence of chain breaks match the earlier run.
`model.turn_off_query_biasing()` disables learned query-dependent sequence
attention in every outer-product layer. With no external sequence weights,
the installed code uses uniform weights over valid MSA rows. Setting
`return_seq_weights=False` only suppresses returned tensors; it does not disable
weighting. Here the model toggle is explicitly switched off, all layer flags
are checked, and returned weights are audited for three small inputs and a
largest 512-row input. `results/WEIGHTING_AUDIT.json` records the checks.

## Execution and restart

Preparation and inference each run as a single grouped job. The GPU runner loads
one model and processes all 396 inputs, writing validated per-input checkpoints.
It first replays the old proper-1gb1 prediction with bias enabled as a technical
control. It then disables bias for every new result. No old biased prediction
is reused as a manuscript-matched output, even if the selected MSA is unchanged.
The three unweighted smoke-test predictions are reused within the new run.
Larger selected workloads run first and provide a memory/timing gate before
continuing. There are no per-input jobs or artificial delays.

Submit from the persistent project root after checking existing job ownership:

```bash
sbatch scripts/prepare.sbatch
# After successful preparation (or using afterok:<preparation_job>):
sbatch scripts/run.sbatch
```

The batch scripts use the sibling `msa_pairformer_pilot/env` and pinned weight
cache. They request bounded CPU, memory, scratch, and GPU resources and use
stage locks. The GPU batch currently performs inference only. For the completed
saved-data analysis, execute on a compute node or on the local workstation:

```bash
python scripts/analyze_full.py --root "$PWD"
python scripts/analyze_rank_tolerant.py --root "$PWD"
python scripts/plot_summary.py --root "$PWD"
python scripts/compare_configurations.py --root "$PWD" --old ../msa_pairformer_full
python scripts/build_inventory.py --root "$PWD"
```

The four metrics are ordinary top-L Jaccard, rank-tolerant top-L Jaccard,
Pearson over all eligible pairs, and Pearson over the union of top-L pairs.
Each is evaluated at minimum sequence separations 6, 12, and 24. L denotes the
homolog-registered overlap, shared between query and donor comparisons. Ties
are retained. Rank tolerance credits a top-L pair when it is in the other
map's top floor(1.5L), using the union of the strict top-L sets as denominator.
The original table, rank-tolerant table, and combined table each have 2,172 rows.

## Provenance

- `provenance/protocol_reconciliation.json`: manuscript passage and hash,
  original and corrected settings, and preserved analysis choices.
- `provenance/submissions.json`: scheduler job IDs.
- `results/RUN.json`: new inference fingerprint, environment, model and weights.
- `results/WEIGHTING_AUDIT.json`: technical replay and uniform-weight audit.
- `results/INFERENCE_COMPLETE.json`: final inference validation.
- `results/analysis/VALIDATED.json` and `results/rank_tolerant/VALIDATED.json`:
  numerical checks and matrix hashes.
- `results/configuration_comparison/`: paired old/new metrics and matrix changes.
- `inventory/`: matrix indexes, file list, and checksums.

Changing depth and weighting together resolves the stated protocol discrepancy;
it does not separately estimate the effect of each change. Historical MSA
Transformer subsets were not recovered, and its registration differs, so
absolute cross-method scores are still not a controlled performance comparison.
