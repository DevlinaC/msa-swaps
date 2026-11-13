# Main structural-score evidence

AF2 records were parsed from recovered raw TM-align logs and include both chain-normalized scores, archived filenames, model and seed. AF3 scores are the original saved table linked to local model paths, not freshly recalculated alignments. Repeated native runs do not enter the swapped table.

Run from the repository root:

```bash
python scripts/summarize_swapped_structures.py
```

This writes one row per swapped prediction with both reference TM-scores (9,050 per AF method), model/seed/sample identifiers, normalization notes, the explicit ≥0.6 threshold and a summary to `data/publication/structure_scores/derived/`. It also preserves a long comparison table with both available AF2 normalizations. No model/seed is paired across AF2 and AF3.

AF2 uses Chain_2, matching the recovered historical extracted scores. The supplied AF3 parser intends reference-normalized Chain_1 scores, but its actual historical invocation and reference/run manifest remain unconfirmed. The saved-score zero-donor conclusion can be checked here; this does not substitute for reproducing TM-align from coordinates or Figure 3's permutation test.
