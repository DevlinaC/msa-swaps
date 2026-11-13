# Remaining sources for the collaboration draft

Updated November 13, 2026. Missing means not recovered in the reviewed local files or targeted source checks; it does not mean the analysis was never performed.

## Ready now

S1–S4 plotting data and scripts; Table S1; GREMLIN inference source; historical ordinary/rank-tolerant overlap functions; corrected 512-row/unweighted Pairformer inputs, matrices, analysis and plots; chimera prediction panels (current S9 A–C). Main structural ensembles and saved scores are already present in the parent project, so there is no need to request every AF run again.

## Ask collaborators when shared access returns

1. **Figure 1 and final artwork:** editable layout and final PNG/SVG/PDF; replace the manuscript-derived README preview with the final master. High-resolution artwork is useful but does not replace underlying numeric sources.
2. **Figure 3:** query-stratified permutation KS code, exact input selection, output, permutations, RNG seed and weighting/exchangeability definition.
3. **S5:** full per-prediction table and calculation/plot code; donor aggregation for the 850 single-sequence/native controls; final medians resolving manuscript/SI differences; rows for named exceptional swaps.
4. **Figures 4–5 and S8:** final run IDs, selected structures, per-model fold/confidence labels, analyzed regions, matching Pairformer matrices and residue mappings, and analysis/plotting code. Many AF3 candidate runs are local; identify intended runs before requesting additional coordinates. Supply missing AF2 case-study data where applicable.
5. **S7:** exact GA95 and GA95-L45Y inputs, output ensembles, rankings/confidence and fold labels. GA88 (2JWS) is not GA95-L45Y.
6. **S9 D–E:** exact expressed construct, raw and processed NMR data, acquisition/processing parameters, assignments and plotting projects. Resolve experimental construct numbering against the 183-residue prediction input.
7. **MSA Transformer historical inference:** production `modules.MSA_TR_af_sampling`, original invocation/environment, selected rows and seeds. `MSA_TR_test.py` uses the first 64 rows and cannot establish the historical 256-row sampling procedure.
8. **Provenance:** actual AF3 TM-align execution/reference manifest, remaining inference configurations and ProteinMPNN source/design/selection records. Original project code now uses MIT; confirm upstream notices, data/asset release terms and the final manuscript citation.

## Available in the author's own project; do not request as missing collaborator data

- Corrected Pairformer swap run: `msa_pairformer_512_unweighted`, now packaged. It supersedes the 256-row/bias-enabled run for this analysis. Matching case-study matrices are a separate need.
- BLAST source location: `/data/chakravartyd2/for_lydia_project/homology_check`. Reconcile the E<0.05 manuscript screen with saved E≤0.001 summaries using raw hits; do not launch BLAST again merely because thresholds differ.
- Main AF ensembles and native/swap score records: the repository now reconstructs a complete model-labelled swapped score table from saved AF2 logs and AF3 table links. Final plotting selections and historical AF3 scoring provenance still need confirmation; do not request all models again.

## Release decisions

Final figure selections, S5 aggregation/medians, Figure 5 category denominators, exact historical MSA Transformer preprocessing, reference-normalized TM-score execution provenance and durable raw-data archive links remain open. This draft should not be described as reproducing the entire manuscript until these are closed.

## Targeted Biowulf source check on November 13

A read-only lookup in the shared `code`, `msatr`, `missing_msat` and `gremlin` folders and the personal project `scripts` folder found `/data/AF_interpretation/msa_swaps/code/MSA_TR_test.py` and the GREMLIN runner, but not `MSA_TR_af_sampling.py`. This was a bounded source-directory search, not an exhaustive scan of other users’ installations or environments. No computation or job submission was performed on the login node.
