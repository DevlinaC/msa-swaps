# Top L Pearson analysis of GREMLIN and MSA Transformer contacts

Item 3 extends the historical fixed-cap GREMLIN separation analysis (preserved [calculation source](../../../provenance/source_scripts/contact_analysis/analyze_gremlin_separation.py)) with Pearson correlations at **|i−j| ≥6, ≥12, and ≥24**. The selected definition is correlation of continuous scores at matching residue pairs in the **union of the two top-L contact sets**, for both GREMLIN and MSA Transformer. The analysis was completed locally on November 2, 2026.

The results support the manuscript's observation that swapped-MSA contact information generally resembles the donor more strongly than the query. That pattern is not universal under every Pearson comparison: MSA Transformer has four query-favoring comparison records across three pairs, and 38 GREMLIN query comparisons cannot supply L jointly eligible pairs. GREMLIN results also depend materially on the existing donor-offset adjustment. These qualifications are part of the result.

## Main results

Each defined query/donor comparison contributes once to its mean. A is the sequence query and B is the MSA donor. Means for A and B use different sample counts where A is undefined; the final column compares only pairs for which both are defined.

| Method | Separation | Mean R versus A | Mean R versus B | Defined A / B | B greater than A among jointly defined pairs |
| --- | ---: | ---: | ---: | ---: | ---: |
| GREMLIN | ≥6 | −0.618221 | 0.935043 | 356 / 362 | 356/356 |
| GREMLIN | ≥12 | −0.585257 | 0.937825 | 352 / 362 | 352/352 |
| GREMLIN | ≥24 | −0.465551 | 0.948687 | 340 / 362 | 340/340 |
| MSA Transformer | ≥6 | −0.593924 | 0.633075 | 362 / 362 | 362/362 |
| MSA Transformer | ≥12 | −0.465961 | 0.719965 | 362 / 362 | 361/362 |
| MSA Transformer | ≥24 | −0.301693 | 0.746040 | 362 / 362 | 359/362 |

GREMLIN donor medians are 0.999308, 0.999031, and 0.998209; MSA Transformer donor medians are 0.658448, 0.756525, and 0.816888. Pair-level values, quartiles, ranges, and length-class summaries are retained in the CSVs. No inferential p-values are reported: residue pairs are dependent, and the 362 swaps repeatedly use the same 34 proteins.

Current publication exports are in [Figure S3](../../../figures/S3/README.md). The older combined Jaccard/Pearson report is preserved in the local-only archive; it is not the current S3 layout.

The combined artwork above is the original analysis export. For current manuscript S3 plotting commands, use the [figure guide](../../../docs/FIGURES.md). Preserved calculation scripts require the original external data layout; they are not portable plotting entry points.

## How the correlation is defined

For each method, swap, reference A/B, and separation cutoff:

1. Establish a shared coordinate system and eligible residue-pair universe. Both inputs must have a finite score for the same residue pair, and that pair must meet the separation cutoff. Unavailable scores are excluded, never replaced with zero.
2. Rank each input's continuous scores within that shared universe. Select its top-L score threshold, retaining all ties at the threshold. Thus a set can contain more than L contacts.
3. Take the union of those two sets. For each pair in that union, retrieve the score from each input, even if the pair was not selected among the other input's top-L scores.
4. Calculate Pearson R between the two score vectors in identical residue-pair order. Independently sorting the vectors would lose pair identity and is not done.

The union has between L and 2L entries when there are no ties; ties may increase its size. The analysis exports the actual pair identities, scores, and set-membership indicators for inspection.

### GREMLIN correspondence and scores

- Scores are the saved continuous `zscore` values. The item 2 z-score ≥1 floor and top-100 limit are **not** applied: item 3 asks for top L among available scored pairs.
- Both native and swapped alignment columns must satisfy gap fraction <0.25. The common universe also excludes pairs absent from either saved score table after GREMLIN's upstream filtering.
- A comparisons use offset zero. B comparisons reuse the Jaccard-selected offset from item 2 for that pair and separation; Pearson R is not used to choose an offset.
- L is the length of the overlap between the original alignment-coordinate intervals after that fixed offset, before removing high-gap positions. It is a sequence-coordinate length, not the number of retained contact pairs or columns.
- Shifted indices outside the overlap are unavailable. The exported `i/j` coordinates are zero-based native-reference coordinates; subtract the recorded offset to recover swapped coordinates.
- Input hashes are checked against item 2 so its offsets cannot silently be applied to altered source tables or alignments.

### MSA Transformer correspondence and scores

- Use the saved continuous contact-score matrices with no probability threshold.
- Center-crop both matrices to their shared dimension L, preserving the historical convention: an odd size difference removes one extra position at the N terminus.
- Use the upper triangle at the requested sequence separation. No offset optimization or GREMLIN gap mask is applied.
- Exported pair coordinates refer to the cropped matrices. Add `crop_start_swapped` or `crop_start_native` from the comparison table to recover each original matrix coordinate.

### Undefined values

A primary top-L correlation is undefined if there are fewer than L shared eligible pairs. L is not silently reduced. Pearson R is also undefined for fewer than two observations or constant/near-constant score vectors. These cases remain explicit rows rather than being assigned zero or removed from the dataset.

All 38 undefined primary values here are GREMLIN query-A comparisons with insufficient shared eligible pairs: 6 at separation 6, 10 at 12, and 22 at 24. Some native/query combinations retain very few common positions after gap filtering, even though each input separately had nonempty contact sets in item 2. All 1,086 GREMLIN donor correlations and all 2,172 MSA Transformer correlations are defined. The full result has **4,344 records: 4,306 defined and 38 undefined**.

## Exceptions and sensitivity checks

### MSA Transformer query-favoring comparisons

| Swap | Separation | R versus A | R versus B |
| --- | ---: | ---: | ---: |
| `1gb1_2mm4_aligned` | ≥12 | 0.204715 | 0.100698 |
| `1gb1_2mm4_aligned` | ≥24 | −0.123186 | −0.179972 |
| `2fs1_2mm4_aligned` | ≥24 | −0.076028 | −0.150217 |
| `2k9h_2mm4_aligned` | ≥24 | −0.315162 | −0.361841 |

All involve donor `2mm4`; the first two pair IDs use repaired inputs, while `2k9h_2mm4_aligned` is outside the eight repaired cases. These observations were retained without excluding cases or changing selection rules. A negative correlation to both references does not establish strong query-like contact information merely because A exceeds B.

### Top L selection affects negative correlations

The union preferentially includes pairs that score highly in either input. When the high-scoring pairs differ, this conditioning can induce negative correlation. Consequently, negative query R values should **not** be presented as evidence of intrinsic biological anticorrelation between folds.

The diagnostic correlation across all shared eligible pairs, without top-L selection, has these means:

| Method | Separation | Mean R versus A | Mean R versus B |
| --- | ---: | ---: | ---: |
| GREMLIN | ≥6 | 0.023444 | 0.979988 |
| GREMLIN | ≥12 | 0.007364 | 0.979762 |
| GREMLIN | ≥24 | −0.000548 | 0.978935 |
| MSA Transformer | ≥6 | 0.018502 | 0.890436 |
| MSA Transformer | ≥12 | 0.005967 | 0.884024 |
| MSA Transformer | ≥24 | 0.009252 | 0.853540 |

Under this diagnostic, donor R exceeds query R for every jointly defined comparison. GREMLIN has 352 jointly defined comparisons at separation 24, because ten query comparisons have no common eligible pairs; all other method/cutoff combinations have 362. These whole-universe coefficients answer a different question from the selected top-L-union statistic and should be labeled separately.

### GREMLIN depends on residue registration

Without the inherited donor offsets, mean GREMLIN donor top-L-union R falls to **0.414696, 0.434067, and 0.475162** at separations 6, 12, and 24. Donor R then exceeds query R in **297/356, 293/352, and 272/340** jointly defined comparisons. The historical offsets therefore materially influence the result. Offset selection used contact overlap, so this is not an independent validation of residue correspondence, even though Pearson itself was not optimized. Center-cropping for MSA Transformer is likewise an assumed coordinate correspondence, not a new sequence alignment.

Native GREMLIN gap-mask alignments have the same provenance limitation as item 2: compatible dimensions and query labels were checked previously, but exact identity of all historical alignment rows remains unconfirmed.

## Files and validation

| File | Contents |
| --- | --- |
| [comparisons.csv](comparisons.csv) | All 4,344 primary comparisons, statuses, L, eligible/selected counts, cutoffs, offsets, and diagnostic correlations |
| [summary.csv](summary.csv) | Pooled and length-class summaries, with separate defined/undefined counts |
| [donor_query_direction.csv](donor_query_direction.csv) | Within-pair A/B directions for primary, whole-universe, and unshifted correlations |
| `selected_residue_pairs.csv.gz` (local-only) | 588,855 selected residue-pair records with both scores and top-L membership flags |
| [input_audit.csv](input_audit.csv) | Dimensions, gap retention, score-table counts, and matrix-symmetry checks |
| [provenance.json](provenance.json) | Input/source hashes, methods, versions, output checksums, limitations, and numerical validation |
| [COMPLETE.json](COMPLETE.json) | Completion marker and provenance checksum |
| `serialized_validation.json` (local-only) | Independent recalculation from exported pair scores, plus output-checksum and link checks |
| [Analysis script](../../../provenance/source_scripts/contact_analysis/analyze_topl_pearson.py) | Reproducible analysis and combined figure |
| Boundary tests (local-only source review) | Pair identity, ties, undefined cases, affine scaling, shifted joins, and center-cropping checks |

Four boundary tests passed. SciPy Pearson coefficients agree with NumPy correlations within `1.23 × 10⁻¹⁵` across primary and diagnostic calculations. Independent recalculation from the exported residue-pair scores reproduces all 4,306 defined primary correlations within `6.67 × 10⁻¹⁶`, with pair identities, separation, membership counts, and output checksums verified. Coverage checks enforce all 362 pairs, two methods, two references, and three cutoffs. Input checks reject invalid matrices, nonfinite scores, inconsistent dimensions/indices, missing files, and changed GREMLIN inputs. The combined figure was visually inspected and its title spacing corrected without changing numerical results.

## Reproduction

From the workspace root, using the same Python environment as item 2:

```bash
python -m unittest discover -s for_github_repo/codes -p 'test_*pearson.py' -v
python for_github_repo/codes/analyze_topl_pearson.py \
  --data-root . \
  --item2-dir for_github_repo/analysis/gremlin_separation \
  --out-dir for_github_repo/analysis/pearson_topl_rerun
```

The script uses one numerical thread and requires a new output directory. Python 3.11.15, NumPy 2.2.6, pandas 2.3.3, SciPy 1.17.1, and Matplotlib 3.11.0 were used. No model inference, downloads, or cluster jobs were required. On Biowulf, run within an allocated compute job.

## Suggested manuscript wording

“We additionally compared continuous contact scores over the union of the top-L residue pairs from each input at sequence separations |i−j| ≥6, ≥12, and ≥24. Donor-associated similarity generally persisted across these settings. Mean donor correlations were 0.935–0.949 for GREMLIN and 0.633–0.746 for MSA Transformer. GREMLIN donor correlations exceeded query correlations in all jointly evaluable comparisons (356, 352, and 340 pairs), while MSA Transformer showed this pattern in 362, 361, and 359 of 362 pairs, respectively. GREMLIN comparisons retained the previously selected global offsets; comparisons lacking L shared eligible pairs were reported as undefined. Negative query correlations under top-L-union selection were accompanied by near-zero query correlations across all eligible pairs, indicating that the negative values depend on contact selection.”

These results address contact-score agreement, not agreement with AlphaFold structural contacts or AlphaFold's internal representations. Item 1 remains pending. The manuscript itself has not been edited.

The preserved original calculation script includes its boundary-test file in provenance collection. That test suite is local-only; raw recalculation with this original source requires the complete source package. The portable S3 plotting command does not need it.
