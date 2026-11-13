# Jaccard spot check

Audit date: 2026-10-30. Existing contact outputs were read locally; no GREMLIN, MSA Transformer, or AlphaFold inference was rerun. Original scripts and scientific result files were not modified.

## Sampling and checks

Seed `20261030` selects 10 target/partner pairs without replacement from each nominal size class (60, 100, 150 residues), using sorted pair IDs and Python's `random.Random`. Two additional corrected cases were randomly selected from repaired pairs not already sampled: `1cqu_2mm4_aligned` and `2ovo_2mm4_aligned`. The repair supplement is reported separately from the 30-pair random sample.

Each pair was compared against A and B using both methods: 128 comparisons total. The audit:

- Imports the bundled GREMLIN v2 and MSA Transformer v3 functions and reads the corrected paths where applicable.
- Reconstructs contact sets with the existing filtering, center-cropping, and GREMLIN donor-offset search.
- Compares recalculated scores and GREMLIN offsets against the saved tables.
- Evaluates standard intersection-over-union manually, with scikit-learn, and as one minus SciPy's Jaccard distance on identical Boolean contact indicators.
- Independently reconstructs MSA Transformer threshold sets and the tolerant numerator using Python sets.
- Requires usable gap masks, valid GREMLIN indices, finite matrices and scores, and nonempty unions. No sampled comparison had an empty contact set. Empty-set conventions are therefore not validated by this sample.

Numerical agreement thresholds were 1e-10 for saved scores and 1e-12 for library comparisons. An absolute difference of 0.05 was used as a descriptive screening flag for tolerant versus standard overlap, not a statistical significance threshold. Repeated targets and equal sampling across unequal-size strata mean the following means should not be interpreted as whole-dataset estimates or independent biological replicates.

## Results

**All 128 recalculated scores exactly match the saved CSV values, and all GREMLIN selected offsets match.** Manual, scikit-learn, and SciPy standard Jaccard agree within 1.12e-16.

Means for the 30 randomly selected pairs:

| Method | Reference | Current score | Standard Jaccard | Mean difference | Comparisons differing by at least 0.05 |
| --- | --- | ---: | ---: | ---: | ---: |
| GREMLIN | A | 0.021908 | 0.021908 | 0 | 0/30 |
| GREMLIN | B | 0.913408 | 0.913408 | 0 | 0/30 |
| MSA Transformer | A | 0.048090 | 0.024472 | 0.023618 | 2/30 |
| MSA Transformer | B | 0.816782 | 0.625030 | 0.191752 | 30/30 |

The largest tolerant-versus-standard difference is 0.272727. For example, `2fs1_1bw6_aligned` versus B is 0.969697 with tolerance and 0.696970 with standard Jaccard. In the two extra repaired cases, the mean MSA Transformer B score is 0.411656 with tolerance versus 0.244598 with standard Jaccard.

The donor-versus-target direction remains the same in all 30 random MSA Transformer pairs: B scores higher than A under both definitions. This supports the direction of the sampled comparison, but does not make the score magnitudes interchangeable or establish a population-wide result.

## Interpretation

The project lead recalls that visual inspection showed contact maps that appeared nearly identical despite lower strict Jaccard scores. Together with concerns about noise and score distributions in GREMLIN versus MSA Transformer, this motivated the tolerant definition as a way to reflect the shared contact patterns more faithfully. This historical rationale was supplied after the numerical spot check; the audit did not independently repeat or validate that visual comparison.

The tolerance credits an exact residue pair selected among the top L contacts in one matrix when it appears among the top 1.5L contacts in the other. It therefore accommodates rank changes around the cutoff, not shifts to neighboring residue pairs. Rank selection alone is invariant to positive uniform rescaling of a matrix, so scaling by itself does not explain a strict-versus-tolerant difference; changes in ranking or selection do. The 1.5L choice remains a methodological parameter rather than a value validated by visual resemblance alone.

No arithmetic mismatch or saved-score reproduction failure was found in this sample. There is no indication that switching between standard-Jaccard libraries would change the results.

The current MSA Transformer score is a **tolerant contact-overlap score**, not ordinary Jaccard. Its difference from standard Jaccard is substantial in the sampled donor comparisons. Figures and tables retaining these values should explicitly label that definition. Using standard Jaccard instead would require regenerating the corresponding results and figures, not merely changing the library call. This audit does not change that scientific choice or expand to all 362 pairs.

Reproduction does not independently validate biological residue correspondence, center-cropping, or the asymmetric GREMLIN offset optimization. Native GREMLIN gap masks use `original_proper_msas_af_sampling/*.fas` candidates: all 34 had matched the GREMLIN logs' alignment dimensions and query labels in the earlier input check, but exact historical identity of all alignment rows is unconfirmed.

## Reproduction and artifacts

This report describes the historical random spot check with seed 20261030, 10 pairs per nominal size class, and two additional repaired pairs. The audit ran with Python 3.11.15, NumPy 2.2.6, pandas 2.3.3, SciPy 1.17.1, and scikit-learn 1.9.0. The original full audit and input checksums remain in the working archive; all per-comparison results and summary values are included here.

The minimal repository includes four illustrative raw-data pairs rather than all 32 audit pairs. Run `python scripts/validate_examples.py` from its root to reproduce those examples and verify their scores against the complete result tables. This validates the packaged calculation; it is not a rerun of the entire random sample.

- [All 128 comparisons](../data/validation/comparisons.csv)
- [Audit summary](../data/validation/summary.csv)
- [Portable contact calculation](../scripts/compute_contacts.py)
- [Example validation](../scripts/validate_examples.py)

The historical filenames in this report identify the audited source versions. Their analysis functions are preserved in `scripts/gremlin_contacts.py` and `scripts/msatr_contacts.py`; fixed cluster paths and batch discovery have been replaced by a portable input manifest. See [file provenance](file_sources.csv).
