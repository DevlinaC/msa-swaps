# Item 1 GREMLIN separation analysis with the revised class cap

This is the current item 1 analysis as of 5 November 2026. It supersedes the
fixed-cap results in `../gremlin_separation/`, which remain intact for provenance.

The original supplied `gremlin_jaccard.py` and `gremlin_jaccard_v2.py` set
`GREMLIN_TOP_N = 100`. This was a downstream contact-selection rule, not a
limit on GREMLIN score generation. The revision raises the cap to **150 only
for the 150-residue class**. The 60- and 100-residue classes retain cap 100.
The z-score ≥1 filter, gap fraction <25%, minimum sequence separations, sorting
and tie handling, and donor offset search −5…+5 are unchanged.

This means **up to 150 eligible contacts**, not necessarily 150 in every input.
Among the 100 native-plus-swapped inputs in the 150-residue class, 80 reach the
cap at separation ≥6, 60 at ≥12, and 60 at ≥24. All remaining inputs use all
eligible contacts passing the original filters.

| Separation | Mean query Jaccard | Mean donor Jaccard | Donor greater than query |
| --- | --- | --- | --- |
| ≥6 | 0.025302 | 0.921258 | 362/362 |
| ≥12 | 0.023002 | 0.915685 | 362/362 |
| ≥24 | 0.017466 | 0.912387 | 362/362 |

All 724 original ≥5 top-100 comparisons were independently replayed; maximum
absolute error was 1.11×10⁻¹⁶. All 1,632 comparisons in the 60- and 100-residue
classes are unchanged. All 2,172 donor/query offsets are identical to the
original separation analysis, so item 2's existing Pearson calculations retain
the same inputs and registration and do not need recomputation.

`contact_selection_audit.csv` records the cap, eligible count, and actual count
for each input and cutoff. `selected_contacts.csv.gz` exports selected pairs
and scores. `changes_from_top100.csv` contains paired old/new outcomes.
`provenance.json` records parameters, input hashes, and output hashes;
`COMPLETE.json` authenticates it. The producing script is
`../../codes/analyze_gremlin_cap150_separation.py`.
