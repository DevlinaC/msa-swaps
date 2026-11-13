# Supporting figure review — November 13, 2026

Reviewed the local `Manuscript/supporting_data/` package against this repository. The detailed figure-folder comparison is retained locally in `provenance/supporting_figure_comparison_20261113.csv` (excluded from Git as a redundant review inventory). The decisions below summarize the review before import. A different hash means different bytes, not necessarily a different scientific result. Source archives, duplicate exports, historical folders and the Word manuscript were not copied wholesale.

Added the supplied SI artwork for S5–S9 and caption transcriptions for S1–S9. The [import manifest](../provenance/supporting_figure_import_20261113.csv) records original and packaged hashes, the SI document hash and any caption-link adaptation. Artwork hashes were checked against the source artwork manifest, and the corresponding embedded media hashes were checked against `SI.docx`. PNGs were visually inspected. PDF and PNG files are unchanged source exports; no new scientific figures or data were calculated.

## Coverage and decisions

| Figure | Review and repository action |
| --- | --- |
| S1–S4 | Existing calculation data are retained. Added source SI captions alongside current generated figures; retained portable repository plotting commands rather than importing alternative path-specific wrappers. Caption terminology and numerical discrepancies below remain explicit. |
| S5 | Added supplied PDF/PNG. No per-prediction calculation table or regeneration script was found in its source folder. S4 is a different calculation and cannot supply these missing results. |
| S6 | Added the older three-panel SI image as reference only. The current four-panel 512-sequence/unweighted package remains under `workflows/pairformer/`; its data already have identical copies in the repo. |
| S7 | Added native PNG showing GA95/GA95-L45Y sequences, two structures and a confidence scale. No model ensembles or fold-count source tables were supplied in this folder. |
| S8 | Added PDF/PNG. Artwork depicts BCCIP and RfaH/NusG structures and summary tables. Exact AF3 model/map assignments remain unresolved; the Figure 4 AF2 and 512-sequence AF3 replay results are not substitutes for this mapping. |
| S9 | Added full SI PDF/PNG including NMR artwork, separate from the existing regenerable A–C prediction panels. Calculation data already have identical copies in the repo. Raw spectra, assignments, experimental construct confirmation and NMR plotting sources remain missing. |

## Issues to reconcile with the manuscript

- **S1:** source review reports within-length-class mean/median TM-scores 0.292/0.291 (181 unordered pairs), versus manuscript 0.32/0.31. The statement that every pair is below 0.6 needs the within-class qualification; the source review reports three cross-class pairs above 0.6.
- **S2:** GREMLIN uses ordinary Jaccard; MSA Transformer uses rank-tolerant Jaccard. The supplied caption title does not distinguish them. Keep the existing 100/100/150 GREMLIN caps.
- **S3:** the implementation uses ordinary Pearson correlation of continuous scores over the top-L union, despite the caption's “rank-tolerant Pearson” wording.
- **S4:** retain the existing metric/mapping definitions and predicted-reference provenance; do not reinterpret this as S5.
- **S5:** manuscript and SI disagree on all single/full-MSA medians. In AF2 GREMLIN, for example, the source review records manuscript 3.61/5.20% versus SI 3.76/5.51%. Exact source tables and aggregation must resolve the discrepancy; no values were changed here.
- **S6:** the older SI image and newer four-panel analysis are distinct. Choosing replacement artwork requires a matching caption, including the additional rank-tolerant panel.
- **S7:** its caption describes ensemble results and panels not shown in the supplied artwork.
- **S8:** its caption describes contact maps and direct representation-to-structure comparisons not displayed in the artwork. Its summary values must be matched to the actual AF3 calculation before claiming reproduction.
- **S9:** the caption describes D and E; the supplied artwork visibly labels D above the two NMR plots but has no separate E label. The experimental methods' 126 + 60 residue description totals 186, while the computational construct is 183 residues. The analyzed CTD mapping does not establish the physical splice boundary.

These are review findings, not corrections to the manuscripts or evidence of new numerical validation. The supplied captions are preserved as reference transcriptions with only their review links adapted. Use each figure README and the current workflow definitions when interpreting them.
