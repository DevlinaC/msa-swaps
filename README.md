# Query sequence biases AlphaFold predictions despite strong MSA-derived structural information

This repository contains the associated files, data and scripts for the paper **“Query sequence biases AlphaFold predictions despite strong MSA-derived structural information”** by Madeleine F. Clore, Lydia C. Schwartz, Devlina Chakravarty, Suchetan Dontha, Joseph F. Thole, G. Marius Clore, and Lauren L. Porter.

> **Release status — collaboration draft (November 13, 2026).** Reproduction coverage is partial, and a GitHub checkout alone does not contain all required inputs. See [figure coverage](docs/FIGURES.md), [publication metadata and pending release items](docs/PUBLICATION_METADATA.md), and the [quick start](#quick-start-reproduce-packaged-figures).

## Overview

The study examines how AlphaFold2 (AF2) and AlphaFold3 (AF3) respond when a protein's query sequence and its multiple sequence alignment (MSA) supply conflicting structural information. **A denotes the query protein; B denotes the MSA donor.**

| Input condition | Query sequence | MSA | Scope per AF method |
| --- | --- | --- | --- |
| Native input | A | A | 34 proteins × 25 predictions |
| Swapped input | A | B | 362 directed swaps × 25 predictions = 9,050 |
| Single sequence | A | None | 34 proteins × 25 predictions |

The dataset contains 14, 10 and 10 proteins in nominal 60-, 100- and 150-residue classes. Swaps are made within each class. Contact predictions retain donor-associated information, while AF predictions often do not realize that contact pattern, supporting a strong influence of the query sequence in these experiments.

![Figure 1. Native and swapped MSA inputs and their associated contact and structure predictions](figures/Figure_1/Figure_1.png)

*Figure 1 from the supplied PowerPoint PDF export. The [figure folder](figures/Figure_1/README.md) contains the PDF, editable source and image provenance.*

## Figures and supporting data

Browse the linked figure guides for artwork, source data, scripts and provenance. **Supplied artwork** preserves the manuscript display; **reproducible plots** are generated from packaged data. The [detailed figure guide](docs/FIGURES.md) gives commands and remaining limitations.

### Main figures

| Figure | Included files and current coverage |
| --- | --- |
| [Figure 1 — study overview](figures/Figure_1/README.md) | PDF, PNG preview and editable PowerPoint; illustrative inference is not newly reproduced. |
| [Figure 2 — query/donor similarity](figures/Figure_2/README.md) | Structural-score and contact-overlap tables, original plotting scripts and PDF/PNG panels. Reproduces the supplied historical contact selection. |
| [Figure 3 — swapped versus single-sequence predictions](figures/Figure_3/README.md) | Figure 3A plots, underlying tables and plotting code. Pooled distributions reproduce; the AF2 statistical grouping correction remains pending. |
| [Figure 4 — structure/contact-map case studies](figures/Figure_4/README.md) | Full PDF/PNG and editable artwork; eight matrices, eight selected AF2 structures, four calculation scripts and verified CSV/PNG/SVG results. All four numerical comparisons reproduce. |
| [Figure 5 — mutant ensembles](figures/Figure_5/README.md) | 100 AF2 models, four A3M inputs, confidence records and configurations. Model provenance is validated; fold assignments and final percentage plots remain pending. |

The related [RfaH/NusG 512-sequence analysis](workflows/case_studies/README.md) includes four matrices, four representative AF3 structures, comparison code and a reproduced results table. Exact correspondence to final S8 panels remains unresolved.

### Supporting figures

| Figure | Included files and current coverage |
| --- | --- |
| [S1 — protein folds and structural similarity](figures/S1/README.md) | Structure renders, TM-score matrices, plotting code and combined/individual exports. Manuscript summary statistics need reconciliation. |
| [S2 — contact overlap](figures/S2/README.md) | Reproducible GREMLIN ordinary Jaccard and MSA Transformer rank-tolerant overlap plots at three sequence separations. |
| [S3 — contact-score correlations](figures/S3/README.md) | Tables and reproducible Pearson-correlation plots over top-L unions; supplied caption terminology needs reconciliation. |
| [S4 — MSA contacts versus predicted structures](figures/S4/README.md) | Tables, mapping/selection audits, calculation sources and reproducible precision/Jaccard panels. |
| [S5 — realization of donor-MSA contacts](figures/S5/README.md) | Supplied SI PDF/PNG and caption. Numerical source tables and control aggregation remain missing. |
| [S6 — MSA Pairformer comparisons](figures/S6/README.md) | Current four-panel 512-sequence/unweighted analysis, matrices, tables and plotting code; older three-panel SI artwork retained separately. |
| [S7 — GA95 and GA95-L45Y](figures/S7/README.md) | Supplied SI PNG and caption. Model ensembles and fold-count sources remain missing. |
| [S8 — AF3 BCCIP and RfaH/NusG comparisons](figures/S8/README.md) | Supplied SI PDF/PNG and caption; exact model/map identities and final panel correspondence remain unresolved. |
| [S9 — NusG–RfaH chimera](figures/S9/README.md) | Reproducible prediction panels A–C and full supplied SI artwork including NMR displays. Experimental source data and construct confirmation remain missing. |

[Table S1 and reference data](data/publication/Table_S1/) include target annotations, sequences and reference coordinates. Supplied SI captions are retained for reference; unresolved caption/artwork differences are described in the [supporting-figure review](docs/SUPPORTING_DATA_REVIEW_20261113.md).

Main structural comparisons use experimental references and a fold-recovery threshold of **TM-score ≥0.6**. S4 instead compares MSA-derived contacts with contacts measured in **rank-1 AF predictions**. The [contact methods](docs/contact_methods.md) describe contact selection, ordinary and rank-tolerant Jaccard, Pearson correlations, and the S4 comparison definitions.

The current Pairformer analysis uses up to **512 sequences with query weighting disabled**. See the [Pairformer workflow](workflows/pairformer/REPOSITORY_NOTES.md) for its inputs, matrices and dependencies. The [structural-score table](data/publication/structure_scores/README.md) contains both query and donor TM-scores for each swapped prediction, with normalization and provenance notes.

## Quick start: reproduce packaged figures

Use Python 3.11 from the repository root. These commands use the inputs included in the Git package:

```bash
python -m pip install -r requirements-structure.txt
python scripts/validate_collaborator_import.py
python scripts/reproduce_case_studies.py
python scripts/plot_main_figures.py
python scripts/reproduce_figures.py --figure all
```

They validate the imported models and confidence records, reproduce five case-study metric tables, plot supplied Figures 2/3A, and regenerate S1–S4, current S6 and S9 prediction panels A–C. They use saved results, not fresh inference. The case-study and main-figure runners require a fresh `--out-dir` on repeat runs; supplementary plotting writes into the documented figure folders.

Figure 2 retains historical contact selection. Figure 3 statistical grouping, Figure 5 fold assignments and the remaining scientific gaps are documented in the [figure guide](docs/FIGURES.md). Supplied full-layout artwork is distinct from automatically regenerated panels.

The following additional commands require the **complete local package**, including records excluded from Git; approved archive download links are not yet available:

```bash
python scripts/validate_publication.py
python scripts/summarize_swapped_structures.py
```

See the [data guide](data/README.md), [upload policy](docs/GIT_DATA_POLICY.md) and [pipeline guide](docs/pipeline.md) for the distinction between packaged plotting, raw-data reanalysis and model inference.

![Figure S2. Contact overlap at three sequence separations](figures/S2/Figure_S2_contact_overlap.png)

## Repository contents

```text
scripts/                 Figure generation, score summaries and validation
workflows/               GREMLIN, MSA Transformer, Pairformer and chimera workflows
data/publication/        Scientific tables, target annotations and coordinates
figures/                 Figure exports
provenance/              Input checksums and preserved source scripts
docs/                    Analysis methods and reproduction instructions
```

Historical example calculations and data are retained locally and excluded from the GitHub draft. They use their original definitions; see the [contact methods](docs/contact_methods.md) before comparing them with the current analyses.

## NIH Biowulf

Run inference, substantive analysis and environment installation on an allocated compute node, with persistent files under `/data`. See [Biowulf guidance](docs/BIOWULF.md) for resource allocation and execution notes.

## Software attribution

This project uses [GREMLIN](https://github.com/sokrypton/GREMLIN), [MSA Transformer / ESM](https://github.com/facebookresearch/esm), [MSA Pairformer](https://github.com/yoakiyama/MSA_Pairformer), [ColabFold](https://github.com/sokrypton/ColabFold) and [AlphaFold 3](https://github.com/google-deepmind/alphafold3). Model weights are obtained separately from their upstream providers.

## Citation and data availability

Provisional citation, based on the November 12, 2026 manuscript draft:

Clore, M. F., Schwartz, L. C., Chakravarty, D., Dontha, S., Thole, J. F., Clore, G. M., and Porter, L. L. (2026). *Query sequence biases AlphaFold predictions despite strong MSA-derived structural information*. Manuscript in preparation.

Download the [BibTeX citation](CITATION.bib). Publication status, venue and DOI remain to be confirmed; this is not a citation to a published article. Permanent data-archive links and data/figure release terms are also pending. The [publication metadata record](docs/PUBLICATION_METADATA.md) documents the source and remaining decisions. Current local data coverage is described in the [data guide](data/README.md).

## Development disclosure

OpenAI Codex was used as a development assistant for analysis and figure-generation code, code review, path handling, refactoring, documentation, repository organization, and validation. This assistance extended the project authors' existing scripts and workflows. The project authors and maintainers remain responsible for reviewing the code, verifying the analyses and figures, and the scientific interpretation of the results.

## License

Original project scripts and general-purpose code are released under the [MIT License](LICENSE).

This license does not automatically apply to the manuscript, figure compositions, scientific datasets, protein structure files, model weights, or third-party code and assets. Those materials retain their applicable terms and attribution requirements. Upstream-derived code retains its original license and notices; the project MIT License does not replace them. See [licensing scope](docs/LICENSING.md) for details.
