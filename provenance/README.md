# Provenance records

`import_manifest_20261113.json` records byte-for-byte source imports before packaging edits. `protected_data_sha256.json` locks scientific input files and saved validation records; narrative README files are excluded because numbering and links were updated. New plot provenance identifies adapted scripts and input hashes. `packaging_environment.json` records the actual local environment used for validation and plotting.

Source edits during packaging: new root documentation and figure commands; path-only adapters for S1–S4; S7→S9 labels/filenames in active chimera exporters and documentation; updated chimera export timestamps; additive structural-score reconciliation and validation. During the initial packaging step, no model inference, structural alignment or scientific contact calculation was rerun. Subsequent collaborator case-study integration replayed the supplied contact calculations, as documented below. Original source snapshots remain under `source_scripts/`, with imported chimera READMEs under that workflow’s provenance folder.

## JSON files eligible for Git

JSON, JSONL and JSON5 remain ignored by default. The explicit exceptions in `.gitignore` cover only the following 15 JSON files (756,984 bytes total at review on November 13, 2026):

| Files | Why included |
| --- | --- |
| `protected_data_sha256.json` in this folder | Required by `scripts/validate_publication.py`; preserves the full-package integrity checks. |
| `data/publication/S2/GREMLIN/{provenance,COMPLETE}.json` | Required to authenticate the S2 GREMLIN plotting table. |
| `data/publication/S2/MSA_Transformer/{provenance,COMPLETE}.json` | Required to authenticate the S2 MSA Transformer plotting table. |
| `data/publication/S3/{provenance,COMPLETE}.json` | Required to authenticate the S3 plotting table. The shared S2/S3 script loads all three datasets. |
| `workflows/pairformer/model_config.json` and `results/INFERENCE_COMPLETE.json` | Required by publication validation to check the selected run settings and completion count. |
| `workflows/pairformer/figure_palette.json`, `results/analysis/VALIDATED.json` and `results/rank_tolerant/VALIDATED.json` | Required by the Pairformer figure generator. |
| `workflows/chimera/data/models.json` and `summary.json` | Required for S9 panels and model-count validation. |
| `figures/Figure_1/provenance.json` | Records the manuscript-derived preview's source and checksums. |

All other JSON files remain excluded, including per-target inference records, bulk manifests, historical metadata and generated export reports. The environment and import manifests described above are currently local-only.

This is the minimum direct-input set for the advertised plotting/validation entry points, plus Figure 1 attribution; it is **not a complete dependency package for full validation or reanalysis**. The protected checksum manifest references 424 JSON files, many of them per-target inference inputs, as well as other excluded data. Full publication validation still requires the complete local package or a future versioned archive. No checks are skipped, and no scientific files or hashes were rewritten to accommodate upload exclusions.

Some preserved records contain original local execution paths. These are provenance labels, not portable input locations. New optional JSON files require an explicit review and exception before inclusion.

## Verified collaborator additions

`collaborator_import_20261113.csv` records byte hashes for 148 unchanged imports from the shared folder and selected archive members. `collaborator_derived_sha256.csv` locks the five newly derived mutant CSVs. `collaborator_case_validation_20261113.csv` records the five successful numerical table comparisons; `collaborator_case_environment_20261113.txt` records the replay environment and script/manifest hashes. Run `scripts/validate_collaborator_import.py` to verify the imported bytes and mutant sequence/confidence correspondence. These CSV/text records are eligible for Git without expanding the JSON allowlist.

`figures_2_3_sources.csv` authenticates the supplied main-figure plotting inputs/scripts and records duplicate source locations. `figures_2_3_populations.csv`, `figures_2_3_exports.csv` and `figures_2_3_environment.txt` record population checks, frozen plot hashes and rendering versions. These records concern plotting only, not corrected Figure 3 statistics.

`figure4_exports.csv` records the unchanged user-supplied PDF, its PNG render, eight table images and the Figure 4 validation record. `figure4_validation_20261113.csv` records all four passing comparisons; the shared `collaborator_case_environment_20261113.txt` captures the replay versions and runner/input-manifest hashes. Original matrices, PDBs, scripts and reference tables remain covered by the collaborator import manifest.

`figure5_source_review_20261113.csv` records the direct original-ZIP review of the four mutant ensembles, including configuration source hashes and model/confidence counts. All source bytes and score values match the existing packaged data; no fold classifications were inferred.

The local-only, Git-ignored `supporting_figure_comparison_20261113.csv` inventories source figure-folder files against repository hashes before import. `supporting_figure_import_20261113.csv` records the selected SI reference artwork and caption imports, source document hash, and original/packaged file hashes. See the [review](../docs/SUPPORTING_DATA_REVIEW_20261113.md) for scope and unresolved differences.

The [scientific upload policy](../docs/GIT_DATA_POLICY.md) also documents exact CSV/PDB/CIF exceptions and retained structural-reference files. Required S3/S4 plotting tables are included; the full protected-source audit still requires excluded local records.

`archive_manifest_20261113.csv` records redundant or superseded files moved to the ignored local archive, their unchanged hashes, and canonical replacements. Original protected-source checksums remain unchanged; full local validation resolves archived paths through this manifest. Archived files are not Git dependencies for figure plotting.
