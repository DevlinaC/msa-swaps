# Provenance records

`import_manifest_20261113.json` records byte-for-byte source imports before packaging edits. `protected_data_sha256.json` locks scientific input files and saved validation records; narrative README files are excluded because numbering and links were updated. New plot provenance identifies adapted scripts and input hashes. `packaging_environment.json` records the actual local environment used for validation and plotting.

Source edits during packaging: new root documentation and figure commands; path-only adapters for S1–S4; S7→S9 labels/filenames in active chimera exporters and documentation; updated chimera export timestamps; additive structural-score reconciliation and validation. No model inference, structural alignment or scientific contact calculation was rerun. Original source snapshots remain under `source_scripts/`, with imported chimera READMEs under that workflow’s provenance folder.

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
