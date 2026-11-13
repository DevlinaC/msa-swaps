# Creating the collaboration repository

The local draft has not been uploaded. Suggested repository: **DevlinaC/msa-swaps**, initially private while the authors revise it. Confirm the repository name/visibility before executing creation commands.

After reviewing the files and confirming that name and visibility, run from the draft root:

```bash
git init -b main
git add .
git status --short
# Review the staged files before committing; editor state and outputs/ are ignored.
git commit -m "Assemble manuscript analysis and figure workflows"
gh repo create DevlinaC/msa-swaps --private --source=. --remote=origin --push
```

These commands assume authenticated GitHub CLI access to DevlinaC. After creation, invite Lauren's GitHub account through repository Settings → Collaborators with write access, then share the actual repository URL. A private repository link alone does not grant access. No invitation or email has been sent.

Original project code is covered by the [MIT License](../LICENSE). Before public release, finalize data/asset release terms, preserve applicable upstream notices, and complete remaining scientific sources, raw-data archive links and the confirmed publication citation. A [provisional manuscript citation and release checklist](PUBLICATION_METADATA.md) are available. See [licensing scope](LICENSING.md). The figure index identifies which results are already reproducible. This collaboration draft can be reviewed before every publication component is complete.

## Temporary exclusions before deciding what to share

The current `.gitignore` excludes local environments/caches, JSON/JSONL/JSON5 files except the [15 reviewed JSON exceptions](../provenance/README.md#json-files-eligible-for-git), full contact-matrix CSVs, compressed CSVs and explicitly listed CSVs of at least 1 MiB. Small summary/metadata CSVs remain eligible. Git does not support size-based ignore rules, so the large-file list is a snapshot and must be reviewed when new files are added. Files are retained locally.

These exclusions include scientific inputs, configuration and provenance needed by the documented reproduction commands. A GitHub upload with these exclusions **will not reproduce all figures by itself**. Before publication, agree which files to track and which to supply in a versioned data archive; then update the README and data links accordingly.

To inspect an exclusion, use `git check-ignore -v path/to/file`. To intentionally include an approved file, remove its explicit exclusion or add a narrow exception after the relevant rule (for example `!workflows/pairformer/model_config.json`). `git add -f path/to/file` is also possible after review. Ignore rules do not remove already tracked files; inspect the index if using these rules in an existing repository. This local draft has not yet been initialized as a Git repository.

## Historical material kept locally

The ignore rules also exclude the October documentation/plots archive, legacy `data/targets.csv`, `data/results/`, `data/examples/` and `data/validation/`, the old story image/poster, November 12 planning snapshots, imported chimera README snapshots and superseded S7 figure exports. These files remain on disk. Historical demo commands require those local-only inputs.

Current manuscript figure exports, preserved scientific calculation sources and editable structure assets remain eligible under these new rules, including chimera PyMOL sessions whose original filenames contain S7. Existing scientific-data exclusions above still apply.
