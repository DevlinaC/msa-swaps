# Collaborating on this manuscript repository

See [collaborator notes and repository status](docs/COLLABORATOR_NOTES.md) for the manuscript version, outstanding files, temporary exclusions and release decisions.

Start with `docs/FIGURES.md` for current numbering and coverage. The old local chimera folder was called Figure_S7; it is S9 in the current manuscript. S7 now refers to GA95/GA95-L45Y.

For each figure addition, provide:

1. Source table with target/donor, condition, AF method, model/seed/sample and inclusion/exclusion fields where applicable.
2. Calculation and plotting code, exact command, environment versions and required inputs.
3. Final PNG plus editable SVG/PDF or layout project; for structure cartoons, include selected coordinates, residue selections and PyMOL scene/session.
4. A short README stating metric definitions, ranking/selection rules, panel mapping and any manual assembly steps.

Do not silently replace an analysis with a different threshold, MSA depth, weighting setting or residue mapping. Write new results under a distinct version, document the change and update the figure index. Preserve original model and table identifiers. Keep credentials, personal editor state and private sharing links out of Git.

For numerical changes, validate counts, unique keys, selection/mapping and a representative independent calculation. For plotting-only changes, confirm the input tables remain unchanged and inspect the regenerated output.

Full raw AF ensembles and other bulky datasets should go in an agreed versioned data archive; put the archive manifest, checksums and permanent links here. Original project code uses the [MIT License](LICENSE). Preserve upstream licenses and notices, and document separate data/asset terms as described in [licensing scope](docs/LICENSING.md).

Suggested workflow: create the repository under DevlinaC (proposed name `msa-swaps`), keep it private during co-author revision, invite Lauren's GitHub account with write access, and use a branch/pull request for scientific changes. No repository or collaborator invitation has been created by this local update.
