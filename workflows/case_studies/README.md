# Figure 4 and RfaH/NusG subsample comparisons

The collaborator delivery of November 13, 2026 supplies four Figure 4 comparisons and one RfaH/NusG subsampling comparison. Original scripts, coordinate files, contact matrices and expected metric tables are preserved byte-for-byte under `sources/`. Their source paths and hashes are recorded in [the import manifest](../../provenance/collaborator_import_20261113.csv).

## Replay

From the repository root, using Python 3.11:

```bash
python -m pip install -r requirements-structure.txt
python scripts/validate_collaborator_import.py
python scripts/reproduce_case_studies.py
```

The runner authenticates imported files, copies them to temporary working directories, bounds numeric library threads to one, runs the supplied scripts unchanged, and checks all five metric tables. It requires identical row order, column names and text labels; numeric cells use relative tolerance 1e-10 and absolute tolerance 1e-12. Expected CSV outputs are removed from the temporary copies before execution so an unwritten output cannot pass. Results, logs, validation and environment information are written to `outputs/case_studies/`; original sources remain unchanged. Select `--case figure4` or `--case rfah_nusg_512` for individual packages. For a repeat run, choose a fresh `--out-dir`.

On Biowulf, execute on an allocated compute node and set `TMPDIR` to requested local scratch. This is a small validation/replay task; do not launch one production job per comparison.

## Scientific definitions

These scripts compare **structures with structures** and **MSA-derived maps with MSA-derived maps**. They do not calculate direct MSA-to-structure agreement. That separate comparison is covered by manuscript S4.

| Package | Inputs and interpretation |
| --- | --- |
| GB98 T25V/T25E | 56 residues; structure ordinary Jaccard/binary Pearson versus MSA rank-tolerant Jaccard/continuous Pearson. The original T25V matrix filename is `T52V.txt`; it is retained exactly. |
| RfaH/NusG | Pairs with at least one residue in CTD 110–162 and sequence separation ≥6; top L uses 53 residues. The Figure 4 script uses `RfaH_unbiased_predicted_cb_contacts.txt`. |
| BCCIP alpha/beta | Maps of size 322/314 are registered through 166 shared nongap positions. Metrics describe that registered subset, not every residue of the full structures. |
| DZZB/mk2h | Structural intra- and interchain contacts remain separate channels; the supplied predictor maps are 43×43. |
| RfaH/NusG 512 | Two matrices and two representative AF3 structures per target; 12 comparison rows. This package does not establish five evaluated contact-map replicates per target. |

The Figure 4 representatives are AF2 structures; the subsample representatives are AF3. Original inference settings, selected alignment rows and exact final S8 assignments remain incomplete. Successful downstream replay does not establish full inference reproducibility. Current Figure 4 editable layouts and PyMOL assets are in [the figure source folder](../../figures/Figure_4/README.md).

The generated table images retain the supplied display conventions. In particular, the GB98 image contains numeric cells without standalone row/column labels; use its labelled CSV and the definitions above when interpreting it outside the combined layout.
