# Figure 5 AF2 mutant inputs and model records

This package contains the four mutant ensembles supplied on November 13, 2026: GB98 T25V-I17T, T25E-K28I, T25E-W43E and RfaH L142K-H152V. Each has 25 AF2 predictions (five model networks × five seeds), for 100 models total. No inference or fold classification was rerun.

| File or directory | Contents |
| --- | --- |
| `inputs/`, `inputs.csv` | Four unchanged A3Ms, query sequences, source tags, archive members and checksums. Renamed A3M paths identify constructs without changing their contents. |
| `models/`, `models.csv` | 100 unchanged PDBs, model/seed/rank identities, source members, model/score hashes and mean pLDDT across all residues. |
| `residue_plddt.csv` | All 8,250 pLDDT values copied from the corresponding score JSONs, indexed by model and 1-based query position. |
| `run_config.csv` | Original archive configuration keys and serialized values. Values are recorded, not adopted as current runtime recommendations. |
| `source_archives.csv` | Archive filenames, sizes and SHA-256 hashes; original ZIPs remain outside the repository. |

Source score JSONs contained only the `plddt` array. Their values are retained in the CSV; their original hashes and member names are in `models.csv`. The mean is the arithmetic mean over every supplied residue; it is not a CTD-specific confidence filter. Coordinates and A3Ms retain original bytes. Run `python scripts/validate_collaborator_import.py` from the repository root to check source hashes, model/seed coverage, query/coordinate sequence correspondence and all residue confidence values against PDB B-factors (tolerance 0.011 for field rounding).

The RfaH archive uses the internal tag `RfaH_MSA_94a41`; the package maps it to `RfaH_L142K_H152V`. The 162-residue input has K142 and V152. The GB98 inputs have the stated mutant residues. These checks establish the supplied sequences, not the complete provenance of their parental design or alignment selection.

**Fold assignments and inclusion in fold percentages remain unassigned/not evaluated.** Do not derive final Figure 5 percentages from rank, confidence or filenames. Final model-level fold labels, analyzed regions, confidence thresholds and denominators are still needed. These four mutant ensembles do not supply the 20 feature-selected ensembles for Figure 5B or the missing AF3 counterpart classifications.

The supplied bar-plot script hard-codes percentages; it has not been promoted as an active numerical analysis. PyMOL sessions, logs, original ZIPs and individual JSON files from the mutant archives were not copied wholesale.

## Direct archive review

The [November 13 source review](../../../provenance/figure5_source_review_20261113.csv) rechecked the original ZIP hashes, all 100 PDB and four A3M byte streams, every score JSON hash and all 8,250 confidence values. Both original configuration objects match the CSV representation exactly (32 keys each). All four ensembles passed. No extra per-model JSON copies are required; the CSV preserves all supplied score content. The source archives are needed only to repeat this archive comparison, not to run the packaged validator.
