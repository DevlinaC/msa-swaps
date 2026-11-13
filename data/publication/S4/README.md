# Item 3 secondary GREMLIN cap revision

The original supplied GREMLIN Jaccard code used `GREMLIN_TOP_N = 100`, which
was also used for item 3's secondary z-score ≥1 comparison. This revision
raises the secondary cap to **150 only for the 150-residue class**, retaining
the z-score threshold, gap filtering, coordinate mapping, and deterministic
tie breaking. The 60- and 100-residue classes remain capped at 100.
The revised secondary selection is named `class_cap_z1` in the tables.

The primary `topL` analysis is unchanged: it already selected 150 contacts in
all 1,680 GREMLIN structural comparisons for the 150-residue class. All 13,440
primary rows across both contact methods are identical to the original results.
The other 5,040 secondary rows also retain their numeric results.

The revision reuses the validated score/distance records exported for the
primary top-150 contacts. Before increasing the cap, the script reconstructed
and verified all 1,680 affected original top-100 selections. It then retained
up to 150 pairs with z-score ≥1. Actual revised counts range from 13 to 150;
retaining the confidence filter does not guarantee 150 eligible pairs.

All 20,160 precision values remain defined. Independent replay of 1,759,566
exported selected pairs reproduced precision to 1.11×10⁻¹⁶, and 6,516 checks
confirmed identical selected pairs across matched structural references.
Input compatibility sensitivity analyses were regenerated with the same
three flagged AF2 native inputs as before.

The complete revision and validation ran together on Biowulf cn0147 as job
34846045, finishing in 32 seconds (4 GB requested; 298 MB reported usage).
No AlphaFold/contact inference or duplicate structure calculation was needed.
The original validated outputs remain in `../full-34533285/`.

`provenance.json` and `COMPLETE.json` record the revision and its source.
`VALIDATED.json` authenticates the independent validation outputs.
`cap_revision_audit.csv` records old/new counts and precision for affected rows.
The revision script is `../../../codes/revise_structure_contact_cap.py`.
