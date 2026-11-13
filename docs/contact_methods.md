# Contact definitions and comparison key

A denotes the query sequence and B the MSA donor. Sequence separation |i−j| ≥6, ≥12 or ≥24 is distinct from a spatial contact-distance threshold.

## S2

GREMLIN uses ordinary set Jaccard, gap fraction <0.25, z-score ≥1 and caps of 100/100/150 for the nominal 60/100/150-residue classes. The original `gremlin_jaccard_v2.py` fixed `top_n=100` for every class; that historical rule is preserved in the legacy examples. S2 changes the largest class cap only and retains the historical donor-offset search. It is not a top-L analysis.

MSA Transformer uses rank-tolerant overlap of top-L contacts. Cutoff ties are retained. The historical comparison center-crops to a common length; this differs from the homolog registration in S4 and Pairformer. For strict sets S,T and expanded top-floor(1.5L) sets Sw,Tw:

- Ordinary Jaccard = |S ∩ T| / |S ∪ T|.
- Rank-tolerant overlap = |(S ∩ Tw) ∪ (T ∩ Sw)| / |S ∪ T|.

Rank tolerance matches the same residue pair at a wider score-rank cutoff; it does not match neighboring positions. Scores using different definitions are not directly interchangeable.

## S3

Pearson R compares continuous scores at the same residue pairs over the union of the two top-L contact sets, for each swapped-versus-native comparison. It is not a correlation of sorted score lists and is not “rank-tolerant Pearson.” Undefined correlations are excluded with recorded counts. Selecting a union of high-scoring pairs can induce negative query correlations; these do not alone establish biological anticorrelation.

## S4

Primary top-L contacts retain cutoff ties, with L the mapped comparison length. GREMLIN uses its eligible gap-filtered score universe without the S2 z-score floor. The cap-based/z-score-filtered calculations in the same data folder are separate sensitivity results. The plot explicitly selects `selection == topL`.

Structural contacts are Cβ–Cβ distances ≤8 Å, with Cα for glycine, filtered by sequence separation. Only mapped eligible residue pairs with observable distances enter the comparison. Donor registration is inferred from shared homolog rows, independently of contact overlap. AF2 rank 1 uses pLDDT; AF3 uses ranking_score among the sampled predictions.

| Color | Contacts inferred by GREMLIN/MSA Transformer from | Structural contacts in rank-1 AF prediction from |
| --- | --- | --- |
| Light gray | Query A sequence + query A MSA | Query A sequence + query A MSA |
| Red | Query A sequence + donor B MSA | Query A sequence + query A MSA |
| Gold | Query A sequence + donor B MSA | Donor B sequence + donor B MSA |
| Dark gray | Query A sequence + donor B MSA | Query A sequence + donor B MSA |

Every S4 reference structure is a prediction, not an experimental structure. Dark gray asks whether the AF output realizes contacts inferred from the same swapped input used to generate it. Each native baseline has 34 observations; every swapped distribution has 362.

Precision = shared contacts / selected eligible MSA-derived contacts. Jaccard = shared contacts / union of selected MSA-derived and structural contacts in the eligible universe. S4A shows precision; S4B uses ordinary Jaccard for both methods.

## Pairformer

The current run uses 512 rows maximum, deterministic diversity selection, no HHfilter and query weighting disabled (uniform sequence weights). All 396 selected inputs and matrices are provided. Homolog-registered common positions, top-L cutoff ties and the same ≥6/12/24 thresholds are used for ordinary Jaccard, rank-tolerant overlap, all-pair Pearson and top-L-union Pearson. The current run and historical 256-row/bias-enabled run changed both depth and weighting; their differences cannot be attributed to either factor alone across the full dataset.

## Interpretation

S2/S3 and Pairformer measure similarity between inferred contact maps. S4 connects those contacts to AF-predicted geometry. Donor-like contacts persist but often are poorly realized in swapped-input AF outputs. These external readouts do not reveal AF's internal weighting of MSA information.
