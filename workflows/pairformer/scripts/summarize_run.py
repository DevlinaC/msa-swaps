#!/usr/bin/env python3
"""Write a concise run summary from validated outputs, without inferred values."""
import csv,json
from pathlib import Path

r=Path(__file__).resolve().parent.parent
with (r/'results/analysis/summary_with_rank_tolerance.csv').open() as f: summary=list(csv.DictReader(f))
lookup={(x['separation'],x['reference'],x['metric']):float(x['mean']) for x in summary}
with (r/'results/analysis/comparisons_with_rank_tolerance.csv').open() as f: comparisons=list(csv.DictReader(f))
metrics=['ordinary_jaccard','tolerant_jaccard','pearson_all_eligible','pearson_union_top_L']
labels=['Ordinary top-L Jaccard','Rank-tolerant top-L Jaccard','Pearson, all eligible pairs','Pearson, union of top-L pairs']
report=json.loads((r/'results/configuration_comparison/VALIDATED.json').read_text())
audit=json.loads((r/'results/WEIGHTING_AUDIT.json').read_text())
inference=json.loads((r/'results/INFERENCE_COMPLETE.json').read_text())
caption=(
 'MSA-Pairformer contact-map similarity for 362 swapped MSAs. Predictions used '
 'up to 512 sequences per MSA with query weighting disabled (uniform sequence weights). '
 'The query sequence was retained, and deeper alignments underwent deterministic '
 'diversity selection without HHfilter. Each swapped map was compared with its '
 'proper query (red) and donor (gold) reference maps at matched residue positions. '
 'Panels show (A) ordinary top-L Jaccard, (B) rank-tolerant top-L Jaccard, '
 '(C) Pearson correlation over all eligible residue pairs, and (D) Pearson '
 'correlation over the union of the two top-L sets, at minimum sequence separations '
 'of 6, 12, and 24 residues. L is the length of the common registered overlap; '
 'rank tolerance credits a top-L pair found within the other map’s top floor(1.5L) '
 'pairs, using the strict top-L union as denominator. Cutoff ties are retained. '
 'Each point represents one swap (362 per group); jitter is horizontal only, '
 'and black bars mark medians.'
)
(r/'results/analysis/figure_legend.txt').write_text(caption+'\n')
lines=['# MSA-Pairformer manuscript matched swap analysis','',
'**Current manuscript configuration: up to 512 sequences per MSA, query weighting disabled.**',
'This replaces the earlier 256-sequence, query-biased swap analysis for manuscript use.',
'The earlier matrices and figures remain in `../msa_pairformer_full/` for provenance and comparison.','',
'## What was corrected','',
'We reran all **34 proper and 362 swapped MSAs** with the same pinned pretrained Cβ',
'predictor, explicitly calling `turn_off_query_biasing()`. Every outer-product layer',
'was checked to have query attention disabled. Returned sequence weights were',
'checked against uniform weights for three small inputs and a largest 512-row input.',
'No biased pilot matrix was reused as a new result. An old biased control was replayed',
f"before switching off bias (maximum absolute difference {audit['control_biased_max_abs_difference']:.6g}).",'',
'We retained the query first and used deterministic diversity selection to a maximum',
'of 512 sequences, without HHfilter. Shallower alignments retained all rows.',
f"Selected depth increased for **{report['depth_increased']}/396 inputs**; every prior selected",
'row is contained in the new subset. Original source MSAs, corrected donor-2mm4 inputs,',
'and homolog-based residue registration are unchanged. The five random RfaH/NusG',
'subsamples described in the manuscript are a separate case-study protocol.','',
'## Results','',
'Means over 362 swaps; each cell is **query / donor** similarity.','',
'| Metric | ≥6 | ≥12 | ≥24 |','| --- | ---: | ---: | ---: |']
for metric,label in zip(metrics,labels):
 lines.append('| '+label+' | '+' | '.join(f"{lookup[s,'query',metric]:.3f} / **{lookup[s,'donor',metric]:.3f}**" for s in ['6','12','24'])+' |')
lines+=['','Number of swaps with donor similarity greater than query similarity:','','| Metric | ≥6 | ≥12 | ≥24 |','| --- | ---: | ---: | ---: |']
for metric,label in zip(metrics,labels):
 counts=[]
 for sep in ['6','12','24']:
  sets={ref:{x['id']:float(x[metric]) for x in comparisons if x['separation']==sep and x['reference']==ref} for ref in ['query','donor']}
  counts.append(sum(sets['donor'][i]>v for i,v in sets['query'].items()))
 lines.append('| '+label+' | '+' | '.join(f'{n}/362' for n in counts)+' |')
lines+=['','Rank-tolerant Jaccard credits a top-L pair if it occurs within the other map’s',
'top floor(1.5L) pairs, with the strict top-L union as denominator. L is the common',
'homolog-registered overlap. All comparisons use the upper triangle at the specified',
'minimum sequence separation; score ties are retained.','',
'All 396 matrices and 2,172 comparison records passed validation. The numerical',
'analysis was independently replayed, including checks against the historical',
'rank-tolerant metric implementation. These comparisons measure similarity between',
'predicted maps, not structural accuracy or AlphaFold’s internal representations.',
'The swaps reuse proteins and are not independent biological replicates.','',
'## Change from the earlier run','',
'| Donor metric | ≥6 old → new | ≥12 old → new | ≥24 old → new |','| --- | ---: | ---: | ---: |']
for metric,label in zip(metrics,labels):
 cells=[]
 for sep in ['6','12','24']:
  x=next(x for x in report['summary'] if str(x['separation'])==sep and x['reference']=='donor' and x['metric']==metric)
  cells.append(f"{x['old_mean']:.3f} → {x['new_mean']:.3f}")
 lines.append('| '+label+' | '+' | '.join(cells)+' |')
lines+=['',f"**{report['matrices_changed']}/396 matrices changed.** Depth and weighting changed together;",
'we cannot attribute the score differences separately to those two factors.',
'Paired values and per-matrix differences are in `results/configuration_comparison/`.','',
'For donor 2mm4 specifically, all 13 swapped inputs and the proper donor reference',
'retain the exact same 18-row alignments in both runs. Its mean donor ordinary',
'Jaccard at ≥6 increases from 0.172 to 0.390, and rank-tolerant Jaccard from',
'0.275 to 0.542. Thus query weighting contributed to that earlier low tail even',
'without a depth change; the causal effect of shallow MSA depth itself remains untested.','',
'## Figures and data','',
'![Current four-panel analysis](results/analysis/contact_comparison.png)','',
'**Figure legend.** '+caption,'',
'Configuration details are stated in this legend rather than in the plot titles.','',
'- [Figure legend text](results/analysis/figure_legend.txt)',
'- [Main figure PDF](results/analysis/contact_comparison.pdf) · [Editable SVG](results/analysis/contact_comparison.svg)',
'- [Upper row PDF](results/analysis/panels/upper_row_jaccard.pdf) · [Lower row PDF](results/analysis/panels/lower_row_pearson.pdf)',
'- Individual panels A–D, PNG/PDF/SVG: `results/analysis/panels/`.',
'- [All four metrics](results/analysis/comparisons_with_rank_tolerance.csv) · [Summary](results/analysis/summary_with_rank_tolerance.csv)',
'- [Paired configuration comparison](results/configuration_comparison/summary.csv)',
'- [Collaborator guide](COLLABORATOR_GUIDE.md) · [Pipeline details](PIPELINE_DETAILS.md)',
'- [Swapped matrix index](inventory/swapped_matrices.tsv) · [All matrix index](inventory/all_matrices.tsv)',
'- [Weighting audit](results/WEIGHTING_AUDIT.json) · [Protocol reconciliation](provenance/protocol_reconciliation.json)','',
'Biowulf: `/data/chakravartyd2/for_lydia_project/msa_pairformer_512_unweighted/`.',
f"All {inference['new']} new matrices were generated in GPU job {inference['job_id']}.",
'Regenerate figures from saved tables with `python scripts/plot_summary.py`.',
'All scripts are under `scripts/`; reading the matrices does not require a GPU.','']
(r/'README.md').write_text('\n'.join(lines))
print('Wrote README from validated summaries.')
