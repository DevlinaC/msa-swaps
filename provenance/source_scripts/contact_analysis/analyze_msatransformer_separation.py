#!/usr/bin/env python3
"""Extend manuscript Figure 2B's rank-tolerant Jaccard to three separations.

Uses saved contact matrices; does not rerun contact inference or AlphaFold.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/lydia_report_mpl')
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'for_github_repo/analysis/msatransformer_separation'
FIG = ROOT / 'reports/figure2_msatransformer_jaccard'
SEPARATIONS = (6, 12, 24)
FIXED = {'1bor', '1bw6', '1cqu', '1gb1', '1vex', '2fs1', '2nc2', '2ovo'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(a, b, sep):
    length = min(len(a), len(b))
    ca, cb = (len(a)-length+1)//2, (len(b)-length+1)//2
    i, j = np.triu_indices(length, k=sep)
    x, y = a[ca+i, ca+j], b[cb+i, cb+j]
    wide = int(1.5*length)
    assert len(x) >= wide
    tx, ty = np.sort(x)[-length], np.sort(y)[-length]
    wx, wy = np.sort(x)[-wide], np.sort(y)[-wide]
    A, B, AW, BW = x >= tx, y >= ty, x >= wx, y >= wy
    credited = (A & BW) | (B & AW)
    union = A | B
    # Independently verify the Boolean arithmetic using residue-pair sets.
    sets = [set(zip(i[m].tolist(), j[m].tolist())) for m in (A, B, AW, BW)]
    sa, sb, saw, sbw = sets
    assert len((sa & sbw) | (sb & saw)) == int(credited.sum())
    assert len(sa | sb) == int(union.sum())
    return dict(L=length, wide_k=wide, crop_swapped=ca, crop_native=cb,
                n_eligible=len(x), n_top_swapped=int(A.sum()), n_top_native=int(B.sum()),
                n_wide_swapped=int(AW.sum()), n_wide_native=int(BW.sum()),
                n_credited=int(credited.sum()), n_union=int(union.sum()),
                n_intersection=int((A & B).sum()),
                tolerant_jaccard=float(credited.sum()/union.sum()),
                strict_jaccard=float((A & B).sum()/union.sum()))


def main():
    assert not OUT.exists(), f'Output already exists: {OUT}'
    sources = {}
    prior = json.loads((ROOT/'for_github_repo/analysis/pearson_topl/provenance.json').read_text())

    def track(path):
        key = str(path.relative_to(ROOT))
        sources[key] = dict(sha256=sha(path), bytes=path.stat().st_size)
        return path

    def matrix(path):
        path = track(path)
        key = str(path.relative_to(ROOT))
        assert sources[key]['sha256'] == prior['source_files'][key]['sha256'], key
        a = np.loadtxt(path, delimiter=',')
        assert a.ndim == 2 and a.shape[0] == a.shape[1] and np.isfinite(a).all()
        assert np.max(abs(a-a.T)) <= 1e-6
        return a

    targets = pd.read_csv(track(ROOT/'for_github_repo/msa-swaps-github/data/targets.csv'))
    classes = dict(zip(targets.pdb_id.str.lower(), targets.length_class))
    saved = pd.read_csv(track(ROOT/'for_github_repo/analysis_csvs/jaccard_results_msatransformer.csv'))
    expected = {(a,b) for a in classes for b in classes if a != b and classes[a] == classes[b]}
    assert len(saved) == 362 and set(zip(saved.Protein_A, saved.Protein_B)) == expected
    legacy_path = track(ROOT/'for_github_repo/TM_Jaccard_codes/Jaccard_code/msa_transformer_jaccard_v3.py')
    spec = importlib.util.spec_from_file_location('legacy_msatr', legacy_path)
    legacy = importlib.util.module_from_spec(spec); spec.loader.exec_module(legacy)
    native = {p:matrix(ROOT/f'original_proper_msas_af_sampling/{p}/full_{p}.csv') for p in sorted(classes)}
    records, baseline = [], []
    for row in saved.sort_values('Aligned_File').itertuples(index=False):
        a, b, pair = row.Protein_A, row.Protein_B, row.Aligned_File
        fixed = b == '2mm4' and a in FIXED
        path = ROOT/(f'2mm4_fixes/{pair}/full_{pair}.csv' if fixed else
                     f'original_junk_msas_af_sampling/{a}/{pair}/full_{pair}.csv')
        swapped = matrix(path)
        for ref, p in (('A',a), ('B',b)):
            for sep in SEPARATIONS:
                result = compare(swapped, native[p], sep)
                records.append(dict(pair=pair, protein_a=a, protein_b=b, reference=ref,
                                    length_class=classes[a], separation=sep, corrected=fixed, **result))
                if sep == 6:
                    saved_value = getattr(row, f'Jaccard_vs_Protein{ref}')
                    original = legacy.tolerant_jaccard(swapped, native[p])
                    error = max(abs(result['tolerant_jaccard']-saved_value), abs(result['tolerant_jaccard']-original))
                    assert error <= 1e-12, (pair, ref, error)
                    baseline.append(dict(pair=pair, reference=ref, saved=saved_value,
                                         reproduced=result['tolerant_jaccard'], max_abs_error=error))
    frame = pd.DataFrame(records)
    assert len(frame) == 2172 and not frame.duplicated(['pair','reference','separation']).any()
    assert frame.tolerant_jaccard.between(0,1).all()
    # Confirm strict top-L membership agrees with the saved Pearson selection.
    pearson = pd.read_csv(track(ROOT/'for_github_repo/analysis/pearson_topl/comparisons.csv'))
    pearson = pearson[pearson.method == 'MSA Transformer']
    joined = frame.merge(pearson, on=['pair','reference','separation'], suffixes=('_new','_prior'), validate='one_to_one')
    assert len(joined) == 2172
    for key in ('L','n_top_swapped','n_top_native','n_union','n_intersection'):
        assert (joined[key+'_new'] == joined[key+'_prior']).all(), key
    summary = []
    for (sep, ref), group in frame.groupby(['separation','reference']):
        for cls in ('all',60,100,150):
            sub = group if cls == 'all' else group[group.length_class == cls]
            v = sub.tolerant_jaccard
            summary.append(dict(separation=sep, reference=ref, length_class=cls, n=len(sub),
                                mean=v.mean(), median=v.median(), q25=v.quantile(.25),
                                q75=v.quantile(.75), min=v.min(), max=v.max()))
    OUT.mkdir(parents=True); FIG.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT/'comparisons.csv',index=False)
    pd.DataFrame(summary).to_csv(OUT/'summary.csv',index=False)
    pd.DataFrame(baseline).to_csv(OUT/'baseline_validation.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig, ax = plt.subplots(figsize=(4.6,3.8))
    values = [frame.loc[(frame.reference=='B')&(frame.separation==s),'tolerant_jaccard'].to_numpy() for s in SEPARATIONS]
    parts = ax.violinplot(values, positions=[1,2,3], widths=.72, showmedians=True, showextrema=True, bw_method='scott')
    for body in parts['bodies']:
        body.set_facecolor('#DAA520'); body.set_edgecolor('#333333'); body.set_alpha(1); body.set_linewidth(.7)
    for key in ('cmedians','cmins','cmaxes','cbars'):
        parts[key].set_color('#333333');parts[key].set_linewidth(.9)
    ax.set_xticks([1,2,3], ['≥6','≥12','≥24'])
    ax.set_xlabel('Minimum sequence separation |i − j|')
    ax.set_ylabel('Rank-tolerant Jaccard')
    ax.set_ylim(-.02,1.08);ax.set_yticks(np.arange(0,1.01,.2))
    ax.set_title('MSA Transformer donor similarity',fontsize=12,pad=13)
    for x in (1,2,3):ax.text(x,1.035,'n = 362',ha='center',fontsize=9,color='#555555')
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',color='#e6e6e6',linewidth=.5);ax.set_axisbelow(True)
    fig.tight_layout()
    for ext in ('png','pdf','svg'):
        fig.savefig(FIG/f'msatransformer_jaccard_donor_separations.{ext}',dpi=400)
    plt.close(fig)
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), script_sha256=sha(Path(__file__)),
        metric='rank-tolerant top-L versus top-floor(1.5L) Jaccard, as manuscript Figure 2B',
        parameters=dict(separations=SEPARATIONS, center_crop=True, cutoff_ties='retain', offset_search=False),
        validation=dict(baseline_comparisons=len(baseline),max_baseline_error=max(r['max_abs_error'] for r in baseline),
                        independent_set_checks=2172,prior_pearson_membership_checks=2172,matrices_verified_against_prior=396),
        source_files=sources, outputs={f.name:sha(f) for f in OUT.glob('*.csv')},
        figures={f.name:sha(f) for f in FIG.iterdir() if f.suffix in ('.png','.pdf','.svg')})
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (OUT/'COMPLETE.json').write_text(json.dumps({'provenance_sha256':sha(OUT/'provenance.json')},indent=2)+'\n')
    print(pd.DataFrame(summary).query("length_class == 'all'").to_string(index=False))
    print('Validated',len(baseline),'baseline comparisons; saved panel to',FIG)


if __name__ == '__main__':
    main()
