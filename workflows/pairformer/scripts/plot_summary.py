#!/usr/bin/env python3
"""Replot validated results with the manuscript palette; do not rerun analysis."""
import argparse
import csv
import hashlib
import io
import json
import os
import sys
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/msa_pairformer_plot_mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/msa_pairformer_plot_cache')
sys.path.insert(0, str(Path(__file__).resolve().parent/'pilot_lib'))
from common import require, sha256, atomic_bytes, write_json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from statistics import median


def plot(root):
    dest=root/'results/analysis'
    ranked=root/'results/rank_tolerant'
    table=dest/'comparisons.csv'; validation=dest/'VALIDATED.json'
    sources=(table,validation,dest/'summary.csv',dest/'paired_direction.csv',
             ranked/'comparisons.csv',ranked/'summary.csv',ranked/'paired_direction.csv',ranked/'VALIDATED.json')
    original_hashes={str(p.relative_to(root)):sha256(p) for p in sources}
    palette=json.loads((root/'figure_palette.json').read_text())['palette']
    with table.open() as f:
        rows=list(csv.DictReader(f))
    require(len(rows)==2172, 'Expected the validated full-run comparison table')
    rank_validation=json.loads((ranked/'VALIDATED.json').read_text())
    require(rank_validation['status']=='PASS','Rank-tolerant validation failed')
    for name in ('comparisons.csv','summary.csv','paired_direction.csv'):
        require(sha256(ranked/name)==rank_validation['output_hashes'][name], 'Rank-tolerant source hash mismatch')
    with (ranked/'comparisons.csv').open() as f:
        rank_rows=list(csv.DictReader(f))
    key=lambda row:(row['id'],row['separation'],row['reference'])
    rank_lookup={key(row):row for row in rank_rows}
    require(len(rank_rows)==len(rank_lookup)==2172,'Duplicate or missing rank-tolerant comparisons')
    require({key(row) for row in rows}==set(rank_lookup),'Comparison keys differ')
    for row in rows:
        rank_row=rank_lookup[key(row)]
        for field in ('L','ordinary_jaccard','top_swap','top_reference','query','donor','source_depth','selected_depth'):
            require(row[field]==rank_row[field], f'Comparison mismatch: {key(row)} {field}')
        row['tolerant_jaccard']=rank_row['tolerant_jaccard']
        for field in ('wide_k','wide_swap','wide_reference','credited','union'):
            row['rank_'+field]=rank_row[field]

    def write_table(path, records):
        buf=io.StringIO(); writer=csv.DictWriter(buf,fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records); atomic_bytes(path,buf.getvalue().encode())

    write_table(dest/'comparisons_with_rank_tolerance.csv',rows)
    with (dest/'summary.csv').open() as f:
        summary=list(csv.DictReader(f))
    with (ranked/'summary.csv').open() as f:
        for row in csv.DictReader(f):
            if row['metric']=='tolerant_jaccard':
                summary.append({field:row[field] for field in summary[0]})
    require(len(summary)==24,'Expected 24 summary rows')
    write_table(dest/'summary_with_rank_tolerance.csv',summary)
    plt.rcParams.update({'font.size':10, 'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none'})
    metrics=('ordinary_jaccard','tolerant_jaccard','pearson_all_eligible','pearson_union_top_L')
    titles=('A  Ordinary top-L Jaccard','B  Rank-tolerant top-L Jaccard',
            'C  Pearson: all eligible pairs','D  Pearson: union of top-L pairs')
    counts={}
    def draw_panel(ax,index):
        metric,title=metrics[index],titles[index]
        for sep_index,sep in enumerate((6,12,24)):
            for ref_index,ref in enumerate(('query','donor')):
                selected=[r for r in rows if int(r['separation'])==sep and
                          r['reference']==ref and r['status']=='ok' and r[metric]!='']
                values=[float(r[metric]) for r in selected]
                counts[f'{metric}/{sep}/{ref}']=len(values)
                require(len(values)==362,'Expected all 362 swaps in every group')
                position=1+3*sep_index+ref_index
                # Stable horizontal jitter only; identical swap offsets across panels.
                jitter=[.56*(int(hashlib.sha256(r['id'].encode()).hexdigest()[:8],16)/0xffffffff-.5)
                        for r in selected]
                ax.scatter([position+x for x in jitter],values,s=10,alpha=.42,
                           color=palette[ref],edgecolors='none',zorder=2)
                ax.plot([position-.31,position+.31],[median(values)]*2,
                        color='black',linewidth=1.7,zorder=3)
        ax.set_xticks([1.5,4.5,7.5],['≥6','≥12','≥24'])
        ax.set_xlabel('Minimum sequence separation')
        ax.set_title(title)
        ax.set_xlim(.4,8.6)
        ax.set_ylim((-.025,1.025) if metric in ('ordinary_jaccard','tolerant_jaccard') else (-1.04,1.04))
        ax.grid(axis='y',alpha=.2)
        ax.set_axisbelow(True)
    figures={}
    def render(indices,shape,size,stem):
        fig,grid=plt.subplots(*shape,figsize=size,constrained_layout=True,squeeze=False)
        axes=grid.ravel()
        for ax,index in zip(axes,indices):
            draw_panel(ax,index)
        axes[0].legend(handles=[Patch(facecolor=palette['query'],label='Query A'),
                               Patch(facecolor=palette['donor'],label='Donor B')],
                       loc='upper center',bbox_to_anchor=(.5,1.22),ncol=2,borderaxespad=0)
        for suffix in ('png','pdf','svg'):
            buffer=io.BytesIO()
            fig.savefig(buffer,format=suffix,dpi=300)
            path=dest/(stem+'.'+suffix); path.parent.mkdir(parents=True,exist_ok=True)
            atomic_bytes(path,buffer.getvalue()); figures[str(path.relative_to(dest))]=sha256(path)
        plt.close(fig)

    render((0,1,2,3),(2,2),(11,8),'contact_comparison')
    render((0,1),(1,2),(11,4.5),'panels/upper_row_jaccard')
    render((2,3),(1,2),(11,4.5),'panels/lower_row_pearson')
    for index,name in enumerate(('panel_A_ordinary_jaccard','panel_B_rank_tolerant_jaccard',
                                 'panel_C_pearson_all_pairs','panel_D_pearson_top_L_union')):
        render((index,),(1,1),(5.5,4.5),'panels/'+name)
    for name,expected in original_hashes.items():
        require(sha256(root/name)==expected, f'Numerical source changed: {name}')
    write_json(dest/'figure_provenance.json',dict(
        palette=palette,palette_source=json.loads((root/'figure_palette.json').read_text()),
        source_hashes=original_hashes,figures=figures,group_counts=counts,
        plotting_code_sha256=sha256(Path(__file__)),
        revision='Four-panel summary plus separate upper/lower rows and individual panels; original numerical results and validation records unchanged',
        export_formats='PNG at 300 dpi; vector PDF with embedded fonts; SVG with editable text and vector points.',
        combined_tables={name:sha256(dest/name) for name in ('comparisons_with_rank_tolerance.csv','summary_with_rank_tolerance.csv')},
        joined_records=len(rows),metrics=metrics,
        inference_configuration=dict(max_depth=512,query_biasing=False,sequence_weighting='uniform'),
        configuration_annotation='Written figure legend: results/analysis/figure_legend.txt; no configuration title on plots.',
        point_plot='All 362 swaps per group; translucent points with deterministic horizontal jitter; black bars mark medians. No vertical jitter or omitted outliers.',
        jitter='SHA256 of swap ID; first 32 bits mapped uniformly to horizontal offsets in [-0.28, 0.28].'))
    print('Updated PNG/PDF/SVG. Numerical table and validation hashes unchanged.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    plot(p.parse_args().root.resolve())
