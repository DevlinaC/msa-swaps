import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

af2_path = Path("AF2_TMs.csv")
af3_path = Path("tm_align_full_results.csv")
gremlin_path = Path("jaccard_results_gremlin.csv")
msatr_path = Path("jaccard_results_msatransformer.csv")

af2 = pd.read_csv(af2_path)
af2_x = af2["AF2_A"].astype(float).dropna()

tm = pd.read_csv(af3_path)
junk = tm[tm["Job_Type"].astype(str).str.lower().eq("junk")].copy()
junk["QueryPDB"] = junk["Job_Label"].str.split("_").str[0]
af3_x = junk.loc[
    junk["Compared_Protein"].astype(str) == junk["QueryPDB"].astype(str),
    "TM_Score"
].astype(float).dropna()

gremlin = pd.read_csv(gremlin_path)
gremlin_x = gremlin["Jaccard_vs_ProteinA"].astype(float).dropna()

msatr = pd.read_csv(msatr_path)
msatr_x = msatr["Jaccard_vs_ProteinA"].astype(float).dropna()

datasets = [af2_x.values, af3_x.values, gremlin_x.values, msatr_x.values]
positions = [1, 2, 3, 4]

fig, ax = plt.subplots(figsize=(9.2, 6.8))

# Similar-to-query region: 90% transparent = alpha 0.10
ax.axhspan(0.6, 1.02, facecolor="#ae0639", alpha=0.10, zorder=0)

ax.text(
    0.965, 0.965, "Similar to Query",
    transform=ax.transAxes,
    ha="right", va="top",
    fontsize=18,
    fontweight="bold",
    color="#ae0639",
    alpha=1.0
)

parts = ax.violinplot(
    datasets,
    positions=positions,
    widths=0.72,
    showmeans=False,
    showmedians=True,
    showextrema=True,
    bw_method="scott"
)

for body in parts["bodies"]:
    body.set_facecolor("0.75")
    body.set_edgecolor("black")
    body.set_alpha(1.0)
    body.set_zorder(2)

for key in ["cmedians", "cbars", "cmins", "cmaxes"]:
    parts[key].set_color("black")
    parts[key].set_linewidth(1.5)
    parts[key].set_zorder(3)

ax.set_xticks(positions)
ax.set_xticklabels(
    [
        "AlphaFold2\nTM-score\n(n=9,050)",
        "AlphaFold3\nTM-score\n(n=9,050)",
        "GREMLIN\nJaccard score\n(n=362)",
        "MSA Transformer\nJaccard score\n(n=362)",
    ],
    fontsize=14
)

ax.set_ylabel("Query similarity", fontsize=16)
ax.set_ylim(0, 1.02)
ax.tick_params(axis="y", labelsize=14)
ax.grid(axis="y", linestyle="--", alpha=0.22, zorder=1)

fig.tight_layout()
fig.savefig("figure2_query_similarity_violins_bg90pct_transparency.png", dpi=300, bbox_inches="tight")
fig.savefig("figure2_query_similarity_violins_bg90pct_transparency.pdf", bbox_inches="tight")
plt.show()
