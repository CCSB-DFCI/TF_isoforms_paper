# coding: utf-8

# In[1]:


from ast import literal_eval
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import sys

sys.path.append("../../..")
from data_loading import *
from isoform_pairwise_metrics import *
from plotting import (
    y1h_pdi_per_tf_gene_plot,
    m1h_activation_per_tf_gene_plot,
    COLOR_PURPLE,
)
from data_loading import load_annotated_TFiso1_collection, load_y1h_pdi_data


# ## functions

# In[2]:


def vals_per_bait(row, kfit_vs, colname, alt_suffix, ascending):
    kmers = []
    seq = row.seq
    seq_len = row.seq_len
    for i in range(0, seq_len - 8):
        kmer = seq[i : i + 8]
        kmers.append(kmer)

    sub = kfit_vs[kfit_vs["seq"].isin(kmers)]
    sub["abs"] = np.abs(sub[colname])
    sub = sub.sort_values(by="abs", ascending=ascending)
    largest_kmer = sub["seq"].iloc[0]
    largest_val = sub[colname].iloc[0]
    xval = sub["affinityEstimate_ref"].iloc[0]
    yval = sub["affinityEstimate_%s" % alt_suffix].iloc[0]

    rgb = sub["%s_rgb" % colname].iloc[0]
    return "%s_%s_%s_%s_%s" % (largest_kmer, largest_val, xval, yval, rgb)


# In[3]:


def get_rgb(row, colname, m):
    x = row[colname]
    rgb = m.to_rgba(x)
    return rgb


# In[4]:


from matplotlib import pyplot as plt
from matplotlib import patches

PAPER_PRESET = {
    "style": "ticks",
    "font": "Helvetica",
    "context": "paper",
    "rc": {
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "axes.linewidth": 0.5,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "axes.edgecolor": "black",
        "xtick.major.pad": 3.0,
        "ytick.major.pad": 3.0,
    },
}
PAPER_FONTSIZE = 10

sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# ## 1. import baits and y1h results

# In[5]:


tfs = load_annotated_TFiso1_collection()


# In[6]:


y1h_baits_f = "../../../../data/internal/Y1H_DNA_baits.fa"


# In[7]:


ids = []
seqs = []

for record in SeqIO.parse(y1h_baits_f, "fasta"):
    ids.append(record.description)
    seqs.append(str(record.seq))

dna = pd.DataFrame()
dna["id"] = ids
dna["seq"] = seqs
dna["seq_len"] = dna.seq.str.len()
dna["id_upper"] = dna["id"].str.upper()
dna


# In[8]:


y1h = load_y1h_pdi_data(add_missing_data=True)
y1h


# ## CREB1

# In[26]:


kfit_dat = pd.read_table("../../../../data/processed/pbms/CREB1kfit_dat.csv", sep=",")
da_dat = pd.read_table("../../../../data/processed/pbms/CREB1da_dat.csv", sep=",")


# In[27]:


kfit_ref = kfit_dat[kfit_dat["cname"] == "CREB1-ref"]
kfit_2 = kfit_dat[kfit_dat["cname"] == "CREB1-alt"]


# In[28]:


kfit_vs = kfit_ref.merge(kfit_2, on=["seq"], suffixes=("_ref", "_alt"))
kfit_vs.head()


# In[29]:


da_2 = da_dat[da_dat["cname"] == "CREB1-alt"]


# In[30]:


kfit_vs = kfit_vs.merge(da_2[["seq", "contrastQ", "contrastQ_cut"]], on="seq")


# In[31]:


kfit_vs.contains_any_motif_ref.value_counts()


# In[32]:


## define the other non TBX5 canonical k-mers
def define_other_kmers(row):
    if row.contains_any_motif_ref == "CREB1 k-mer":
        return "* CREB1 k-mer"
    else:
        return "other k-mer"


kfit_vs["contains_any_motif"] = kfit_vs.apply(define_other_kmers, axis=1)
kfit_vs.contains_any_motif.value_counts()


# In[37]:


kfit_vs.affinityEstimate_ref.max()


# In[40]:


fig, ax = plt.subplots(1, 1, figsize=(2, 2))

markers = [",", "."]
titles = ["CREB1-alt"]

for i, suffix in enumerate(["alt"]):
    for k, motif in enumerate(["other k-mer", "* CREB1 k-mer"]):
        for j, qval in enumerate(
            ["(0.1,1]", "(0.01,0.1]", "(0.001,0.01]", "[0,0.001]"]
        ):
            sub = kfit_vs[
                (kfit_vs["contrastQ_cut"] == qval)
                & (kfit_vs["contains_any_motif"] == motif)
            ]
            xs = sub["affinityEstimate_ref"]
            ys = sub["affinityEstimate_%s" % suffix]

            color = sns.color_palette("Spectral_r", n_colors=4)[j]
            marker = markers[k]

            if marker == "o":
                ax.scatter(
                    xs,
                    ys,
                    15,
                    marker=marker,
                    edgecolors="black",
                    facecolors=color,
                    alpha=1,
                    linewidth=0.5,
                )
            elif marker == ",":
                ax.scatter(
                    xs,
                    ys,
                    1,
                    marker=".",
                    edgecolors=color,
                    facecolors="none",
                    alpha=0.5,
                )
            else:
                ax.scatter(
                    xs,
                    ys,
                    30,
                    marker=marker,
                    edgecolors=color,
                    facecolors="white",
                    alpha=1,
                )

    ax.set_xlim((6.8, 11.6))
    ax.set_ylim((6.8, 11.6))
    ax.plot(
        [6.8, 11.6],
        [6.8, 11.6],
        color="black",
        linestyle="dashed",
        linewidth=1,
        zorder=1,
    )
    ax.set_xticks([7, 8, 9, 10, 11])
    ax.set_yticks([7, 8, 9, 10, 11])
    ax.set_xlabel("reference PBM affinity")
    ax.set_title(titles[i])

ax.set_ylabel("alternative PBM affinity")
fig.savefig(
    "../../../../figures/CREB1_isoforms_pbm_scatter.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# In[43]:


# add colors to kfit_vs
cmap_name = "Spectral"
colname = "contrastDifference_alt"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(vmin=kfit_vs[colname].min(), vmax=kfit_vs[colname].max())
m = cm.ScalarMappable(norm=norm, cmap=cmap)

kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)


# In[56]:


# add colors to kfit_vs
cmap_name = "Spectral"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(
    vmin=kfit_vs[["contrastDifference_alt"]].min(),
    vmax=kfit_vs[["contrastDifference_alt"]].max(),
)
m = cm.ScalarMappable(norm=norm, cmap=cmap)

colname = "contrastDifference_alt"
kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)


# In[57]:


creb1_y1h = (
    y1h.loc[y1h["gene_symbol"] == "CREB1", y1h.columns[1:]]
    .copy()
    .set_index("clone_acc")
)
creb1_y1h = creb1_y1h.loc[:, creb1_y1h.any(axis=0)]


# In[58]:


creb1_baits = list(creb1_y1h.columns)
creb1_baits = dna[dna["id_upper"].isin(creb1_y1h)]


# In[59]:


creb1_baits["tmp_alt"] = creb1_baits.apply(
    vals_per_bait,
    axis=1,
    kfit_vs=kfit_vs,
    colname="contrastDifference_alt",
    alt_suffix="alt",
    ascending=False,
)
creb1_baits["val_kmer_alt"] = (
    creb1_baits["tmp_alt"].str.split("_", expand=True)[0].astype(str)
)
creb1_baits["val_diff_alt"] = (
    creb1_baits["tmp_alt"].str.split("_", expand=True)[1].astype(float)
)
creb1_baits["xval_alt"] = (
    creb1_baits["tmp_alt"].str.split("_", expand=True)[2].astype(float)
)
creb1_baits["yval_alt"] = (
    creb1_baits["tmp_alt"].str.split("_", expand=True)[3].astype(float)
)
creb1_baits["rgb_alt"] = (
    creb1_baits["tmp_alt"].str.split("_", expand=True)[4].astype(str)
)
creb1_baits["rgb_alt"] = creb1_baits["rgb_alt"].apply(literal_eval)


# In[60]:


colors = creb1_baits[["id_upper", "rgb_alt"]].set_index("id_upper").T
colors.index = ["CREB1-1"]
colors.loc["CREB1-2"] = [(0, 0, 0, 1)] * len(creb1_baits)
colors = colors.loc[["CREB1-2", "CREB1-1"]]


# In[61]:


annot = creb1_baits[["id_upper", "val_diff_alt"]].set_index("id_upper").T
annot.index = ["CREB1-1"]
annot.loc["CREB1-2"] = ["NA"] * len(creb1_baits)
annot = annot.loc[["CREB1-2", "CREB1-1"]]


# In[62]:


fig, ax = plt.subplots(1, 1, figsize=(2, 2))
y1h_pdi_per_tf_gene_plot(
    "CREB1", data=y1h, ax=ax, iso_order=["CREB1-2", "CREB1-1"], bait_colors=colors
)

plt.colorbar(
    m,
    ax=ax,
    orientation="horizontal",
    shrink=0.75,
    label="largest ∆ PBM affinity across 8-mers in bait",
)
plt.savefig(
    "../../../../figures/CREB1_y1h_with_pbm.pdf", bbox_inches="tight", dpi="figure"
)


# In[63]:


annot


# In[64]:


colors
