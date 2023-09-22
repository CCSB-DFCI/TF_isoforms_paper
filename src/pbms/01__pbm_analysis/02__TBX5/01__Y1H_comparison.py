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


# ## TBX5

# In[9]:


kfit_dat = pd.read_table("../../../../data/processed/pbms/TBX5kfit_dat.csv", sep=",")
da_dat = pd.read_table("../../../../data/processed/pbms/TBX5da_dat.csv", sep=",")


# In[10]:


kfit_ref = kfit_dat[kfit_dat["cname"] == "TBX5C05-REF"]
kfit_3 = kfit_dat[kfit_dat["cname"] == "TBX5A05"]
kfit_2 = kfit_dat[kfit_dat["cname"] == "TBX5B05"]


# In[11]:


kfit_vs = kfit_ref.merge(kfit_2, on=["seq"], suffixes=("_ref", ""))
kfit_vs = kfit_vs.merge(kfit_3, on=["seq"], suffixes=("_iso2", "_iso3"))
kfit_vs.head()


# In[12]:


da_3 = da_dat[da_dat["cname"] == "TBX5A05"]
da_2 = da_dat[da_dat["cname"] == "TBX5B05"]


# In[13]:


kfit_vs = kfit_vs.merge(da_2[["seq", "contrastQ", "contrastQ_cut"]], on="seq")
kfit_vs = kfit_vs.merge(
    da_3[["seq", "contrastQ", "contrastQ_cut"]], on="seq", suffixes=("_iso2", "_iso3")
)


# In[14]:


kfit_vs.contains_any_motif_ref.value_counts()


# In[15]:


## define the other non TBX5 canonical k-mers
def define_other_kmers(row):
    if row.contains_any_motif_ref == "* TBX5 k-mer":
        return row.contains_any_motif_ref
    elif "ACGTGT" in row.seq or "ACACGT" in row.seq:
        return "* ACGTGT k-mer"
    else:
        return "other k-mer"


kfit_vs["contains_any_motif"] = kfit_vs.apply(define_other_kmers, axis=1)
kfit_vs.contains_any_motif.value_counts()


# In[16]:


kfit_vs.affinityEstimate_iso3.min()


# In[17]:


fig, axarr = plt.subplots(1, 2, figsize=(4, 2))

markers = [",", ".", "o"]
titles = ["TBX5-2", "TBX5-3"]

for i, suffix in enumerate(["iso2", "iso3"]):
    ax = axarr[i]

    for k, motif in enumerate(["other k-mer", "* TBX5 k-mer", "* ACGTGT k-mer"]):
        for j, qval in enumerate(
            ["(0.1,1]", "(0.01,0.1]", "(0.001,0.01]", "[0,0.001]"]
        ):
            sub = kfit_vs[
                (kfit_vs["contrastQ_cut_%s" % suffix] == qval)
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

    ax.set_xlim((9.2, 13))
    ax.set_ylim((9.2, 13))
    ax.plot(
        [9.2, 13], [9.2, 13], color="black", linestyle="dashed", linewidth=1, zorder=1
    )
    ax.set_xticks([10, 11, 12, 13])
    ax.set_yticks([10, 11, 12, 13])
    ax.set_xlabel("reference PBM affinity")
    ax.set_title(titles[i])

axarr[0].set_ylabel("alternative PBM affinity")
fig.savefig(
    "../../../../figures/TBX5_isoforms_pbm_scatter.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# In[18]:


# add colors to kfit_vs
cmap_name = "Spectral"
colname = "contrastDifference_iso2"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(vmin=kfit_vs[colname].min(), vmax=kfit_vs[colname].max())
m = cm.ScalarMappable(norm=norm, cmap=cmap)

kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)


# In[19]:


kfit_vs[["contrastDifference_iso2", "contrastDifference_iso3"]].min().min()


# In[20]:


# add colors to kfit_vs
cmap_name = "Spectral"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(
    vmin=kfit_vs[["contrastDifference_iso2", "contrastDifference_iso3"]].min().min(),
    vmax=kfit_vs[["contrastDifference_iso2", "contrastDifference_iso3"]].max().max(),
)
m = cm.ScalarMappable(norm=norm, cmap=cmap)

colname = "contrastDifference_iso2"
kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)

colname = "contrastDifference_iso3"
kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)


# In[21]:


tbx5_y1h = y1h.loc[y1h["tf"] == "TBX5", y1h.columns[1:]].copy().set_index("unique_acc")
tbx5_y1h = tbx5_y1h.loc[:, tbx5_y1h.any(axis=0)]


# In[22]:


tbx5_baits = list(tbx5_y1h.columns)
tbx5_baits = dna[dna["id_upper"].isin(tbx5_baits)]


# In[23]:


tbx5_baits["tmp_iso2"] = tbx5_baits.apply(
    vals_per_bait,
    axis=1,
    kfit_vs=kfit_vs,
    colname="contrastDifference_iso2",
    alt_suffix="iso2",
    ascending=False,
)
tbx5_baits["val_kmer_iso2"] = (
    tbx5_baits["tmp_iso2"].str.split("_", expand=True)[0].astype(str)
)
tbx5_baits["val_diff_iso2"] = (
    tbx5_baits["tmp_iso2"].str.split("_", expand=True)[1].astype(float)
)
tbx5_baits["xval_iso2"] = (
    tbx5_baits["tmp_iso2"].str.split("_", expand=True)[2].astype(float)
)
tbx5_baits["yval_iso2"] = (
    tbx5_baits["tmp_iso2"].str.split("_", expand=True)[3].astype(float)
)
tbx5_baits["rgb_iso2"] = (
    tbx5_baits["tmp_iso2"].str.split("_", expand=True)[4].astype(str)
)
tbx5_baits["rgb_iso2"] = tbx5_baits["rgb_iso2"].apply(literal_eval)


# In[24]:


tbx5_baits["tmp_iso3"] = tbx5_baits.apply(
    vals_per_bait,
    axis=1,
    kfit_vs=kfit_vs,
    colname="contrastDifference_iso3",
    alt_suffix="iso3",
    ascending=False,
)
tbx5_baits["val_kmer_iso3"] = (
    tbx5_baits["tmp_iso3"].str.split("_", expand=True)[0].astype(str)
)
tbx5_baits["val_diff_iso3"] = (
    tbx5_baits["tmp_iso3"].str.split("_", expand=True)[1].astype(float)
)
tbx5_baits["xval_iso3"] = (
    tbx5_baits["tmp_iso3"].str.split("_", expand=True)[2].astype(float)
)
tbx5_baits["yval_iso3"] = (
    tbx5_baits["tmp_iso3"].str.split("_", expand=True)[3].astype(float)
)
tbx5_baits["rgb_iso3"] = (
    tbx5_baits["tmp_iso3"].str.split("_", expand=True)[4].astype(str)
)
tbx5_baits["rgb_iso3"] = tbx5_baits["rgb_iso3"].apply(literal_eval)


# In[25]:


colors = tbx5_baits[["id_upper", "rgb_iso2", "rgb_iso3"]].set_index("id_upper").T
colors.index = ["TBX5-2", "TBX5-3"]
colors.loc["TBX5-1"] = [(0, 0, 0, 1)] * len(tbx5_baits)
colors = colors.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]


# In[26]:


annot = (
    tbx5_baits[["id_upper", "val_diff_iso2", "val_diff_iso3"]].set_index("id_upper").T
)
annot.index = ["TBX5-2", "TBX5-3"]
annot.loc["TBX5-1"] = ["NA"] * len(tbx5_baits)
annot = annot.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]


# In[34]:


fig, ax = plt.subplots(1, 1, figsize=(4, 2))
y1h_pdi_per_tf_gene_plot(
    "TBX5",
    data=y1h,
    ax=ax,
    iso_order=["TBX5-1", "TBX5-2", "TBX5-3"],
    bait_colors=colors,
)

plt.colorbar(
    m,
    ax=ax,
    orientation="horizontal",
    shrink=0.75,
    label="largest ∆ PBM affinity across 8-mers in bait",
)
plt.savefig(
    "../../../../figures/TBX5_y1h_with_pbm.pdf", bbox_inches="tight", dpi="figure"
)


# In[35]:


sns.color_palette("Spectral_r", n_colors=4)
