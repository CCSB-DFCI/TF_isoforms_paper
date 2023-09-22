# coding: utf-8

# In[1]:


from ast import literal_eval
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import patches
import matplotlib.cm as cm
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import sys

# import utils
sys.path.append("../")

from data_loading import *
from isoform_pairwise_metrics import *
from plotting import (
    y1h_pdi_per_tf_gene_plot,
    m1h_activation_per_tf_gene_plot,
    COLOR_PURPLE,
)
from data_loading import load_annotated_TFiso1_collection, load_y1h_pdi_data


# In[2]:


tfs = load_annotated_TFiso1_collection()


# In[3]:


y1h_baits_f = "../../data/internal/Y1H_DNA_baits.fa"


# ## 1. import baits and y1h results

# In[4]:


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


# In[5]:


y1h = load_y1h_pdi_data(add_missing_data=True)
y1h


# ## CREB1

# In[6]:


kfit_dat = pd.read_table("../../data/processed/pbms/CREB1kfit_dat.csv", sep=",")
da_dat = pd.read_table("../../data/processed/pbms/CREB1da_dat.csv", sep=",")


# In[7]:


kfit_ref = kfit_dat[kfit_dat["cname"] == "CREB1-ref"]
kfit_alt = kfit_dat[kfit_dat["cname"] == "CREB1-alt"]


# In[8]:


kfit_vs = kfit_ref.merge(kfit_alt, on=["seq"], suffixes=("_ref", "_alt"))


# In[9]:


da_alt = da_dat[da_dat["cname"] == "CREB1-alt"]


# In[10]:


kfit_vs = kfit_vs.merge(
    da_alt[
        ["seq", "contrastQ", "contrastQ_cut", "contrastAverage", "contrastDifference"]
    ],
    on="seq",
)


# In[11]:


kfit_vs.contains_any_motif_ref.value_counts()


# In[12]:


kfit_vs.affinityEstimate_ref.max()


# In[13]:


creb1_y1h = (
    y1h.loc[y1h["tf"] == "CREB1", y1h.columns[1:]].copy().set_index("unique_acc")
)
creb1_y1h = creb1_y1h.loc[:, creb1_y1h.any(axis=0)]
creb1_y1h


# In[14]:


creb1_baits = list(creb1_y1h.columns)
creb1_baits = dna[dna["id_upper"].isin(creb1_baits)]
creb1_baits


# In[15]:


kfit_vs.columns


# In[16]:


fig, ax = plt.subplots(1, 1, figsize=(2, 2))

markers = [",", "."]
titles = ["CREB1-1"]

for k, motif in enumerate(["*other k-mer", "CREB1 k-mer"]):
    for j, qval in enumerate(["(0.1,1]", "(0.01,0.1]", "(0.001,0.01]", "[0,0.001]"]):
        sub = kfit_vs[
            (kfit_vs["contrastQ_cut"] == qval)
            & (kfit_vs["contains_any_motif_ref"] == motif)
        ]
        xs = sub["affinityEstimate_ref"]
        ys = sub["affinityEstimate_alt"]
        ts = sub["contrastQ"]

        color = sns.color_palette("Spectral_r", n_colors=4)[j]
        marker = markers[k]

        if marker == "o":
            ax.scatter(
                xs,
                ys,
                30,
                marker=".",
                edgecolors="black",
                facecolors=color,
                alpha=1,
                linewidth=0.5,
                zorder=10,
            )
        elif marker == ",":
            ax.scatter(
                xs,
                ys,
                30,
                marker=".",
                edgecolors=color,
                facecolors=color,
                alpha=0.5,
                zorder=10,
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
                linewidth=0.5,
                zorder=10,
            )

ax.set_xlim((6.8, 11.6))
ax.set_ylim((6.8, 11.6))
ax.plot(
    [6.8, 11.6], [6.8, 11.6], color="black", linestyle="dashed", linewidth=1, zorder=1
)
# ax.set_xticks([6, 8, 10, 12])
# ax.set_yticks([6, 8, 10, 12])
ax.set_xlabel("CREB1 reference")
ax.set_ylabel("CREB1 alternative")
ax.set_title("CREB1-1")

# # annotate
# for i, row in creb1_baits.iterrows():
#     ax.annotate(row.id_upper, xy=(row.xval, row.yval),
#                 xytext=(5, -5), textcoords='offset points', arrowprops = dict(arrowstyle="-"),
#                 ha="left", va="top", fontsize=7,
#                 bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))

fig.savefig(
    "../../figures/CREB1_isoforms_pbm_scatter.pdf", dpi="figure", bbox_inches="tight"
)


# In[17]:


# add colors to kfit_vs
cmap_name = "Spectral"
colname = "contrastDifference_alt"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(vmin=kfit_vs[colname].min(), vmax=kfit_vs[colname].max())
m = cm.ScalarMappable(norm=norm, cmap=cmap)


def get_rgb(row, colname, m):
    x = row[colname]
    rgb = m.to_rgba(x)
    return rgb


kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)
kfit_vs.sample(5)


# In[18]:


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


# In[19]:


creb1_baits["tmp"] = creb1_baits.apply(
    vals_per_bait,
    axis=1,
    kfit_vs=kfit_vs,
    colname="contrastDifference_alt",
    alt_suffix="alt",
    ascending=False,
)
creb1_baits["val_kmer"] = creb1_baits["tmp"].str.split("_", expand=True)[0].astype(str)
creb1_baits["val_diff"] = (
    creb1_baits["tmp"].str.split("_", expand=True)[1].astype(float)
)
creb1_baits["xval"] = creb1_baits["tmp"].str.split("_", expand=True)[2].astype(float)
creb1_baits["yval"] = creb1_baits["tmp"].str.split("_", expand=True)[3].astype(float)
creb1_baits["rgb"] = creb1_baits["tmp"].str.split("_", expand=True)[4].astype(str)
creb1_baits["rgb"] = creb1_baits["rgb"].apply(literal_eval)
creb1_baits


# In[20]:


colors = creb1_baits[["id_upper", "rgb"]].set_index("id_upper").T
colors.index = ["CREB1-1"]
colors.loc["CREB1-2"] = [(0, 0, 0, 1)] * len(creb1_baits)
colors = colors.loc[["CREB1-2", "CREB1-1"]]
colors


# In[21]:


annot = creb1_baits[["id_upper", "val_diff"]].set_index("id_upper").T
annot.index = ["CREB1-1"]
annot.loc["CREB1-2"] = ["NA"] * len(creb1_baits)
annot = annot.loc[["CREB1-2", "CREB1-1"]]
annot


# In[22]:


fig, ax = plt.subplots(1, 1, figsize=(2, 2.5))
y1h_pdi_per_tf_gene_plot(
    "CREB1",
    data=y1h,
    ax=ax,
    iso_order=["CREB1-2", "CREB1-1"],
    bait_colors=colors,
    bait_annot=annot,
)

plt.colorbar(
    m,
    ax=ax,
    orientation="horizontal",
    shrink=0.75,
    label="largest ∆ PBM affinity\nacross 8-mers in bait",
)
plt.savefig("../../figures/CREB1_y1h_with_pbm.pdf", bbox_inches="tight", dpi="figure")


# In[23]:


def y1h_pdi_per_tf_gene_plot(
    gene_name,
    data,
    ax=None,
    min_n_isoforms=1,
    min_n_partners=1,
    iso_order=None,
    bait_colors=None,
    bait_annot=None,
):
    tf = (
        data.loc[data["tf"] == gene_name, data.columns[1:]]
        .copy()
        .set_index("unique_acc")
    )
    tf.index = tf.index.map(isoform_display_name)
    tf = tf.loc[:, tf.any(axis=0)]
    if ax is None:
        ax = plt.gca()
    if tf.shape[0] < min_n_isoforms or tf.shape[1] < min_n_partners:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No PDI data available",
            ha="center",
            va="center",
            fontsize=30,
            fontweight="bold",
            color="grey",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    if iso_order is None:
        tf = tf
    else:
        tf = tf.loc[iso_order, :]
    binary_profile_matrix(
        tf,
        ax=ax,
        column_label_rotation=90,
        bait_colors=bait_colors,
        bait_annot=bait_annot,
    )
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )


# ## TBX5

# In[24]:


kfit_dat = pd.read_table("../../data/processed/pbms/TBX5kfit_dat.csv", sep=",")
da_dat = pd.read_table("../../data/processed/pbms/TBX5da_dat.csv", sep=",")


# In[25]:


kfit_ref = kfit_dat[kfit_dat["cname"] == "TBX5C05-REF"]
kfit_3 = kfit_dat[kfit_dat["cname"] == "TBX5A05"]
kfit_2 = kfit_dat[kfit_dat["cname"] == "TBX5B05"]


# In[26]:


kfit_vs = kfit_ref.merge(kfit_2, on=["seq"], suffixes=("_ref", ""))
kfit_vs = kfit_vs.merge(kfit_3, on=["seq"], suffixes=("_iso2", "_iso3"))
kfit_vs.head()


# In[27]:


da_3 = da_dat[da_dat["cname"] == "TBX5A05"]
da_2 = da_dat[da_dat["cname"] == "TBX5B05"]


# In[28]:


kfit_vs = kfit_vs.merge(da_2[["seq", "contrastQ", "contrastQ_cut"]], on="seq")
kfit_vs = kfit_vs.merge(
    da_3[["seq", "contrastQ", "contrastQ_cut"]], on="seq", suffixes=("_iso2", "_iso3")
)


# In[29]:


kfit_vs.contains_any_motif_ref.value_counts()


# In[30]:


kfit_vs.affinityEstimate_iso3.min()


# In[31]:


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
                & (kfit_vs["contains_any_motif_ref"] == motif)
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
    "../../figures/TBX5_isoforms_pbm_scatter.pdf", dpi="figure", bbox_inches="tight"
)


# In[32]:


# add colors to kfit_vs
cmap_name = "Spectral"
colname = "contrastDifference_iso2"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(vmin=kfit_vs[colname].min(), vmax=kfit_vs[colname].max())
m = cm.ScalarMappable(norm=norm, cmap=cmap)

kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)
kfit_vs.sample(5)


# In[33]:


kfit_vs[["contrastDifference_iso2", "contrastDifference_iso3"]].min().min()


# In[34]:


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
kfit_vs.sample(5)


# In[35]:


tbx5_y1h = y1h.loc[y1h["tf"] == "TBX5", y1h.columns[1:]].copy().set_index("unique_acc")
tbx5_y1h = tbx5_y1h.loc[:, tbx5_y1h.any(axis=0)]
tbx5_y1h


# In[36]:


tbx5_baits = list(tbx5_y1h.columns)
tbx5_baits = dna[dna["id_upper"].isin(tbx5_baits)]
tbx5_baits


# In[37]:


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


# In[38]:


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


# In[39]:


colors = tbx5_baits[["id_upper", "rgb_iso2", "rgb_iso3"]].set_index("id_upper").T
colors.index = ["TBX5-2", "TBX5-3"]
colors.loc["TBX5-1"] = [(0, 0, 0, 1)] * len(tbx5_baits)
colors = colors.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
colors


# In[40]:


annot = (
    tbx5_baits[["id_upper", "val_diff_iso2", "val_diff_iso3"]].set_index("id_upper").T
)
annot.index = ["TBX5-2", "TBX5-3"]
annot.loc["TBX5-1"] = ["NA"] * len(tbx5_baits)
annot = annot.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
annot


# In[41]:


def binary_profile_matrix(
    data,
    ax=None,
    box_size=0.7,
    fill_color="black",
    border_color="black",
    column_label_rotation=40,
    bait_colors=None,
    bait_annot=None,
):
    """Used for edgotyping: displays binary results with a grid of boxes

    Empty box for negative, filled box for positive and
    missing box for missing values.

    (Copied over from ccsblib)

    Args:
        data (pandas.DataFrame): boolean values
        ax (matplotlib.axes.Axes): axes to draw on
        box_size (float): area of the boxes between 0 and 1
        fill_color (str): color of filled sqaure
        border_color (str): color of outside of square
        column_label_roation: angle in degrees of top labels
    Examples:
        Display results of some dummy edgotyping data:
        .. plot::
            :context: close-figs
            >>> import pandas as pd
            >>> df = pd.DataFrame(index=['Allele A', 'Allele B', 'Allele C'],
            ...                   columns=['Partner W', 'Partner X', 'Parner Y', 'Partner Z'],
            ...                   data=[[True] * 4,
            ...                         [True, False, True, False],
            ...                         [False, False, np.nan, False]])
            >>> binary_profile_matrix(df)
        You can switch the rows and columns by transposing the input DataFrame:
        .. plot::
            :context: close-figs
             >>> binary_profile_matrix(df.T,
             ...                            fill_color='grey',
             ...                            border_color='black',
             ...                            box_size=0.6,
             ...                            column_label_rotation=90)
    """
    if box_size > 1.0 or box_size < 0.0:
        raise ValueError("box_size must be between 0-1")
    if ax is None:
        ax = plt.gca()
    ax.set_aspect("equal")
    positives = [
        (i, j)
        for i in range(data.shape[1])
        for j in range(data.shape[0])
        if pd.notnull(data.iloc[j, i]) and data.iloc[j, i] == 1
    ]
    negatives = [
        (i, j)
        for i in range(data.shape[1])
        for j in range(data.shape[0])
        if pd.notnull(data.iloc[j, i]) and data.iloc[j, i] == 0
    ]
    for x, y in negatives:
        if bait_colors is None:
            neg_fill = False
            facecolor = None
        else:
            neg_fill = True
            facecolor = bait_colors.iloc[y, x]

        r = patches.Rectangle(
            (x - box_size / 2, y - box_size / 2),
            box_size,
            box_size,
            fill=neg_fill,
            facecolor=facecolor,
            edgecolor=border_color,
            linewidth=1,
        )
        ax.add_patch(r)

        if bait_annot is not None:
            annot = bait_annot.iloc[y, x]
            ax.text(
                x,
                y,
                "{:.2f}".format(annot),
                fontsize=6,
                color="black",
                ha="center",
                va="center",
            )

    for x, y in positives:
        if bait_colors is None:
            facecolor = fill_color
            linewidth = 1
        else:
            facecolor = bait_colors.iloc[y, x]
            linewidth = 3

        r = patches.Rectangle(
            (x - box_size / 2, y - box_size / 2),
            box_size,
            box_size,
            fill=True,
            facecolor=facecolor,
            edgecolor=border_color,
            linewidth=linewidth,
        )
        ax.add_patch(r)

        if bait_annot is not None:
            annot = bait_annot.iloc[y, x]
            if annot != "NA":
                ax.text(
                    x,
                    y,
                    "{:.2f}".format(annot),
                    fontsize=6,
                    color="black",
                    ha="center",
                    va="center",
                )

    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.tick_top()
    ax.set_xticks(range(data.shape[1]))
    ax.xaxis.set_tick_params(length=0)
    ax.xaxis.set_ticklabels(data.columns, rotation=column_label_rotation, ha="center")
    ax.set_yticks(range(data.shape[0]))
    ax.yaxis.set_tick_params(length=0)
    ax.set_yticklabels(data.index)
    ax.set_ylim((-0.5, data.shape[0] - 0.5))
    ax.set_xlim((-0.5, data.shape[1] - 0.5))
    ax.invert_yaxis()


def isoform_display_name(s):
    """Convert clone accession ID to display friendly format"""
    return s.split("|")[0] + "-" + s.split("|")[1].split("/")[0]


# In[42]:


fig, ax = plt.subplots(1, 1, figsize=(5, 5))
y1h_pdi_per_tf_gene_plot(
    "TBX5",
    data=y1h,
    ax=ax,
    iso_order=["TBX5-1", "TBX5-2", "TBX5-3"],
    bait_colors=colors,
    bait_annot=annot,
)

plt.colorbar(
    m,
    ax=ax,
    orientation="horizontal",
    shrink=0.75,
    label="largest ∆ PBM affinity across 8-mers in bait",
)
plt.savefig("../../figures/TBX5_y1h_with_pbm.pdf", bbox_inches="tight", dpi="figure")


# In[43]:


np.mean([-0.72, -0.38, -0.88, -0.60, -0.62, -0.16])


# In[44]:


np.mean([-0.36, -0.88, -0.72, -0.53, -0.72, -0.51, -0.13])


# In[45]:


def vals_per_bait_last50(row, kfit_vs, colname, alt_suffix, ascending):
    kmers = []
    seq = row.seq
    seq_len = row.seq_len
    for i in range(seq_len - 100, seq_len - 8):
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


# In[46]:


tbx5_baits["tmp_iso2_l5"] = tbx5_baits.apply(
    vals_per_bait_last50,
    axis=1,
    kfit_vs=kfit_vs,
    colname="contrastDifference_iso2",
    alt_suffix="iso2",
    ascending=False,
)
tbx5_baits["val_kmer_iso2_l5"] = (
    tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[0].astype(str)
)
tbx5_baits["val_diff_iso2_l5"] = (
    tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[1].astype(float)
)
tbx5_baits["xval_iso2_l5"] = (
    tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[2].astype(float)
)
tbx5_baits["yval_iso2_l5"] = (
    tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[3].astype(float)
)
tbx5_baits["rgb_iso2_l5"] = (
    tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[4].astype(str)
)
tbx5_baits["rgb_iso2_l5"] = tbx5_baits["rgb_iso2_l5"].apply(literal_eval)


# In[47]:


tbx5_baits["tmp_iso3_l5"] = tbx5_baits.apply(
    vals_per_bait_last50,
    axis=1,
    kfit_vs=kfit_vs,
    colname="contrastDifference_iso3",
    alt_suffix="iso3",
    ascending=False,
)
tbx5_baits["val_kmer_iso3_l5"] = (
    tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[0].astype(str)
)
tbx5_baits["val_diff_iso3_l5"] = (
    tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[1].astype(float)
)
tbx5_baits["xval_iso3_l5"] = (
    tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[2].astype(float)
)
tbx5_baits["yval_iso3_l5"] = (
    tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[3].astype(float)
)
tbx5_baits["rgb_iso3_l5"] = (
    tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[4].astype(str)
)
tbx5_baits["rgb_iso3_l5"] = tbx5_baits["rgb_iso3_l5"].apply(literal_eval)


# In[48]:


colors = tbx5_baits[["id_upper", "rgb_iso2_l5", "rgb_iso3_l5"]].set_index("id_upper").T
colors.index = ["TBX5-2", "TBX5-3"]
colors.loc["TBX5-1"] = [(0, 0, 0, 1)] * len(tbx5_baits)
colors = colors.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
colors


# In[49]:


annot = (
    tbx5_baits[["id_upper", "val_diff_iso2_l5", "val_diff_iso3_l5"]]
    .set_index("id_upper")
    .T
)
annot.index = ["TBX5-2", "TBX5-3"]
annot.loc["TBX5-1"] = ["NA"] * len(tbx5_baits)
annot = annot.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
annot


# In[50]:


fig, ax = plt.subplots(1, 1, figsize=(5, 5))
y1h_pdi_per_tf_gene_plot(
    "TBX5",
    data=y1h,
    ax=ax,
    iso_order=["TBX5-1", "TBX5-2", "TBX5-3"],
    bait_colors=colors,
    bait_annot=annot,
)

plt.colorbar(
    m,
    ax=ax,
    orientation="horizontal",
    shrink=0.75,
    label="largest ∆ PBM affinity across 8-mers in bait",
)
plt.savefig(
    "../../figures/TBX5_y1h_with_pbm_last50.pdf", bbox_inches="tight", dpi="figure"
)
