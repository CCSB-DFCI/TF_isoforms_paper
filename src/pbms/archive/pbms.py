#!/usr/bin/env python
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
from plotting import y1h_pdi_per_tf_gene_plot, m1h_activation_per_tf_gene_plot, COLOR_PURPLE
from data_loading import (load_annotated_TFiso1_collection, 
                          load_y1h_pdi_data, 
                          load_developmental_tissue_expression_remapped,
                          load_gtex_remapped)


# In[2]:


PAPER_PRESET = {"style": "ticks", "font": "Helvetica", "context": "paper", 
                "rc": {"font.size":7,"axes.titlesize":7,
                       "axes.labelsize":7, 'axes.linewidth':0.5,
                       "legend.fontsize":6, "xtick.labelsize":6,
                       "ytick.labelsize":6, "xtick.major.size": 3.0,
                       "ytick.major.size": 3.0, "axes.edgecolor": "black",
                       "xtick.major.pad": 3.0, "ytick.major.pad": 3.0}}
PAPER_FONTSIZE = 7


# In[3]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[4]:


np.random.seed(2023)


# In[5]:


tfs = load_annotated_TFiso1_collection()


# In[6]:


y1h_baits_f = "../../data/internal/Y1H_DNA_baits.fa"


# ## 1. import baits and y1h results

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


#y1h = load_y1h_pdi_data(add_missing_data=True) # erroring out
y1h = load_y1h_pdi_data()
y1h


# ## CREB1

# In[9]:


kfit_dat = pd.read_table("../../data/processed/pbms/CREB1kfit_dat.csv", sep=",")
da_dat = pd.read_table("../../data/processed/pbms/CREB1da_dat.csv", sep=",")


# In[10]:


kfit_ref = kfit_dat[kfit_dat["cname"] == "CREB1-ref"]
kfit_alt = kfit_dat[kfit_dat["cname"] == "CREB1-alt"]


# In[11]:


kfit_vs = kfit_ref.merge(kfit_alt, on=["seq"], suffixes=("_ref", "_alt"))


# In[12]:


da_alt = da_dat[da_dat["cname"] == "CREB1-alt"]


# In[13]:


kfit_vs = kfit_vs.merge(da_alt[["seq", "contrastQ", "contrastQ_cut", "contrastAverage", "contrastDifference"]], on="seq")


# In[14]:


kfit_vs.contains_any_motif_ref.value_counts()


# In[15]:


kfit_vs.affinityEstimate_ref.max()


# In[16]:


creb1_y1h = (y1h.loc[y1h["gene_symbol"] == "CREB1", y1h.columns[1:]].copy().set_index("clone_acc"))
creb1_y1h = creb1_y1h.loc[:, creb1_y1h.any(axis=0)]
creb1_y1h


# In[17]:


creb1_baits = list(creb1_y1h.columns)
creb1_baits = dna[dna["id_upper"].isin(creb1_baits)]
creb1_baits


# In[18]:


kfit_vs.columns


# In[19]:


fig, ax = plt.subplots(1, 1, figsize=(1.2, 1.2))

markers = [",", "."]
titles = ["CREB1-1"]

for k, motif in enumerate(["*other k-mer", "CREB1 k-mer"]):
    for j, qval in enumerate(["(0.1,1]", "(0.01,0.1]", "(0.001,0.01]", "[0,0.001]"]):

        sub = kfit_vs[(kfit_vs["contrastQ_cut"] == qval) & 
                      (kfit_vs["contains_any_motif_ref"] == motif)]
        xs = sub["affinityEstimate_ref"]
        ys = sub["affinityEstimate_alt"]
        ts = sub["contrastQ"]
        
        color = sns.color_palette("Spectral_r", n_colors=4)[j]
        marker = markers[k]

        if marker == "o":
            ax.scatter(xs, ys, 30, marker=".", edgecolors="black", facecolors=color, alpha=1, linewidth=0.5,
                       zorder=10)
        elif marker == ",":
            ax.scatter(xs, ys, 30, marker=".", edgecolors=color, facecolors=color, alpha=0.5,
                       zorder=10)
        else:
            ax.scatter(xs, ys, 30, marker=marker, edgecolors=color, facecolors='white', alpha=1, linewidth=0.5,
                       zorder=10)

ax.set_xlim((6.8, 11.6))
ax.set_xticks([7, 8, 9, 10, 11])
ax.set_ylim((6.8, 11.6))
ax.set_yticks([7, 8, 9, 10, 11])
ax.plot([6.8, 11.6], [6.8, 11.6], color="black", linestyle="dashed", linewidth=1, zorder=1)
# ax.set_xticks([6, 8, 10, 12])
# ax.set_yticks([6, 8, 10, 12])
ax.set_xlabel("CREB1 reference")
ax.set_ylabel("CREB1 alternative")
ax.set_title("CREB1-1")

for loc in ['top', 'right']:
    ax.spines[loc].set_visible(False)
    
fig.savefig("../../figures/fig3/CREB1_isoforms_pbm_scatter.pdf", dpi="figure", bbox_inches="tight")


# In[20]:


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


# In[21]:


def vals_per_bait(row, kfit_vs, colname, alt_suffix, ascending):
    
    kmers = []
    seq = row.seq
    seq_len = row.seq_len
    for i in range(0, seq_len-8):
        kmer = seq[i:i+8]
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


# In[22]:


creb1_baits["tmp"] = creb1_baits.apply(vals_per_bait, axis=1, 
                                       kfit_vs=kfit_vs, colname="contrastDifference_alt", alt_suffix="alt",
                                       ascending=False)
creb1_baits["val_kmer"] = creb1_baits["tmp"].str.split("_", expand=True)[0].astype(str)
creb1_baits["val_diff"] = creb1_baits["tmp"].str.split("_", expand=True)[1].astype(float)
creb1_baits["xval"] = creb1_baits["tmp"].str.split("_", expand=True)[2].astype(float)
creb1_baits["yval"] = creb1_baits["tmp"].str.split("_", expand=True)[3].astype(float)
creb1_baits["rgb"] = creb1_baits["tmp"].str.split("_", expand=True)[4].astype(str)
creb1_baits["rgb"] = creb1_baits["rgb"].apply(literal_eval)
creb1_baits


# In[23]:


colors = creb1_baits[["id_upper", "rgb"]].set_index("id_upper").T
colors.index = ["CREB1-1"]
colors.loc["CREB1-2"] = [(0, 0, 0, 1)] * len(creb1_baits)
colors = colors.loc[["CREB1-2", "CREB1-1"]]
colors


# In[24]:


annot = creb1_baits[["id_upper", "val_diff"]].set_index("id_upper").T
annot.index = ["CREB1-1"]
annot.loc["CREB1-2"] = ["NA"] * len(creb1_baits)
annot = annot.loc[["CREB1-2", "CREB1-1"]]
annot


# In[25]:


fig, ax = plt.subplots(1, 1, figsize=(2, 2))
y1h_pdi_per_tf_gene_plot("CREB1", data=y1h, ax=ax, 
                         iso_order=["CREB1-2", "CREB1-1"], bait_colors=colors, bait_annot=annot)

plt.colorbar(m, ax=ax, orientation="horizontal", shrink=0.75, label="largest ∆ PBM affinity\nacross 8-mers in bait")
plt.savefig('../../figures/fig3/CREB1_y1h_with_pbm.pdf', bbox_inches='tight', dpi="figure")


# ## TBX5

# In[26]:


kfit_dat = pd.read_table("../../data/processed/pbms/TBX5kfit_dat.csv", sep=",")
da_dat = pd.read_table("../../data/processed/pbms/TBX5da_dat.csv", sep=",")


# In[27]:


kfit_ref = kfit_dat[kfit_dat["cname"] == "TBX5C05-REF"]
kfit_3 = kfit_dat[kfit_dat["cname"] == "TBX5A05"]
kfit_2 = kfit_dat[kfit_dat["cname"] == "TBX5B05"]


# In[28]:


kfit_vs = kfit_ref.merge(kfit_2, on=["seq"], suffixes=("_ref", ""))
kfit_vs = kfit_vs.merge(kfit_3, on=["seq"], suffixes=("_iso2", "_iso3"))
kfit_vs.head()


# In[29]:


da_3 = da_dat[da_dat["cname"] == "TBX5A05"]
da_2 = da_dat[da_dat["cname"] == "TBX5B05"]


# In[30]:


kfit_vs = kfit_vs.merge(da_2[["seq", "contrastQ", "contrastQ_cut"]], on="seq")
kfit_vs = kfit_vs.merge(da_3[["seq", "contrastQ", "contrastQ_cut"]], on="seq", suffixes=("_iso2", "_iso3"))


# In[31]:


kfit_vs.contains_any_motif_ref.value_counts()


# In[32]:


def add_motif_category(row, motif_name, motif, motif_rc):
    if motif in row.seq or motif_rc in row.seq:
        return motif_name
    else:
        return row.contains_any_motif_ref

kfit_vs["contains_any_motif_ref"] = kfit_vs.apply(add_motif_category, motif_name="* ACGTGT k-mer",
                                                  motif="ACGTGT", motif_rc="ACACGT", axis=1)
kfit_vs.contains_any_motif_ref.value_counts()


# In[33]:


fig, axarr = plt.subplots(1, 2, figsize=(3, 1.5))

markers = [",", ".", "o"]
titles = ["TBX5-2", "TBX5-3"]

for i, suffix in enumerate(["iso2", "iso3"]):
    
    ax = axarr[i]

    for k, motif in enumerate(["other k-mer", "* TBX5 k-mer", "* ACGTGT k-mer"]):
        for j, qval in enumerate(["(0.1,1]", "(0.01,0.1]", "(0.001,0.01]", "[0,0.001]"]):
            
            sub = kfit_vs[(kfit_vs["contrastQ_cut_%s" % suffix] == qval) & 
                          (kfit_vs["contains_any_motif_ref"] == motif)]
            xs = sub["affinityEstimate_ref"]
            ys = sub["affinityEstimate_%s" % suffix]
            
            color = sns.color_palette("Spectral_r", n_colors=4)[j]
            marker = markers[k]
            
            if marker == "o":
                ax.scatter(xs, ys, 15, marker=marker, edgecolors="black", facecolors=color, alpha=1, linewidth=0.5)
            elif marker == ",":
                ax.scatter(xs, ys, 1, marker=".", edgecolors=color, facecolors='none', alpha=0.5)
            else:
                ax.scatter(xs, ys, 30, marker=marker, edgecolors=color, facecolors='white', alpha=1)
                
    
    ax.set_xlim((9.2, 13))
    ax.set_ylim((9.2, 13))
    ax.plot([9.2, 13], [9.2, 13], color="black", linestyle="dashed", linewidth=1, zorder=1)
    ax.set_xticks([10, 11, 12, 13])
    ax.set_yticks([10, 11, 12, 13])
    ax.set_xlabel("reference PBM affinity")
    ax.set_title(titles[i])
    
axarr[0].set_ylabel("alternative PBM affinity")

for loc in ['top', 'right']:
    axarr[0].spines[loc].set_visible(False)
    axarr[1].spines[loc].set_visible(False)
    
fig.savefig("../../figures/fig3/TBX5_isoforms_pbm_scatter.pdf", dpi="figure", bbox_inches="tight")


# In[34]:


# add colors to kfit_vs
cmap_name = "Spectral"
colname = "contrastDifference_iso2"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(vmin=kfit_vs[colname].min(), vmax=kfit_vs[colname].max())
m = cm.ScalarMappable(norm=norm, cmap=cmap)

kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)
kfit_vs.sample(5)


# In[35]:


kfit_vs[["contrastDifference_iso2", "contrastDifference_iso3"]].min().min()


# In[36]:


# add colors to kfit_vs
cmap_name = "Spectral"

cmap = cm.get_cmap(cmap_name)
norm = mpl.colors.Normalize(vmin=kfit_vs[["contrastDifference_iso2", "contrastDifference_iso3"]].min().min(), 
                            vmax=kfit_vs[["contrastDifference_iso2", "contrastDifference_iso3"]].max().max())
m = cm.ScalarMappable(norm=norm, cmap=cmap)

colname = "contrastDifference_iso2"
kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)

colname = "contrastDifference_iso3"
kfit_vs["%s_rgb" % colname] = kfit_vs.apply(get_rgb, colname=colname, m=m, axis=1)
kfit_vs.sample(5)


# In[37]:


tbx5_y1h = (y1h.loc[y1h["gene_symbol"] == "TBX5", y1h.columns[1:]].copy().set_index("clone_acc"))
tbx5_y1h = tbx5_y1h.loc[:, tbx5_y1h.any(axis=0)]
tbx5_y1h


# In[38]:


tbx5_baits = list(tbx5_y1h.columns)
tbx5_baits = dna[dna["id_upper"].isin(tbx5_baits)]
tbx5_baits


# In[39]:


tbx5_baits["tmp_iso2"] = tbx5_baits.apply(vals_per_bait, axis=1, 
                                       kfit_vs=kfit_vs, colname="contrastDifference_iso2", alt_suffix="iso2",
                                       ascending=False)
tbx5_baits["val_kmer_iso2"] = tbx5_baits["tmp_iso2"].str.split("_", expand=True)[0].astype(str)
tbx5_baits["val_diff_iso2"] = tbx5_baits["tmp_iso2"].str.split("_", expand=True)[1].astype(float)
tbx5_baits["xval_iso2"] = tbx5_baits["tmp_iso2"].str.split("_", expand=True)[2].astype(float)
tbx5_baits["yval_iso2"] = tbx5_baits["tmp_iso2"].str.split("_", expand=True)[3].astype(float)
tbx5_baits["rgb_iso2"] = tbx5_baits["tmp_iso2"].str.split("_", expand=True)[4].astype(str)
tbx5_baits["rgb_iso2"] = tbx5_baits["rgb_iso2"].apply(literal_eval)


# In[40]:


tbx5_baits["tmp_iso3"] = tbx5_baits.apply(vals_per_bait, axis=1, 
                                       kfit_vs=kfit_vs, colname="contrastDifference_iso3", alt_suffix="iso3",
                                       ascending=False)
tbx5_baits["val_kmer_iso3"] = tbx5_baits["tmp_iso3"].str.split("_", expand=True)[0].astype(str)
tbx5_baits["val_diff_iso3"] = tbx5_baits["tmp_iso3"].str.split("_", expand=True)[1].astype(float)
tbx5_baits["xval_iso3"] = tbx5_baits["tmp_iso3"].str.split("_", expand=True)[2].astype(float)
tbx5_baits["yval_iso3"] = tbx5_baits["tmp_iso3"].str.split("_", expand=True)[3].astype(float)
tbx5_baits["rgb_iso3"] = tbx5_baits["tmp_iso3"].str.split("_", expand=True)[4].astype(str)
tbx5_baits["rgb_iso3"] = tbx5_baits["rgb_iso3"].apply(literal_eval)


# In[41]:


colors = tbx5_baits[["id_upper", "rgb_iso2", "rgb_iso3"]].set_index("id_upper").T
colors.index = ["TBX5-2", "TBX5-3"]
colors.loc["TBX5-1"] = [(0, 0, 0, 1)] * len(tbx5_baits)
colors = colors.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
colors


# In[42]:


annot = tbx5_baits[["id_upper", "val_diff_iso2", "val_diff_iso3"]].set_index("id_upper").T
annot.index = ["TBX5-2", "TBX5-3"]
annot.loc["TBX5-1"] = ["NA"] * len(tbx5_baits)
annot = annot.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
annot


# In[43]:


fig, ax = plt.subplots(1, 1, figsize=(5, 5))
y1h_pdi_per_tf_gene_plot("TBX5", data=y1h, ax=ax, 
                         iso_order=["TBX5-1", "TBX5-2", "TBX5-3"], bait_colors=colors, bait_annot=annot)

plt.colorbar(m, ax=ax, orientation="horizontal", shrink=0.75, label="largest ∆ PBM affinity across 8-mers in bait")
plt.savefig('../../figures/fig3/TBX5_y1h_with_pbm.pdf', bbox_inches='tight', dpi="figure")


# In[44]:


np.mean([-0.72, -0.38, -0.88, -0.60, -0.62, -0.16])


# In[45]:


np.mean([-0.36, -0.88, -0.72, -0.53, -0.72, -0.51, -0.13])


# In[46]:


def vals_per_bait_last50(row, kfit_vs, colname, alt_suffix, ascending):
    
    kmers = []
    seq = row.seq
    seq_len = row.seq_len
    for i in range(seq_len-100, seq_len-8):
        kmer = seq[i:i+8]
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


# In[47]:


tbx5_baits["tmp_iso2_l5"] = tbx5_baits.apply(vals_per_bait_last50, axis=1, 
                                       kfit_vs=kfit_vs, colname="contrastDifference_iso2", alt_suffix="iso2",
                                       ascending=False)
tbx5_baits["val_kmer_iso2_l5"] = tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[0].astype(str)
tbx5_baits["val_diff_iso2_l5"] = tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[1].astype(float)
tbx5_baits["xval_iso2_l5"] = tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[2].astype(float)
tbx5_baits["yval_iso2_l5"] = tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[3].astype(float)
tbx5_baits["rgb_iso2_l5"] = tbx5_baits["tmp_iso2_l5"].str.split("_", expand=True)[4].astype(str)
tbx5_baits["rgb_iso2_l5"] = tbx5_baits["rgb_iso2_l5"].apply(literal_eval)


# In[48]:


tbx5_baits["tmp_iso3_l5"] = tbx5_baits.apply(vals_per_bait_last50, axis=1, 
                                       kfit_vs=kfit_vs, colname="contrastDifference_iso3", alt_suffix="iso3",
                                       ascending=False)
tbx5_baits["val_kmer_iso3_l5"] = tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[0].astype(str)
tbx5_baits["val_diff_iso3_l5"] = tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[1].astype(float)
tbx5_baits["xval_iso3_l5"] = tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[2].astype(float)
tbx5_baits["yval_iso3_l5"] = tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[3].astype(float)
tbx5_baits["rgb_iso3_l5"] = tbx5_baits["tmp_iso3_l5"].str.split("_", expand=True)[4].astype(str)
tbx5_baits["rgb_iso3_l5"] = tbx5_baits["rgb_iso3_l5"].apply(literal_eval)


# In[49]:


colors = tbx5_baits[["id_upper", "rgb_iso2_l5", "rgb_iso3_l5"]].set_index("id_upper").T
colors.index = ["TBX5-2", "TBX5-3"]
colors.loc["TBX5-1"] = [(0, 0, 0, 1)] * len(tbx5_baits)
colors = colors.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
colors


# In[50]:


annot = tbx5_baits[["id_upper", "val_diff_iso2_l5", "val_diff_iso3_l5"]].set_index("id_upper").T
annot.index = ["TBX5-2", "TBX5-3"]
annot.loc["TBX5-1"] = ["NA"] * len(tbx5_baits)
annot = annot.loc[["TBX5-1", "TBX5-2", "TBX5-3"]]
annot


# In[51]:


fig, ax = plt.subplots(1, 1, figsize=(5, 5))
y1h_pdi_per_tf_gene_plot("TBX5", data=y1h, ax=ax, 
                         iso_order=["TBX5-1", "TBX5-2", "TBX5-3"], bait_colors=colors, bait_annot=annot)

plt.colorbar(m, ax=ax, orientation="horizontal", shrink=0.75, label="largest ∆ PBM affinity across 8-mers in bait")
plt.savefig('../../figures/TBX5_y1h_with_pbm_last50.pdf', bbox_inches='tight', dpi="figure")


# ### expression levels of TBX5

# In[52]:


def developmental_tissue_expression_plot(gene_name, figsize, ylim, means, cols, fig_suffix):
    locs = [x for x in list(means.index) if x.split("|")[0] == gene_name]
    
    # include isos that aren't cloned
    locs = list(set(locs + [x for x in list(means.index) if x.split(" ")[1][:-4] == gene_name]))
    
    n_isos = len(means.loc[locs])
    palette = sns.color_palette("Set2")
    fig, axes = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(figsize)
    ### bar chart ###
    (means.loc[locs, cols]
          .T
          .plot.bar(ax=axes[0],
                    legend=False,
                    width=0.7,
                    color=list(palette)))
    ### percentages ###
    raw_means = 2 ** means.loc[locs, cols] - 1.
    (raw_means.div(raw_means.sum(axis=0))
              .T.plot.bar(ax=axes[1], 
                          stacked=True,
                          legend=False,
                          color=list(palette)))
    axes[0].set_ylabel('log2(tpm + 1)\n')
    axes[0].set_ylim(ylim)
    axes[1].set_ylabel('percent')
    axes[1].set_yticklabels(['{:.0%}'.format(t) for t in axes[1].get_yticks()])
    axes[1].legend(loc='lower left', bbox_to_anchor=(1, 0), frameon=False)
    axes[0].axhline(y=1, color='black', linewidth=0.5, linestyle="dashed")
    
    axes[0].spines['top'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    
    plt.subplots_adjust(hspace=0.25)
    plt.savefig('../../figures/fig3/expression_' + gene_name + '_' + fig_suffix + '.pdf',
                bbox_inches='tight')


# In[53]:


df_dev, metadata_dev, genes_dev = load_developmental_tissue_expression_remapped()

rename_dev_stage = {'8 week post conception,embryo': '08',
'11 week post conception,late embryo': '11',
'embryo,7 week post conception': '07',
'infant': 'infant',
'10 week post conception,late embryo': '10',
'young adult': 'young adult',
'13 week post conception,late embryo': '13',
'16 week post conception,late embryo': '16',
'4 week post conception,embryo': '04',
'neonate': 'neonate',
'19 week post conception,late embryo': '19',
'9 week post conception,late embryo': '09',
'adolescent': 'adolescent',
'5 week post conception,embryo': '05',
'embryo,6 week post conception': '06',
'12 week post conception,late embryo': '12',
'18 week post conception,late embryo': '18',
'toddler': 'toddler',
'elderly': 'elderly',
'middle adult': 'adult',
'school age child': 'child'}

metadata_dev['dev_stage'] = metadata_dev['Developmental_Stage'].map(rename_dev_stage)
means_dev = (df_dev.groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']), axis=1)
           .mean())


# In[54]:


df_gtex, metadata_gtex, genes_gtex = load_gtex_remapped()

exclusion_list_gtex = {'Cells - Leukemia cell line (CML)',
                       'Cells - EBV-transformed lymphocytes',
                       'Cells - Cultured fibroblasts'}

df_gtex = df_gtex.loc[:, ~df_gtex.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex = metadata_gtex.loc[~metadata_gtex['body_site'].isin(exclusion_list_gtex), :]

means_gtex = df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean()


# In[55]:


heart_cols = [x for x in means_dev.columns if "heart" in x]
ovary_cols = [x for x in means_dev.columns if "ovary" in x]
developmental_tissue_expression_plot("TBX5", (3.5, 1.75), (0, 12), means_dev, heart_cols, 
                                     "means_dev_heart")


# In[56]:


heart_cols = [x for x in means_gtex.columns if "Heart" in x]
developmental_tissue_expression_plot("TBX5", (0.5, 1.75), (0, 12), means_gtex, heart_cols, 
                                     "means_gtex_heart")


# In[ ]:




