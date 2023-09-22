#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings

warnings.filterwarnings("ignore")


# In[2]:


import matplotlib as mpl
import matplotlib.pyplot as plt
import met_brewer
import pandas as pd
import numpy as np
import seaborn as sns
import sys
import upsetplot

import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smt

from Bio.Seq import Seq
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
from scipy.stats import pearsonr

import plotting
from plotting import PAPER_PRESET, PAPER_FONTSIZE, annotate_pval, mimic_r_boxplot


get_ipython().run_line_magic("matplotlib", "inline")
get_ipython().run_line_magic("config", "InlineBackend.figure_format = 'svg'")
mpl.rcParams["figure.autolayout"] = False


# In[3]:


from data_loading import (
    load_annotated_6k_collection,
    load_annotated_gencode_tfs,
    load_valid_isoform_clones,
    load_developmental_tissue_expression_gencode,
)


# In[4]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[5]:


pal = {
    "ref": sns.color_palette("Set2")[0],
    "rewire": sns.color_palette("Set2")[2],
    "DN": sns.color_palette("Set2")[1],
    "NA": "lightgray",
    "likely": "darkgray",
}


# ## variables

# In[6]:


dn_data_f = "../data/processed/DN_cats_Joung.tsv"


# In[7]:


brca_cnts_f = (
    "../data/processed/Nathans_analysis/Breast_cancer/isoCounts.BreastCancer.txt"
)
brca_tx_f = (
    "../data/processed/Nathans_analysis/Breast_cancer/transcript.BreastCancer.txt"
)
pam50_f = "../data/processed/Nathans_analysis/Breast_cancer/groups.PAM50.txt"
bulk_f = (
    "../data/processed/Nathans_analysis/Breast_cancer/groups.BreastCancer_ratios.txt"
)


# ## 1. import data

# In[8]:


dn_data = pd.read_table(dn_data_f)
dn_data.head()


# In[9]:


skiprows = list(range(96320, 96387)) + list(range(99680, 99687))
brca = pd.read_table(brca_tx_f, sep="\t", skiprows=skiprows)
brca.shape


# In[10]:


pam50_samps = pd.read_table(pam50_f, header=None)
pam50_samps.columns = ["file", "samp_type_id", "samp_type"]
pam50_samps["tcga_id"] = pam50_samps["file"].str.split(".", expand=True)[0]
pam50_samps.samp_type.value_counts()


# In[11]:


bulk_samps = pd.read_table(bulk_f, header=None)
bulk_samps.columns = ["file", "samp_type_id", "samp_type"]
bulk_samps["tcga_id"] = bulk_samps["file"].str.split(".", expand=True)[0]
bulk_samps.samp_type.value_counts()


# ## 2. map sample types

# In[12]:


brca_samps = list(bulk_samps[bulk_samps["samp_type"] != "controls"]["tcga_id"])
print("# breast cancer samples: %s" % len(brca_samps))


# same ctrls in both
ctrl_samps = list(pam50_samps[pam50_samps["samp_type"] == "controls"]["tcga_id"])
print("# control samples: %s" % len(ctrl_samps))

luma_samps = list(pam50_samps[pam50_samps["samp_type"] == "Luminal A"]["tcga_id"])
print("# Luminal A samples: %s" % len(luma_samps))

lumb_samps = list(pam50_samps[pam50_samps["samp_type"] == "Luminal B"]["tcga_id"])
print("# Luminal B samples: %s" % len(lumb_samps))

tn_samps = list(pam50_samps[pam50_samps["samp_type"] == "Basal-like"]["tcga_id"])
print("# Basal-like samples: %s" % len(tn_samps))

her2_samps = list(pam50_samps[pam50_samps["samp_type"] == "HER2-enriched"]["tcga_id"])
print("# HER2-enriched samples: %s" % len(her2_samps))

norm_samps = list(pam50_samps[pam50_samps["samp_type"] == "Normal-like"]["tcga_id"])
print("# Normal-like samples: %s" % len(norm_samps))


# In[13]:


# one brca samp is weirdly missing, remove
brca_samps = [x for x in brca_samps if x in brca.columns]
len(brca_samps)


# In[14]:


## patient is the 3rd value in the barcode
## source: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
bulk_samps["patient_id"] = bulk_samps["tcga_id"].str.split("-", expand=True)[2]
pam50_samps["patient_id"] = pam50_samps["tcga_id"].str.split("-", expand=True)[2]


# In[15]:


tcga_samps = bulk_samps.merge(
    pam50_samps, on=["tcga_id", "patient_id"], how="outer", suffixes=("_brca", "_pam50")
)
print(len(tcga_samps))


# In[16]:


tcga_ctrls = tcga_samps[
    (tcga_samps["samp_type_brca"] == "controls")
    | (tcga_samps["samp_type_pam50"] == "controls")
]
len(tcga_ctrls)


# In[17]:


tcga_tumors = tcga_samps[
    (tcga_samps["samp_type_brca"] != "controls")
    | (tcga_samps["samp_type_pam50"] != "controls")
]
len(tcga_tumors)


# In[18]:


tcga_paired = tcga_ctrls.merge(
    tcga_tumors, on=["patient_id"], suffixes=("_ctrl", "_tumor")
)
print(len(tcga_paired))


# ## 3. aggregate TF iso expression across transcripts

# In[19]:


tfs = load_annotated_6k_collection()


# In[20]:


tf_id_map = pd.DataFrame()
gene_names = []
iso_ids = []
clone_accs = []
enst_ids = []

for tf in tfs:
    db = tfs[tf]
    for iso in db.isoforms:
        enst_id = iso.ensembl_transcript_ids
        try:
            clone_acc = iso.clone_acc
        except:
            clone_acc = "none"

        if enst_id is None:
            gene_names.append(tf)
            iso_ids.append(iso.name)
            clone_accs.append(clone_acc)
            enst_ids.append("none")
        else:
            for enst_id in iso.ensembl_transcript_ids:
                gene_names.append(tf)
                iso_ids.append(iso.name)
                clone_accs.append(clone_acc)
                enst_ids.append(enst_id)

tf_id_map["gene_name"] = gene_names
tf_id_map["iso_id"] = iso_ids
tf_id_map["clone_acc"] = clone_accs
tf_id_map["enst_id"] = enst_ids
print(len(tf_id_map))
tf_id_map.sample(5)


# In[21]:


def merge_id(row):
    if row.enst_id == "none":
        return row.clone_acc
    else:
        return row.enst_id


tf_id_map["merge_id"] = tf_id_map.apply(merge_id, axis=1)
tf_id_map.sample(5)


# In[22]:


dd = tf_id_map[["iso_id", "gene_name"]].drop_duplicates()
print(len(dd))
gene_dict = {row.iso_id: row.gene_name for i, row in dd.iterrows()}


# In[23]:


brca_cols = [x for x in brca.columns if x != "UID"]
len(brca_cols)


# In[24]:


brca = brca.merge(tf_id_map, left_on="UID", right_on="merge_id")
len(brca)


# In[25]:


brca_isos = brca.groupby("iso_id")[brca_cols].agg("sum").reset_index()
len(brca_isos)


# In[26]:


brca_isos.head()


# ## 4. calculate isoform ratios

# In[27]:


# calculate isoform ratios, set anything w gene-level exp <= 1 to nan
brca_genes = pd.Series(
    index=brca_isos.iso_id, data=brca_isos.iso_id.map(gene_dict).values
)

brca_idx = brca_isos.set_index("iso_id", inplace=False)
brca_idx = brca_idx[brca_cols]
brca_gene_sum = brca_idx.groupby(brca_genes).transform("sum")

f_brca = brca_idx / brca_gene_sum
f_brca_nan = f_brca * (brca_gene_sum >= 1).applymap(
    lambda x: {False: np.nan, True: 1}[x]
)


# ## 5. calculate median expression + ratio in each PAM50 type

# In[28]:


tcga_paired_ctrls = list(tcga_paired["tcga_id_ctrl"].unique())
tcga_paired_tumors = list(tcga_paired["tcga_id_tumor"].unique())


# In[29]:


brca_isos["med_brca_tpm"] = brca_isos[brca_samps].median(axis=1)
brca_isos["med_ctrl_tpm"] = brca_isos[ctrl_samps].median(axis=1)
brca_isos["med_luma_tpm"] = brca_isos[luma_samps].median(axis=1)
brca_isos["med_lumb_tpm"] = brca_isos[lumb_samps].median(axis=1)
brca_isos["med_tn_tpm"] = brca_isos[tn_samps].median(axis=1)
brca_isos["med_her2_tpm"] = brca_isos[her2_samps].median(axis=1)
brca_isos["med_norm_tpm"] = brca_isos[norm_samps].median(axis=1)
brca_isos["med_paired-brca_tpm"] = brca_isos[tcga_paired_tumors].median(axis=1)
brca_isos["med_paired-ctrls_tpm"] = brca_isos[tcga_paired_ctrls].median(axis=1)


# In[30]:


f_brca_nan["med_brca_rationan"] = f_brca_nan[brca_samps].median(axis=1)
f_brca_nan["med_ctrl_rationan"] = f_brca_nan[ctrl_samps].median(axis=1)
f_brca_nan["med_luma_rationan"] = f_brca_nan[luma_samps].median(axis=1)
f_brca_nan["med_lumb_rationan"] = f_brca_nan[lumb_samps].median(axis=1)
f_brca_nan["med_tn_rationan"] = f_brca_nan[tn_samps].median(axis=1)
f_brca_nan["med_her2_rationan"] = f_brca_nan[her2_samps].median(axis=1)
f_brca_nan["med_norm_rationan"] = f_brca_nan[norm_samps].median(axis=1)
f_brca_nan["med_paired-brca_rationan"] = f_brca_nan[tcga_paired_tumors].median(axis=1)
f_brca_nan["med_paired-ctrls_rationan"] = f_brca_nan[tcga_paired_ctrls].median(axis=1)


# In[31]:


f_brca["med_brca_ratio"] = f_brca[brca_samps].median(axis=1)
f_brca["med_ctrl_ratio"] = f_brca[ctrl_samps].median(axis=1)
f_brca["med_luma_ratio"] = f_brca[luma_samps].median(axis=1)
f_brca["med_lumb_ratio"] = f_brca[lumb_samps].median(axis=1)
f_brca["med_tn_ratio"] = f_brca[tn_samps].median(axis=1)
f_brca["med_her2_ratio"] = f_brca[her2_samps].median(axis=1)
f_brca["med_norm_ratio"] = f_brca[norm_samps].median(axis=1)
f_brca["med_paired-brca_ratio"] = f_brca[tcga_paired_tumors].median(axis=1)
f_brca["med_paired-ctrls_ratio"] = f_brca[tcga_paired_ctrls].median(axis=1)


# ## 7. calculate expression/ratio change per isoform across paired samps

# In[32]:


paired_ctrl_samps = list(tcga_paired["tcga_id_ctrl"])
print(len(paired_ctrl_samps))
paired_tumor_samps = list(tcga_paired["tcga_id_tumor"])
print(len(paired_tumor_samps))


# In[33]:


## calculate p-value using wilcoxon
def paired_pval(row, ctrl_cols, tumor_cols):
    x = row[ctrl_cols]
    x = [x for x in x if not pd.isnull(x)]
    y = row[tumor_cols]
    y = [y for y in y if not pd.isnull(y)]

    try:
        stat, p = wilcoxon(x, y)
        return p
    except:
        return np.nan


## calculate p-value using wilcoxon
def paired_stat(row, ctrl_cols, tumor_cols):
    x = row[ctrl_cols]
    y = row[tumor_cols]

    try:
        stat, p = wilcoxon(x, y)
        return stat
    except:
        return np.nan


f_brca["wilcox_pval"] = f_brca.apply(
    paired_pval, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1
)
f_brca["wilcox_stat"] = f_brca.apply(
    paired_stat, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1
)
print(len(f_brca))

f_brca_filt = f_brca[~pd.isnull(f_brca["wilcox_pval"])]
print(len(f_brca_filt))

f_brca_filt["wilcox_padj"] = smt.multipletests(
    list(f_brca_filt["wilcox_pval"]), alpha=0.05, method="fdr_bh"
)[1]

f_brca_nan["wilcox_pval"] = f_brca_nan.apply(
    paired_pval, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1
)
f_brca_nan["wilcox_stat"] = f_brca_nan.apply(
    paired_stat, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1
)
print(len(f_brca_nan))

f_brca_nan_filt = f_brca_nan[~pd.isnull(f_brca_nan["wilcox_pval"])]
print(len(f_brca_nan_filt))

f_brca_nan_filt["wilcox_padj"] = smt.multipletests(
    list(f_brca_nan_filt["wilcox_pval"]), alpha=0.05, method="fdr_bh"
)[1]


# In[34]:


for i, row in tcga_paired.iterrows():
    f_brca_filt["paired-diff_%s_ratio" % (i + 1)] = (
        f_brca_filt[row.tcga_id_tumor] - f_brca_filt[row.tcga_id_ctrl]
    )
    f_brca_nan_filt["paired-diff_%s_rationan" % (i + 1)] = f_brca_nan_filt[
        row.tcga_id_tumor
    ].fillna(0) - f_brca_nan[row.tcga_id_ctrl].fillna(0)


# In[35]:


paired_ratio_cols = [x for x in f_brca_filt.columns if "paired-diff_" in x]
paired_rationan_cols = [x for x in f_brca_nan_filt.columns if "paired-diff_" in x]


# In[36]:


f_brca_filt["med_paired-diff_ratio"] = f_brca_filt[paired_ratio_cols].median(axis=1)
f_brca_nan_filt["med_paired-diff_rationan"] = f_brca_nan_filt[
    paired_rationan_cols
].median(axis=1)


# ## 6. merge median expression w/ DN cats

# In[37]:


f_brca_filt = f_brca_filt.merge(tf_id_map, on="iso_id")
f_brca_nan_filt = f_brca_nan_filt.merge(tf_id_map, on="iso_id")


# In[38]:


f_brca_nan_med_cols = (
    ["clone_acc"]
    + [x for x in f_brca_nan_filt.columns if "med_" in x]
    + [x for x in f_brca_nan_filt.columns if "wilcox" in x]
)


# In[39]:


f_brca_med_cols = (
    ["clone_acc"]
    + [x for x in f_brca_filt.columns if "med_" in x]
    + [x for x in f_brca_filt.columns if "wilcox" in x]
)


# In[40]:


dn_data_exp = dn_data.merge(
    f_brca_nan_filt[f_brca_nan_med_cols], left_on="tf1p0_id", right_on="clone_acc"
)
dn_data_exp = dn_data_exp.merge(
    f_brca_filt[f_brca_med_cols],
    left_on="tf1p0_id",
    right_on="clone_acc",
    suffixes=("_nan", ""),
).drop_duplicates()
print(len(dn_data_exp))


# In[41]:


dn_data_exp["neglog_padj"] = -np.log10(dn_data_exp["wilcox_padj"])
dn_data_exp["neglog_padj_nan"] = -np.log10(dn_data_exp["wilcox_padj_nan"])


# ## 7. plots

# In[42]:


len(
    dn_data_exp[
        (dn_data_exp["wilcox_padj"] < 0.05)
        & (dn_data_exp["dn_cat"].isin(["ref", "rewire", "DN"]))
    ]
)


# In[43]:


len(dn_data_exp[(dn_data_exp["wilcox_padj"] < 0.05) & (dn_data_exp["dn_cat"] == "ref")])


# In[44]:


len(dn_data_exp[(dn_data_exp["wilcox_padj"] < 0.05) & (dn_data_exp["dn_cat"] == "DN")])


# In[45]:


len(
    dn_data_exp[
        (dn_data_exp["wilcox_padj"] < 0.05) & (dn_data_exp["dn_cat"] == "rewire")
    ]
)


# In[46]:


len(dn_data_exp[dn_data_exp["dn_cat"].isin(["ref", "rewire", "DN"])])


# In[47]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(
    data=dn_data_exp[dn_data_exp["dn_cat"].isin(["ref", "rewire", "DN"])],
    x="med_paired-diff_ratio",
    y="neglog_padj",
    hue="dn_cat",
    palette=pal,
    linewidth=0.25,
    edgecolor="black",
    alpha=0.8,
    zorder=10,
    **{"s": 8}
)

ax.set_xlabel("median isoform difference (tumor - normal)")
ax.set_ylabel("-log10(Wilcoxon adjusted p-value)")
ax.set_title("TF isoforms in breast cancer\n")

ax.set_xlim((-0.25, 0.12))
ax.set_ylim((-0.25, 14))
ax.axhline(y=-np.log10(0.05), linestyle="dashed", color="black", linewidth=0.5)
ax.axvline(x=0, linestyle="dashed", color="black", linewidth=0.5)

ax.get_legend().remove()

fig.savefig("../figures/BRCA_Volcano.pdf", dpi="figure", bbox_inches="tight")


# In[48]:


dn_data_exp[dn_data_exp["gene_name"] == "KLF7"][
    ["tf1p0_id", "dn_cat", "med_paired-diff_ratio", "neglog_padj"]
]


# In[50]:


dn_data_exp[dn_data_exp["dn_cat"].isin(["ref", "rewire", "DN"])][
    [
        "gene_name",
        "family",
        "tf1p0_id",
        "dn_cat",
        "med_paired-diff_ratio",
        "neglog_padj",
    ]
].sort_values(by="neglog_padj", ascending=False).head(20)


# In[51]:


brca_isos = brca_isos.merge(
    tf_id_map[["iso_id", "gene_name"]], on="iso_id"
).drop_duplicates()
len(brca_isos)


# In[52]:


brca_isos_paired = brca_isos[
    ["gene_name", "iso_id"] + tcga_paired_ctrls + tcga_paired_tumors
]
new_ctrl_cols = ["normal - %s" % (i + 1) for i, x in enumerate(tcga_paired_ctrls)]
new_tumor_cols = ["tumor - %s" % (i + 1) for i, x in enumerate(tcga_paired_tumors)]
brca_isos_paired.columns = ["gene_name", "iso_id"] + new_ctrl_cols + new_tumor_cols


# In[53]:


def brca_expression_plot(
    gene_name, figsize, ylim, df, cols, fig_suffix, ctrls_line, tumor_line
):
    df_sub = df[df["gene_name"] == gene_name]
    df_sub.set_index("iso_id", inplace=True)
    df_sub = df_sub[cols].drop_duplicates()
    # print(df_sub.head())
    n_isos = len(df_sub)
    palette = sns.color_palette("husl", as_cmap=False, n_colors=n_isos)
    fig, axes = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(figsize)
    ### bar chart ###
    (df_sub.T.plot.bar(ax=axes[0], legend=False, width=0.7, color=list(palette)))
    ### percentages ###
    (
        df_sub.div(df_sub.sum(axis=0)).T.plot.bar(
            ax=axes[1], stacked=True, legend=False, color=list(palette)
        )
    )
    axes[0].set_yscale("symlog")
    axes[0].set_ylabel("tpm")
    # axes[0].set_ylim(ylim)
    axes[1].set_ylabel("percent")
    axes[1].set_yticklabels(["{:.0%}".format(t) for t in axes[1].get_yticks()])
    axes[1].legend(loc="lower left", bbox_to_anchor=(1, 0))
    axes[0].axhline(y=1, color="black", linewidth=0.5, linestyle="dashed")

    # add medians
    axes[1].plot(
        ctrls_line[0], ctrls_line[1], color="black", linewidth=0.5, linestyle="dashed"
    )
    axes[1].plot(
        tumor_line[0], tumor_line[1], color="black", linewidth=0.5, linestyle="dashed"
    )

    plt.subplots_adjust(hspace=0.25)
    plt.savefig(
        "../figures/brca_" + gene_name + "_" + fig_suffix + ".pdf", bbox_inches="tight"
    )


# In[54]:


dn_data_exp[dn_data_exp["gene_name"] == "PKNOX1"][
    [
        "gene_name",
        "tf1p0_id",
        "dn_cat",
        "med_paired-brca_ratio",
        "med_paired-ctrls_ratio",
        "med_paired-diff_ratio",
    ]
]


# In[55]:


cols = new_ctrl_cols[0:35] + new_tumor_cols[0:35]
brca_expression_plot(
    "PKNOX1",
    (9, 3),
    (0, 6),
    brca_isos_paired,
    cols,
    "paired",
    ([0, 34.5], [0.62, 0.62]),
    ([34.5, 70], [0.55, 0.55]),
)


# In[56]:


f_brca_paired = f_brca_filt[
    ["gene_name", "iso_id"] + tcga_paired_ctrls + tcga_paired_tumors
]
new_ctrl_cols = ["normal - %s" % (i + 1) for i, x in enumerate(tcga_paired_ctrls)]
new_tumor_cols = ["tumor - %s" % (i + 1) for i, x in enumerate(tcga_paired_tumors)]
f_brca_paired.columns = ["gene_name", "iso_id"] + new_ctrl_cols + new_tumor_cols


# In[57]:


f_brca_paired_melt = pd.melt(f_brca_paired, id_vars=["gene_name", "iso_id"])
f_brca_paired_melt["samp"] = f_brca_paired_melt["variable"].str.split(" ", expand=True)[
    0
]


# In[58]:


dn_data_exp = dn_data_exp.merge(
    tf_id_map[["gene_name", "clone_acc", "iso_id", "merge_id"]],
    left_on=["gene_name", "tf1p0_id"],
    right_on=["gene_name", "clone_acc"],
).drop_duplicates()
print(len(dn_data_exp))
dn_data_exp[dn_data_exp["gene_name"] == "PKNOX1"]


# In[59]:


tmp = f_brca_paired_melt[f_brca_paired_melt["gene_name"] == "PKNOX1"]

fig = plt.figure(figsize=(2, 2))

ax = sns.boxplot(
    data=tmp,
    x="iso_id",
    y="value",
    hue="samp",
    fliersize=0,
    palette={"normal": "gray", "tumor": sns.color_palette("Set2")[3]},
)
mimic_r_boxplot(ax)

sns.swarmplot(
    data=tmp,
    x="iso_id",
    y="value",
    hue="samp",
    palette={"normal": "gray", "tumor": sns.color_palette("Set2")[3]},
    ax=ax,
    size=1,
    edgecolor="black",
    linewidth=0.5,
    alpha=0.5,
    split=True,
)

# annotate w p-vals
ys = [0.87, 0.7, 0.4, 0.4]
for i, iso in enumerate(tmp.iso_id.unique()):
    print(iso)
    padj = dn_data_exp[dn_data_exp["iso_id"] == iso]["wilcox_padj"].iloc[0]

    annotate_pval(ax, i - 0.2, i + 0.2, ys[i], 0, ys[i], padj, PAPER_FONTSIZE)

ax.set_xlabel("")
ax.set_ylabel("isoform ratio")
ax.set_title("PKNOX1 isoforms in breast cancer")
ax.get_legend().remove()
ax.set_ylim((-0.05, 1))

fig.savefig("../figures/BRCA_PKNOX1_boxplot.pdf", dpi="figure", bbox_inches="tight")


# ## read in TCGA PanCan data

# In[60]:


pancan = pd.read_table("../data/external/PanCan_PKNOX1_w_header.txt")
pancan.head()


# In[62]:


cols = [
    "event_id",
    "event_type",
    "event_chr",
    "event_coordinates",
    "alt_region_coordinates",
    "gene_name",
]
samp_cols = [x.split(".")[0] for x in pancan.columns if x.startswith("TCGA")]
new_cols = cols + samp_cols
pancan.columns = new_cols
pancan.head()


# In[66]:


for i, row in tcga_paired.iterrows():
    try:
        pancan["paired-diff_%s_ratio" % (i + 1)] = (
            pancan[row.tcga_id_tumor] - pancan[row.tcga_id_ctrl]
        )
        pancan["paired-diff_%s_rationan" % (i + 1)] = pancan[row.tcga_id_tumor].fillna(
            0
        ) - pancan[row.tcga_id_ctrl].fillna(0)
    except KeyError:
        print("missing ctrl %s or tumor %s" % (row.tcga_id_ctrl, row.tcga_id_tumor))


# In[75]:


diff_cols = [x for x in pancan.columns if "paired-diff" in x]
pancan_m = pd.melt(pancan[diff_cols])
pancan_m["type"] = pancan_m["variable"].apply(lambda row: row.split("_")[-1])
pancan_m


# In[79]:


sns.distplot(pancan_m[pancan_m["type"] == "ratio"]["value"])


# In[80]:


pancan_m[pancan_m["type"] == "ratio"]["value"].median()


# In[81]:


pancan_m[pancan_m["type"] == "ratio"]["value"].mean()


# In[ ]:
