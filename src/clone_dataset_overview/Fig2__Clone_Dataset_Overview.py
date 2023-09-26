# coding: utf-8

# # Figure 2: overview of TFIso1.0 clone collection/dataset

# In[245]:


import matplotlib as mpl
import met_brewer
import pandas as pd
import numpy as np
import seaborn as sns
import sys

from matplotlib import pyplot as plt
from scipy import stats

# import utils
sys.path.append("../")

from data_loading import (
    load_valid_isoform_clones,
    load_y2h_isoform_data,
    load_y1h_pdi_data,
    load_m1h_activation_data,
    load_annotated_gencode_tfs,
    load_annotated_TFiso1_collection,
    load_developmental_tissue_expression_remapped,
    load_gtex_remapped,
    load_tf_families,
)

from plotting import (
    PAPER_PRESET,
    PAPER_FONTSIZE,
    nice_boxplot,
    nice_violinplot,
    mimic_r_boxplot,
)


get_ipython().run_line_magic("matplotlib", "inline")
get_ipython().run_line_magic("config", "InlineBackend.figure_format = 'svg'")
mpl.rcParams["figure.autolayout"] = False


# In[2]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[3]:


np.random.seed(2023)


# In[22]:


rename_dev_stage = {
    "8 week post conception,embryo": "08",
    "11 week post conception,late embryo": "11",
    "embryo,7 week post conception": "07",
    "infant": "infant",
    "10 week post conception,late embryo": "10",
    "young adult": "young adult",
    "13 week post conception,late embryo": "13",
    "16 week post conception,late embryo": "16",
    "4 week post conception,embryo": "04",
    "neonate": "neonate",
    "19 week post conception,late embryo": "19",
    "9 week post conception,late embryo": "09",
    "adolescent": "adolescent",
    "5 week post conception,embryo": "05",
    "embryo,6 week post conception": "06",
    "12 week post conception,late embryo": "12",
    "18 week post conception,late embryo": "18",
    "toddler": "toddler",
    "elderly": "elderly",
    "middle adult": "adult",
    "school age child": "child",
}


# ## 1. load clone collection and gencode

# In[4]:


genc_tfs = load_annotated_gencode_tfs()


# In[5]:


clone_tfs = load_annotated_TFiso1_collection()


# ## 2. splicing category figures
#
# ### first: GENCODE

# In[6]:


alt_n = 0
alt_c = 0
alt_int = 0
alt_5ss = 0
alt_3ss = 0
exon_sk = 0
mut_ex = 0
intron_ret = 0
tot = 0

for tf in genc_tfs.keys():
    gene = genc_tfs[tf]
    ref = gene.reference_isoform.name
    alts = [x.name for x in gene.alternative_isoforms]
    for alt in alts:
        splicing_cats = gene.splicing_categories(ref, alt)

        if splicing_cats["alternative N-terminal"]:
            alt_n += 1
        if splicing_cats["alternative C-terminal"]:
            alt_c += 1
        if splicing_cats["alternative internal exon"]:
            alt_int += 1
        if splicing_cats["alternative 5' splice site"]:
            alt_5ss += 1
        if splicing_cats["alternative 3' splice site"]:
            alt_3ss += 1
        if splicing_cats["exon skipping"]:
            exon_sk += 1
        if splicing_cats["mutually exclusive exons"]:
            mut_ex += 1
        if splicing_cats["intron retention"]:
            intron_ret += 1

        tot += 1

genc_df = pd.DataFrame.from_dict(
    {
        "alt. N-terminal": [alt_n],
        "alt. C-terminal": [alt_c],
        "alt. internal exon": [alt_int],
        "alt. 5' splice site": [alt_5ss],
        "alt. 3' splice site": [alt_3ss],
        "exon skipping": [exon_sk],
        "mutually exclusive exons": [mut_ex],
        "intron retention": [intron_ret],
        "total": tot,
    }
)
genc_df.index = ["gencode"]
genc_df


# ### next: TFIso1.0

# In[7]:


alt_n = 0
alt_c = 0
alt_int = 0
alt_5ss = 0
alt_3ss = 0
exon_sk = 0
mut_ex = 0
intron_ret = 0
tot = 0

for tf in clone_tfs.keys():
    gene = clone_tfs[tf]
    ref = gene.cloned_reference_isoform.name
    alts = [x.name for x in gene.cloned_isoforms if x.name != ref]
    for alt in alts:
        splicing_cats = gene.splicing_categories(ref, alt)

        if splicing_cats["alternative N-terminal"]:
            alt_n += 1
        if splicing_cats["alternative C-terminal"]:
            alt_c += 1
        if splicing_cats["alternative internal exon"]:
            alt_int += 1
        if splicing_cats["alternative 5' splice site"]:
            alt_5ss += 1
        if splicing_cats["alternative 3' splice site"]:
            alt_3ss += 1
        if splicing_cats["exon skipping"]:
            exon_sk += 1
        if splicing_cats["mutually exclusive exons"]:
            mut_ex += 1
        if splicing_cats["intron retention"]:
            intron_ret += 1

        tot += 1

clone_df = pd.DataFrame.from_dict(
    {
        "alt. N-terminal": [alt_n],
        "alt. C-terminal": [alt_c],
        "alt. internal exon": [alt_int],
        "alt. 5' splice site": [alt_5ss],
        "alt. 3' splice site": [alt_3ss],
        "exon skipping": [exon_sk],
        "mutually exclusive exons": [mut_ex],
        "intron retention": [intron_ret],
        "total": tot,
    }
)
clone_df.index = ["TFIso1.0"]
clone_df


# ### then: novel isoforms

# In[8]:


alt_n = 0
alt_c = 0
alt_int = 0
alt_5ss = 0
alt_3ss = 0
exon_sk = 0
mut_ex = 0
intron_ret = 0
tot = 0

for tf in clone_tfs.keys():
    gene = clone_tfs[tf]
    ref = gene.cloned_reference_isoform.name
    novels = [x.name for x in gene.cloned_isoforms if x.is_novel_isoform()]
    for novel in novels:
        splicing_cats = gene.splicing_categories(ref, novel)

        if splicing_cats["alternative N-terminal"]:
            alt_n += 1
        if splicing_cats["alternative C-terminal"]:
            alt_c += 1
        if splicing_cats["alternative internal exon"]:
            alt_int += 1
        if splicing_cats["alternative 5' splice site"]:
            alt_5ss += 1
        if splicing_cats["alternative 3' splice site"]:
            alt_3ss += 1
        if splicing_cats["exon skipping"]:
            exon_sk += 1
        if splicing_cats["mutually exclusive exons"]:
            mut_ex += 1
        if splicing_cats["intron retention"]:
            intron_ret += 1

        tot += 1

novel_df = pd.DataFrame.from_dict(
    {
        "alt. N-terminal": [alt_n],
        "alt. C-terminal": [alt_c],
        "alt. internal exon": [alt_int],
        "alt. 5' splice site": [alt_5ss],
        "alt. 3' splice site": [alt_3ss],
        "exon skipping": [exon_sk],
        "mutually exclusive exons": [mut_ex],
        "intron retention": [intron_ret],
        "total": tot,
    }
)
novel_df.index = ["TFIso1.0 - novel"]
novel_df


# ### make plot

# In[9]:


splicing = genc_df.append(clone_df).append(novel_df)
splicing_tot = splicing["total"]
splicing = splicing.drop("total", axis=1)
splicing_perc = splicing.divide(splicing_tot, axis="rows").reset_index()
splicing_perc


# In[10]:


splicing_perc_melt = pd.melt(splicing_perc, id_vars="index")


# In[11]:


colors = met_brewer.met_brew(name="VanGogh2")
sns.palplot(colors)


# In[15]:


fig = plt.figure(figsize=(3.2, 1.75))

ax = sns.barplot(
    data=splicing_perc_melt,
    x="variable",
    y="value",
    hue="index",
    palette={
        "gencode": colors[2],
        "TFIso1.0": colors[4],
        "TFIso1.0 - novel": colors[7],
    },
)
ax.set_xlabel("")
ax.set_xticklabels(
    list(splicing_perc_melt["variable"].unique()), ha="right", va="top", rotation=30
)
ax.set_ylabel("% of alternative isoforms")

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))

fig.savefig("../../figures/fig2/splicing_cats.pdf", dpi="figure", bbox_inches="tight")


# ## 3. plot expression of novel isoforms

# In[18]:


status_map = {}

# only loop through clone collection
for tf in clone_tfs.keys():
    gene = clone_tfs[tf]

    try:
        annot_ref = gene.reference_isoform.name
    except:
        annot_ref = "none"

    try:
        annot_alt = gene.alternative_isoforms
    except:
        annot_alt = []

    for iso in gene.cloned_isoforms:
        if iso.name == annot_ref:
            status_map[iso.clone_acc] = {"gene_name": tf, "status": "ref"}
        elif iso.is_novel_isoform():
            status_map[iso.clone_acc] = {"gene_name": tf, "status": "novel"}
        else:
            status_map[iso.clone_acc] = {"gene_name": tf, "status": "alt"}


# In[19]:


status_map = pd.DataFrame.from_dict(status_map, orient="index")
status_map


# In[23]:


df_dev, metadata_dev, genes_dev = load_developmental_tissue_expression_remapped()

rename_dev_stage = {
    "8 week post conception,embryo": "08",
    "11 week post conception,late embryo": "11",
    "embryo,7 week post conception": "07",
    "infant": "infant",
    "10 week post conception,late embryo": "10",
    "young adult": "young adult",
    "13 week post conception,late embryo": "13",
    "16 week post conception,late embryo": "16",
    "4 week post conception,embryo": "04",
    "neonate": "neonate",
    "19 week post conception,late embryo": "19",
    "9 week post conception,late embryo": "09",
    "adolescent": "adolescent",
    "5 week post conception,embryo": "05",
    "embryo,6 week post conception": "06",
    "12 week post conception,late embryo": "12",
    "18 week post conception,late embryo": "18",
    "toddler": "toddler",
    "elderly": "elderly",
    "middle adult": "adult",
    "school age child": "child",
}

metadata_dev["dev_stage"] = metadata_dev["Developmental_Stage"].map(rename_dev_stage)
means_dev = df_dev.groupby(
    df_dev.columns.map(metadata_dev["organism_part"] + " " + metadata_dev["dev_stage"]),
    axis=1,
).mean()


# In[25]:


df_gtex, metadata_gtex, genes_gtex = load_gtex_remapped()

exclusion_list_gtex = {
    "Cells - Leukemia cell line (CML)",
    "Cells - EBV-transformed lymphocytes",
    "Cells - Cultured fibroblasts",
}

df_gtex = df_gtex.loc[
    :, ~df_gtex.columns.map(metadata_gtex["body_site"]).isin(exclusion_list_gtex)
]
metadata_gtex = metadata_gtex.loc[
    ~metadata_gtex["body_site"].isin(exclusion_list_gtex), :
]

means_gtex = df_gtex.groupby(
    df_gtex.columns.map(metadata_gtex["body_site"]), axis=1
).mean()


# In[27]:


metadata_gtex_dummy = pd.read_table(
    "../../data/processed/metadata_gtex_dummy.csv", sep=",", index_col=0
)


# In[28]:


# use same downsample as fig1
means_gtex_downsample = df_gtex.groupby(
    df_gtex.columns.map(metadata_gtex_dummy["body_site"]), axis=1
).mean()


# In[29]:


means_dev["median"] = means_dev.median(axis=1)
means_dev["max"] = means_dev.max(axis=1)

means_gtex_downsample["median"] = means_gtex_downsample.median(axis=1)
means_gtex_downsample["max"] = means_gtex_downsample.max(axis=1)


# In[30]:


dev_mm = means_dev[["median", "max"]].reset_index()
gtex_ds_mm = means_gtex_downsample[["median", "max"]].reset_index()


# In[31]:


dev_mm["clone_acc"] = dev_mm["UID"].str.split(" ", expand=True)[0]
gtex_ds_mm["clone_acc"] = gtex_ds_mm["UID"].str.split(" ", expand=True)[0]
mm = dev_mm[dev_mm["clone_acc"] != "noclone"].merge(
    gtex_ds_mm[gtex_ds_mm["clone_acc"] != "noclone"],
    on="clone_acc",
    suffixes=("_dev", "_gtex_ds"),
)


# In[32]:


status_map = status_map.reset_index()
status_map["clone_acc"] = status_map["index"].str.split(" ", expand=True)[0]


# In[33]:


exp_nov = status_map.merge(mm, on="clone_acc")
exp_nov_melt = pd.melt(
    exp_nov,
    id_vars=["index", "gene_name", "status", "clone_acc"],
    value_vars=["median_dev", "max_dev", "median_gtex_ds", "max_gtex_ds"],
)
exp_nov_melt["measurement"] = exp_nov_melt["variable"].str.split("_", expand=True)[0]


# In[34]:


colors = met_brewer.met_brew(name="Monet")
sns.palplot(colors)


# In[36]:


fig = plt.figure(figsize=(3, 2))

exp_nov_melt["value_log2"] = np.log2(exp_nov_melt["value"] + 1)
ax = sns.boxplot(
    data=exp_nov_melt[exp_nov_melt["variable"].str.contains("dev")],
    x="status",
    y="value",
    hue="measurement",
    palette={"median": colors[7], "max": colors[6]},
    flierprops={"marker": "o"},
    fliersize=4,
    notch=True,
)

mimic_r_boxplot(ax)

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))

ax.set_xlabel("clone category")
ax.set_ylabel("isoform expression (tpm)")

fig.savefig(
    "../../figures/fig2/novel_isos.dev_expr_boxplot.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# In[278]:


fig = plt.figure(figsize=(1.3, 1.3))

ax = sns.boxplot(
    data=exp_nov_melt[exp_nov_melt["variable"].str.contains("dev")],
    x="status",
    y="value_log2",
    hue="measurement",
    palette={"median": colors[7], "max": colors[6]},
    flierprops={"marker": "o"},
    fliersize=4,
    notch=True,
)

mimic_r_boxplot(ax)

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))

ax.set_xlabel("clone category")
ax.set_ylabel("isoform expression (tpm)")

ticks = [0, 1, 5, 10, 20, 30]
ticklabels = [0, 1, 5, 10, 20, 30]
ax.set_yticks([np.log2(y + 1) for y in ticks])
ax.set_yticklabels(ticklabels)
ax.tick_params(axis="y", labelsize=fontsize - 2)
plt.title("Developmental RNA-seq")

fig.savefig(
    "../../figures/fig2/novel_isos.dev_expr_boxplot.log2.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# In[279]:


fig = plt.figure(figsize=(1.3, 1.3))

ax = sns.boxplot(
    data=exp_nov_melt[exp_nov_melt["variable"].str.contains("gtex_ds")],
    x="status",
    y="value_log2",
    hue="measurement",
    palette={"median": colors[7], "max": colors[6]},
    flierprops={"marker": "o"},
    fliersize=4,
    notch=True,
)

mimic_r_boxplot(ax)

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))

ax.set_xlabel("clone category")
ax.set_ylabel("isoform expression (tpm)")

ticks = [0, 1, 5, 10, 20, 30, 50]
ticklabels = [0, 1, 5, 10, 20, 30, 50]
ax.set_yticks([np.log2(y + 1) for y in ticks])
ax.set_yticklabels(ticklabels)
ax.tick_params(axis="y", labelsize=fontsize - 2)
plt.title("GTEx (down-sampled)")

fig.savefig(
    "../../figures/fig2/novel_isos.gtex_ds_expr_boxplot.log2.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# In[59]:


fig = plt.figure(figsize=(1.7, 1.7))

ax = sns.boxplot(
    data=exp_nov_melt[exp_nov_melt["variable"].str.contains("gtex_ds")],
    x="status",
    y="value",
    hue="measurement",
    palette={"median": colors[7], "max": colors[6]},
    flierprops={"marker": "o"},
    fliersize=4,
    notch=True,
)

mimic_r_boxplot(ax)

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))

ax.set_xlabel("clone category")
ax.set_ylabel("isoform expression (tpm)")


plt.title("GTEx (down-sampled)")

fig.savefig(
    "../../figures/fig2/novel_isos.gtex_ds_expr_boxplot.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# ## 4. plot representation of families across assays

# In[176]:


fam = load_tf_families()


# ### gencode

# In[177]:


len(genc_tfs)


# In[178]:


genc_df = {k: genc_tfs[k].GENCODE_isoforms for k in genc_tfs.keys()}
genc_df = {k: [v.name for v in values] for k, values in genc_df.items()}
genc_df = [(k, v) for k, sublist in genc_df.items() for v in sublist]
genc_df = pd.DataFrame(genc_df, columns=["gene", "isoform"])
genc_df["family"] = genc_df["gene"].map(fam)
genc_df.sample(5)


# In[179]:


leave_separate = [
    "C2H2 ZF",
    "Homeodomain",
    "bHLH",
    "Nuclear receptor",
    "bZIP",
    "Forkhead",
    "Ets",
]


def rename_family(row):
    if row.family in leave_separate:
        return row.family
    else:
        return "Other"


genc_df["family_renamed"] = genc_df.apply(rename_family, axis=1)
genc_df.sample(5)


# In[180]:


genc_vc = genc_df.groupby("family_renamed")["isoform"].agg("count").reset_index()
genc_vc["source"] = "GENCODE"
genc_vc


# ### TF Iso 1.0

# In[181]:


len(clone_tfs)


# In[219]:


clone_df = {k: clone_tfs[k].cloned_isoforms for k in clone_tfs.keys()}
clone_df = {k: [v.clone_acc for v in values] for k, values in clone_df.items()}
clone_df = [(k, v) for k, sublist in clone_df.items() for v in sublist]
clone_df = pd.DataFrame(clone_df, columns=["gene", "isoform"])
clone_df["family"] = clone_df["gene"].map(fam)
clone_df.sample(5)


# In[220]:


len(clone_df)


# In[221]:


len(clone_df.gene.unique())


# In[183]:


def rename_family(row):
    if row.family in leave_separate:
        return row.family
    else:
        return "Other"


clone_df["family_renamed"] = clone_df.apply(rename_family, axis=1)


# In[184]:


clone_vc = clone_df.groupby("family_renamed")["isoform"].agg("count").reset_index()
clone_vc["source"] = "TFIso1.0"


# ### Y1H

# In[185]:


y1h = load_y1h_pdi_data()
y1h["family"] = y1h["gene_symbol"].map(fam)
y1h["family_renamed"] = y1h.apply(rename_family, axis=1)
y1h.sample(5)


# In[207]:


len(y1h)


# In[208]:


len(y1h.tf.unique())


# In[186]:


baits = [
    x
    for x in y1h.columns
    if x not in ["gene_symbol", "clone_acc", "family", "family_renamed"]
]
y1h["any_true"] = y1h[baits].sum(axis=1)
y1h


# In[187]:


y1h_vc = y1h.groupby("family_renamed")["clone_acc"].agg("count").reset_index()
y1h_vc.columns = ["family_renamed", "isoform"]
y1h_vc["source"] = "Y1H (all)"


# In[188]:


y1h_any_vc = (
    y1h[y1h["any_true"] > 0]
    .groupby("family_renamed")["clone_acc"]
    .agg("count")
    .reset_index()
)
y1h_any_vc.columns = ["family_renamed", "isoform"]
y1h_any_vc["source"] = "Y1H (≥1 PDI)"


# ### Y2H

# In[189]:


y2h = load_y2h_isoform_data()
y2h["family"] = y2h["ad_gene_symbol"].map(fam)
y2h["family_renamed"] = y2h.apply(rename_family, axis=1)
y2h.sample(5)


# In[209]:


len(y2h.ad_clone_acc.unique())


# In[210]:


len(y2h.ad_gene_symbol.unique())


# In[190]:


y2h_vc = y2h.groupby("family_renamed")["ad_clone_acc"].agg("count").reset_index()
y2h_vc.columns = ["family_renamed", "isoform"]
y2h_vc["source"] = "Y2H (all)"


# In[191]:


y2h_any_vc = (
    y2h[y2h["Y2H_result"] == True]
    .groupby("family_renamed")["ad_clone_acc"]
    .agg("count")
    .reset_index()
)
y2h_any_vc.columns = ["family_renamed", "isoform"]
y2h_any_vc["source"] = "Y2H (≥1 PPI)"


# ### activation

# In[192]:


m1h = load_m1h_activation_data()
m1h["M1H_mean"] = m1h[["M1H_rep1", "M1H_rep2", "M1H_rep3"]].mean(axis=1)
m1h["family"] = m1h["gene"].map(fam)
m1h["family_renamed"] = m1h.apply(rename_family, axis=1)
m1h.sample(5)


# In[211]:


len(m1h.clone_acc.unique())


# In[212]:


len(m1h.gene.unique())


# In[193]:


m1h_vc = m1h.groupby("family_renamed")["clone_acc"].agg("count").reset_index()
m1h_vc.columns = ["family_renamed", "isoform"]
m1h_vc["source"] = "M1H (all)"


# In[198]:


m1h_any_vc = (
    m1h[m1h["M1H_mean"].abs() > 1]
    .groupby("family_renamed")["clone_acc"]
    .agg("count")
    .reset_index()
)
m1h_any_vc.columns = ["family_renamed", "isoform"]
m1h_any_vc["source"] = "M1H (≥2-fold activ.)"


# ### merge for plot

# In[199]:


mrg_vc = genc_vc.append(clone_vc)
mrg_vc = mrg_vc.append(y1h_vc).append(y1h_any_vc)
mrg_vc = mrg_vc.append(y2h_vc).append(y2h_any_vc)
mrg_vc = mrg_vc.append(m1h_vc).append(m1h_any_vc)
mrg_vc


# In[293]:


mrg_piv = pd.pivot_table(
    mrg_vc, values="isoform", columns="source", index="family_renamed"
)
mrg_piv = mrg_piv.fillna(0)
mrg_piv = (mrg_piv / mrg_piv.sum(axis=0)) * 100
mrg_piv = mrg_piv.T
mrg_piv = mrg_piv.reindex(
    [
        "GENCODE",
        "TFIso1.0",
        "Y1H (all)",
        "Y1H (≥1 PDI)",
        "Y2H (all)",
        "Y2H (≥1 PPI)",
        "M1H (all)",
        "M1H (≥2-fold activ.)",
    ]
)
mrg_piv = mrg_piv.reset_index()

mrg_piv = mrg_piv[
    [
        "source",
        "Other",
        "Ets",
        "Forkhead",
        "bZIP",
        "Nuclear receptor",
        "bHLH",
        "Homeodomain",
        "C2H2 ZF",
    ]
]
mrg_piv


# In[294]:


colors = met_brewer.met_brew(name="Hokusai1")
colors.append("lightgrey")
colors = colors[::-1]
# colors[7] = "lightgrey"
sns.palplot(colors)


# In[302]:


ax = mrg_piv.plot.bar(x="source", stacked=True, color=colors, figsize=(2, 2))

ax.set_ylabel("% of isoforms")
ax.set_xlabel("")

plt.legend()
handles, labels = ax.get_legend_handles_labels()
ax.legend(reversed(handles), reversed(labels), loc=2, bbox_to_anchor=(1.01, 1))

plt.savefig("../../figures/fig2/assay_families.detailed.pdf", bbox_inches="tight")


# In[311]:


ax = mrg_piv[
    mrg_piv["source"].isin(
        ["GENCODE", "TFIso1.0", "Y1H (all)", "Y2H (all)", "M1H (all)"]
    )
].plot.bar(x="source", stacked=True, color=colors, figsize=(1.1, 1.5))

ax.set_ylabel("% of isoforms")
ax.set_xlabel("")

plt.legend()
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    reversed(handles),
    reversed(labels),
    loc=2,
    bbox_to_anchor=(1.01, 1),
    borderpad=0.25,
    handlelength=1,
    handletextpad=0.2,
)

plt.savefig("../../figures/fig2/assay_families.pdf", bbox_inches="tight")


# ## 5. plot biochemical activities of TFs across ref/alt/novel

# In[230]:


mane_select_clones = {
    tf.MANE_select_isoform.clone_acc
    for tf in clone_tfs.values()
    if tf.cloned_MANE_select_isoform
}


# In[232]:


iso = load_valid_isoform_clones()
iso["is_longest_isoform"] = iso["clone_acc"].isin(
    iso.sort_values("num_aa", ascending=False)
    .groupby("gene")
    .nth(0)["clone_acc"]
    .values
)
iso["category"] = "alternative"
iso.loc[iso["clone_acc"].isin(mane_select_clones), "category"] = "reference"
iso.loc[iso["is_novel_isoform"], "category"] = "novel"


# In[236]:


len(iso["gene"].unique())


# In[235]:


genes_w_ref = list(iso[iso["category"] == "reference"]["gene"].unique())
len(genes_w_ref)


# In[237]:


# subset iso df to only genes w MANE select isoform
iso_sub = iso[iso["gene"].isin(genes_w_ref)]
len(iso_sub)


# In[241]:


iso_sub["valid_ppi_test"] = iso["clone_acc"].map(
    y2h.groupby("ad_clone_acc").apply(
        lambda rows: (
            (rows["Y2H_result"] == True) | (rows["Y2H_result"] == False)
        ).any()
    )
)
iso_sub["at_least_one_ppi"] = iso["clone_acc"].map(
    y2h.groupby("ad_clone_acc").apply(lambda rows: ((rows["Y2H_result"] == True)).any())
)


y1h = y1h.drop_duplicates("clone_acc")
iso_sub["at_least_one_pdi"] = iso_sub["clone_acc"].map(
    y1h.drop(columns=["gene_symbol"]).set_index("clone_acc").sum(axis=1) > 0
)

iso_sub["at_least_two_fold_activation"] = iso_sub["clone_acc"].map(
    m1h.drop(columns=["gene"]).set_index("clone_acc").mean(axis=1).abs() > 1
)


# In[242]:


iso_sub.category.value_counts()


# In[313]:


colors = met_brewer.met_brew(name="Monet")
sns.palplot(colors)


# In[318]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2.1, h=1.5)
cats = ["reference", "alternative", "novel"]
positives = []
tested = []
for cat in cats:
    positives.append(
        iso_sub.loc[
            iso_sub["valid_ppi_test"] & (iso_sub["category"] == cat), "at_least_one_ppi"
        ].sum()
    )
    tested.append(
        iso_sub.loc[
            iso_sub["valid_ppi_test"] & (iso_sub["category"] == cat), "at_least_one_ppi"
        ]
        .notnull()
        .sum()
    )
for cat in cats:
    positives.append(
        iso_sub.loc[(iso_sub["category"] == cat), "at_least_one_pdi"].sum()
    )
    tested.append(
        iso_sub.loc[(iso_sub["category"] == cat), "at_least_one_pdi"].notnull().sum()
    )
for cat in cats:
    positives.append(
        iso_sub.loc[(iso_sub["category"] == cat), "at_least_two_fold_activation"].sum()
    )
    tested.append(
        iso_sub.loc[(iso_sub["category"] == cat), "at_least_two_fold_activation"]
        .notnull()
        .sum()
    )
for cat in cats:
    tested_iso = (
        iso_sub["valid_ppi_test"]
        & iso_sub["at_least_two_fold_activation"].notnull()
        & iso_sub["at_least_one_pdi"].notnull()
        & (iso_sub["category"] == cat)
    )
    positives.append(
        (
            iso_sub.loc[tested_iso, "at_least_one_ppi"]
            | iso_sub.loc[tested_iso, "at_least_two_fold_activation"]
            | iso_sub.loc[tested_iso, "at_least_one_pdi"]
        ).sum()
    )
    tested.append(tested_iso.sum())

vals = [p / n for p, n in zip(positives, tested)]
# errs = [np.sqrt(((p / n) * (1 - (p / n)) / n)) for p, n in zip(positives, tested)]

pos = np.array(positives)
neg = np.array(tested) - pos
fracs = np.array(vals)
intv = stats.beta.interval(0.6827, pos + 1, neg + 1)
errs = [fracs - intv[0], intv[1] - fracs]
errs[0][pos == 0] = 0.0
errs[1][neg == 0] = 0.0

offset = 0.5
x_pos = (
    [i for i in range(3)]
    + [i + offset for i in range(3, 6)]
    + [i + offset * 2 for i in range(6, 9)]
    + [i + offset * 3 for i in range(9, 12)]
)


ax.bar(x=x_pos, height=vals, color=[colors[0], colors[1], colors[2]] * 3)
ax.errorbar(
    x=x_pos, y=vals, yerr=errs, color="black", fmt="none", linewidth=1, capsize=1
)
ax.set_ylim(0, 1.1)
ax.set_yticks(np.linspace(0, 1, 11))
ax.set_yticks(np.linspace(0, 1, 21), minor=True)
ax.set_yticklabels([f"{y:.0%}" for y in ax.get_yticks()])
for loc in ["top", "bottom", "right"]:
    ax.spines[loc].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_xticks([0, 1, 2])

ax.set_xticklabels(
    ["Reference", "Alternative", "Novel"],
    rotation=45,
    ha="right",
    va="top",
    fontweight="bold",
)
[t.set_color(colors[i]) for i, t in enumerate(ax.xaxis.get_ticklabels())]

ax.set_ylabel("Proportion of isoforms")
ax.text(y=1.025, x=x_pos[1], s="≥ 1 PPI", fontsize=6, va="bottom", ha="center")
ax.text(y=1.025, x=x_pos[4], s="≥ 1 PDI", fontsize=6, va="bottom", ha="center")
ax.text(
    y=1.025,
    x=x_pos[7],
    s="≥ 2-fold\nactivation/\nrepression",
    fontsize=6,
    va="bottom",
    ha="center",
)
ax.text(
    y=1.025,
    x=x_pos[10],
    s="Any one\nof three",
    fontsize=6,
    fontstyle="italic",
    va="bottom",
    ha="center",
)
fig.savefig(
    "../../figures/fig2/at-least-some-assay-result_ref-vs-alt-vs-novel_bar.pdf",
    bbox_inches="tight",
)


# In[263]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2.5, h=1.5)
cats = ["reference", "alternative", "novel"]
positives = []
tested = []
for cat in cats:
    positives.append(
        iso_sub.loc[
            iso_sub["valid_ppi_test"] & (iso_sub["category"] == cat), "at_least_one_ppi"
        ].sum()
    )

for cat in cats:
    positives.append(
        iso_sub.loc[(iso_sub["category"] == cat), "at_least_one_pdi"].sum()
    )

for cat in cats:
    positives.append(
        iso_sub.loc[(iso_sub["category"] == cat), "at_least_two_fold_activation"].sum()
    )

for cat in cats:
    positives.append(
        (
            iso_sub.loc[(iso_sub["category"] == cat), "at_least_one_ppi"].fillna(False)
            | iso_sub.loc[
                (iso_sub["category"] == cat), "at_least_two_fold_activation"
            ].fillna(False)
            | iso_sub.loc[(iso_sub["category"] == cat), "at_least_one_pdi"].fillna(
                False
            )
        ).sum()
    )

tested = [(iso_sub["category"] == cat).sum() for cat in cats] * 4
vals = [p / n for p, n in zip(positives, tested)]

pos = np.array(positives)
neg = np.array(tested) - pos
fracs = np.array(vals)
intv = stats.beta.interval(0.6827, pos + 1, neg + 1)
errs = [fracs - intv[0], intv[1] - fracs]
errs[0][pos == 0] = 0.0
errs[1][neg == 0] = 0.0

offset = 0.5
x_pos = (
    [i for i in range(3)]
    + [i + offset for i in range(3, 6)]
    + [i + offset * 2 for i in range(6, 9)]
    + [i + offset * 3 for i in range(9, 12)]
)
ax.bar(x=x_pos, height=vals, color=[colors[0], colors[1], colors[2]] * 3)
ax.errorbar(
    x=x_pos, y=vals, yerr=errs, color="black", fmt="none", capsize=2, linewidth=1
)
ax.set_ylim(0, 1)
ax.set_yticks(np.linspace(0, 1, 11))
ax.set_yticks(np.linspace(0, 1, 21), minor=True)
ax.set_yticklabels([f"{y:.0%}" for y in ax.get_yticks()])
for loc in ["top", "bottom", "right"]:
    ax.spines[loc].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_xticks([0, 1, 2])

ax.set_xticklabels(
    ["Reference", "Alternative", "Novel"],
    rotation=45,
    ha="right",
    va="top",
    fontweight="bold",
)
[t.set_color(colors[i]) for i, t in enumerate(ax.xaxis.get_ticklabels())]

ax.set_ylabel("Proportion of isoforms")
ax.text(y=1, x=x_pos[1], s="≥ 1 PPI", fontsize=7, va="bottom", ha="center")
ax.text(y=1, x=x_pos[4], s="≥ 1 PDI", fontsize=7, va="bottom", ha="center")
ax.text(
    y=1,
    x=x_pos[7],
    s="≥ 2-fold\nactivation/\nrepression",
    fontsize=7,
    va="bottom",
    ha="center",
)
ax.text(y=1, x=x_pos[10], s="Any one\nof three", fontsize=7, va="bottom", ha="center")
fig.savefig(
    "../../figures/fig2/at-least-some-assay-result_ref-vs-alt-vs-novel_absolute_bar.pdf",
    bbox_inches="tight",
)


# ## 6. network ball

# In[273]:


# table of edges
#    - clone to (edge + clone_id) + to duplicate
# table of nodes
#    - clone to gene
#    - dna vs isoform vs

ppi = load_isoform_and_paralog_y2h_data()
ppi = ppi.loc[
    (ppi["category"] == "tf_isoform_ppis") & (ppi["Y2H_result"] == True),
    ["ad_clone_acc", "ad_gene_symbol", "db_gene_symbol"],
]
ppi = ppi.rename(columns={"ad_clone_acc": "isoform", "db_gene_symbol": "partner"})
ppi["partner"] = ppi["partner"] + "-" + ppi["ad_gene_symbol"]
pdi = pd.read_csv("../../data/internal/a2_juan_pdi_w_unique_isoacc.tsv", sep="\t")
clones = load_valid_isoform_clones()
pdi = pdi.loc[pdi["clone_acc"].isin(clones["clone_acc"]), :]
pdi["partner"] = pdi["bait"] + "-" + pdi["gene_symbol"]
pdi["isoform"] = pdi["clone_acc"]
edges = pd.concat(
    [ppi.loc[:, ["isoform", "partner"]], pdi.loc[:, ["isoform", "partner"]]]
)
edges.to_csv("../../output/edges.tsv", sep="\t", index=False)

clones = clones.rename(columns={"clone_acc": "node_id"})
clones["type"] = "isoform"
dna = pd.DataFrame(data=pdi["partner"].unique(), columns=["node_id"])
dna["type"] = "DNA"
proteins = pd.DataFrame(data=ppi["partner"].unique(), columns=["node_id"])
proteins["type"] = "Protein"
nodes = pd.concat([clones, proteins, dna], sort=True)
nodes.to_csv("../../output/node_table.tsv", sep="\t", index=False)
