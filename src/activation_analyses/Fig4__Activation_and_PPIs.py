#!/usr/bin/env python
# coding: utf-8

# # Fig 4: Activation and PPIs/Dimerization

# In[1]:


import numpy as np
import matplotlib as mpl
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import seaborn as sns
import sys

from statannotations.Annotator import Annotator

# import utils
sys.path.append("../")
sys.path.append("../data_loading")

import plotting
from plotting import (
    mimic_r_boxplot,
    violinplot_reflected,
    y2h_ppi_per_tf_gene_plot,
    y1h_pdi_per_tf_gene_plot,
    m1h_activation_per_tf_gene_plot,
)

from data_loading import (
    load_y2h_isoform_data,
    load_m1h_activation_data,
    load_ppi_partner_categories,
    load_annotated_TFiso1_collection,
    load_human_tf_db,
    load_y1h_pdi_data,
    load_tf_families,
    DIMERIZING_TF_FAMILIES,
    load_full_y2h_data_including_controls,
    load_valid_isoform_clones,
)

from isoform_pairwise_metrics import load_ref_vs_alt_isoforms_table, _add_PPI_columns


# In[2]:


PAPER_PRESET = {
    "style": "ticks",
    "font": "Helvetica",
    "context": "paper",
    "rc": {
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "axes.linewidth": 0.5,
        "legend.fontsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "axes.edgecolor": "black",
        "xtick.major.pad": 3.0,
        "ytick.major.pad": 3.0,
    },
}
PAPER_FONTSIZE = 7


# In[3]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[4]:


np.random.seed(2023)


# ## 1. load TFs, assay data (Y2H, M1H), Pfam domains, and PPI partner categories/cofactors

# In[5]:


tfs = load_annotated_TFiso1_collection()
pairs = load_ref_vs_alt_isoforms_table()

# RORC-1 alt iso is causing an error - filter out here - there's no data for it?
pairs = pairs[pairs["clone_acc_alt"] != "RORC|1/6|05F11"]

pairs["ref_iso"] = pairs["clone_acc_ref"].apply(
    lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
)
pairs["alt_iso"] = pairs["clone_acc_alt"].apply(
    lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
)
pairs["f_disorder_difference"] = pairs.apply(
    lambda x: tfs[x["gene_symbol"]].disordered_fraction_of_different_regions(
        x["ref_iso"], x["alt_iso"]
    ),
    axis=1,
)


# In[6]:


y2h = load_y2h_isoform_data()
m1h = load_m1h_activation_data()
m1h["mean"] = m1h[["M1H_rep1", "M1H_rep2", "M1H_rep3"]].mean(axis=1)
cats = load_ppi_partner_categories()


# In[7]:


df = pd.read_excel(
    "../../data/external/Geiger-et-al_MCP_2012_Supplementary-Table-2.xlsx", skiprows=1
)
hek_avrg = df[["iBAQ HEK293_1", "iBAQ HEK293_2", "iBAQ HEK293_3"]].mean(axis=1)
print((hek_avrg > 0).sum(), "proteins expressed in HEK293 proteome")
hek_expressed_genes = set(
    df.loc[(hek_avrg > 0) & df["Gene Names"].notnull(), "Gene Names"]
    .str.split(";")
    .explode()
    .values
)
all_partners = set(y2h["db_gene_symbol"].unique())
print(
    "of {} PPI partners, {} are expressed in HEK293 cells".format(
        len(all_partners), len(all_partners.intersection(hek_expressed_genes))
    )
)


# In[8]:


# now add Pfam AD/RDs
pfam = pd.read_csv(
    "../../data/external/Pfam-A.clans.tsv",
    sep="\t",
    names=["pfam_accession", "clan", "clan_name", "short_name", "name"],
)


# In[9]:


cof = pd.read_csv(
    "../../data/external/AnimalTFDB3_Homo_sapiens_TF_cofactors.txt", sep="\t"
)
if cof["Symbol"].duplicated().any():
    raise UserWarning("unexpected duplicates")


# ## 2. categorize effector domain changes between alt/ref iso

# In[10]:


dom = pd.concat(
    [g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()]
)

# add activation or repression annotation from effector domain studies
effector_domain_type = {}
for tf in tfs.values():
    for d in tf.cloned_reference_isoform.aa_seq_features:
        if d.category == "effector_domain":
            effector_domain_type[d.accession] = d.name
dom["type"] = dom["accession"].map(effector_domain_type)


# In[11]:


# add activation or repression annotation from Pfam domains directly
pfam_ad = pfam[
    (pfam["name"].str.contains("transcription activation"))
    | (pfam["name"].str.contains("transactivation"))
    | (pfam["short_name"].str.contains("TAD"))
].copy()
pfam_ad["type"] = "AD"

# RD
pfam_rd = pfam[
    (pfam["short_name"].str.contains("NRIP1_repr"))
    | (pfam["name"].str.contains("KRAB"))
].copy()
pfam_rd["type"] = "RD"
pfam_effs = pd.concat([pfam_ad, pfam_rd])


def get_pfam_type(row):
    if not pd.isnull(row["type"]):
        return row["type"]
    else:
        pfam_sub = pfam_effs[pfam_effs["pfam_accession"] == row["accession"]]
        if len(pfam_sub) > 0:
            return pfam_sub["type"].iloc[0]
        else:
            return np.nan


dom["type_incl_pfam"] = dom.apply(get_pfam_type, axis=1)


# In[12]:


# considering Pfam and Effector domains
def fraction_of_effector_domains_removed(row, effector_type):
    ds = dom.loc[
        (dom["alt_iso"] == row["alt_iso"]) & (dom["type_incl_pfam"] == effector_type), :
    ]
    if ds.shape[0] == 0:
        return np.nan
    return ds[["deletion", "frameshift"]].sum().sum() / ds["length"].sum()


def insertion_in_effector_domains(row, effector_type):
    ds = dom.loc[
        (dom["alt_iso"] == row["alt_iso"]) & (dom["type_incl_pfam"] == effector_type), :
    ]
    if ds.shape[0] == 0:
        return np.nan
    return ds["insertion"].sum()


def domain_length(row, effector_type):
    ds = dom.loc[
        (dom["alt_iso"] == row["alt_iso"]) & (dom["type_incl_pfam"] == effector_type), :
    ]
    if ds.shape[0] == 0:
        return np.nan
    return ds["length"].sum()


for effector_type in ["AD", "RD", "Bif"]:
    pairs[
        "fraction_of_{}_domains_removed_incl_pfam".format(effector_type)
    ] = pairs.apply(
        fraction_of_effector_domains_removed, effector_type=effector_type, axis=1
    )
    pairs["insertion_in_{}_domains_incl_pfam".format(effector_type)] = pairs.apply(
        insertion_in_effector_domains, effector_type=effector_type, axis=1
    )
    pairs["length_of_{}_domains_incl_pfam".format(effector_type)] = pairs.apply(
        domain_length, effector_type=effector_type, axis=1
    )


# ## 3. plot number of annotated domains (of various types) across isoforms

# In[13]:


# plot of number of activation domains per ref iso
# fraction of sequnce within effector domains
def count_effector_domains(gene):
    iso = gene.cloned_reference_isoform
    c = 0
    for d in iso.aa_seq_features:
        if d.category == "effector_domain":
            c += 1
    return c


n_effector = [count_effector_domains(tf) for tf in tfs.values()]

fig, ax = plt.subplots(1, 1, figsize=(2, 1))
ax.hist(
    n_effector, range=(-0.25, max(n_effector) + 0.25), bins=(max(n_effector) * 2 + 1)
)
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel("Number of genes ({} total)".format(len(tfs)))
ax.set_xlabel("Effector domains in reference isoform")


# In[14]:


# plot of number of activation domains per ref iso
# fraction of sequnce within effector domains
def count_Soto_effector_domains(gene):
    iso = gene.cloned_reference_isoform
    c = 0
    for d in iso.aa_seq_features:
        if d.category == "effector_domain" and d.accession.startswith("Soto"):
            c += 1
    return c


n_effector_soto = [count_Soto_effector_domains(tf) for tf in tfs.values()]

fig, ax = plt.subplots(1, 1, figsize=(2, 1))
ax.hist(
    n_effector_soto,
    range=(-0.25, max(n_effector) + 0.25),
    bins=(max(n_effector) * 2 + 1),
)
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel("Number of genes ({} total)".format(len(tfs)))
ax.set_xlabel("Effector domains in reference isoform")
ax.set_title("Soto et al. data")


# In[15]:


# plot of number of activation domains per ref iso
# fraction of sequnce within effector domains
def count_Bintu_effector_domains(gene):
    iso = gene.cloned_reference_isoform
    c = 0
    for d in iso.aa_seq_features:
        if d.category == "effector_domain" and not d.accession.startswith("Soto"):
            c += 1
    return c


n_effector_bintu = [count_Bintu_effector_domains(tf) for tf in tfs.values()]

fig, ax = plt.subplots(1, 1, figsize=(2, 1))
ax.hist(
    n_effector_bintu,
    range=(-0.25, max(n_effector) + 0.25),
    bins=(max(n_effector) * 2 + 1),
)
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel("Number of genes ({} total)".format(len(tfs)))
ax.set_xlabel("Effector domains in reference isoform")
ax.set_title("Data from Bintu lab papers")


# In[16]:


counter = {
    "Soto": {"AD": 0, "RD": 0, "Bif": 0},
    "Tycko": {"AD": 0, "RD": 0},
    "DelRosso": {"AD": 0, "RD": 0},
}
for tf in tfs.values():
    has_effector = False
    for domain in tf.cloned_reference_isoform.aa_seq_features:
        if domain.category == "effector_domain":
            counter[domain.accession.split("_")[0]][domain.name] += 1
counter


# In[17]:


fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.bar(
    x=range(6),
    height=[counter[x][y] for x in counter.keys() for y in ["AD", "RD"]],
    color=[sns.color_palette("Set2")[0], sns.color_palette("Set2")[3]] * 3,
)
ax.set_xticks([0.5, 2.5, 4.5])
ax.set_xticklabels(["Soto et al.", "Tycko et al.", "DelRosso et al."])
ax.set_ylabel("Number of effector domains")


# ## 4. summary plot looking at gain/loss of domains and activity

# In[18]:


m1h["gte_2_fold"] = m1h["mean"].abs() >= 1
pairs["m1h_gte_2_fold_at_least_one_iso_per_gene"] = pairs["gene_symbol"].map(
    m1h.groupby("gene_symbol")["gte_2_fold"].any()
)
pairs["abs_activation_fold_change_log2"] = pairs["activation_fold_change_log2"].abs()


# In[19]:


# create a color map of domain length
# sum up lengths of all domains (plot only includes examples w 1 type of domain)
pairs["tot_dom_length_incl_pfam"] = pairs[
    [
        "length_of_AD_domains_incl_pfam",
        "length_of_RD_domains_incl_pfam",
        "length_of_Bif_domains_incl_pfam",
    ]
].sum(axis=1)
t_dom_length_incl_pfam = pairs.loc[:, "tot_dom_length_incl_pfam"].values
t_dom_length_incl_pfam = t_dom_length_incl_pfam[t_dom_length_incl_pfam > 0]

# using min and max makes colors too hard too read - cut off
cmap = sns.color_palette("flare", as_cmap=True)
norm = plt.Normalize(25, 250)
palette_dom_length_incl_pfam = {
    value: cmap(norm(value)) for value in t_dom_length_incl_pfam
}


def re_color(row, palette):
    if row["tot_dom_length_incl_pfam"] == 0:
        color = sns.color_palette("flare")[0]
    else:
        color = palette[row["tot_dom_length_incl_pfam"]]
    return color


pairs["color_dom_length_incl_pfam"] = pairs.apply(
    re_color, axis=1, palette=palette_dom_length_incl_pfam
)


# In[127]:


df = pairs.copy()
df = df.loc[
    df["activation_fold_change_log2"].notnull()
    & df["m1h_gte_2_fold_at_least_one_iso_per_gene"],
    :,
]
palette = palette_dom_length_incl_pfam
hue = "tot_dom_length_incl_pfam"
color = "color_dom_length_incl_pfam"
t = t_dom_length_incl_pfam

gs_kw = dict(width_ratios=[0.4, 0.4, 0.5, 1.1, 1.1, 1.75, 0.85])

fig, axs = plt.subplots(1, 7, sharey=True, gridspec_kw=gs_kw)
fig.set_size_inches(w=8.2, h=2)

point_size = 6


tot_loss_activ = df.loc[
    (df["fraction_of_AD_domains_removed_incl_pfam"] == 1)
    & (
        df["fraction_of_RD_domains_removed_incl_pfam"].isnull()
        | (df["fraction_of_RD_domains_removed_incl_pfam"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed_incl_pfam"].isnull()
        | (df["fraction_of_Bif_domains_removed_incl_pfam"] == 0)
    ),
    :,
]
axs[0].set_title("activation\ndomain")
sns.swarmplot(
    data=tot_loss_activ,
    y="activation_fold_change_log2",
    x="fraction_of_AD_domains_removed_incl_pfam",
    size=point_size,
    clip_on=False,
    ax=axs[0],
    palette=palette,
    hue=hue,
    linewidth=1,
    edgecolor="black",
    alpha=1,
)
axs[0].set_xticks([])
axs[0].set_xlabel("")
axs[0].get_legend().remove()

tbx5_y = df.loc[
    (df["clone_acc_alt"] == "TBX5|3/3|08H01"), "activation_fold_change_log2"
].values[0]
for point in axs[0].collections:
    for x, y in point.get_offsets():
        if np.isclose(tbx5_y, y):
            print("found: %s, %s" % (x, y))
            axs[0].annotate(
                "TBX5-3",
                xy=(x, y),
                xytext=(0, -10),
                textcoords="offset points",
                arrowprops=dict(
                    arrowstyle="-", connectionstyle="arc3,rad=0.3", color="black"
                ),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )

tot_loss_repr = df.loc[
    (df["fraction_of_RD_domains_removed_incl_pfam"] == 1)
    & (
        df["fraction_of_AD_domains_removed_incl_pfam"].isnull()
        | (df["fraction_of_AD_domains_removed_incl_pfam"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed_incl_pfam"].isnull()
        | (df["fraction_of_Bif_domains_removed_incl_pfam"] == 0)
    ),
    :,
]
axs[1].set_title("repression\ndomain")
sns.swarmplot(
    data=tot_loss_repr,
    y="activation_fold_change_log2",
    x="fraction_of_RD_domains_removed_incl_pfam",
    size=point_size,
    clip_on=False,
    ax=axs[1],
    palette=palette,
    hue=hue,
    linewidth=1,
    edgecolor="black",
    alpha=1,
)
axs[1].set_xticks([])
axs[1].set_xlabel("")
axs[1].get_legend().remove()

dlx1_y = df.loc[
    (df["clone_acc_alt"] == "DLX1|2/2|07E09"), "activation_fold_change_log2"
].values[0]
for point in axs[1].collections:
    for x, y in point.get_offsets():
        if np.isclose(dlx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[1].annotate(
                "DLX1-2",
                xy=(x, y),
                xytext=(0, 10),
                textcoords="offset points",
                arrowprops=dict(
                    arrowstyle="-", connectionstyle="arc3,rad=0.3", color="black"
                ),
                ha="center",
                va="bottom",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )


tot_loss_both = df.loc[
    (df["fraction_of_AD_domains_removed_incl_pfam"] == 1)
    & (df["fraction_of_RD_domains_removed_incl_pfam"] == 1),
    :,
]
axs[2].set_title("both activ. &\nrepr. domains")
sns.swarmplot(
    data=tot_loss_both,
    y="activation_fold_change_log2",
    x="fraction_of_RD_domains_removed_incl_pfam",
    size=point_size,
    clip_on=False,
    ax=axs[2],
    palette=palette,
    hue=hue,
    linewidth=1,
    edgecolor="black",
    alpha=1,
)
axs[2].set_xticks([])
axs[2].set_xlabel("")
axs[2].get_legend().remove()


# now partial loss
axs[3].set_title("activation\ndomain")
partial_loss_activ = df.loc[
    (
        df["m1h_gte_2_fold_at_least_one_iso_per_gene"]
        & (df["fraction_of_AD_domains_removed_incl_pfam"] > 0)
        & (df["fraction_of_AD_domains_removed_incl_pfam"] < 1)
        & (
            df["fraction_of_RD_domains_removed_incl_pfam"].isnull()
            | (df["fraction_of_RD_domains_removed_incl_pfam"] == 0)
        )
        & (
            df["fraction_of_Bif_domains_removed_incl_pfam"].isnull()
            | (df["fraction_of_Bif_domains_removed_incl_pfam"] == 0)
        )
    ),
    :,
]
axs[3].scatter(
    partial_loss_activ.loc[:, "fraction_of_AD_domains_removed_incl_pfam"].values,
    partial_loss_activ.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_activ.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[3].set_xlabel("")
axs[3].set_xlim(1, 0)
axs[3].set_xticks([0.99, 0.5, 0.01])
axs[3].set_xticklabels([f"{x:.0%}" for x in axs[3].get_xticks()])


axs[4].set_title("repression\ndomain")
partial_loss_repr = df.loc[
    (
        df["m1h_gte_2_fold_at_least_one_iso_per_gene"]
        & (df["fraction_of_RD_domains_removed_incl_pfam"] > 0)
        & (df["fraction_of_RD_domains_removed_incl_pfam"] < 1)
        & (
            df["fraction_of_AD_domains_removed_incl_pfam"].isnull()
            | (df["fraction_of_AD_domains_removed_incl_pfam"] == 0)
        )
        & (
            df["fraction_of_Bif_domains_removed_incl_pfam"].isnull()
            | (df["fraction_of_Bif_domains_removed_incl_pfam"] == 0)
        )
    ),
    :,
]

axs[4].scatter(
    partial_loss_repr.loc[:, "fraction_of_RD_domains_removed_incl_pfam"].values,
    partial_loss_repr.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_repr.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[4].set_xlabel("")
axs[4].set_xlim(1, 0)
axs[4].set_xticks([0.99, 0.5, 0.01])
axs[4].set_xticklabels([f"{x:.0%}" for x in axs[4].get_xticks()])


all_retained = df.loc[
    (
        (df["fraction_of_AD_domains_removed_incl_pfam"] == 0)
        | (df["fraction_of_RD_domains_removed_incl_pfam"] == 0)
        | (df["fraction_of_Bif_domains_removed_incl_pfam"] == 0)
    )
    & (
        df["fraction_of_AD_domains_removed_incl_pfam"].isnull()
        | (df["fraction_of_AD_domains_removed_incl_pfam"] == 0)
    )
    & (
        df["fraction_of_RD_domains_removed_incl_pfam"].isnull()
        | (df["fraction_of_RD_domains_removed_incl_pfam"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed_incl_pfam"].isnull()
        | (df["fraction_of_Bif_domains_removed_incl_pfam"] == 0)
    ),
    :,
]
axs[5].set_title("All effector domains\nin alt. iso.")
sns.swarmplot(
    data=all_retained,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[5],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    hue=hue,
    palette=palette,
)
axs[5].set_xticks([])
axs[5].set_xlabel("")
axs[5].get_legend().remove()

# annotate pbx1 and rfx3
pbx1_y = df.loc[
    (df["clone_acc_alt"] == "PBX1|2/2|02C05"), "activation_fold_change_log2"
].values[0]
rfx3_y = df.loc[
    (df["clone_acc_alt"] == "RFX3|3/5|08G08"), "activation_fold_change_log2"
].values[0]
rfx4_y = df.loc[
    (df["clone_acc_alt"] == "RFX3|4/5|11D09"), "activation_fold_change_log2"
].values[0]
creb5_y = df.loc[
    (df["clone_acc_alt"] == "CREB5|2/3|08A12"), "activation_fold_change_log2"
].values[0]
for point in axs[5].collections:
    for x, y in point.get_offsets():
        if np.isclose(pbx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[5].annotate(
                "PBX1-2",
                xy=(x, y),
                xytext=(-10, -10),
                textcoords="offset points",
                arrowprops=dict(
                    arrowstyle="-", connectionstyle="arc3,rad=0.3", color="black"
                ),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx3_y, y):
            print("found: %s, %s" % (x, y))
            axs[5].annotate(
                "RFX3-3",
                xy=(x, y),
                xytext=(-8, -7),
                textcoords="offset points",
                arrowprops=dict(
                    arrowstyle="-", connectionstyle="arc3,rad=-0.2", color="black"
                ),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx4_y, y):
            print("found: %s, %s" % (x, y))
            axs[5].annotate(
                "RFX3-4",
                xy=(x, y),
                xytext=(8, -8),
                textcoords="offset points",
                arrowprops=dict(
                    arrowstyle="-", connectionstyle="arc3,rad=0.3", color="black"
                ),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(creb5_y, y):
            print("found: %s, %s" % (x, y))
            axs[5].annotate(
                "CREB5-2",
                xy=(x, y),
                xytext=(-23, -8),
                textcoords="offset points",
                arrowprops=dict(
                    arrowstyle="-", connectionstyle="arc3,rad=0.3", color="black"
                ),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )

# missing stuff
incl = pd.concat(
    [
        tot_loss_activ,
        tot_loss_repr,
        tot_loss_both,
        partial_loss_activ,
        partial_loss_repr,
        all_retained,
    ]
)

no_annot = df.loc[
    (~df.index.isin(incl.index.values))
    & (pd.isnull(df["fraction_of_AD_domains_removed_incl_pfam"]))
    & (pd.isnull(df["fraction_of_RD_domains_removed_incl_pfam"]))
    & (pd.isnull(df["fraction_of_Bif_domains_removed_incl_pfam"]))
]
axs[6].set_title("No annotated\neffector domains")
sns.swarmplot(
    data=no_annot,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[6],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    color=sns.color_palette("flare")[0],
)
axs[6].set_xticks([])
axs[6].set_xlabel("")


# add colorbar
# mirror figure
fig2, axs2 = plt.subplots(1, 7, sharey=True, gridspec_kw=gs_kw)
fig2.set_size_inches(w=8.2, h=2)
map1 = axs2[6].imshow(np.stack([t, t]), cmap="flare", vmin=25, vmax=250)
cbar = fig.colorbar(map1, ax=axs[6], aspect=30, pad=0.2)
cbar.set_ticks([25, 75, 150, 250])
cbar.set_ticklabels(["<=25", "75", "150", ">=250"])
cbar.set_label("# AA in annotated domain", labelpad=0)


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylim(-7.5, 7.5)
    ax.axhline(y=0, color="black", linewidth=1, linestyle="dashed", zorder=1)
for ax in axs[1:]:
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(which="both", length=0)
    ax.set_ylabel("")
axs[0].set_ylabel("log2(activation fold change)")
fig.savefig(
    "../../figures/fig4/activation_vs_domain_removal_incl_pfam_colored_by_dom_length.pdf",
    bbox_inches="tight",
)


# ## 5. pie charts of PPI categories

# In[21]:


len(cats)


# In[22]:


len(cats.gene_symbol_partner.unique())


# In[23]:


cats.category.value_counts()


# In[24]:


y2h_nonan = y2h[~pd.isnull(y2h["Y2H_result"])]
len(y2h_nonan.db_gene_symbol.unique())


# In[25]:


ggi = y2h_nonan[["ad_gene_symbol", "db_gene_symbol"]].drop_duplicates()
ggi


# In[26]:


# limiting df to those that are in the y2h iso data
cats_y2h = cats[cats["gene_symbol_partner"].isin(ggi["db_gene_symbol"])]
len(cats_y2h)


# In[27]:


cats_dupe = (
    cats_y2h.groupby("gene_symbol_partner")["category"].agg("count").reset_index()
)
cats_dupe[cats_dupe["category"] > 1].head()


# In[28]:


# keeping partners that are assigned more than 1 category for now...
cats_y2h[cats_y2h["gene_symbol_partner"] == "BARD1"]


# In[29]:


def categorize_PPI_partner(row):
    if row["category"] == "TF":
        return "TF"
    elif row["category"] == "cofactor":
        return "cofactor"
    elif row["category"] == "signaling":
        return "signaling"
    else:
        return "other"


cats_y2h["partner_category"] = cats_y2h.apply(categorize_PPI_partner, axis=1)
cats_y2h.partner_category.value_counts()


# In[30]:


ys = np.array(
    [
        len(cats_y2h[cats_y2h["partner_category"] == "TF"]),
        len(cats_y2h[cats_y2h["partner_category"] == "cofactor"]),
        len(cats_y2h[cats_y2h["partner_category"] == "signaling"]),
        len(cats_y2h[cats_y2h["partner_category"] == "other"]),
    ]
)
labels = ["TF", "cofactor", "signaling", "other"]
colors = [
    sns.color_palette("Set2")[2],
    sns.color_palette("Set2")[1],
    sns.color_palette("Set2")[5],
    "darkgray",
]

fig, ax = plt.subplots(figsize=(1.2, 1.2), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(
    ys,
    colors=colors,
    labels=labels,
    autopct="%1.0f%%",
    startangle=-45,
    explode=(0.02, 0.2, 0.05, 0.05),
)

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")

fig.savefig(
    "../../figures/fig4/PPIs-gene-level-manual-categories_simplified.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# In[31]:


cofactor_partners = set(
    cats_y2h.loc[cats_y2h["category"] == "cofactor", "gene_symbol_partner"].unique()
)
coactivator_partners = set(
    cats_y2h.loc[
        cats_y2h["cofactor_type"] == "coactivator", "gene_symbol_partner"
    ].unique()
)
corepressor_partners = set(
    cats_y2h.loc[
        cats_y2h["cofactor_type"] == "corepressor", "gene_symbol_partner"
    ].unique()
)
both_partners = set(
    cats_y2h.loc[cats_y2h["cofactor_type"] == "both", "gene_symbol_partner"].unique()
)
signaling_partners = set(
    cats_y2h.loc[cats_y2h["category"] == "signaling", "gene_symbol_partner"].unique()
)
cofactor_animal_db = set(cof["Symbol"].unique())
tf_gene_symbols = set(load_human_tf_db()["HGNC symbol"].values)


# In[32]:


def categorize_cofactors(row):
    if row["partner_category"] == "cofactor":
        if row["gene_symbol_partner"] in both_partners:
            return "both"
        elif row["gene_symbol_partner"] in coactivator_partners:
            return "coactivator"
        elif row["gene_symbol_partner"] in corepressor_partners:
            return "corepressor"
        else:
            return "unknown"
    else:
        return "NA"


cats_y2h["cofactor_type"] = cats_y2h.apply(categorize_cofactors, axis=1)
cats_y2h.cofactor_type.value_counts()


# In[33]:


cofacs = cats_y2h[cats_y2h["partner_category"] == "cofactor"]

ys = np.array(
    [
        len(cofacs[cofacs["cofactor_type"] == "coactivator"]),
        len(cofacs[cofacs["cofactor_type"] == "corepressor"]),
        len(cofacs[cofacs["cofactor_type"] == "both"]),
        len(cofacs[cofacs["cofactor_type"] == "unknown"]),
    ]
)
labels = ["coactivator", "corepressor", "both", "unknown"]
colors = [
    sns.color_palette("Set2")[0],
    sns.color_palette("Set2")[3],
    sns.color_palette("Set2")[4],
    "darkgray",
]

fig, ax = plt.subplots(figsize=(1.2, 1.2), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(
    ys,
    colors=colors,
    labels=labels,
    autopct="%1.0f%%",
    startangle=60,
    explode=(0.05, 0.05, 0.05, 0.1),
)

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")

fig.savefig(
    "../../figures/fig4/PPIs-gene-level-manual-categories_cofactors.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# ## 6. plot the relationship between gain/loss of PPIs and changes in activity
#
# # TO-DO: Luke is updating the PPI list and we will use gain/loss of cofactors not separated by coactivator/corepressor and absolute fc plots

# In[34]:


def add_restricted_ppi_columns(pairs, rows, label):
    pairs_cf = pairs[["clone_acc_ref", "clone_acc_alt"]].copy()
    _add_PPI_columns(df=pairs_cf, y2h=y2h.loc[rows, :])
    return pd.merge(
        pairs,
        pairs_cf,
        how="left",
        on=["clone_acc_ref", "clone_acc_alt"],
        suffixes=("", "_" + label),
    )


# In[35]:


pairs = add_restricted_ppi_columns(
    pairs, rows=y2h["db_gene_symbol"].isin(cofactor_partners), label="cofactors"
)
pairs = add_restricted_ppi_columns(
    pairs,
    rows=y2h["db_gene_symbol"].isin(cofactor_animal_db),
    label="cofactors_animal_db",
)
pairs = add_restricted_ppi_columns(
    pairs,
    rows=(
        y2h["db_gene_symbol"].isin(coactivator_partners)
        & y2h["db_gene_symbol"].isin(hek_expressed_genes)
    ),
    label="coactivators_HEK",
)
pairs = add_restricted_ppi_columns(
    pairs,
    rows=(
        y2h["db_gene_symbol"].isin(coactivator_partners)
        & ~y2h["db_gene_symbol"].isin(hek_expressed_genes)
    ),
    label="coactivators_not_HEK",
)
pairs = add_restricted_ppi_columns(
    pairs, rows=~y2h["db_gene_symbol"].isin(cofactor_partners), label="not_cofactors"
)
pairs = add_restricted_ppi_columns(
    pairs, rows=y2h["db_gene_symbol"].isin(coactivator_partners), label="coactivators"
)
pairs = add_restricted_ppi_columns(
    pairs, rows=y2h["db_gene_symbol"].isin(corepressor_partners), label="corepressors"
)
pairs = add_restricted_ppi_columns(
    pairs,
    rows=(
        y2h["db_gene_symbol"].isin(corepressor_partners)
        & y2h["db_gene_symbol"].isin(hek_expressed_genes)
    ),
    label="corepressors_HEK",
)
pairs = add_restricted_ppi_columns(
    pairs, rows=y2h["db_gene_symbol"].isin(tf_gene_symbols), label="tfs"
)

pairs = add_restricted_ppi_columns(
    pairs,
    rows=(
        y2h["db_gene_symbol"].isin(signaling_partners)
        & y2h["db_gene_symbol"].isin(hek_expressed_genes)
    ),
    label="signaling_HEK",
)


# In[36]:


def bar_activation_vs_ppi(x, y, pairs=pairs, x_label=None, y_label=None, color=None):
    """
    TODO:
        - calculate p-value properly
            - this requires permuting in some smart way
            - one question is whether the genes are the number of independent data points or the isoforms are
            - I think the answer is the isoforms are

    """
    df = pairs.copy()
    if x_label is None:
        x_label = x
    if y_label is None:
        y_label = y
    if color is None:
        color = sns.color_palette("Set2")[1]
    fig, ax = plt.subplots(1, 1, figsize=(1.75, 1.5))

    def bin_delta_ppi(delta_ppi):
        if pd.isnull(delta_ppi):
            return np.nan
        if delta_ppi < 0:
            return "loss"
        elif delta_ppi > 0:
            return "gain"
        elif delta_ppi == 0:
            return "equal"
        else:
            raise ValueError(delta_ppi)

    df[x + "_binned"] = df[x].apply(bin_delta_ppi)
    sns.stripplot(
        data=df,
        x=x + "_binned",
        y=y,
        order=["loss", "equal", "gain"],
        alpha=0.75,
        color=color,
        linewidth=1,
        edgecolor="black",
        ax=ax,
    )
    if False:
        sns.pointplot(
            data=df,
            x=x + "_binned",
            y=y,
            order=["loss", "equal", "gain"],
            alpha=0.5,
            color="black",
            ax=ax,
        )
    if True:
        sns.boxplot(
            data=df,
            x=x + "_binned",
            y=y,
            order=["loss", "equal", "gain"],
            fliersize=0,
            color=color,
            ax=ax,
        )
        mimic_r_boxplot(ax)
    else:
        sns.violinplot(
            data=df,
            x=x + "_binned",
            y=y,
            order=["loss", "equal", "gain"],
            color="lightgrey",
            ax=ax,
        )
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    with_data = df[x].notnull() & df[y].notnull()
    n_pair = with_data.sum()
    n_iso = len(
        set(df.loc[with_data, ["clone_acc_ref", "clone_acc_alt"]].values.flatten())
    )
    n_gene = df.loc[with_data, "gene_symbol"].nunique()
    scc, p_scc = stats.spearmanr(
        df.loc[df[x].notnull() & df[y].notnull(), x].values,
        df.loc[df[x].notnull() & df[y].notnull(), y].values,
    )
    ax.text(
        s=f"{n_pair:d} pairs\n{n_iso:d} isoforms\n{n_gene:d} genes\nSpearman r = {scc:.2f}\np = {p_scc:.2f}",
        x=1.03,
        y=0.95,
        ha="left",
        va="top",
        transform=ax.transAxes,
    )
    # ax.set_ylim(-4, 4) # NOTE cuts outlier TODO add broken axis
    ax.axhline(y=0, color="black", linestyle="dashed", linewidth=1)
    for pos in ["top", "bottom", "right"]:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    fig.savefig(f"../../figures/fig4/{x}-vs-{y}_scatter.pdf", bbox_inches="tight")


# In[37]:


pairs["activation_abs_fold_change"] = pairs["activation_fold_change_log2"].abs()


# In[38]:


bar_activation_vs_ppi(
    pairs=pairs.loc[
        pairs["at_least_one_isoform_in_gene_abs_activation_gte_2fold"] == True, :
    ],
    x="PPI_delta_n",
    y="activation_fold_change_log2",
    x_label="Net difference in number of PPIs",
    y_label="log-2 fold-change in activation",
)


# In[39]:


bar_activation_vs_ppi(
    pairs=pairs.loc[
        pairs["at_least_one_isoform_in_gene_abs_activation_gte_2fold"] == True, :
    ],
    x="PPI_delta_n_corepressors_HEK",
    y="activation_fold_change_log2",
    x_label="Difference in number of corepressor PPIs\nrestricted to those expressed in HEK293 cells",
    y_label="log-2 fold-change in activation",
    color=sns.color_palette("Set2")[3],
)


# In[40]:


bar_activation_vs_ppi(
    pairs=pairs.loc[
        pairs["at_least_one_isoform_in_gene_abs_activation_gte_2fold"] == True, :
    ],
    x="PPI_delta_n_coactivators",
    y="activation_fold_change_log2",
    x_label="Difference in number of coactivator PPIs\nrestricted to those expressed in HEK293 cells",
    y_label="log-2 fold-change in activation",
    color=sns.color_palette("Set2")[0],
)


# In[41]:


bar_activation_vs_ppi(
    pairs=pairs.loc[
        pairs["at_least_one_isoform_in_gene_abs_activation_gte_2fold"] == True, :
    ],
    x="PPI_delta_n_signaling_HEK",
    y="activation_fold_change_log2",
    x_label="Difference in number of signaling PPIs\nrestricted to those expressed in HEK293 cells",
    y_label="log-2 fold-change in activation",
    color=sns.color_palette("Set2")[5],
)


# In[42]:


pairs[pairs["gene_symbol"] == "TCF4"].sort_values(by="activation_fold_change_log2")[
    [
        "clone_acc_ref",
        "clone_acc_alt",
        "at_least_one_isoform_in_gene_abs_activation_gte_2fold",
        "PPI_delta_n_coactivators_HEK",
        "activation_fold_change_log2",
    ]
]


# ## 7. dimerizing TFs: plot conservation of interactions across isoforms

# In[43]:


ppi = load_y2h_isoform_data(require_at_least_one_ppi_per_isoform=True)
ppi.head()


# In[44]:


n_iso_per_ppi = (
    (ppi.groupby(["ad_gene_symbol", "db_gene_symbol"]))["ad_clone_acc"]
    .agg("count")
    .reset_index()
)
n_true_iso_per_ppi = (
    (ppi[ppi["Y2H_result"] == True].groupby(["ad_gene_symbol", "db_gene_symbol"]))[
        "ad_clone_acc"
    ]
    .agg("count")
    .reset_index()
)
n_iso_per_ppi = n_iso_per_ppi.merge(
    n_true_iso_per_ppi, on=["ad_gene_symbol", "db_gene_symbol"], how="left"
)
n_iso_per_ppi["f_iso_positive"] = (
    n_iso_per_ppi["ad_clone_acc_y"] / n_iso_per_ppi["ad_clone_acc_x"]
)
n_iso_per_ppi = n_iso_per_ppi[["ad_gene_symbol", "db_gene_symbol", "f_iso_positive"]]
n_iso_per_ppi


# In[45]:


n_iso_per_ppi = pd.merge(
    n_iso_per_ppi,
    ppi.groupby(["ad_gene_symbol", "db_gene_symbol"])["Y2H_result"]
    .apply(lambda x: x.notnull().sum())
    .rename("n_iso_successfully_tested")
    .reset_index(),
    how="left",
    on=["ad_gene_symbol", "db_gene_symbol"],
)
if n_iso_per_ppi["n_iso_successfully_tested"].isnull().any():
    raise UserWarning("unexpected missing values")
n_iso_per_ppi = n_iso_per_ppi.loc[n_iso_per_ppi["n_iso_successfully_tested"] >= 2, :]


# In[46]:


tf_fam = load_tf_families()


# In[47]:


n_iso_per_ppi["db_is_tf"] = n_iso_per_ppi["db_gene_symbol"].isin(tf_fam)
n_iso_per_ppi["ad_tf_family"] = n_iso_per_ppi["ad_gene_symbol"].map(tf_fam)
n_iso_per_ppi["db_tf_family"] = n_iso_per_ppi["db_gene_symbol"].map(tf_fam)


# In[48]:


def tf_tf_dimer_ppi_catagories(row):
    is_dimer_ad = row["ad_tf_family"] in DIMERIZING_TF_FAMILIES
    if pd.isnull(row["db_tf_family"]):
        if is_dimer_ad:
            return "Dimerizing TF / other"
        else:
            return "Non-dimerizing TF / other"
    else:  # TF-TF PPI
        if is_dimer_ad:
            if row["db_tf_family"] == row["ad_tf_family"]:
                return "Dimerzing TF-TF PPI"
            else:
                return "Outside family TF-TF PPI"
        else:
            if row["db_tf_family"] == row["ad_tf_family"]:
                return "Non-dimer TF-TF PPI"
            else:
                return "Outside family TF-TF PPI"


n_iso_per_ppi["dimer_cat"] = n_iso_per_ppi.apply(tf_tf_dimer_ppi_catagories, axis=1)

n_iso_per_ppi.head()


# In[49]:


# permutation test
cats = [
    "Dimerzing TF-TF PPI",
    "Dimerizing TF / other",
    "Outside family TF-TF PPI",
    "Non-dimer TF-TF PPI",
    "Non-dimerizing TF / other",
]

rnd_vals = {}
obs_val = {}
for cat in cats[1:]:
    rnd_vals[cat] = []
    obs_val[cat] = (
        n_iso_per_ppi.loc[
            n_iso_per_ppi["dimer_cat"] == cats[0], "f_iso_positive"
        ].mean()
        - n_iso_per_ppi.loc[n_iso_per_ppi["dimer_cat"] == cat, "f_iso_positive"].mean()
    )
for _i in range(1000):
    rand_cats = (
        n_iso_per_ppi["dimer_cat"].sample(frac=1).reset_index(drop=True)
    )  # equivalent to shuffle
    for cat in cats[1:]:
        rnd_vals[cat].append(
            n_iso_per_ppi.loc[rand_cats == cats[0], "f_iso_positive"].mean()
            - n_iso_per_ppi.loc[rand_cats == cat, "f_iso_positive"].mean()
        )

pvals = {}
for cat in cats[1:]:
    pvals[cat] = sum(v >= obs_val[cat] for v in rnd_vals[cat]) / len(rnd_vals[cat])
pvals


# In[50]:


# n for each bar
# axis labels etc.
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2, h=2)
violinplot_reflected(
    data=n_iso_per_ppi,
    x="dimer_cat",
    y="f_iso_positive",
    inner=None,
    cut=0,
    color=sns.color_palette("Set2")[0],
    order=cats,
    ax=ax,
)

pairs_to_test = [(cats[0], cat) for cat in cats[1:]]
annotator = Annotator(
    ax=ax,
    pairs=pairs_to_test,
    data=n_iso_per_ppi,
    x="dimer_cat",
    y="f_iso_positive",
    order=cats,
)
annotator.configure(loc="outside")

# TODO: set this programmatically
annotator.annotate_custom_annotations(["****", "****", "****", "****"])

# TODO: set the statistic and error bar calculation
sns.pointplot(
    data=n_iso_per_ppi,
    x="dimer_cat",
    y="f_iso_positive",
    order=cats,
    color="black",
    linestyles="",
    edgecolors="black",
    linewidth=1,
    ax=ax,
)
ax.set_ylabel("Fraction isoforms interacting")
ax.set_xlabel("")

ax.set_xticklabels(cats, rotation=30, ha="right", va="top")

ax.set_ylim(0, 1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

fig.savefig("../../figures/fig4/n_iso_ppi_by_dimer_cat.pdf", bbox_inches="tight")


# ## 8. dimerizing TFs: plot M1H change

# In[51]:


len(pairs)


# In[52]:


len(pairs[~pd.isnull(pairs["activation_abs_fold_change"])])


# In[53]:


pairs[~pd.isnull(pairs["activation_abs_fold_change"])].groupby(
    "is_dimerizing_TF_family"
)["activation_abs_fold_change"].agg("mean")


# In[54]:


# absolute difference
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=1.2, h=2)
sns.boxplot(
    data=pairs,
    x="is_dimerizing_TF_family",
    y="activation_abs_fold_change",
    order=[True, False],
    ax=ax,
    boxprops={"facecolor": "w"},
)  # , fliersize=0)
# sns.stripplot(data=df, x='is_dimerizing_tf', y='m1h_abs_change', order=[True, False], ax=ax, color='.25')
ax.set_xlabel("Is dimerizing TF")
ax.set_ylabel(
    "absolute log2 fold change of M1H activation\nbetween reference and alternative isoforms"
)
ax.xaxis.set_tick_params(rotation=90, length=0)
for loc in ["bottom", "right", "top"]:
    ax.spines[loc].set_visible(False)

# p-value
pval = stats.mannwhitneyu(
    pairs.loc[
        pairs["is_dimerizing_TF_family"]
        & pairs["activation_abs_fold_change"].notnull(),
        "activation_abs_fold_change",
    ].values,
    pairs.loc[
        ~pairs["is_dimerizing_TF_family"]
        & pairs["activation_abs_fold_change"].notnull(),
        "activation_abs_fold_change",
    ].values,
)[1]
print(pval)
annotator = Annotator(
    ax=ax,
    pairs=[(True, False)],
    data=pairs,
    x="is_dimerizing_TF_family",
    y="activation_abs_fold_change",
    order=[True, False],
)
annotator.configure(loc="outside")
# TODO: set this programmatically
annotator.annotate_custom_annotations(["P = {:.3f}".format(pval)])
plt.savefig(
    "../../figures/fig4/M1H-absolute-change_by-dimerizing-TF-family_boxplot.pdf",
    bbox_inches="tight",
)


# ## 9. dimerizing TFs: retained interactions heatmap

# In[55]:


# TF-TF binding
tftf = (
    ppi.loc[ppi["Y2H_result"] == True, ["ad_gene_symbol", "db_gene_symbol"]]
    .drop_duplicates()
    .copy()
)
tftf["ad_dbd"] = tftf["ad_gene_symbol"].map(tf_fam)
tftf["db_dbd"] = tftf["db_gene_symbol"].map(tf_fam)
tftf = tftf.dropna()
tftf.head()


# In[56]:


# TF-TF rewiring
tftf = pd.merge(
    tftf,
    (
        ppi.groupby(["ad_gene_symbol", "db_gene_symbol"])["Y2H_result"].apply(
            lambda x: (x == False).sum()
        )
        / ppi.groupby(["ad_gene_symbol", "db_gene_symbol"])["Y2H_result"].apply(
            lambda x: (x.notnull().sum())
        )
    ).reset_index(),
    how="left",
    on=["ad_gene_symbol", "db_gene_symbol"],
)


# In[57]:


tftf.loc[tftf["db_dbd"] == "bHLH", "ad_gene_symbol"].value_counts()


# In[58]:


DIMERIZING_TF_FAMILIES


# In[59]:


fams = tf_fam.value_counts().index
fams = list(
    filter(lambda x: x in tftf["ad_dbd"].unique() or x in tftf["db_dbd"].unique(), fams)
)
len(fams)


# In[60]:


num_pairs = [
    ((tftf["ad_dbd"] == x) & (tftf["db_dbd"] == y)).sum() for x in fams for y in fams
]
num_pairs = np.reshape(num_pairs, (-1, len(fams)))
num_pairs = pd.DataFrame(num_pairs, index=fams, columns=fams)
num_pairs.head()


# In[61]:


rewiring = [
    tftf.loc[(tftf["ad_dbd"] == x) & (tftf["db_dbd"] == y), "Y2H_result"].mean()
    for x in fams
    for y in fams
]
rewiring = np.reshape(rewiring, (-1, len(fams)))
rewiring = pd.DataFrame(rewiring, index=fams, columns=fams)
rewiring.head()


# In[62]:


# limit to fams w 10 DB
num_pairs_sum = num_pairs.sum(axis=0)
filt_fams = num_pairs_sum[num_pairs_sum >= 3]
filt_fams


# In[63]:


rewiring_filt = rewiring.loc[list(filt_fams.index), list(filt_fams.index)]
rewiring_filt


# In[64]:


fig = plt.figure(figsize=(2.2, 2))
g = sns.heatmap(rewiring_filt, cmap="viridis_r", cbar_kws={"label": "rewiring score"})

# highlight the squares corresponding to dimerizing pairs
for i, fam in enumerate(list(filt_fams.index)):
    if fam in DIMERIZING_TF_FAMILIES:
        g.add_patch(Rectangle((i, i), 1, 1, fill=False, edgecolor="black", lw=1))

g.set_xlabel("TF isoform family")
g.set_ylabel("Y2H partner TF family")
fig.savefig(
    "../../figures/fig4/Dimerizing_PPI_heatmap.small.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# In[65]:


fig = plt.figure(figsize=(4.5, 4.5))
g = sns.heatmap(
    rewiring,
    cmap="viridis_r",
    cbar_kws={"label": "rewiring score", "aspect": 40},
    annot=num_pairs,
    annot_kws={"fontsize": 6},
)

# highlight the squares corresponding to dimerizing pairs
for i, fam in enumerate(list(rewiring.index)):
    if fam in DIMERIZING_TF_FAMILIES:
        g.add_patch(Rectangle((i, i), 1, 1, fill=False, edgecolor="black", lw=1))

g.set_xlabel("TF isoform family")
g.set_ylabel("Y2H partner TF family")
fig.savefig(
    "../../figures/fig4/Dimerizing_PPI_heatmap.all.pdf",
    dpi="figure",
    bbox_inches="tight",
)


# ## 10. isoform example vignettes

# In[66]:


# reload data since we edited dfs above
y2h = load_full_y2h_data_including_controls()
y1h = load_y1h_pdi_data(add_missing_data=True)
m1h = load_m1h_activation_data(add_missing_data=True)
isoforms = load_valid_isoform_clones()
y2h = y2h.loc[y2h["ad_clone_acc"].isin(isoforms["clone_acc"]).values, :]
y1h = y1h.loc[y1h["clone_acc"].isin(isoforms["clone_acc"]).values, :]
m1h = m1h.loc[m1h["clone_acc"].isin(isoforms["clone_acc"].values), :]

tfs = load_annotated_TFiso1_collection()


# ### RFX3

# In[67]:


gene_name = "RFX3"


# In[68]:


fig, ax = plt.subplots(1, 1, figsize=(2, 0.8))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(0, 6))
plt.savefig(
    "../../figures/fig4/{}_m1h-profile.pdf".format(gene_name), bbox_inches="tight"
)


# In[69]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig(
    "../../figures/fig4/{}_exon_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# In[70]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig(
    "../../figures/fig4/{}_protein_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# In[71]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig(
    "../../figures/fig4/{}_protein_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# ### PBX1

# In[72]:


gene_name = "PBX1"


# In[73]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig(
    "../../figures/fig4/{}_y2h-profile.pdf".format(gene_name), bbox_inches="tight"
)


# In[74]:


fig, ax = plt.subplots(1, 1, figsize=(2, 0.7))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(0, 2.5))
plt.savefig(
    "../../figures/fig4/{}_m1h-profile.pdf".format(gene_name), bbox_inches="tight"
)


# In[75]:


fig, ax = plt.subplots(figsize=(5, 2))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig(
    "../../figures/fig4/{}_exon_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# In[76]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig(
    "../../figures/fig4/{}_protein_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# ### CREB5

# In[77]:


gene_name = "CREB5"


# In[78]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(0.8, 0.8))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig(
    "../../figures/fig4/{}_y2h-profile.pdf".format(gene_name), bbox_inches="tight"
)


# In[79]:


fig, ax = plt.subplots(1, 1, figsize=(2, 0.6))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(-2.2, 2.2))
plt.savefig(
    "../../figures/fig4/{}_m1h-profile.pdf".format(gene_name), bbox_inches="tight"
)


# In[80]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig(
    "../../figures/fig4/{}_exon_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# In[81]:


fig, ax = plt.subplots(figsize=(5, 0.7))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig(
    "../../figures/fig4/{}_protein_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# In[125]:


gene_name = "DLX1"
fig, ax = plt.subplots(figsize=(5, 0.7))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig(
    "../../figures/fig4/{}_protein_diagram.pdf".format(gene_name),
    bbox_inches="tight",
    dpi="figure",
)


# In[ ]:
