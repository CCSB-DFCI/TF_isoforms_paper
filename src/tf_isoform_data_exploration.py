# coding: utf-8

# # TF isoform data
#
#
# ## TODO
#
# - check valid clones only in all data
# - get sequence information
# - remove partners that didn't test positive with any isoform
# - remove single isoforms
#     - but need to add first to paralogs
# - Look into effect of number of PPIs per TF
# - combine categories to get paralog data

# In[1]:


import os
from itertools import combinations

import numpy as np
from scipy import stats
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

import ccsblib
from plotting import validation_plot

from data_loading import (
    load_full_y2h_data_including_controls,
    load_y1h_pdi_data,
    load_m1h_activation_data,
    load_valid_isoform_clones,
    load_seq_comparison_data,
)
from isoform_pairwise_metrics import pairs_of_isoforms_comparison_table

get_ipython().run_line_magic("matplotlib", "inline")

y2h = load_full_y2h_data_including_controls()
y1h = load_y1h_pdi_data()
# m1h = load_m1h_activation_data()
isoforms = load_valid_isoform_clones()
idents = load_seq_comparison_data()

# tfs = tfs.drop(columns=['tpm_stdev'])
# tfs = tfs.set_index(['gene_symbol', 'isoacc', 'tiss'])
# tfs = tfs.unstack(level='tiss')
# tfs.columns = tfs.columns.get_level_values(1)

# tfs.to_csv('tf_isoform_tissue_tpms.tsv', sep='\t')


# In[2]:


y2h.head()


# In[3]:


print("Numbers for the isoform experiment (not paralogs or controls):")
print(
    y2h.loc[(y2h["category"] == "tf_isoform_ppis"), "ad_clone_acc"].nunique(),
    "isoforms attempted to test for PPIs",
)
print(
    y2h.loc[
        (y2h["category"] == "tf_isoform_ppis") & (y2h["Y2H_result"].notnull()),
        "ad_clone_acc",
    ].nunique(),
    "isoforms successfully tested for at least one PPI",
)
print(
    y2h.loc[
        (y2h["category"] == "tf_isoform_ppis") & (y2h["Y2H_result"] == True),
        "ad_clone_acc",
    ].nunique(),
    "isoforms with at least one positive PPI",
)
print(
    ((y2h["category"] == "tf_isoform_ppis") & (y2h["Y2H_result"] == True)).sum(),
    "positive PPIs, involving",
    y2h.loc[
        (y2h["category"] == "tf_isoform_ppis") & (y2h["Y2H_result"] == True),
        "db_gene_symbol",
    ].nunique(),
    "protein partners",
)


# In[4]:


yang = pd.read_excel(
    "../data/external/Yang_et_al_Cell_2014_Table_S2.xlsx", sheet_name="2B-Isoform PPIs"
)


# In[5]:


# isoforms with >= 2 tested and at least one positive PPI
# count pairwise combinations and genes
yang["2_or_more_iso_per_gene"] = yang["Gene_Symbol"].map(
    (yang.groupby("Gene_Symbol")["Isoform_ID"].nunique() > 1)
)


# In[6]:


yang["2_or_more_partners_per_iso"] = yang["Isoform_ID"].map(
    yang["Isoform_ID"].value_counts() >= 2
)


# In[7]:


yang["2_or_more_iso_with_2_or_more_partners_per_gene"] = (
    yang["Gene_Symbol"]
    .map(
        yang.loc[yang["2_or_more_partners_per_iso"], :]
        .groupby("Gene_Symbol")["Isoform_ID"]
        .nunique()
        >= 2
    )
    .fillna(False)
)


# In[8]:


yang.loc[
    yang["2_or_more_iso_with_2_or_more_partners_per_gene"], "Gene_Symbol"
].nunique()


# In[9]:


import itertools


def count_pairs_per_gene(df):
    pairs_count = 0
    for iso_a, iso_b in itertools.combinations(df["Isoform_ID"].unique(), 2):
        if (
            (
                (df["Isoform_ID"] == iso_a) & (df["Interaction_Found"] == "positive")
            ).any()
            and (
                (df["Isoform_ID"] == iso_b) & (df["Interaction_Found"] == "positive")
            ).any()
            and len(
                set(df.loc[df["Isoform_ID"] == iso_a, "Interactor_ID"]).intersection(
                    set(df.loc[df["Isoform_ID"] == iso_b, "Interactor_ID"])
                )
            )
            >= 2
        ):
            pairs_count += 1
    return pairs_count


iso_pairs_per_gene = yang.groupby("Gene_Symbol").apply(count_pairs_per_gene)
print(iso_pairs_per_gene.sum(), (iso_pairs_per_gene >= 1).sum())


# In[10]:


yang.loc[yang["2_or_more_iso_with_2_or_more_partners_per_gene"], :]


# In[11]:


isoforms["gene"].nunique()


# In[12]:


isoforms.shape


# In[13]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(3.5, 2.5)
isoforms["gene"].value_counts().plot.hist(ax=ax, range=(0.75, 8.25), bins=15)
ax.set_ylabel("Number of TF genes")
ax.set_xlabel("Number of isoforms per gene")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xticks(range(1, 9))
plt.savefig("../figures/isoform_clones_per_tf_gene.pdf", bbox_inches="tight")


# In[14]:


from data_loading import load_tf_families

all_tfs = load_tf_families()
all_tfs.shape[0]


# In[15]:


iso_y1h_pos = set(y1h.loc[y1h[y1h.columns[2:]].any(axis=1), "clone_acc"].values)


# In[16]:


iso_gte_1_pos_ppi_iso_data_only = set(
    y2h.loc[
        (y2h["category"] == "tf_isoform_ppis") & (y2h["Y2H_result"] == True),
        "ad_clone_acc",
    ].unique()
)
iso_gte_1_pos_ppi_all_data = set(
    y2h.loc[(y2h["Y2H_result"] == True), "ad_clone_acc"].unique()
)
pos_in_only_non_isoform_ppi_data = iso_gte_1_pos_ppi_all_data.difference(
    iso_gte_1_pos_ppi_iso_data_only
)
print(
    len(pos_in_only_non_isoform_ppi_data),
    "isoforms that are Y2H positive in a dataset outside the isoform",
)
print(pos_in_only_non_isoform_ppi_data)
y2h.loc[
    (y2h["ad_clone_acc"].isin(pos_in_only_non_isoform_ppi_data))
    & (y2h["Y2H_result"] == True),
    :,
]


# In[17]:


y2h = load_full_y2h_data_including_controls()
# restict to TF isoform data (i.e. not paralogs etc.)
ppi = y2h.loc[
    (y2h["category"] == "tf_isoform_ppis"),
    ["category", "ad_clone_acc", "ad_gene_symbol", "db_gene_symbol", "Y2H_result"],
].copy()
# at least one positive per PPI partner
ppi = ppi.loc[
    ppi.groupby(["ad_gene_symbol", "db_gene_symbol"])["Y2H_result"].transform(
        lambda row: (row == True).any()
    ),
    :,
]
# at least one successfully tested PPI per isoform
ppi = ppi.loc[
    ppi.groupby("ad_clone_acc")["Y2H_result"].transform(lambda x: (x.notnull().any())),
    :,
]
# at least two partners per isoform
ppi = ppi.loc[
    ppi.groupby("ad_gene_symbol")["ad_clone_acc"].transform(lambda x: x.nunique() >= 2),
    :,
]
y1h = load_y1h_pdi_data()
m1h = load_m1h_activation_data()
isoforms = load_valid_isoform_clones()
iso_pairs = pairs_of_isoforms_comparison_table(
    isoforms=isoforms, y2h=ppi, y1h=y1h, m1h=m1h
)
iso_pairs["both_iso_y2h_pos"] = iso_pairs["clone_acc_a"].isin(
    iso_gte_1_pos_ppi_iso_data_only
) & iso_pairs["clone_acc_b"].isin(iso_gte_1_pos_ppi_iso_data_only)
iso_pairs["both_iso_y2h_pos_all_data"] = iso_pairs["clone_acc_a"].isin(
    iso_gte_1_pos_ppi_all_data
) & iso_pairs["clone_acc_b"].isin(iso_gte_1_pos_ppi_all_data)


# In[18]:


iso_pairs["both_iso_y2h_or_y1h_pos"] = iso_pairs["clone_acc_a"].isin(
    iso_gte_1_pos_ppi_all_data.union(iso_y1h_pos)
) & iso_pairs["clone_acc_b"].isin(iso_gte_1_pos_ppi_all_data.union(iso_y1h_pos))


# In[19]:


iso_pairs["isoform_a"] = iso_pairs["clone_acc_a"].apply(
    lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
)
iso_pairs["isoform_b"] = iso_pairs["clone_acc_b"].apply(
    lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
)


# In[20]:


iso_pairs.loc[
    (iso_pairs["ppi_n_tested"] >= 2) & iso_pairs["both_iso_y2h_pos_all_data"], :
].shape


# In[21]:


genes_nonzero_pair = set(
    iso_pairs.loc[
        (iso_pairs["ppi_n_tested"] >= 2) & iso_pairs["both_iso_y2h_pos_all_data"],
        "tf_gene_symbol",
    ].unique()
)
len(genes_nonzero_pair)


# In[22]:


iso_pairs.loc[
    (iso_pairs["ppi_n_tested"] >= 2) & iso_pairs["both_iso_y2h_or_y1h_pos"], :
].shape


# In[23]:


# Y1H rescue
genes_y1h_rescue = set(
    iso_pairs.loc[
        (iso_pairs["ppi_n_tested"] >= 2) & iso_pairs["both_iso_y2h_or_y1h_pos"],
        "tf_gene_symbol",
    ].unique()
)
print(len(genes_y1h_rescue))
print(genes_y1h_rescue.difference(genes_nonzero_pair))


# In[24]:


rescue_pairs = iso_pairs.loc[
    (iso_pairs["ppi_n_tested"] >= 2)
    & (iso_pairs["both_iso_y2h_or_y1h_pos"] ^ iso_pairs["both_iso_y2h_pos_all_data"]),
    :,
]
rescue_iso = set(
    rescue_pairs[["clone_acc_a", "clone_acc_b"]].values.flatten()
).difference(iso_gte_1_pos_ppi_all_data)
print(len(rescue_iso))
print(rescue_iso)
rescue_pairs


# In[25]:


(
    iso_pairs.loc[
        (iso_pairs["ppi_n_tested"] >= 2) & iso_pairs["both_iso_y2h_pos_all_data"],
        [
            "tf_gene_symbol",
            "isoform_a",
            "isoform_b",
            "ppi_n_tested",
            "ppi_n_shared",
            "ppi_n_min",
            "ppi_jaccard",
        ],
    ].to_csv("../output/non_zero_isoform_pairs.tsv", index=False, sep="\t")
)


# In[26]:


iso_pairs.loc[
    (iso_pairs["ppi_n_tested"] >= 2)
    & ~iso_pairs["both_iso_y2h_pos"]
    & iso_pairs["both_iso_y2h_pos_all_data"],
    :,
]


# In[27]:


iso_pairs.loc[
    (iso_pairs["ppi_n_tested"] >= 2)
    & iso_pairs["both_iso_y2h_pos"]
    & (iso_pairs["both_iso_y2h_pos_all_data"] ^ iso_pairs["both_iso_y2h_pos"]),
    :,
]


# In[28]:


ppi


# In[29]:


iso_pairs["pdi_n_max"] = iso_pairs["pdi_n_tested"] - (
    iso_pairs["pdi_n_min"] - iso_pairs["pdi_n_shared"]
)


# In[30]:


iso_pairs.loc[
    iso_pairs["pdi_jaccard"].notnull(),
    [
        "tf_gene_symbol",
        "clone_acc_a",
        "clone_acc_b",
        "pdi_jaccard",
        "pdi_n_min",
        "pdi_n_max",
    ],
].to_csv("../output/PDI_isoform_pairs.tsv", sep="\t", index=False)


# In[31]:


iso_pairs.shape


# In[32]:


(isoforms.groupby("gene").size() >= 2).sum()


# In[33]:


iso_pairs.notnull().sum()


# In[34]:


fig, ax = plt.subplots(1, 1)
ax.hist(iso_pairs["ppi_jaccard"], range=(0, 1), bins=10)
ax.set_ylabel("Number of pairs of isoforms")
ax.set_xlabel("PPI: Jaccard similarity")
plt.savefig("../figures/ppi_jaccard_dist.pdf")


# In[35]:


fig, ax = plt.subplots(1, 1)
ax.hist(iso_pairs["aa_seq_pct_identity"], range=(0, 100), bins=25)
ax.set_ylabel("Number of pairs of isoforms")
ax.set_xlabel("AA sequence % identity")
plt.savefig("../figures/aa_seq_pct_id_hist.pdf")


# In[36]:


# for col in iso_pairs.columns[4:]:
#    plt.hist(iso_pairs[col], bins=30)
#    plt.xlabel(col)
#    plt.ylabel('count')
#    plt.show()


# In[37]:


iso_pairs.head()


# ### Explore PDI/PPI profile differences and sequence similarity
#
# Look into:
# - PPI and PDI Jaccard distribution
#   - Faceted by degree
# - PPI versus PDI Jaccard/Simpson
#   - Faceted by number of interactors
# - PPI/PDI Jaccard versus sequence similarity
#
# TODO
# - Compare randomly selected isoforms from different genes (from paralog data)

# In[38]:


#  Jaccard/Simpson for PPI/PDI, histograms
def make_jaccard_simpson_hist(int_type, iso_pairs=iso_pairs, main_title="", nb=30):
    # int_type - ppi or pdi
    # nb - number of bins
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 2.5))
    fig.suptitle(main_title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.87])  # adjust to allow for subtitle

    ax1.hist(iso_pairs[int_type + "_jaccard"], bins=nb)
    ax1.set_title(int_type.upper() + " Jaccard")
    ax1.set_xlabel("Jaccard index")
    ax1.set_ylabel("Number of isoform pairs")

    ax2.hist(iso_pairs[int_type + "_simpson"], bins=nb)
    ax2.set_title(int_type.upper() + " Simpson")
    ax2.set_xlabel("Simpson")


make_jaccard_simpson_hist("ppi", iso_pairs, "All")
make_jaccard_simpson_hist("pdi", iso_pairs, "All")
make_jaccard_simpson_hist(
    "ppi", iso_pairs.loc[iso_pairs["ppi_n_min"] >= 1, :], main_title="Min. degree 1"
)
make_jaccard_simpson_hist(
    "pdi", iso_pairs.loc[iso_pairs["pdi_n_min"] >= 1, :], main_title="Min. degree 1"
)


# In[39]:


# Jaccard vs Simpson for PPI and PDI
def make_jaccard_simpson_plot(df):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

    alpha = 0.3
    ax1.scatter(df.ppi_jaccard, df.ppi_simpson, alpha=alpha)
    ax1.set_title("PPI")
    ax1.set_aspect("equal")
    ax2.scatter(df.pdi_jaccard, df.pdi_simpson, alpha=alpha)
    ax2.set_title("PDI")
    ax2.set_aspect("equal")

    fig.tight_layout()
    fig.text(0.5, 0.05, "Jaccard index", ha="center")
    fig.text(-0.02, 0.5, "Simpson", va="center", rotation="vertical")


make_jaccard_simpson_plot(
    iso_pairs
)  # only plotted for cases  with 1+ interactions for both iso. of a pair


# In[40]:


# Is there a relationship between degree and interaction profile similarity?
# For now use number of pdi/ppi tested as estimate of degree.
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True, figsize=(12, 6))
alpha = 0.3

ax1.scatter(iso_pairs.ppi_n_tested, iso_pairs.ppi_jaccard, alpha=alpha)
ax1.set_xlabel("Number of PPIs tested for TF in isoform pair")
ax1.set_ylabel("Jaccard index for isoform pair")

ax2.scatter(iso_pairs.pdi_n_tested, iso_pairs.pdi_jaccard, alpha=alpha)
ax2.set_xlabel("Number of PDIs tested for TF in isoform pair")
ax2.set_ylabel("Jaccard index for isoform pair")

fig.tight_layout()

# Luke question - how to determine if there is no bias between degree and Jaccard?
stats.spearmanr(iso_pairs.ppi_n_tested, iso_pairs.ppi_jaccard)
stats.spearmanr(iso_pairs.pdi_n_tested, iso_pairs.pdi_jaccard)


# In[41]:


iso_pairs.sort_values(["ppi_jaccard", "aa_seq_pct_identity"], ascending=[True, False])


# In[42]:


x = iso_pairs.loc[iso_pairs["ppi_jaccard"].notnull(), "aa_seq_pct_identity"].values
y = iso_pairs.loc[iso_pairs["ppi_jaccard"].notnull(), "ppi_jaccard"].values
fig, ax = plt.subplots(1, 1)
ax.scatter(x, y, alpha=0.3)
ax.set_ylabel("PPI Jaccard similarity")
ax.set_xlabel("AA sequence % identity")
print(stats.spearmanr(x, y))
plt.savefig("../figures/ppi_jaccard_vs_aa_id.pdf", bbox_inches="tight")


# In[43]:


x = iso_pairs.loc[iso_pairs["pdi_jaccard"].notnull(), "aa_seq_pct_identity"].values
y = iso_pairs.loc[iso_pairs["pdi_jaccard"].notnull(), "pdi_jaccard"].values
fig, ax = plt.subplots(1, 1)
ax.scatter(x, y, alpha=0.3)
ax.set_ylabel("PDI Jaccard similarity")
ax.set_xlabel("AA sequence % identity")
print(stats.spearmanr(x, y))
plt.savefig("../figures/pdi_jaccard_vs_aa_id.pdf", bbox_inches="tight")


# In[44]:


x = iso_pairs.loc[
    iso_pairs["activation_fold_change"].notnull(), "aa_seq_pct_identity"
].values
y = iso_pairs.loc[
    iso_pairs["activation_fold_change"].notnull(), "activation_fold_change"
].values
fig, ax = plt.subplots(1, 1)
ax.scatter(x, y, alpha=0.3)
ax.set_ylabel("M1H fold change activation")
ax.set_xlabel("AA sequence % identity")
print(stats.spearmanr(x, y))
plt.savefig("../figures/m1h_fold_change_vs_aa_id.pdf", bbox_inches="tight")


# In[45]:


# Jaccard/Simpson versus similarity for PPI and PDI
def make_profile_diff_vs_seq_plot(df):
    fig, axs = plt.subplots(1, 4, sharex=True, sharey=True, figsize=(15, 3))
    fig.tight_layout()
    (ax1, ax2, ax3, ax4) = axs

    for ax in axs:
        ax.set_xlabel("AA seq. identity (%)")

    alpha = 0.3

    ax1.scatter(df.aa_seq_pct_id, df.ppi_jaccard, alpha=alpha)
    ax1.set_title("PPI")
    ax1.set_ylabel("Jaccard index")
    ax2.scatter(df.aa_seq_pct_id, df.ppi_simpson, alpha=alpha)
    ax2.set_title("PPI")
    ax2.set_ylabel("Simpson")
    ax3.scatter(df.aa_seq_pct_id, df.pdi_jaccard, alpha=alpha)
    ax3.set_title("PDI")
    ax3.set_ylabel("Jaccard index")
    ax4.scatter(df.aa_seq_pct_id, df.pdi_simpson, alpha=alpha)
    ax4.set_title("PDI")
    ax4.set_ylabel("Simpson")


make_profile_diff_vs_seq_plot(iso_pairs)


# In[46]:


# Jaccard/Simpson versus similarity for M1H
fig, ax1 = plt.subplots(1, sharex=True, sharey=True)
fig.tight_layout()
alpha = 0.3
ax1.scatter(iso_pairs.aa_seq_pct_id, iso_pairs.activation_fold_change, alpha=alpha)
ax1.set_xlabel("AA seq. identity (%)")
ax1.set_ylabel("M1H activation fold change")
plt.savefig("../figures/M1H_fold_change_vs_aa_id.pdf", bbox_inches="tight")


# In[47]:


# get isoform-specific y1h and y2h degree
deg = y1h.loc[:, y1h.columns[2:]].sum(axis=1).rename("y1h_degree")
y1h_nd = y1h[["gene_symbol", "clone_acc"]].join(deg, how="left")
y2h_nd = (
    y2h.loc[y2h.category == "tf_isoform_ppis", :]
    .groupby("ad_clone_acc")["ad_clone_acc"]
    .count()
    .rename("y2h_degree")
)
int_nd = y1h_nd.join(y2h_nd, how="outer", on="clone_acc")
int_nd.head()


# In[48]:


# plot the degree distributions
# degere of y1h
plt.hist(int_nd.y1h_degree.dropna(), bins=30)
plt.xlabel("Degree of TF")
plt.ylabel("Count")
plt.show()
plt.hist(int_nd.y2h_degree.dropna(), bins=30)
plt.xlabel("Degree of TF")
plt.ylabel("Count")


# In[49]:


# compare pdi versus ppi degree
fig, (ax1, ax2) = plt.subplots(1, 2)
alpha = 0.3

ax1.scatter(int_nd.y1h_degree, int_nd.y2h_degree, alpha=alpha)
ax1.set_aspect("equal")
ax1.set_xlabel("PDI degree")
ax1.set_ylabel("PPI degree")

ax2.scatter(int_nd.y1h_degree, int_nd.y2h_degree, alpha=alpha)
ax2.set_aspect("equal")
ax2.set_xlabel("PDI degree")
ax2.set_ylabel("PPI degree")
ax2.set_title("Zoomed in")
ax2.set_xlim(0, 20)
ax2.set_xticks(range(0, 20, 2))
ax2.set_ylim(0, 20)
ax2.set_yticks(range(0, 20, 2))
# ax2.set_yscale('log')
# ax2.set_xscale('log')


# In[50]:


x = "ppi_jaccard"
y = "activation_fold_change"
xy = iso_pairs.loc[iso_pairs[x].notnull() & iso_pairs[y].notnull(), :]
plt.scatter(xy[x], xy[y])
plt.xlabel("PPI Jaccard index")
plt.ylabel("Activation fold change")
print(stats.spearmanr(xy[x], xy[y]))
plt.savefig("../figures/ppi_vs_activation.pdf", bbox_inches="tight")


# In[51]:


x = "pdi_jaccard"
y = "activation_fold_change"
xy = iso_pairs.loc[iso_pairs[x].notnull() & iso_pairs[y].notnull(), :]
print(xy.shape)
plt.scatter(xy[x], xy[y])
plt.xlabel("PDI Jaccard index")
plt.ylabel("Activation fold change")
print(stats.spearmanr(xy[x], xy[y]))
plt.savefig("../figures/pdi_vs_activation.pdf", bbox_inches="tight")


# In[52]:


x = "pdi_jaccard"
y = "ppi_jaccard"
xy = iso_pairs.loc[iso_pairs[x].notnull() & iso_pairs[y].notnull(), :]
print(xy.shape)
plt.scatter(xy[x], xy[y], alpha=0.3)
plt.xlabel("PDI Jaccard index")
plt.ylabel("PPI Jaccard index")
stats.spearmanr(xy[x], xy[y])


# In[53]:


x = "pdi_jaccard"
y = "ppi_jaccard"
print(iso_pairs["ppi_n_min"] >= 1)
print(iso_pairs[x].notnull())


# In[54]:


xy = iso_pairs.loc[
    iso_pairs[x].notnull()
    & iso_pairs[y].notnull()
    & (iso_pairs["ppi_n_min"] >= 1)
    & (iso_pairs["pdi_n_min"] >= 1),
    :,
]
print(xy.shape)
plt.scatter(xy[x], xy[y], alpha=0.3)
plt.xlabel("PDI Jaccard index")
plt.ylabel("PPI Jaccard index")
stats.spearmanr(xy[x], xy[y])


# In[55]:


# Look at direction of activation change? I.e. is the isoform with less binding partners
# the one with activation closer to 0?


# In[56]:


iso_pairs.sort_values("ppi_n_min_diff", ascending=False).head(20)


# In[57]:


# Check if this result is robust against requiring at least
# one interaction partner for both genes.
# To see if the effect is driven by non-functional isoforms.
x = "ppi_jaccard"
y = "activation_fold_change"
xy = iso_pairs.loc[
    iso_pairs[x].notnull() & iso_pairs[y].notnull() & (iso_pairs["ppi_n_min"] >= 1), :
]
print(xy.shape)
stats.spearmanr(xy[x], xy[y])


# In[58]:


xy.sort_values("activation_fold_change", ascending=False).head()


# In[59]:


# overlaps

a = set(iso_pairs.loc[iso_pairs["ppi_jaccard"].notnull(), "tf_gene_symbol"].unique())
b = set(iso_pairs.loc[iso_pairs["pdi_jaccard"].notnull(), "tf_gene_symbol"].unique())
c = set(
    iso_pairs.loc[
        iso_pairs["activation_fold_change"].notnull(), "tf_gene_symbol"
    ].unique()
)
from matplotlib_venn import venn3

venn3([a, b, c], set_labels=["Y2H PPI", "Y1H PDI", "M1H activation"])
plt.savefig("../figures/tf_gene_data_integration_venn.pdf")


# In[61]:


ppi["Y2H_result"].value_counts()


# In[62]:


ppi["Y2H_result"].notnull().sum() / ppi.shape[0]


# In[63]:


def remake_yang_et_al_fig_2b(values, color="navy"):
    ax.plot(
        range(1, values.shape[0] + 1),
        sorted(values),
        marker="D",
        markersize=3,
        color=color,
        clip_on=False,
    )
    ax.set_xlabel("Isoform pairs", fontsize=14)
    ax.set_ylabel("Interaction profile dissimilarity\n(Jaccard distance)", fontsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels(["{:g}".format(t) for t in ax.get_yticks()], fontsize=12)
    ax.set_ylim(ax.get_ylim()[0], 1)


fig, ax = plt.subplots(1, 1)
values = iso_pairs.loc[
    iso_pairs["ppi_jaccard"].notnull() & (iso_pairs["ppi_n_min"] >= 1), "ppi_jaccard"
].values
values = 1.0 - values
remake_yang_et_al_fig_2b(values)
plt.savefig("../figures/remake_yang_et_al_fig4B.pdf", bbox_inches="tight")


# In[64]:


fig, ax = plt.subplots(1, 1)
values = iso_pairs.loc[
    iso_pairs["pdi_jaccard"].notnull() & (iso_pairs["pdi_n_min"] >= 1), "pdi_jaccard"
].values
values = 1.0 - values
remake_yang_et_al_fig_2b(values)
plt.savefig("../figures/remake_yang_et_al_fig4B_PDI.pdf", bbox_inches="tight")


# In[65]:


# combine PPI and PDI
iso_pairs["ppi_pdi_jaccard"] = (
    iso_pairs["ppi_n_shared"] + iso_pairs["pdi_n_shared"]
) / (iso_pairs["ppi_n_tested"] + iso_pairs["pdi_n_tested"])
fig, ax = plt.subplots(1, 1)
values = iso_pairs.loc[
    iso_pairs["ppi_pdi_jaccard"].notnull()
    & ((iso_pairs["ppi_n_min"] >= 1) | (iso_pairs["pdi_n_min"] >= 1)),
    "ppi_pdi_jaccard",
].values
values = 1.0 - values
remake_yang_et_al_fig_2b(values)
plt.savefig(
    "../figures/remake_yang_et_al_fig4B_combine_PDI_PPI.pdf", bbox_inches="tight"
)


# In[66]:


# all on one plot
fig, ax = plt.subplots(1, 1)
values_ppi = (
    1
    - iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) & (iso_pairs["pdi_n_min"] >= 1)),
        "ppi_jaccard",
    ].values
)
values_pdi = (
    1
    - iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) & (iso_pairs["pdi_n_min"] >= 1)),
        "pdi_jaccard",
    ].values
)
values_ppi_pdi = (
    1
    - iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) & (iso_pairs["pdi_n_min"] >= 1)),
        "ppi_pdi_jaccard",
    ].values
)
remake_yang_et_al_fig_2b(values_ppi, color="C1")
remake_yang_et_al_fig_2b(values_pdi, color="C2")
remake_yang_et_al_fig_2b(values_ppi_pdi, color="C0")
ax.text(ax.get_xlim()[1], values_ppi.max(), "PPI", va="center", fontsize=14, color="C1")
ax.text(ax.get_xlim()[1], values_pdi.max(), "PDI", va="center", fontsize=14, color="C2")
ax.text(
    ax.get_xlim()[1],
    values_ppi_pdi.max(),
    "PPI+PDI",
    va="center",
    fontsize=14,
    color="C0",
)
plt.savefig("../figures/remake_yang_et_al_fig4B_all_one_plot.pdf", bbox_inches="tight")


# In[67]:


fig, ax = plt.subplots(1, 1)
remake_yang_et_al_fig_2b(values_pdi, color="C2")
ax.text(ax.get_xlim()[1], values_pdi.max(), "PDI", va="center", fontsize=14, color="C2")
plt.savefig("../figures/remake_yang_et_al_fig4B_PDI.pdf", bbox_inches="tight")


# In[68]:


# TEMP: make with (2020-07-30)
ys = (
    1
    - iso_pairs.loc[
        (iso_pairs["ppi_n_tested"] >= 2) & iso_pairs["both_iso_y2h_pos"], "ppi_jaccard"
    ].values
)


fig, ax = plt.subplots(1, 1)
remake_yang_et_al_fig_2b(ys, color="C0")
ax.text(ax.get_xlim()[1], ys.max(), "PPI", va="center", fontsize=14, color="C0")
plt.savefig("../figures/ppi_jaccard_20200730.pdf", bbox_inches="tight")


# In[69]:


fig, ax = plt.subplots(1, 1)
remake_yang_et_al_fig_2b(values_ppi, color="C1")
ax.text(ax.get_xlim()[1], values_ppi.max(), "PPI", va="center", fontsize=14, color="C1")
plt.savefig("../figures/remake_yang_et_al_fig4B_PPI.pdf", bbox_inches="tight")


# In[70]:


fig, ax = plt.subplots(1, 1)
remake_yang_et_al_fig_2b(values_ppi_pdi, color="C0")
ax.text(
    ax.get_xlim()[1],
    values_ppi_pdi.max(),
    "PPI+PDI",
    va="center",
    fontsize=14,
    color="C0",
)
plt.savefig("../figures/remake_yang_et_al_fig4B_both.pdf", bbox_inches="tight")


# In[71]:


# product
fig, ax = plt.subplots(1, 1)
values_ppi = (
    1
    - iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) & (iso_pairs["pdi_n_min"] >= 1)),
        "ppi_jaccard",
    ].values
)
values_pdi = (
    1
    - iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) & (iso_pairs["pdi_n_min"] >= 1)),
        "pdi_jaccard",
    ].values
)
values_ppi_pdi = 1 - (
    iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) & (iso_pairs["pdi_n_min"] >= 1)),
        "ppi_jaccard",
    ].values
    * iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) & (iso_pairs["pdi_n_min"] >= 1)),
        "pdi_jaccard",
    ].values
)
remake_yang_et_al_fig_2b(values_ppi, color="C1")
remake_yang_et_al_fig_2b(values_pdi, color="C2")
remake_yang_et_al_fig_2b(values_ppi_pdi, color="C0")
ax.text(ax.get_xlim()[1], values_ppi.max(), "PPI", va="center", fontsize=14, color="C1")
ax.text(
    ax.get_xlim()[1],
    values_pdi.max() - 0.03,
    "PDI",
    va="center",
    fontsize=14,
    color="C2",
)
ax.text(
    ax.get_xlim()[1],
    values_ppi_pdi.max() + 0.03,
    "PPI+PDI",
    va="center",
    fontsize=14,
    color="C0",
)
plt.savefig(
    "../figures/remake_yang_et_al_fig4B_all_one_plot_product.pdf", bbox_inches="tight"
)


# In[72]:


fig, ax = plt.subplots(1, 1)
remake_yang_et_al_fig_2b(values_ppi_pdi, color="C0")
ax.text(
    ax.get_xlim()[1],
    values_ppi_pdi.max(),
    "PPI+PDI",
    va="center",
    fontsize=14,
    color="C0",
)
plt.savefig(
    "../figures/remake_yang_et_al_fig4B_pdi_ppi_product.pdf", bbox_inches="tight"
)


# In[73]:


values_ppi


# In[74]:


values_ppi_pdi


# In[75]:


plt.hist(values_ppi_pdi - values_ppi)


# In[76]:


stats.ttest_rel(values_ppi, values_ppi_pdi)


# In[77]:


stats.wilcoxon(values_ppi, values_ppi_pdi)


# In[78]:


stats.ks_2samp(values_ppi, values_ppi_pdi)


# In[79]:


def remake_yang_et_al_fig_2b_pct_x(values, color="navy"):
    ax.plot(
        np.linspace(0, 100, values.shape[0] + 1)[1:],
        sorted(values),
        marker="D",
        markersize=3,
        color=color,
        clip_on=False,
    )
    ax.set_xlabel("Percentage Isoform pairs", fontsize=14)
    ax.set_ylabel("Interaction profile dissimilarity\n(Jaccard distance)", fontsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels(["{:g}".format(t) for t in ax.get_yticks()], fontsize=12)
    ax.set_ylim(ax.get_ylim()[0], 1)


fig, ax = plt.subplots(1, 1)
values_ppi = (
    1
    - iso_pairs.loc[
        iso_pairs["ppi_jaccard"].notnull() & ((iso_pairs["ppi_n_min"] >= 1)),
        "ppi_jaccard",
    ].values
)
values_pdi = (
    1
    - iso_pairs.loc[
        iso_pairs["pdi_jaccard"].notnull() & ((iso_pairs["pdi_n_min"] >= 1)),
        "pdi_jaccard",
    ].values
)
values_ppi_pdi = (
    1
    - iso_pairs.loc[
        iso_pairs["ppi_pdi_jaccard"].notnull()
        & ((iso_pairs["ppi_n_min"] >= 1) | (iso_pairs["pdi_n_min"] >= 1)),
        "ppi_pdi_jaccard",
    ].values
)
remake_yang_et_al_fig_2b_pct_x(values_ppi, color="C1")
remake_yang_et_al_fig_2b_pct_x(values_pdi, color="C2")
remake_yang_et_al_fig_2b_pct_x(values_ppi_pdi, color="C0")
ax.text(ax.get_xlim()[1], 0.93, "PPI", va="center", fontsize=14, color="C1")
ax.text(ax.get_xlim()[1], 0.86, "PDI", va="center", fontsize=14, color="C2")
ax.text(ax.get_xlim()[1], 1, "PPI+PDI", va="center", fontsize=14, color="C0")
plt.savefig(
    "../figures/remake_yang_et_al_fig4B_all_one_plot_pct_isoforms.pdf",
    bbox_inches="tight",
)


# In[80]:


stats.ks_2samp(values_ppi, values_ppi_pdi)


# In[81]:


stats.ks_2samp(values_pdi, values_ppi_pdi)


# In[82]:


iso_pairs.loc[(iso_pairs["pdi_jaccard"] == 0) & (iso_pairs["pdi_n_min"] >= 1), :]
