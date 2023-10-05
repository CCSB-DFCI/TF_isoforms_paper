# coding: utf-8

# ## TODO
#
# - compare different families
#     - by how much splicing difference across development
#     - metric to compare pie charts
#         - for the highest expressed isoform, what proportion does it compose, filter for > 1 TPM at gene level, take the IQR?
# - make better plots
# - identify examples?
# - take a look by hand at the isoforms that don't show up in GTEx or dev-tissue
#     - could it be a short read mapping problem
# - look at the ORFcaptureseq long read data

# In[1]:


from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys

# import utils
sys.path.append("../")

from data_loading import (
    load_valid_isoform_clones,
    load_annotated_TFiso1_collection,
    load_developmental_tissue_expression,
)

pd.set_option("display.max_columns", 100)


# In[2]:


df, metadata, genes = load_developmental_tissue_expression()


# In[3]:


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
metadata["dev_stage"] = metadata["Developmental_Stage"].map(rename_dev_stage)
means = df.groupby(
    df.columns.map(metadata["organism_part"] + " " + metadata["dev_stage"]), axis=1
).mean()


# In[7]:


means = df.groupby(
    [df.columns.map(metadata["organism_part"]), df.columns.map(metadata["dev_stage"])],
    axis=1,
).mean()


# In[165]:


rfx4["dev_stage"].unique()


# In[172]:


# RFX4 example


dev_order = [
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "16",
    "18",
    "19",
    "neonate",
    "infant",
    "toddler",
    "child",
    "young adult",
    "adolescent" "adult",
    "elderly",
]

rfx4 = df.loc[df.index.str.startswith("RFX4|")].T.copy()
rfx4 = 2**rfx4 - 1
for column in rfx4.columns:
    rfx4[column] = rfx4[column] / rfx4.sum(axis=1)
rfx4["tissue"] = rfx4.index.map(metadata["organism_part"])
rfx4["dev_stage"] = rfx4.index.map(metadata["dev_stage"])
rfx4 = rfx4.loc[rfx4["tissue"].isin(["hindbrain", "forebrain"]), :]
fig, axs = plt.subplots(1, 2)
fig.set_size_inches(h=3, w=12)

axs[0].set_title("Forebrain")
sns.swarmplot(
    data=rfx4.loc[rfx4["tissue"] == "forebrain", :],
    y="RFX4|1/3|10C04 RFX4-203",
    x="dev_stage",
    order=[
        x
        for x in dev_order
        if x in rfx4.loc[rfx4["tissue"] == "forebrain", "dev_stage"].unique()
    ],
    ax=axs[0],
)
sns.pointplot(
    data=rfx4.loc[rfx4["tissue"] == "forebrain", :],
    y="RFX4|1/3|10C04 RFX4-203",
    x="dev_stage",
    order=[
        x
        for x in dev_order
        if x in rfx4.loc[rfx4["tissue"] == "forebrain", "dev_stage"].unique()
    ],
    ax=axs[0],
)
axs[1].set_title("Hindbrain")
sns.swarmplot(
    data=rfx4.loc[rfx4["tissue"] == "hindbrain", :],
    y="RFX4|1/3|10C04 RFX4-203",
    x="dev_stage",
    order=[
        x
        for x in dev_order
        if x in rfx4.loc[rfx4["tissue"] == "hindbrain", "dev_stage"].unique()
    ],
    ax=axs[1],
)
sns.pointplot(
    data=rfx4.loc[rfx4["tissue"] == "hindbrain", :],
    y="RFX4|1/3|10C04 RFX4-203",
    x="dev_stage",
    order=[
        x
        for x in dev_order
        if x in rfx4.loc[rfx4["tissue"] == "hindbrain", "dev_stage"].unique()
    ],
    ax=axs[1],
)
for ax in axs:
    ax.xaxis.set_ticklabels(ax.get_xticklabels(), rotation=90)
plt.savefig("../figures/RFX4_brain_dev_isoform_fraction.pdf", bbox_inches="tight")


# In[128]:


from scipy import stats

# take the max IQR per tissue
most_expressed_isoform_per_gene = means.mean(axis=1).groupby(genes).idxmax()
raw_means = 2**means - 1
f = raw_means / raw_means.groupby(genes).transform("sum")
f = f.loc[f.index.isin(most_expressed_isoform_per_gene.values), :]
f.index = f.index.map({v: k for k, v in most_expressed_isoform_per_gene.iteritems()})
f = f * (means.groupby(genes).sum() >= 1).applymap(lambda x: 1 if x else np.nan)
iqr = (
    f.T.reset_index()
    .groupby("level_0")
    .agg(lambda x: stats.iqr(x, nan_policy="omit"))
    .T
)  # .loc['forebrain', :])  # TEST
# .max())


# In[95]:


iqr.sort_values(ascending=False).head()


# In[129]:


# groupby tf family and plot
from data_loading import load_tf_families

tf_fam = load_tf_families()
# iqr = iqr.to_frame(name='iqr_isoform_ratio')
iqr["family"] = iqr.index.map(tf_fam)


# In[79]:


big_families = (
    iqr["family"].value_counts()[iqr["family"].value_counts() >= 5].index.values
)


# In[110]:


# In[131]:


import seaborn as sns

tissues = iqr.columns[:-1]
fig, axs = plt.subplots(3, 3, sharex=True)
fig.set_size_inches(12, 12)
for tissue, ax in zip(tissues, axs.flatten()):
    sns.stripplot(
        data=iqr.loc[iqr["family"].isin(big_families)], x="family", y=tissue, ax=ax
    )
    sns.boxplot(
        data=iqr.loc[iqr["family"].isin(big_families)], x="family", y=tissue, ax=ax
    )
    ax.set_title(tissue)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(0, 0.75)

axs[-1, 0].xaxis.set_ticklabels(ax.get_xticklabels(), rotation=90)
axs[1, 0].set_ylabel("IQR of primary isoform ratio across development")
axs[2, 1].axis("off")
axs[2, 2].axis("off")
plt.savefig(
    "../figures/IQR-isoform-ratio_across-development_by-tissue.pdf", bbox_inches="tight"
)


# In[5]:


def developmental_tissue_expression_plot(gene_name, means=means):
    fig, axes = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(48, 6)
    ### bar chart ###
    (means.loc[genes == gene_name, :].T.plot.bar(ax=axes[0], legend=False, width=0.7))
    ### percentages ###
    raw_means = 2 ** means.loc[genes == gene_name] - 1.0
    (
        raw_means.div(raw_means.sum(axis=0)).T.plot.bar(
            ax=axes[1], stacked=True, legend=False
        )
    )
    axes[0].set_ylabel("Mean log2(TPM + 1)")
    axes[1].set_ylabel("Fraction for each isoform")
    axes[1].set_yticklabels(["{:.0%}".format(t) for t in axes[1].get_yticks()])
    axes[1].legend(loc="lower left", bbox_to_anchor=(1, 0))
    axes[0].axhline(y=1, color="grey", linewidth=0.5)
    plt.subplots_adjust(hspace=0.05)
    fig_dir = Path(
        "/Users/lukelambourne/Dropbox (Partners HealthCare)/TF_isoforms/TF_gene_summary_pages/media"
    )
    plt.savefig(
        fig_dir / (gene_name + "_developmental-tissue-expression.svg"),
        bbox_inches="tight",
    )


# 1 try and make wider
# 2 try interactive...
developmental_tissue_expression_plot("ATF2")


# In[17]:


import altair as alt
import hvplot.pandas


def interactive_developmental_tissue_expression_plot(gene_name, means=means):
    """

    work in progress....
    - choose a backend: altair, plotly, bokeh/holoviews

    """
    # fig, axes = plt.subplots(2, 1, sharex=True)
    # fig.set_size_inches(48, 6)
    ### bar chart ###
    fig = means.loc[genes == gene_name, :].T.plot.bar(
        backend="plotly",
        barmode="group",
        width=3000,
        height=500,
        labels={"index": "", "value": "Mean log2(TPM + 1)"},
    )
    fig.update_layout(plot_bgcolor="white")
    fig.add_hline(y=1, line_width=0.5, line_color="grey")
    # save
    ### percentages ###
    """
    raw_means = 2 ** means.loc[genes == gene_name] - 1.
    (raw_means.div(raw_means.sum(axis=0))
              .T.plot.bar(ax=axes[1], 
                          stacked=True,
                          legend=False))
    axes[1].set_ylabel('Fraction for each isoform')
    axes[1].set_yticklabels(['{:.0%}'.format(t) for t in axes[1].get_yticks()])
    """
    fig_dir = Path(
        "/Users/lukelambourne/Dropbox (Partners HealthCare)/TF_isoforms/TF_gene_summary_pages/media"
    )
    # plt.savefig(fig_dir / (gene_name + '_developmental-tissue-expression.svg'),
    #            bbox_inches='tight')
    fig.write_html(
        fig_dir / (gene_name + "_developmental-tissue-expression_interactive.html")
    )
    return fig


interactive_developmental_tissue_expression_plot("ATF2")


# In[34]:


# NOTE: will take a long time
for tf in genes.unique():
    developmental_tissue_expression_plot(tf)
    plt.close(plt.gcf())


# In[4]:


# check isoform IDs relative to GTEx

from data_loading import load_gtex_remapped

gtex, metadata_gtex, genes_gtex = load_gtex_remapped()


# In[5]:


mean_gtex = gtex.groupby(gtex.columns.map(metadata_gtex["body_site"]), axis=1).mean()


# In[134]:


print(
    "{:.0%} isoforms ≥ 1 TPM in any GTEx tissue".format(
        (mean_gtex >= 1).any(axis=1).sum() / mean_gtex.shape[0]
    )
)
print(
    "{:.0%} isoforms ≥ 1 TPM in any developmental tissue".format(
        (means >= 1).any(axis=1).sum() / means.shape[0]
    )
)

print(
    "{:.0%} cloned isoforms ≥ 1 TPM in any GTEx tissue".format(
        (mean_gtex.loc[~mean_gtex.index.str.startswith("noclone"), :] >= 1)
        .any(axis=1)
        .sum()
        / (~mean_gtex.index.str.startswith("noclone")).sum()
    )
)
print(
    "{:.0%} cloned isoforms ≥ 1 TPM in any developmental tissue".format(
        (means.loc[~means.index.str.startswith("noclone"), :] >= 1).any(axis=1).sum()
        / (~means.index.str.startswith("noclone")).sum()
    )
)

print(
    "{:.0%} novel isoforms ≥ 1 TPM in any GTEx tissue".format(
        (mean_gtex.loc[mean_gtex.index.str.endswith("nomatch"), :] >= 1)
        .any(axis=1)
        .sum()
        / (mean_gtex.index.str.endswith("nomatch")).sum()
    )
)
print(
    "{:.0%} novel isoforms ≥ 1 TPM in any developmental tissue".format(
        (means.loc[means.index.str.endswith("nomatch"), :] >= 1).any(axis=1).sum()
        / (means.index.str.endswith("nomatch")).sum()
    )
)


# In[135]:


(~means.index.str.startswith("noclone")).sum()


# In[31]:


means.index.str.endswith("nomatch").sum()


# In[136]:


not_in_gtex = set(
    mean_gtex.loc[
        (~mean_gtex.index.str.startswith("noclone")) & (mean_gtex < 1).all(axis=1), :
    ].index.values
)
in_dev = set(
    means.loc[
        (~means.index.str.startswith("noclone")) & (means >= 1).any(axis=1), :
    ].index.values
)
print(len(not_in_gtex))
print(len(in_dev.intersection(not_in_gtex)))


# In[137]:


not_in_gtex.difference(in_dev)


# In[24]:


in_dev.intersection(not_in_gtex)


# In[19]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=24, h=4)
(
    (
        means.loc[
            means.index.isin(not_in_gtex) & ~means.index.str.startswith("noclone"), :
        ]
        >= 1
    )
    .sum()
    .plot.bar(ax=ax)
)
ax.set_ylabel("Number of clones, not seen in GTEx, with TPM ≥ 1")
plt.savefig(
    "../figures/n-clones_not-in-GTEx_in-dev_by-tissue-and-dev-stage.pdf",
    bbox_inches="tight",
)


# In[20]:


# get the biggest effect (amount of expression and fraction of isoform)
means.loc[
    means.index.isin(not_in_gtex) & ~means.index.str.startswith("noclone"), :
].idxmax().unique()


# ## Examples
#
# - MAZ - weird that the two cloned isoforms are not there in GTEx but are there in all samples in this data
# - ATF3 - no particular pattern, seems weird again
# - NFIA - also weird, high everywhere, not there in GTEx
# - FOS - also weird
# - FOSB - only in toddler testis
# - ZBTB18 - good example of weird, two isoforms with only one showing up in GTEx
# - ARNT2 - totally crazy
# - ZNF207 - interesting but more tissue specific splicing than development specific

# In[21]:


# look for big change in percent of isoform + high overall expression levels
raw_means = 2**means - 1.0
pct = (raw_means / raw_means.groupby(genes).transform("sum")) * 100
pct.head()


# In[22]:


# goes from less than 10% to more than 50% and has TPM > 1 in both?
tissue = "liver"
cols = [c for c in pct.columns if c.startswith(tissue)]
big_change = ((pct.loc[:, cols] < 10) & (means.loc[:, cols] >= 2)).any(axis=1) & (
    (pct.loc[:, cols] > 50) & (means.loc[:, cols] >= 2)
).any(axis=1)
print(tissue, big_change.sum())
print(genes[big_change[big_change].index].unique())


# In[73]:


big_change[big_change].index


# In[72]:


for tf in genes.unique():
    interactive_developmental_tissue_expression_plot(tf)


# - NR2F2 looks interesting
# - RFX4 seems to have brain specific developmental splicing
# - HHEX in liver

# In[51]:


interactive_developmental_tissue_expression_plot("ZBTB18")


# In[105]:


# improve the plots
# match the IDs
# compare to GTEx for adult tissues
# compare variation across families (nuclear receptors have more varation across development?)
pd.set_option("display.max_columns", 200)
means.sort_index().head()


# In[100]:


mean_gtex.columns


# In[31]:


# for conditions with at least 1 TPM gene level expression
# what's the fraction of each isoform
# take the max
# remove the most expressed isoform
fig, axs = plt.subplots(1, 2, sharey=True)
fig.set_size_inches(w=6, h=2)
most_expressed_isoform_per_gene = means.mean(axis=1).groupby(genes).idxmax()
raw_means = 2**means - 1
f = raw_means / raw_means.groupby(genes).transform("sum")
f = f * (means.groupby(genes).transform("sum") >= 1).applymap(
    lambda x: 1 if x else np.nan
)
f = f.loc[~f.index.isin(most_expressed_isoform_per_gene.values), :]
f = f.max(axis=1)
f.plot.hist(bins=50, range=(0, 1), ax=axs[1])
axs[0].set_xlabel("Fraction of gene expression\nmaximum across conditions")
axs[1].set_xlabel("Fraction of gene expression\nmaximum across conditions")
axs[0].set_ylabel("Alternative isoforms")
axs[1].set_ylabel("")
axs[1].set_title("Developmental dataset")


most_expressed_isoform_per_gene = mean_gtex.mean(axis=1).groupby(genes_gtex).idxmax()
raw_means = 2**mean_gtex - 1
f = raw_means / raw_means.groupby(genes_gtex).transform("sum")
f = f * (mean_gtex.groupby(genes_gtex).transform("sum") >= 1).applymap(
    lambda x: 1 if x else np.nan
)
f = f.loc[~f.index.isin(most_expressed_isoform_per_gene.values), :]
f = f.max(axis=1)
f.plot.hist(bins=50, range=(0, 1), ax=axs[0])
axs[0].set_title("GTEx")
fig.savefig(
    "../figures/fraction-gene-expression_alternative-isoforms.pdf", bbox_inches="tight"
)
# take max across either


# In[23]:


f.sort_values(ascending=False).head()


# In[23]:


# split by cloned, novel, not cloned, sum genes
paired_tissues = [
    ("liver adult", "Liver"),
    ("heart young adult", "Heart - Atrial Appendage"),
    ("testis adult", "Testis"),
]
for cm_tissue, gtex_tissue in paired_tissues:
    paired = pd.merge(
        means.loc[:, [cm_tissue]],
        mean_gtex.loc[:, [gtex_tissue]],
        how="inner",
        left_index=True,
        right_index=True,
    )
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(4, 4)
    paired.plot.scatter(x=gtex_tissue, y=cm_tissue, alpha=0.2, ax=ax)
    upper = paired.max().max()
    ax.set_ylim(0, upper + 0.5)
    ax.set_xlim(0, upper + 0.5)
    ax.set_xlabel("GTEx – Mean log2(TPM + 1)")
    ax.set_ylabel("Developmental dataset – Mean log2(TPM + 1)")
    r = paired.corr().loc[cm_tissue, gtex_tissue]
    ax.set_title("{} – R^2 = {:.2f}".format(gtex_tissue, r**2))
    plt.savefig(
        "../figures/GTEx-vs-dev_{}.pdf".format(gtex_tissue), bbox_inches="tight"
    )


# In[24]:


# split by cloned, novel, not cloned, sum genes
paired_tissues = [
    ("liver adult", "Liver"),
    ("heart young adult", "Heart - Atrial Appendage"),
    ("testis adult", "Testis"),
]
for cm_tissue, gtex_tissue in paired_tissues:
    fig, axs = plt.subplots(1, 3)
    fig.set_size_inches(w=12, h=4)
    paired = pd.merge(
        means.loc[:, [cm_tissue]],
        mean_gtex.loc[:, [gtex_tissue]],
        how="inner",
        left_index=True,
        right_index=True,
    )
    # cloned existing
    (
        paired.loc[
            ~paired.index.str.startswith("noclone")
            & ~paired.index.str.endswith("nomatch"),
            :,
        ].plot.scatter(x=gtex_tissue, y=cm_tissue, alpha=0.2, ax=axs[0])
    )
    r = (
        paired.loc[
            ~paired.index.str.startswith("noclone")
            & ~paired.index.str.endswith("nomatch"),
            :,
        ]
        .corr()
        .loc[cm_tissue, gtex_tissue]
    )
    axs[0].set_title(
        "{} – cloned GENCODE isoforms – R^2 = {:.2f}".format(
            gtex_tissue.split()[0], r**2
        ),
        fontsize=10,
    )
    # non-cloned
    (
        paired.loc[paired.index.str.startswith("noclone"), :].plot.scatter(
            x=gtex_tissue, y=cm_tissue, alpha=0.2, ax=axs[1]
        )
    )
    r = (
        paired.loc[paired.index.str.startswith("noclone"), :]
        .corr()
        .loc[cm_tissue, gtex_tissue]
    )
    axs[1].set_title(
        "{} – uncloned GENCODE isoforms – R^2 = {:.2f}".format(
            gtex_tissue.split()[0], r**2
        ),
        fontsize=10,
    )

    # novel
    (
        paired.loc[paired.index.str.endswith("nomatch"), :].plot.scatter(
            x=gtex_tissue, y=cm_tissue, alpha=0.2, ax=axs[2]
        )
    )
    r = (
        paired.loc[paired.index.str.endswith("nomatch"), :]
        .corr()
        .loc[cm_tissue, gtex_tissue]
    )
    axs[2].set_title(
        "{} – novel isoforms – R^2 = {:.2f}".format(gtex_tissue.split()[0], r**2),
        fontsize=10,
    )

    upper = paired.max().max()
    for ax in axs:
        ax.set_ylim(0, upper + 0.5)
        ax.set_xlim(0, upper + 0.5)
        ax.set_xlabel("GTEx – Mean log2(TPM + 1)")
        ax.set_ylabel("Developmental dataset – Mean log2(TPM + 1)")
    plt.savefig(
        "../figures/GTEx-vs-dev_{}_split.pdf".format(gtex_tissue), bbox_inches="tight"
    )


# In[25]:


raw_means_gtex = 2**mean_gtex - 1.0


# In[26]:


# log val 1, is 1
# log(a * b) = loga + logb
# log(a+b) = ???


# split by cloned, novel, not cloned, sum genes
paired_tissues = [
    ("liver adult", "Liver"),
    ("heart young adult", "Heart - Atrial Appendage"),
    ("testis adult", "Testis"),
]
for cm_tissue, gtex_tissue in paired_tissues:
    paired = pd.merge(
        (raw_means.groupby(genes).sum() + 1).apply(np.log2).loc[:, [cm_tissue]],
        (raw_means_gtex.groupby(genes).sum() + 1).apply(np.log2).loc[:, [gtex_tissue]],
        how="inner",
        left_index=True,
        right_index=True,
    )
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(4, 4)
    paired.plot.scatter(x=gtex_tissue, y=cm_tissue, alpha=0.2, ax=ax)
    upper = paired.max().max()
    ax.set_ylim(0, upper + 0.5)
    ax.set_xlim(0, upper + 0.5)
    ax.set_xlabel("GTEx – Mean log2(TPM + 1)")
    ax.set_ylabel("Developmental dataset – Mean log2(TPM + 1)")
    r = paired.corr().loc[cm_tissue, gtex_tissue]
    ax.set_title("{} – R^2 = {:.2f}".format(gtex_tissue, r**2))
    plt.savefig(
        "../figures/GTEx-vs-dev_{}_gene-level.pdf".format(gtex_tissue),
        bbox_inches="tight",
    )


# In[ ]:


# per gene metric to quantify varaition in splicing across development
# use to compare families
# can make a scatter with gene expression level metric to control for that...
