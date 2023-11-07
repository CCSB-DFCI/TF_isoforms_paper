# coding: utf-8

# In[1]:


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import pandas as pd
from scipy import stats
import seaborn as sns

from data_loading import (
    load_m1h_activation_data,
    load_annotated_TFiso1_collection,
    load_full_y2h_data_including_controls,
    load_valid_isoform_clones,
    load_ref_vs_alt_isoforms_table,
)


# In[2]:


m1h = load_m1h_activation_data()
tfs = load_annotated_TFiso1_collection()
m1h["mean"] = m1h[["M1H_rep1", "M1H_rep2", "M1H_rep3"]].mean(axis=1)

df = load_ref_vs_alt_isoforms_table()

# RORC-1 alt iso is causing an error - filter out here - there's no data for it?
df = df[df["clone_acc_alt"] != "RORC|1/6|05F11"]

df["ref_iso"] = df["clone_acc_ref"].apply(
    lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
)
df["alt_iso"] = df["clone_acc_alt"].apply(
    lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
)
df["f_disorder_difference"] = df.apply(
    lambda x: tfs[x["gene_symbol"]].disordered_fraction_of_different_regions(
        x["ref_iso"], x["alt_iso"]
    ),
    axis=1,
)


def disorder_changes_category(f):
    if pd.isnull(f):
        return np.nan
    elif f == 0:
        return "fully ordered"
    elif f == 1:
        return "fully disordered"
    elif f > 0 and f < 1:
        return "partially disordered"
    else:
        raise ValueError("Invalid fraction")


df["f_disorder_difference_cat"] = df["f_disorder_difference"].apply(
    disorder_changes_category
)
m1h["gte_2_fold"] = m1h["mean"].abs() >= 1
df["m1h_gte_2_fold_at_least_one_iso_per_gene"] = df["gene_symbol"].map(
    m1h.groupby("gene_symbol")["gte_2_fold"].any()
)
df["abs_activation_fold_change_log2"] = df["activation_fold_change_log2"].abs()

dom = pd.concat(
    [g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()]
)

# add activation or repression
effector_domain_type = {}
for tf in tfs.values():
    for d in tf.cloned_reference_isoform.aa_seq_features:
        if d.category == "effector_domain":
            effector_domain_type[d.accession] = d.name
dom["type"] = dom["accession"].map(effector_domain_type)


# considering EFFECTOR DOMAINS ONLY
def fraction_of_effector_domains_removed(row, effector_type):
    ds = dom.loc[(dom["alt_iso"] == row["alt_iso"]) & (dom["type"] == effector_type), :]
    if ds.shape[0] == 0:
        return np.nan
    return ds[["deletion", "frameshift"]].sum().sum() / ds["length"].sum()


def insertion_in_effector_domains(row, effector_type):
    ds = dom.loc[(dom["alt_iso"] == row["alt_iso"]) & (dom["type"] == effector_type), :]
    if ds.shape[0] == 0:
        return np.nan
    return ds["insertion"].sum()


def domain_length(row, effector_type):
    ds = dom.loc[(dom["alt_iso"] == row["alt_iso"]) & (dom["type"] == effector_type), :]
    if ds.shape[0] == 0:
        return np.nan
    return ds["length"].sum()


for effector_type in ["AD", "RD", "Bif"]:
    df["fraction_of_{}_domains_removed".format(effector_type)] = df.apply(
        fraction_of_effector_domains_removed, effector_type=effector_type, axis=1
    )
    df["insertion_in_{}_domains".format(effector_type)] = df.apply(
        insertion_in_effector_domains, effector_type=effector_type, axis=1
    )
    df["length_of_{}_domains".format(effector_type)] = df.apply(
        domain_length, effector_type=effector_type, axis=1
    )


# In[3]:


# now add Pfam AD/RDs
pfam = pd.read_csv(
    "../data/external/Pfam-A.clans.tsv",
    sep="\t",
    names=["pfam_accession", "clan", "clan_name", "short_name", "name"],
)
# AD
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


# get pfam type
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
    df["fraction_of_{}_domains_removed_incl_pfam".format(effector_type)] = df.apply(
        fraction_of_effector_domains_removed, effector_type=effector_type, axis=1
    )
    df["insertion_in_{}_domains_incl_pfam".format(effector_type)] = df.apply(
        insertion_in_effector_domains, effector_type=effector_type, axis=1
    )
    df["length_of_{}_domains_incl_pfam".format(effector_type)] = df.apply(
        domain_length, effector_type=effector_type, axis=1
    )


# In[4]:


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

fig, ax = plt.subplots(1, 1)
ax.hist(
    n_effector, range=(-0.25, max(n_effector) + 0.25), bins=(max(n_effector) * 2 + 1)
)
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel("Number of genes ({} total)".format(len(tfs)))
ax.set_xlabel("Effector domains in reference isoform")


# In[5]:


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

fig, ax = plt.subplots(1, 1)
ax.hist(
    n_effector_soto,
    range=(-0.25, max(n_effector) + 0.25),
    bins=(max(n_effector) * 2 + 1),
)
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel("Number of genes ({} total)".format(len(tfs)))
ax.set_xlabel("Effector domains in reference isoform")
ax.set_title("Soto et al. data")


# In[6]:


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

fig, ax = plt.subplots(1, 1)
ax.hist(
    n_effector_bintu,
    range=(-0.25, max(n_effector) + 0.25),
    bins=(max(n_effector) * 2 + 1),
)
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel("Number of genes ({} total)".format(len(tfs)))
ax.set_xlabel("Effector domains in reference isoform")
ax.set_title("Data from Bintu lab papers")


# In[7]:


# Summary stats
counter = {"Soto": 0, "Tycko": 0, "DelRosso": 0}
genes_with_effector = 0
genes_with_soto_effector = 0
for tf in tfs.values():
    has_effector = False
    has_soto_effector = False
    for domain in tf.cloned_reference_isoform.aa_seq_features:
        if domain.category == "effector_domain":
            has_effector = True
            if domain.accession.startswith("Soto"):
                has_soto_effector = True
            counter[domain.accession.split("_")[0]] += 1
    if has_effector:
        genes_with_effector += 1
    if has_soto_effector:
        genes_with_soto_effector += 1
print(counter)
print(genes_with_effector, genes_with_soto_effector, len(tfs))


# In[8]:


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


# In[9]:


fig, ax = plt.subplots(1, 1)
ax.bar(
    x=range(6),
    height=[counter[x][y] for x in counter.keys() for y in ["AD", "RD"]],
    color=["green", "red"] * 3,
)
ax.set_xticks([0.5, 2.5, 4.5])
ax.set_xticklabels(["Soto et al.", "Tycko et al.", "DelRosso et al."])
ax.set_ylabel("Number of effector domains")


# In[10]:


# length distribution of Soto vs Bintu lab domains
len_soto = []
len_bintu = []
for tf in tfs.values():
    for domain in tf.cloned_reference_isoform.aa_seq_features:
        if domain.category != "effector_domain":
            continue
        if domain.accession.startswith("Soto"):
            len_soto.append(len(dom))
        else:
            len_bintu.append(len(dom))

fig, axs = plt.subplots(ncols=1, nrows=2)
fig.set_size_inches(h=5, w=5)
axs[0].hist(len_soto, range=(1, 350), bins=35)
axs[1].hist(len_bintu, range=(1, 350), bins=35)
axs[1].set_xlabel("Size of effector domain (#aa)")
for ax in axs:
    ax.set_ylabel("Number of domains")
axs[0].set_title("Soto et al.")
axs[1].set_title("Tycko et al. + DelRosso et al.")


# In[11]:


# compare activation and repression


# In[12]:


# look at tf with 13
list(tfs.values())[n_effector.index(13)].protein_diagram()


# In[13]:


tfs["TP63"]["TP63-1"].aa_seq_features[5].accession


# In[14]:


print("Effector domain types:")
dom["type"].value_counts()


# In[15]:


print("Effector domain types including Pfam:")
dom["type_incl_pfam"].value_counts()


# In[16]:


print("Number of different types of effector domains per gene")
dom.groupby("gene_symbol")["type"].nunique().value_counts().sort_index()


# In[17]:


print("Number of different types of effector domains per gene including Pfam")
dom.groupby("gene_symbol")["type_incl_pfam"].nunique().value_counts()


# In[18]:


(
    df.loc[
        df["activation_fold_change_log2"].notnull(),
        [
            "fraction_of_AD_domains_removed",
            "fraction_of_RD_domains_removed",
            "fraction_of_Bif_domains_removed",
        ],
    ]
    .notnull()
    .groupby(
        [
            "fraction_of_AD_domains_removed",
            "fraction_of_RD_domains_removed",
            "fraction_of_Bif_domains_removed",
        ]
    )
    .size()
)


# In[19]:


# color map
t_disorder = df.loc[:, "f_disorder_difference"].values
norm = plt.Normalize(np.nanmin(t_disorder), np.nanmax(t_disorder))
cmap = sns.color_palette("flare", as_cmap=True)
palette_disorder = {value: cmap(norm(value)) for value in t_disorder}


def re_color(row, palette):
    if pd.isnull(row["f_disorder_difference"]):
        color = palette[0]
    else:
        color = palette[row["f_disorder_difference"]]
    return color


df["color_disorder"] = df.apply(re_color, axis=1, palette=palette_disorder)

# set na f_disorder_difference to 0 for now
df["f_disorder_difference"].fillna(0, inplace=True)


# In[20]:


# color map - length instead
# sum up lengths of all domains (plot only includes examples w 1 type of domain)
df["tot_dom_length"] = df[
    ["length_of_AD_domains", "length_of_RD_domains", "length_of_Bif_domains"]
].sum(axis=1)
t_dom_length = df.loc[:, "tot_dom_length"].values
t_dom_length = t_dom_length[t_dom_length > 0]

# using min and max makes colors too hard too read - cut off
norm = plt.Normalize(25, 250)
palette_dom_length = {value: cmap(norm(value)) for value in t_dom_length}


def re_color2(row, palette):
    if row["tot_dom_length"] == 0:
        color = sns.color_palette("flare")[0]
    else:
        color = palette[row["tot_dom_length"]]
    return color


df["color_dom_length"] = df.apply(re_color2, axis=1, palette=palette_dom_length)


# In[21]:


# color map - length instead
# sum up lengths of all domains (plot only includes examples w 1 type of domain)
df["tot_dom_length_incl_pfam"] = df[
    [
        "length_of_AD_domains_incl_pfam",
        "length_of_RD_domains_incl_pfam",
        "length_of_Bif_domains_incl_pfam",
    ]
].sum(axis=1)
t_dom_length_incl_pfam = df.loc[:, "tot_dom_length_incl_pfam"].values
t_dom_length_incl_pfam = t_dom_length_incl_pfam[t_dom_length_incl_pfam > 0]

# using min and max makes colors too hard too read - cut off
norm = plt.Normalize(25, 250)
palette_dom_length_incl_pfam = {
    value: cmap(norm(value)) for value in t_dom_length_incl_pfam
}


def re_color3(row, palette):
    if row["tot_dom_length_incl_pfam"] == 0:
        color = sns.color_palette("flare")[0]
    else:
        color = palette[row["tot_dom_length_incl_pfam"]]
    return color


df["color_dom_length_incl_pfam"] = df.apply(
    re_color3, axis=1, palette=palette_dom_length_incl_pfam
)
df.sample(5)


# In[22]:


df = df.loc[
    df["activation_fold_change_log2"].notnull()
    & df["m1h_gte_2_fold_at_least_one_iso_per_gene"],
    :,
]
palette = palette_disorder
hue = "f_disorder_difference"
color = "color_disorder"
t = t_disorder

gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])

fig, axs = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig.set_size_inches(w=8.2, h=2)

point_size = 6


tot_loss_activ = df.loc[
    (df["fraction_of_AD_domains_removed"] == 1)
    & (
        df["fraction_of_RD_domains_removed"].isnull()
        | (df["fraction_of_RD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed"].isnull()
        | (df["fraction_of_Bif_domains_removed"] == 0)
    ),
    :,
]
axs[0].set_title("activation\ndomain", fontsize=10)
sns.swarmplot(
    data=tot_loss_activ,
    y="activation_fold_change_log2",
    x="fraction_of_AD_domains_removed",
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

tot_loss_repr = df.loc[
    (df["fraction_of_RD_domains_removed"] == 1)
    & (
        df["fraction_of_AD_domains_removed"].isnull()
        | (df["fraction_of_AD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed"].isnull()
        | (df["fraction_of_Bif_domains_removed"] == 0)
    ),
    :,
]
axs[1].set_title("repression\ndomain", fontsize=10)
sns.swarmplot(
    data=tot_loss_repr,
    y="activation_fold_change_log2",
    x="fraction_of_RD_domains_removed",
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


tot_loss_both = df.loc[
    (df["fraction_of_AD_domains_removed"] == 1)
    & (df["fraction_of_RD_domains_removed"] == 1),
    :,
]
# axs[2].set_title('both activ. &\nrepr. domains',
# fontsize=10)
# sns.swarmplot(data=tot_loss_both,
#               y='activation_fold_change_log2',
#               x='fraction_of_RD_domains_removed',
#               size=point_size,
#             clip_on=False,
#               ax=axs[2],
#               palette=palette,
#               hue='f_disorder_difference',
#                linewidth=1,
#                edgecolor="black",
#               alpha=1)
# axs[2].set_xticks([])
# axs[2].set_xlabel('')
# axs[2].get_legend().remove()


# now partial loss
axs[2].set_title("activation\ndomain", fontsize=10)
partial_loss_activ = df.loc[
    (
        df["m1h_gte_2_fold_at_least_one_iso_per_gene"]
        & (df["fraction_of_AD_domains_removed"] > 0)
        & (df["fraction_of_AD_domains_removed"] < 1)
        & (
            df["fraction_of_RD_domains_removed"].isnull()
            | (df["fraction_of_RD_domains_removed"] == 0)
        )
        & (
            df["fraction_of_Bif_domains_removed"].isnull()
            | (df["fraction_of_Bif_domains_removed"] == 0)
        )
    ),
    :,
]
axs[2].scatter(
    partial_loss_activ.loc[:, "fraction_of_AD_domains_removed"].values,
    partial_loss_activ.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_activ.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[2].set_xlabel("")
axs[2].set_xlim(1, 0)
axs[2].set_xticks([0.99, 0.5, 0.01])
axs[2].set_xticklabels([f"{x:.0%}" for x in axs[2].get_xticks()])


axs[3].set_title("repression\ndomain", fontsize=10)
partial_loss_repr = df.loc[
    (
        df["m1h_gte_2_fold_at_least_one_iso_per_gene"]
        & (df["fraction_of_RD_domains_removed"] > 0)
        & (df["fraction_of_RD_domains_removed"] < 1)
        & (
            df["fraction_of_AD_domains_removed"].isnull()
            | (df["fraction_of_AD_domains_removed"] == 0)
        )
        & (
            df["fraction_of_Bif_domains_removed"].isnull()
            | (df["fraction_of_Bif_domains_removed"] == 0)
        )
    ),
    :,
]

axs[3].scatter(
    partial_loss_repr.loc[:, "fraction_of_RD_domains_removed"].values,
    partial_loss_repr.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_repr.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[3].set_xlabel("")
axs[3].set_xlim(1, 0)
axs[3].set_xticks([0.99, 0.5, 0.01])
axs[3].set_xticklabels([f"{x:.0%}" for x in axs[3].get_xticks()])


all_retained = df.loc[
    (
        (df["fraction_of_AD_domains_removed"] == 0)
        | (df["fraction_of_RD_domains_removed"] == 0)
        | (df["fraction_of_Bif_domains_removed"] == 0)
    )
    & (
        df["fraction_of_AD_domains_removed"].isnull()
        | (df["fraction_of_AD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_RD_domains_removed"].isnull()
        | (df["fraction_of_RD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed"].isnull()
        | (df["fraction_of_Bif_domains_removed"] == 0)
    ),
    :,
]
axs[4].set_title("All effector domains\nin alt. iso.", fontsize=10)
sns.swarmplot(
    data=all_retained,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[4],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    hue=hue,
    palette=palette,
)
axs[4].set_xticks([])
axs[4].set_xlabel("")
axs[4].get_legend().remove()

# annotate pbx1
pbx1_y = df.loc[
    (df["clone_acc_alt"] == "PBX1|2/2|02C05"), "activation_fold_change_log2"
].values[0]
for point in axs[4].collections:
    for x, y in point.get_offsets():
        if np.isclose(pbx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "PBX1-2",
                xy=(x, y),
                xytext=(1, -20),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.3"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )

# missing stuff
incl = (
    tot_loss_activ.append(tot_loss_repr)
    .append(tot_loss_both)
    .append(partial_loss_activ)
    .append(partial_loss_repr)
    .append(all_retained)
)
no_annot = df.loc[
    (~df.index.isin(incl.index.values))
    & (pd.isnull(df["fraction_of_AD_domains_removed"]))
    & (pd.isnull(df["fraction_of_RD_domains_removed"]))
    & (pd.isnull(df["fraction_of_Bif_domains_removed"]))
]
axs[5].set_title("No annotated\neffector domains", fontsize=10)
sns.swarmplot(
    data=no_annot,
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

# annotate RFX3-3
rfx3_y = df.loc[
    (df["clone_acc_alt"] == "RFX3|3/5|08G08"), "activation_fold_change_log2"
].values[0]
rfx4_y = df.loc[
    (df["clone_acc_alt"] == "RFX3|4/5|11D09"), "activation_fold_change_log2"
].values[0]
for point in axs[5].collections:
    for x, y in point.get_offsets():
        if np.isclose(rfx3_y, y):
            print("found: %s, %s" % (x, y))
            axs[5].annotate(
                "RFX3-3",
                xy=(x, y),
                xytext=(1, -20),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0.7"),
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
                xytext=(-5, -10),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.7"),
                ha="right",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )


# add colorbar
# mirror figure
gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])
fig2, axs2 = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig2.set_size_inches(w=8.2, h=2)
map1 = axs2[5].imshow(np.stack([t, t]), cmap="flare")
fig.colorbar(
    map1, ax=axs[5], aspect=60, label="% alt. iso. seq. diff.\nin disordered regions"
)


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylim(-7.5, 7.5)
    ax.axhline(y=0, color="black", linewidth=1, linestyle="dashed")
for ax in axs[1:]:
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(which="both", length=0)
    ax.set_ylabel("")
axs[0].set_ylabel("log2(activation fold change)")
fig.savefig(
    "../figures/activation_vs_domain_removal_colored_by_disorder.pdf",
    bbox_inches="tight",
)


# In[23]:


df = df.loc[
    df["activation_fold_change_log2"].notnull()
    & df["m1h_gte_2_fold_at_least_one_iso_per_gene"],
    :,
]
palette = palette_disorder
hue = "f_disorder_difference"
color = "color_disorder"
t = t_disorder

gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])

fig, axs = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
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
axs[0].set_title("activation\ndomain", fontsize=10)
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
axs[1].set_title("repression\ndomain", fontsize=10)
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


tot_loss_both = df.loc[
    (df["fraction_of_AD_domains_removed_incl_pfam"] == 1)
    & (df["fraction_of_RD_domains_removed_incl_pfam"] == 1),
    :,
]
# axs[2].set_title('both activ. &\nrepr. domains',
# fontsize=10)
# sns.swarmplot(data=tot_loss_both,
#               y='activation_fold_change_log2',
#               x='fraction_of_RD_domains_removed',
#               size=point_size,
#             clip_on=False,
#               ax=axs[2],
#               palette=palette,
#               hue='f_disorder_difference',
#                linewidth=1,
#                edgecolor="black",
#               alpha=1)
# axs[2].set_xticks([])
# axs[2].set_xlabel('')
# axs[2].get_legend().remove()


# now partial loss
axs[2].set_title("activation\ndomain", fontsize=10)
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
axs[2].scatter(
    partial_loss_activ.loc[:, "fraction_of_AD_domains_removed_incl_pfam"].values,
    partial_loss_activ.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_activ.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[2].set_xlabel("")
axs[2].set_xlim(1, 0)
axs[2].set_xticks([0.99, 0.5, 0.01])
axs[2].set_xticklabels([f"{x:.0%}" for x in axs[2].get_xticks()])


axs[3].set_title("repression\ndomain", fontsize=10)
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

axs[3].scatter(
    partial_loss_repr.loc[:, "fraction_of_RD_domains_removed_incl_pfam"].values,
    partial_loss_repr.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_repr.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[3].set_xlabel("")
axs[3].set_xlim(1, 0)
axs[3].set_xticks([0.99, 0.5, 0.01])
axs[3].set_xticklabels([f"{x:.0%}" for x in axs[3].get_xticks()])


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
axs[4].set_title("All effector domains\nin alt. iso.", fontsize=10)
sns.swarmplot(
    data=all_retained,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[4],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    hue=hue,
    palette=palette,
)
axs[4].set_xticks([])
axs[4].set_xlabel("")
axs[4].get_legend().remove()

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
for point in axs[4].collections:
    for x, y in point.get_offsets():
        if np.isclose(pbx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "PBX1-2",
                xy=(x, y),
                xytext=(-10, -20),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.3"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx3_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "RFX3-3",
                xy=(x, y),
                xytext=(-8, -7),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.2"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx4_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "RFX3-4",
                xy=(x, y),
                xytext=(8, -8),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0.3"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )

# missing stuff
incl = (
    tot_loss_activ.append(tot_loss_repr)
    .append(tot_loss_both)
    .append(partial_loss_activ)
    .append(partial_loss_repr)
    .append(all_retained)
)
no_annot = df.loc[
    (~df.index.isin(incl.index.values))
    & (pd.isnull(df["fraction_of_AD_domains_removed_incl_pfam"]))
    & (pd.isnull(df["fraction_of_RD_domains_removed_incl_pfam"]))
    & (pd.isnull(df["fraction_of_Bif_domains_removed_incl_pfam"]))
]
axs[5].set_title("No annotated\neffector domains", fontsize=10)
sns.swarmplot(
    data=no_annot,
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


# add colorbar
# mirror figure
gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])
fig2, axs2 = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig2.set_size_inches(w=8.2, h=2)
map1 = axs2[5].imshow(np.stack([t, t]), cmap="flare")
fig.colorbar(
    map1, ax=axs[5], aspect=60, label="% alt. iso. seq. diff.\nin disordered regions"
)


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylim(-7.5, 7.5)
    ax.axhline(y=0, color="black", linewidth=1, linestyle="dashed")
for ax in axs[1:]:
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(which="both", length=0)
    ax.set_ylabel("")
axs[0].set_ylabel("log2(activation fold change)")
fig.savefig(
    "../figures/activation_vs_domain_removal_incl_pfam_colored_by_disorder.pdf",
    bbox_inches="tight",
)


# In[24]:


tot_loss_repr[
    [
        "gene_symbol",
        "clone_acc_ref",
        "clone_acc_alt",
        "f_disorder_difference",
        "color_disorder",
        "activation_fold_change_log2",
    ]
]


# In[25]:


df = df.loc[
    df["activation_fold_change_log2"].notnull()
    & df["m1h_gte_2_fold_at_least_one_iso_per_gene"],
    :,
]
palette = palette_dom_length
hue = "tot_dom_length"
color = "color_dom_length"
t = t_dom_length

gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])

fig, axs = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig.set_size_inches(w=8.2, h=2)

point_size = 6


tot_loss_activ = df.loc[
    (df["fraction_of_AD_domains_removed"] == 1)
    & (
        df["fraction_of_RD_domains_removed"].isnull()
        | (df["fraction_of_RD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed"].isnull()
        | (df["fraction_of_Bif_domains_removed"] == 0)
    ),
    :,
]
axs[0].set_title("activation\ndomain", fontsize=10)
sns.swarmplot(
    data=tot_loss_activ,
    y="activation_fold_change_log2",
    x="fraction_of_AD_domains_removed",
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

tot_loss_repr = df.loc[
    (df["fraction_of_RD_domains_removed"] == 1)
    & (
        df["fraction_of_AD_domains_removed"].isnull()
        | (df["fraction_of_AD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed"].isnull()
        | (df["fraction_of_Bif_domains_removed"] == 0)
    ),
    :,
]
axs[1].set_title("repression\ndomain", fontsize=10)
sns.swarmplot(
    data=tot_loss_repr,
    y="activation_fold_change_log2",
    x="fraction_of_RD_domains_removed",
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


tot_loss_both = df.loc[
    (df["fraction_of_AD_domains_removed"] == 1)
    & (df["fraction_of_RD_domains_removed"] == 1),
    :,
]
# axs[2].set_title('both activ. &\nrepr. domains',
# fontsize=10)
# sns.swarmplot(data=tot_loss_both,
#               y='activation_fold_change_log2',
#               x='fraction_of_RD_domains_removed',
#               size=point_size,
#             clip_on=False,
#               ax=axs[2],
#               palette=palette,
#               hue='f_disorder_difference',
#                linewidth=1,
#                edgecolor="black",
#               alpha=1)
# axs[2].set_xticks([])
# axs[2].set_xlabel('')
# axs[2].get_legend().remove()


# now partial loss
axs[2].set_title("activation\ndomain", fontsize=10)
partial_loss_activ = df.loc[
    (
        df["m1h_gte_2_fold_at_least_one_iso_per_gene"]
        & (df["fraction_of_AD_domains_removed"] > 0)
        & (df["fraction_of_AD_domains_removed"] < 1)
        & (
            df["fraction_of_RD_domains_removed"].isnull()
            | (df["fraction_of_RD_domains_removed"] == 0)
        )
        & (
            df["fraction_of_Bif_domains_removed"].isnull()
            | (df["fraction_of_Bif_domains_removed"] == 0)
        )
    ),
    :,
]
axs[2].scatter(
    partial_loss_activ.loc[:, "fraction_of_AD_domains_removed"].values,
    partial_loss_activ.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_activ.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[2].set_xlabel("")
axs[2].set_xlim(1, 0)
axs[2].set_xticks([0.99, 0.5, 0.01])
axs[2].set_xticklabels([f"{x:.0%}" for x in axs[2].get_xticks()])


axs[3].set_title("repression\ndomain", fontsize=10)
partial_loss_repr = df.loc[
    (
        df["m1h_gte_2_fold_at_least_one_iso_per_gene"]
        & (df["fraction_of_RD_domains_removed"] > 0)
        & (df["fraction_of_RD_domains_removed"] < 1)
        & (
            df["fraction_of_AD_domains_removed"].isnull()
            | (df["fraction_of_AD_domains_removed"] == 0)
        )
        & (
            df["fraction_of_Bif_domains_removed"].isnull()
            | (df["fraction_of_Bif_domains_removed"] == 0)
        )
    ),
    :,
]

axs[3].scatter(
    partial_loss_repr.loc[:, "fraction_of_RD_domains_removed"].values,
    partial_loss_repr.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_repr.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[3].set_xlabel("")
axs[3].set_xlim(1, 0)
axs[3].set_xticks([0.99, 0.5, 0.01])
axs[3].set_xticklabels([f"{x:.0%}" for x in axs[3].get_xticks()])


all_retained = df.loc[
    (
        (df["fraction_of_AD_domains_removed"] == 0)
        | (df["fraction_of_RD_domains_removed"] == 0)
        | (df["fraction_of_Bif_domains_removed"] == 0)
    )
    & (
        df["fraction_of_AD_domains_removed"].isnull()
        | (df["fraction_of_AD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_RD_domains_removed"].isnull()
        | (df["fraction_of_RD_domains_removed"] == 0)
    )
    & (
        df["fraction_of_Bif_domains_removed"].isnull()
        | (df["fraction_of_Bif_domains_removed"] == 0)
    ),
    :,
]
axs[4].set_title("All effector domains\nin alt. iso.", fontsize=10)
sns.swarmplot(
    data=all_retained,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[4],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    hue=hue,
    palette=palette,
)
axs[4].set_xticks([])
axs[4].set_xlabel("")
axs[4].get_legend().remove()

# annotate pbx1
pbx1_y = df.loc[
    (df["clone_acc_alt"] == "PBX1|2/2|02C05"), "activation_fold_change_log2"
].values[0]
for point in axs[4].collections:
    for x, y in point.get_offsets():
        if np.isclose(pbx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "PBX1-2",
                xy=(x, y),
                xytext=(1, -20),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.3"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )

# missing stuff
incl = (
    tot_loss_activ.append(tot_loss_repr)
    .append(tot_loss_both)
    .append(partial_loss_activ)
    .append(partial_loss_repr)
    .append(all_retained)
)
no_annot = df.loc[
    (~df.index.isin(incl.index.values))
    & (pd.isnull(df["fraction_of_AD_domains_removed"]))
    & (pd.isnull(df["fraction_of_RD_domains_removed"]))
    & (pd.isnull(df["fraction_of_Bif_domains_removed"]))
]
axs[5].set_title("No annotated\neffector domains", fontsize=10)
sns.swarmplot(
    data=no_annot,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[5],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    color=sns.color_palette("flare")[0],
)
axs[5].set_xticks([])
axs[5].set_xlabel("")

# annotate RFX3-3
rfx3_y = df.loc[
    (df["clone_acc_alt"] == "RFX3|3/5|08G08"), "activation_fold_change_log2"
].values[0]
rfx4_y = df.loc[
    (df["clone_acc_alt"] == "RFX3|4/5|11D09"), "activation_fold_change_log2"
].values[0]
for point in axs[5].collections:
    for x, y in point.get_offsets():
        if np.isclose(rfx3_y, y):
            print("found: %s, %s" % (x, y))
            axs[5].annotate(
                "RFX3-3",
                xy=(x, y),
                xytext=(1, -20),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0.7"),
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
                xytext=(-5, -10),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.7"),
                ha="right",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )

# add colorbar
# mirror figure
gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])
fig2, axs2 = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig2.set_size_inches(w=8.2, h=2)
map1 = axs2[4].imshow(np.stack([t, t]), cmap="flare", vmin=25, vmax=250)
cbar = fig.colorbar(map1, ax=axs[4], aspect=60, label="# AA in annotated domain")
cbar.set_ticks([25, 75, 150, 250])
cbar.set_ticklabels(["<=25", "75", "150", ">=250"])


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylim(-7.5, 7.5)
    ax.axhline(y=0, color="black", linewidth=1, linestyle="dashed")
for ax in axs[1:]:
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(which="both", length=0)
    ax.set_ylabel("")
axs[0].set_ylabel("log2(activation fold change)")
fig.savefig(
    "../figures/activation_vs_domain_removal_colored_by_dom_length.pdf",
    bbox_inches="tight",
)


# In[26]:


df = df.loc[
    df["activation_fold_change_log2"].notnull()
    & df["m1h_gte_2_fold_at_least_one_iso_per_gene"],
    :,
]
palette = palette_dom_length_incl_pfam
hue = "tot_dom_length_incl_pfam"
color = "color_dom_length_incl_pfam"
t = t_dom_length_incl_pfam

gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])

fig, axs = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
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
axs[0].set_title("activation\ndomain", fontsize=10)
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
axs[1].set_title("repression\ndomain", fontsize=10)
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


tot_loss_both = df.loc[
    (df["fraction_of_AD_domains_removed_incl_pfam"] == 1)
    & (df["fraction_of_RD_domains_removed_incl_pfam"] == 1),
    :,
]
# axs[2].set_title('both activ. &\nrepr. domains',
# fontsize=10)
# sns.swarmplot(data=tot_loss_both,
#               y='activation_fold_change_log2',
#               x='fraction_of_RD_domains_removed',
#               size=point_size,
#             clip_on=False,
#               ax=axs[2],
#               palette=palette,
#               hue='f_disorder_difference',
#                linewidth=1,
#                edgecolor="black",
#               alpha=1)
# axs[2].set_xticks([])
# axs[2].set_xlabel('')
# axs[2].get_legend().remove()


# now partial loss
axs[2].set_title("activation\ndomain", fontsize=10)
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
axs[2].scatter(
    partial_loss_activ.loc[:, "fraction_of_AD_domains_removed_incl_pfam"].values,
    partial_loss_activ.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_activ.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[2].set_xlabel("")
axs[2].set_xlim(1, 0)
axs[2].set_xticks([0.99, 0.5, 0.01])
axs[2].set_xticklabels([f"{x:.0%}" for x in axs[2].get_xticks()])


axs[3].set_title("repression\ndomain", fontsize=10)
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

axs[3].scatter(
    partial_loss_repr.loc[:, "fraction_of_RD_domains_removed_incl_pfam"].values,
    partial_loss_repr.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_repr.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[3].set_xlabel("")
axs[3].set_xlim(1, 0)
axs[3].set_xticks([0.99, 0.5, 0.01])
axs[3].set_xticklabels([f"{x:.0%}" for x in axs[3].get_xticks()])


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
axs[4].set_title("All effector domains\nin alt. iso.", fontsize=10)
sns.swarmplot(
    data=all_retained,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[4],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    hue=hue,
    palette=palette,
)
axs[4].set_xticks([])
axs[4].set_xlabel("")
axs[4].get_legend().remove()

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
for point in axs[4].collections:
    for x, y in point.get_offsets():
        if np.isclose(pbx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "PBX1-2",
                xy=(x, y),
                xytext=(-10, -20),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.3"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx3_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "RFX3-3",
                xy=(x, y),
                xytext=(-8, -7),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.2"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx4_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "RFX3-4",
                xy=(x, y),
                xytext=(8, -8),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0.3"),
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
axs[5].set_title("No annotated\neffector domains", fontsize=10)
sns.swarmplot(
    data=no_annot,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[5],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    color=sns.color_palette("flare")[0],
)
axs[5].set_xticks([])
axs[5].set_xlabel("")


# add colorbar
# mirror figure
gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])
fig2, axs2 = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig2.set_size_inches(w=8.2, h=2)
map1 = axs2[4].imshow(np.stack([t, t]), cmap="flare", vmin=25, vmax=250)
cbar = fig.colorbar(map1, ax=axs[4], aspect=60, label="# AA in annotated domain")
cbar.set_ticks([25, 75, 150, 250])
cbar.set_ticklabels(["<=25", "75", "150", ">=250"])


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylim(-7.5, 7.5)
    ax.axhline(y=0, color="black", linewidth=1, linestyle="dashed")
for ax in axs[1:]:
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(which="both", length=0)
    ax.set_ylabel("")
axs[0].set_ylabel("log2(activation fold change)")
fig.savefig(
    "../figures/activation_vs_domain_removal_incl_pfam_colored_by_dom_length.pdf",
    bbox_inches="tight",
)


# In[27]:


df = df.loc[
    df["activation_fold_change_log2"].notnull()
    & df["m1h_gte_2_fold_at_least_one_iso_per_gene"],
    :,
]
palette = palette_dom_length_incl_pfam
hue = "tot_dom_length_incl_pfam"
color = "color_dom_length_incl_pfam"
t = t_dom_length_incl_pfam

gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])

fig, axs = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig.set_size_inches(w=8.2, h=2)

point_size = 6


tot_loss_activ = df.loc[(df["fraction_of_AD_domains_removed_incl_pfam"] == 1), :]
axs[0].set_title("activation\ndomain", fontsize=10)
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
axs[1].set_title("repression\ndomain", fontsize=10)
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


tot_loss_both = df.loc[
    (df["fraction_of_AD_domains_removed_incl_pfam"] == 1)
    & (df["fraction_of_RD_domains_removed_incl_pfam"] == 1),
    :,
]
# axs[2].set_title('both activ. &\nrepr. domains',
# fontsize=10)
# sns.swarmplot(data=tot_loss_both,
#               y='activation_fold_change_log2',
#               x='fraction_of_RD_domains_removed',
#               size=point_size,
#             clip_on=False,
#               ax=axs[2],
#               palette=palette,
#               hue='f_disorder_difference',
#                linewidth=1,
#                edgecolor="black",
#               alpha=1)
# axs[2].set_xticks([])
# axs[2].set_xlabel('')
# axs[2].get_legend().remove()


# now partial loss
axs[2].set_title("activation\ndomain", fontsize=10)
partial_loss_activ = df.loc[
    (
        df["m1h_gte_2_fold_at_least_one_iso_per_gene"]
        & (df["fraction_of_AD_domains_removed_incl_pfam"] > 0)
        & (df["fraction_of_AD_domains_removed_incl_pfam"] < 1)
    ),
    :,
]
axs[2].scatter(
    partial_loss_activ.loc[:, "fraction_of_AD_domains_removed_incl_pfam"].values,
    partial_loss_activ.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_activ.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[2].set_xlabel("")
axs[2].set_xlim(1, 0)
axs[2].set_xticks([0.99, 0.5, 0.01])
axs[2].set_xticklabels([f"{x:.0%}" for x in axs[2].get_xticks()])


axs[3].set_title("repression\ndomain", fontsize=10)
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

axs[3].scatter(
    partial_loss_repr.loc[:, "fraction_of_RD_domains_removed_incl_pfam"].values,
    partial_loss_repr.loc[:, "activation_fold_change_log2"].values,
    alpha=1,
    s=point_size**2,
    c=partial_loss_repr.loc[:, color].values,
    linewidth=1,
    edgecolor="black",
    clip_on=False,
)
axs[3].set_xlabel("")
axs[3].set_xlim(1, 0)
axs[3].set_xticks([0.99, 0.5, 0.01])
axs[3].set_xticklabels([f"{x:.0%}" for x in axs[3].get_xticks()])


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
axs[4].set_title("All effector domains\nin alt. iso.", fontsize=10)
sns.swarmplot(
    data=all_retained,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[4],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    hue=hue,
    palette=palette,
)
axs[4].set_xticks([])
axs[4].set_xlabel("")
axs[4].get_legend().remove()

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
for point in axs[4].collections:
    for x, y in point.get_offsets():
        if np.isclose(pbx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "PBX1-2",
                xy=(x, y),
                xytext=(-10, -20),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.3"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx3_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "RFX3-3",
                xy=(x, y),
                xytext=(-8, -7),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.2"),
                ha="center",
                va="top",
                fontsize=7,
                bbox=dict(boxstyle="square,pad=0", fc="none", ec="none"),
            )
        if np.isclose(rfx4_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate(
                "RFX3-4",
                xy=(x, y),
                xytext=(8, -8),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0.3"),
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
axs[5].set_title("No annotated\neffector domains", fontsize=10)
sns.swarmplot(
    data=no_annot,
    y="activation_fold_change_log2",
    x="m1h_gte_2_fold_at_least_one_iso_per_gene",
    size=point_size,
    clip_on=False,
    ax=axs[5],
    linewidth=1,
    edgecolor="black",
    alpha=1,
    color=sns.color_palette("flare")[0],
)
axs[5].set_xticks([])
axs[5].set_xlabel("")


# add colorbar
# mirror figure
gs_kw = dict(width_ratios=[0.5, 0.5, 1.2, 1.2, 1.2, 1.6])
fig2, axs2 = plt.subplots(1, 6, sharey=True, gridspec_kw=gs_kw)
fig2.set_size_inches(w=8.2, h=2)
map1 = axs2[4].imshow(np.stack([t, t]), cmap="flare", vmin=25, vmax=250)
cbar = fig.colorbar(map1, ax=axs[4], aspect=60, label="# AA in annotated domain")
cbar.set_ticks([25, 75, 150, 250])
cbar.set_ticklabels(["<=25", "75", "150", ">=250"])


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylim(-7.5, 7.5)
    ax.axhline(y=0, color="black", linewidth=1, linestyle="dashed")
for ax in axs[1:]:
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(which="both", length=0)
    ax.set_ylabel("")
axs[0].set_ylabel("log2(activation fold change)")
fig.savefig(
    "../figures/activation_vs_domain_removal_incl_pfam_colored_by_dom_length_lossen-requirements.pdf",
    bbox_inches="tight",
)


# In[28]:


tfs["TBX5"].protein_diagram()


# In[29]:


df[df["gene_symbol"] == "PBX1"][
    [
        "gene_symbol",
        "clone_acc_ref",
        "clone_acc_alt",
        "fraction_of_AD_domains_removed",
        "fraction_of_RD_domains_removed",
        "fraction_of_Bif_domains_removed",
        "m1h_gte_2_fold_at_least_one_iso_per_gene",
        "activation_fold_change_log2",
    ]
]


# In[30]:


df.plot.scatter(x="f_disorder_difference", y="abs_activation_fold_change_log2")


# In[31]:


x_var = "f_disorder_difference"
y_var = "activation_fold_change_log2"
x, y = df.loc[df[x_var].notnull() & df[y_var].notnull(), [x_var, y_var]].values.T
print(stats.pearsonr(x, y))
print(stats.spearmanr(x, y))


# In[32]:


x_var = "f_disorder_difference"
y_var = "abs_activation_fold_change_log2"
x, y = df.loc[df[x_var].notnull() & df[y_var].notnull(), [x_var, y_var]].values.T
print(stats.pearsonr(x, y))
print(stats.spearmanr(x, y))


# In[33]:


sns.boxplot(data=df, x="f_disorder_difference_cat", y="activation_fold_change_log2")


# In[34]:


sns.boxplot(data=df, x="f_disorder_difference_cat", y="abs_activation_fold_change_log2")


# In[35]:


sns.boxplot(
    data=df.loc[(df["m1h_gte_2_fold_at_least_one_iso_per_gene"] == True), :],
    x="f_disorder_difference_cat",
    y="abs_activation_fold_change_log2",
)


# In[36]:


stats.mannwhitneyu(
    df.loc[
        (df["m1h_gte_2_fold_at_least_one_iso_per_gene"] == True)
        & (df["f_disorder_difference_cat"] == "fully ordered"),
        "abs_activation_fold_change_log2",
    ].values,
    df.loc[
        (df["m1h_gte_2_fold_at_least_one_iso_per_gene"] == True)
        & (df["f_disorder_difference_cat"] == "fully disordered"),
        "abs_activation_fold_change_log2",
    ].values,
)


# Interestingly, I don't see a difference between disordered and ordered changes and
# activation levels. Check for confounding factors like the size of the changes.

# In[37]:


sns.boxplot(
    data=df,
    x="f_disorder_difference_cat",
    y="aa_seq_pct_identity",
    order=["fully disordered", "fully ordered", "partially disordered"],
)


# In[38]:


df.plot.scatter(x="f_disorder_difference", y="abs_activation_fold_change_log2")


# In[39]:


# x axis as size of change (or aa %id), y axis as act diff, split by dis and ordered


# In[40]:


# check examples
(
    df.loc[df["f_disorder_difference_cat"] == "fully ordered", :].sort_values(
        "abs_activation_fold_change_log2", ascending=False
    )
).head()


# In[41]:


df["f_disorder_difference_cat"].value_counts()


# In[42]:


# TODO move PPI stuff to different notebook
sns.swarmplot(
    data=df.loc[(df["n_positive_PPI_ref"] > 0) & (df["n_positive_PPI_alt"] > 0)],
    x="f_disorder_difference_cat",
    y="PPI_jaccard",
)
sns.boxplot(
    data=df.loc[(df["n_positive_PPI_ref"] > 0) & (df["n_positive_PPI_alt"] > 0)],
    x="f_disorder_difference_cat",
    y="PPI_jaccard",
)


# In[43]:


x = df.loc[
    (df["n_positive_PPI_ref"] > 0)
    & (df["n_positive_PPI_alt"] > 0)
    & (df["f_disorder_difference_cat"] == "fully ordered"),
    "PPI_jaccard",
].values
y = df.loc[
    (df["n_positive_PPI_ref"] > 0)
    & (df["n_positive_PPI_alt"] > 0)
    & (df["f_disorder_difference_cat"] == "fully disordered"),
    "PPI_jaccard",
].values
stats.mannwhitneyu(x, y)


# In[44]:


(
    df.loc[
        (df["n_positive_PPI_ref"] > 0) & (df["n_positive_PPI_alt"] > 0), :
    ].plot.scatter(x="f_disorder_difference", y="PPI_jaccard")
)


# In[45]:


x_var = "f_disorder_difference"
y_var = "PPI_jaccard"
x, y = df.loc[
    (df["n_positive_PPI_ref"] > 0)
    & (df["n_positive_PPI_alt"] > 0)
    & df[x_var].notnull()
    & df[y_var].notnull(),
    [x_var, y_var],
].values.T
print(stats.pearsonr(x, y))
print(stats.spearmanr(x, y))


# In[46]:


x_var = "f_disorder_difference"
y_var = "PPI_jaccard"
x, y = df.loc[
    (df["n_positive_PPI_ref"] > 0)
    & (df["n_positive_PPI_alt"] > 0)
    & (df["f_disorder_difference"] > 0)
    & (df["f_disorder_difference"] < 1)
    & df[x_var].notnull()
    & df[y_var].notnull(),
    [x_var, y_var],
].values.T
print(stats.pearsonr(x, y))
print(stats.spearmanr(x, y))
