# coding: utf-8

# TODO
# - disorder from alphafold structures
# - protein sections matching PDI changes
# - count charges AA

# In[1]:


import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

from data_loading import *
from isoform_pairwise_metrics import *
from plotting import (
    y1h_pdi_per_tf_gene_plot,
    m1h_activation_per_tf_gene_plot,
    COLOR_PURPLE,
)
from data_loading import load_annotated_TFiso1_collection, load_y1h_pdi_data


# In[2]:


disorder = pd.read_csv(
    "../data/processed/TFiso1_disorder-and-ss_from-alphafold.tsv", sep="\t"
)
tfs = load_annotated_TFiso1_collection()

# TODO: fix missing data
clones_with_disorder_data = set(disorder["clone_name"].unique())
for tf in tfs.values():
    for iso in tf.cloned_isoforms:
        if iso.name not in clones_with_disorder_data:
            print("missing disorder data for {}".format(iso.name))
            continue
        iso.disorder = disorder.loc[
            disorder["clone_name"] == iso.name, "is_disordered"
        ].values

# sanity check
for tf in tfs.values():
    for iso in tf.cloned_isoforms:
        if hasattr(iso, "disorder"):
            if len(iso.disorder) != len(iso.aa_seq):
                raise UserWarning(
                    "inconsistent amino acid sequence and disordered residues data for {}".format(
                        iso.name
                    )
                )


def disordered_fraction_of_different_regions(gene, ref_iso_name, alt_iso_name):
    algn = gene.pairwise_changes_relative_to_reference(ref_iso_name, alt_iso_name)
    if not hasattr(gene[ref_iso_name], "disorder") or not hasattr(
        gene[alt_iso_name], "disorder"
    ):
        return np.nan
    ref_iter = iter(gene[ref_iso_name].disorder)
    alt_iter = iter(gene[alt_iso_name].disorder)
    merged_disorder = []
    for pos in algn:
        if pos == "I":
            merged_disorder.append(next(alt_iter))
        elif pos == "D":
            merged_disorder.append(next(ref_iter))
        else:
            merged_disorder.append(next(ref_iter))
            next(alt_iter)

    return np.mean(
        [
            is_disordered
            for pos, is_disordered in zip(algn, merged_disorder)
            if pos != "M"
        ]
    )


disordered_fraction_of_different_regions(tfs["CREB1"], "CREB1-2", "CREB1-1")


# In[3]:


# TODO move to isolib
def n_aa_change_from_feature(
    gene, ref_iso_name, alt_iso_name, domain_start, domain_end
):
    algn = gene.pairwise_changes_relative_to_reference(ref_iso_name, alt_iso_name)

    def _coords_transform_aa_seq_to_alignment(i, alignment):
        if i > len(alignment.replace("I", "")):
            raise ValueError("position is not in isoform AA sequence")
        aa_seq_indices = [
            "" if c == "I" else len(alignment[:j].replace("I", ""))
            for j, c in enumerate(alignment)
        ]
        return aa_seq_indices.index(i)

    start = _coords_transform_aa_seq_to_alignment(domain_start, algn)
    end = _coords_transform_aa_seq_to_alignment(domain_end - 1, algn) + 1

    if not all(x == "M" for x in algn[start:end]):
        return 0  # change is within the domain

    big_number = 9999999999999999999999999
    c_dist = big_number
    n_dist = big_number
    for i, l in enumerate(reversed(algn[:start])):
        if l != "M":
            c_dist = i + 1
            break
    for i, l in enumerate(algn[end:]):
        if l != "M":
            n_dist = i + 1
            break
    if c_dist == big_number and n_dist == big_number:
        raise UserWarning("problem calculating distance")
    return min([c_dist, n_dist])


def n_aa_to_all_features(self, ref_iso_name):
    results = []
    ref_iso = self._iso_dict[ref_iso_name]
    row = {"gene": self.name, "ref_iso": ref_iso_name}
    for aa_feature in ref_iso.aa_seq_features:
        for alt_iso_name, alt_iso in self._iso_dict.items():
            if alt_iso_name == ref_iso_name:
                continue
            row.update(
                {
                    "alt_iso": alt_iso_name,
                    "accession": aa_feature.accession,
                    "category": aa_feature.category,
                    "start_in_ref_iso": aa_feature.start,
                    "end_in_ref_iso": aa_feature.end,
                    "length": aa_feature.end - aa_feature.start,
                }
            )
            row.update(
                {
                    "n_aa_change_to_domain": n_aa_change_from_feature(
                        self,
                        ref_iso_name,
                        alt_iso_name,
                        aa_feature.start,
                        aa_feature.end,
                    )
                }
            )
            results.append(row.copy())
    results = pd.DataFrame(results)
    return results


gene = tfs["TBX5"]
ref_iso_name = "TBX5-1"
alt_iso_name = "TBX5-3"

# for dom in gene[ref_iso_name].aa_seq_features:
#    print(dom)
#    n = n_aa_change_from_feature(gene, ref_iso_name, alt_iso_name, dom.start, dom.end)
#    print(n)

# get for all domains
# get DBDs
# find min over DBDs

dist = pd.concat(
    [n_aa_to_all_features(g, g.cloned_reference_isoform.name) for g in tfs.values()]
)
dist["is_DBD"] = dist["accession"].isin(load_dbd_accessions())


# In[4]:


y1h = load_y1h_pdi_data()
y1h = y1h.drop_duplicates()  # TODO: why is this here?
n_pdi = y1h.drop(columns="gene_symbol").set_index("clone_acc").sum(axis=1)


# In[5]:


tfs = load_annotated_TFiso1_collection()

df = pd.concat(
    [g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()]
)
df["is_DBD"] = df["accession"].isin(load_dbd_accessions())
df["is_DBD_flank"] = df["accession"].str.endswith("_flank_N") | df[
    "accession"
].str.endswith("_flank_C")
df_new = (
    df.loc[df["is_DBD"], :]
    .groupby(["gene", "ref_iso", "alt_iso"])[["deletion", "frameshift"]]
    .sum()
    .sum(axis=1)
    / df.loc[df["is_DBD"], :].groupby(["gene", "ref_iso", "alt_iso"])["length"].sum()
).to_frame(name="dbd_fraction")

df_new["dbd_insertion_n_aa"] = (
    df.loc[df["is_DBD"], :].groupby(["gene", "ref_iso", "alt_iso"])["insertion"].sum()
)

df_new["dbd_n_aa_to_change"] = (
    dist.loc[dist["is_DBD"], :]
    .groupby(["gene", "ref_iso", "alt_iso"])["n_aa_change_to_domain"]
    .min()
)

# flank affected
df_new["dbd_flank_affected"] = (
    df.loc[df["is_DBD_flank"], :]
    .groupby(["gene", "ref_iso", "alt_iso"])[["deletion", "insertion", "frameshift"]]
    .sum()
    .sum(axis=1)
    > 0
)
df = df_new.reset_index()
df["dbd_pct_lost"] = df["dbd_fraction"] * 100.0


def dbd_affected_categories(pct_lost):
    if pct_lost < 0:
        raise ValueError("negative percent value")
    elif pct_lost == 0:
        return "Full DBD in\nalternative isoform"
    elif pct_lost >= 100:
        return "Full loss\nof DBD"
    else:
        return "Partial loss\nof DBD"


df["dbd_affected"] = df["dbd_pct_lost"].apply(dbd_affected_categories)
df["dbd_or_flank_affected"] = df["dbd_affected"]
df.loc[
    (df["dbd_affected"] == "Full DBD in\nalternative isoform")
    & df["dbd_flank_affected"],
    "dbd_or_flank_affected",
] = "DBD flank affected"

isoforms = load_valid_isoform_clones()
y1h = load_y1h_pdi_data()
y1h = y1h.drop_duplicates()  # TODO: why is this here?
n_pdi = y1h.drop(columns="gene_symbol").set_index("clone_acc").sum(axis=1)
n_pdi.index = n_pdi.index.map(
    lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
)


# map each isoform to change in PDI vs reference
def delta_pdi(row):
    iso_acc = row["alt_iso"]
    ref_acc = row["ref_iso"]
    if iso_acc == ref_acc:
        return np.nan
    n_ref = n_pdi.get(ref_acc, np.nan)
    n_iso = n_pdi.get(iso_acc, np.nan)
    if n_ref == 0:
        return np.nan
    return (n_iso - n_ref) / n_ref


df["delta_pdi"] = df.apply(delta_pdi, axis=1)
df = df.dropna(subset=["delta_pdi"])

df["tf_family"] = df["gene"].map(lambda x: tfs[x].tf_family)
df["delta_pdi_trunc"] = df["delta_pdi"].clip(upper=1)

if (
    ((df["dbd_fraction"] > 0) | (df["dbd_insertion_n_aa"] > 0))
    & (df["dbd_n_aa_to_change"] > 0)
).any():
    raise UserWarning("something wrong with calculations")
if (
    (df["dbd_fraction"] == 0)
    & (df["dbd_insertion_n_aa"] == 0)
    & (df["dbd_n_aa_to_change"] == 0)
).any():
    raise UserWarning("something wrong with calculations")


# In[6]:


df.loc[df["dbd_insertion_n_aa"] > 0]


# In[7]:


# count
print(len(tfs), "TF genes")
print(
    sum([len(tf.isoforms[0].aa_seq_features) > 0 for tf in tfs.values()]),
    "TF genes with at least one Pfam domain in longest cloned isoform",
)
print(sum([len(tf.isoforms[0].dna_binding_domains) > 0 for tf in tfs.values()]))
tfs_no_dbd = {
    k: v
    for k, v in tfs.items()
    if len(v.isoforms[0].dna_binding_domains) == 0
    and len(v.isoforms[0].aa_seq_features) > 0
}


# In[8]:


fig, ax = plt.subplots(1, 1)
ax.scatter(df["dbd_pct_lost"].values, df["delta_pdi"].values, alpha=0.5)
ax.axhline(0, linestyle="--", color="grey")
ax.set_ylabel("Delta PDIs in alt. iso / # PDIs in ref")
ax.set_xlabel("% DBD removed in alternative isoform")
plt.savefig("../figures/DBD_change_vs_PDI_scatter.pdf", bbox_inches="tight")


# In[9]:


df.sort_values("delta_pdi", ascending=False).head()


# In[10]:


fig, axes = plt.subplots(1, 2)
axes[0].hist(
    df.loc[df["dbd_pct_lost"] == 0, "delta_pdi"].values * 100,
    range=(-100, 500),
    bins=6 * 10,
)
axes[1].hist(
    df.loc[df["dbd_pct_lost"] > 0, "delta_pdi"].values * 100,
    range=(-100, 500),
    bins=6 * 10,
)
for ax in axes:
    ax.axvline(0, linestyle="--", color="grey")
    ax.set_xlim(-110, 510)


# In[11]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(3.2, 3)
df["delta_pdi_trunc"] = df["delta_pdi"].clip(upper=1)
sns.swarmplot(data=df, y="delta_pdi_trunc", x="dbd_affected", size=3.2, ax=ax, alpha=1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_xlabel("")
ax.set_yticks([-1, 0, 1])
ax.set_yticks(np.linspace(-1, 1, 9), minor=True)
ax.set_yticklabels(["-100%", "0", "+≥100%"])
ax.invert_xaxis()
ax.set_ylabel("Change in number of PDI\nin alternative isoform")
plt.savefig("../figures/DBD_change_vs_PDI_truncated_inverted.pdf", bbox_inches="tight")


# In[12]:


df["dbd_or_flank_affected"].value_counts().index.values


# In[13]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=3.8, h=3)
df["delta_pdi_trunc"] = df["delta_pdi"].clip(upper=1)
sns.swarmplot(
    data=df,
    y="delta_pdi_trunc",
    x="dbd_or_flank_affected",
    size=3.2,
    order=[
        "Full DBD in\nalternative isoform",
        "Partial or full\nloss of DBD",
        "DBD flank affected",
    ],
    ax=ax,
    alpha=1,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_xlabel("")
ax.set_yticks([-1, 0, 1])
ax.set_yticks(np.linspace(-1, 1, 9), minor=True)
ax.set_yticklabels(["-100%", "0", "+≥100%"])
ax.set_ylabel("Change in number of PDI\nin alternative isoform")
plt.savefig(
    "../figures/DBD_or_flank_change_vs_PDI_truncated_inverted.pdf", bbox_inches="tight"
)


# In[14]:


df["tf_family_merged"] = df["tf_family"].map(
    lambda x: x
    if x in ["C2H2 ZF", "bHLH", "Homeodomain", "Nuclear receptor"]
    else "other"
)


# In[15]:


# fraction of 0 PDI alternative isoforms by family
df["full_PDI_loss"] = df["delta_pdi_trunc"] == -1
tf_families = ["C2H2 ZF", "bHLH", "Homeodomain", "Nuclear receptor", "other"]
fs = [
    df.loc[df["tf_family_merged"] == fam, "full_PDI_loss"].mean() for fam in tf_families
]
ns = [(df["tf_family_merged"] == fam).sum() for fam in tf_families]
es = [np.sqrt((f * (1 - f)) / n) for f, n in zip(fs, ns)]

fig, ax = plt.subplots(1, 1)
ax.bar(x=[fam + f"\n(N = {n})" for fam, n in zip(tf_families, ns)], height=fs, yerr=es)
ax.set_ylim(0, 1)
ax.set_ylabel("Fraction of full PDI loss alternative isoforms")
fig.savefig(
    "../figures/full-PDI-loss-alternative-isoforms-by-family.pdf", bbox_inches="tight"
)


# In[16]:


# TODO: move to data_loading.py
dis = pd.read_csv(
    "../data/processed/TFiso1_disorder-and-ss_from-alphafold.tsv", sep="\t"
)
n_aa = dis.groupby("clone_name").size().rename("n_aa").to_frame()
n_aa["n_aa_disordered"] = dis.groupby("clone_name")["is_disordered"].sum()
n_aa["n_aa_ordered"] = n_aa["n_aa"] - n_aa["n_aa_disordered"]
for c in n_aa.columns:
    df[f"delta_{c}"] = df["ref_iso"].map(n_aa[c]) - df["alt_iso"].map(n_aa[c])
    df[f"abs_delta_{c}"] = df[f"delta_{c}"].abs()


# In[17]:


df["f_disorder_delta_aa"] = df["abs_delta_n_aa_disordered"] / (
    df["abs_delta_n_aa_disordered"] + df["abs_delta_n_aa_ordered"]
)


# In[18]:


df["pdi_affected"] = df["delta_pdi"] != 0
# for isoforms with full DBD + flanks, split by PDI affected vs not
# and plot changes in sequence, disorder and charged residues
fig, axs = plt.subplots(1, 2)
sns.swarmplot(
    data=df.loc[df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform", :],
    x="pdi_affected",
    y="abs_delta_n_aa_disordered",
    ax=axs[0],
)
sns.swarmplot(
    data=df.loc[df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform", :],
    x="pdi_affected",
    y="abs_delta_n_aa_ordered",
    ax=axs[1],
)


# In[19]:


fig, ax = plt.subplots(1, 1)
sns.swarmplot(
    data=df.loc[df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform", :],
    x="pdi_affected",
    y="f_disorder_delta_aa",
    ax=ax,
)
sns.pointplot(
    data=df.loc[df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform", :],
    x="pdi_affected",
    y="f_disorder_delta_aa",
    ax=ax,
)


# In[20]:


# charged_aas = ("K", "R", "H", "D", "E", "C", "Y")


def delta_n_K_or_R(row):
    aa_ref = tfs[row["gene"]][row["ref_iso"]].aa_seq
    aa_alt = tfs[row["gene"]][row["alt_iso"]].aa_seq
    return (aa_alt.count("K") + aa_alt.count("R")) - (
        aa_ref.count("K") + aa_ref.count("R")
    )


df["delta_n_K_or_R"] = df.apply(delta_n_K_or_R, axis=1)


# In[21]:


# for alternative isoforms, containing the DBD, that are associated with a change in
# DNA binding vs staying the same, is there a positively charged region in the different
# amino acid sequence?

# three categories: lose, same, gain
# variable is max postive charged residue count in 10aa sliding window,


# NOTE I've updated this since the copy and paste from domain/motif notebook
def isoform_specific_regions(gene, subset=None):
    """The name is a bit misleading because it's not specific to one isoform but just
       not common to all isoforms.

    Returns: dict(frozenset: list(str)): isoform IDs and list of contiguous AA sequences
                                         that map to them only

    """
    algn = gene.genomic_alignment_of_aa_seqs(subset=subset)
    subset_prev = None
    isr = {}
    len_algn = len(list(algn.values())[0])
    for i in range(len_algn):
        subset = frozenset(
            {k for k, v in algn.items() if v[i] != "-"}
        )  # the isoforms that have an aa at that genomic position
        if subset_prev is None:
            if (len(subset) < len(algn)) and (len(subset) > 0):
                start = i
                subset_prev = subset
        else:
            if subset != subset_prev:
                if (len(subset_prev) < len(algn)) and (len(subset_prev) > 0):
                    subseq = (algn[list(subset_prev)[0]][start:i], start, i)
                    isr[subset_prev] = isr.get(subset_prev, []) + [subseq]
                start = i
                subset_prev = subset
            elif i == (len_algn - 1):
                if (len(subset_prev) < len(algn)) and (len(subset_prev) > 0):
                    subseq = (algn[list(subset_prev)[0]][start:], start, i + 1)
                    isr[subset_prev] = isr.get(subset_prev, []) + [subseq]
                start = i
                subset_prev = subset
    merged = {}
    for iso_subset, subseqs in isr.items():
        merged[iso_subset] = []
        prev_end = np.inf
        prev_subseq = ""
        for subseq, start, end in subseqs:
            if start <= prev_end + 2:
                prev_subseq += subseq
                prev_end = end
            else:
                if prev_subseq != "":
                    merged[iso_subset].append(prev_subseq)
                prev_subseq = subseq
                prev_end = end
        merged[iso_subset].append(prev_subseq)
    merged
    return merged


def count_K_or_R_sliding_window(row, window_size=6):
    # get different amino acid sequence regions
    specific_aa_regions = isoform_specific_regions(
        tfs[row["gene"]], subset=[row["ref_iso"], row["alt_iso"]]
    )
    aa_seq_regions_in_ref_not_in_alt = specific_aa_regions.get(
        frozenset([row["ref_iso"]]), []
    )
    aa_seq_regions_in_alt_not_in_ref = specific_aa_regions.get(
        frozenset([row["alt_iso"]]), []
    )
    max_count = 0
    for aa_seq in aa_seq_regions_in_ref_not_in_alt:
        for i in range(len(aa_seq) - window_size):
            count = aa_seq[i : i + window_size].count("K") + aa_seq[
                i : i + window_size
            ].count("R")
            if count > max_count:
                max_count = count
    return max_count


df["count_K_or_R_sliding_window"] = df.apply(count_K_or_R_sliding_window, axis=1)


# In[22]:


df["count_K_or_R_sliding_window"].value_counts()


# In[23]:


df["dbd_or_flank_affected"].value_counts()


# In[24]:


# only isoforms with DBD
# three categories: lose, same, gain
# variable is max postive charged residue count in 10aa sliding window,
a = df.loc[
    (df["dbd_affected"] == "Full DBD in\nalternative isoform") & (df["delta_pdi"] < 0),
    "count_K_or_R_sliding_window",
].values
b = df.loc[
    (df["dbd_affected"] == "Full DBD in\nalternative isoform") & (df["delta_pdi"] == 0),
    "count_K_or_R_sliding_window",
].values
print(a)
print(b)
print(a.mean())
print(b.mean())
print(stats.mannwhitneyu(a, b))


# In[25]:


isoform_specific_regions(tfs["KLF7"], subset=["KLF7-1", "KLF7-4"])


# In[26]:


y, x = df.loc[
    df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform",
    ["delta_pdi_trunc", "delta_n_K_or_R"],
].values.T

fig, ax = plt.subplots(1, 1)
ax.set_title("r = {:.2f}, p = {:.2f}".format(*stats.pearsonr(x, y)))
ax.scatter(x=x, y=y)
ax.set_xlabel("reference vs alternative (Arg + Lys) count difference")
ax.set_yticks([-1, 0, 1])
ax.set_yticks(np.linspace(-1, 1, 9), minor=True)
ax.set_yticklabels(["-100%", "0", "+≥100%"])
ax.set_ylabel("Change in number of PDI\nin alternative isoform")
fig.savefig(
    "../figures/PDI-change-vs-Arg-Lys-count-diff_full-DBD-and-flanks-proteins_scatter.pdf",
    bbox_inches="tight",
)


# In[27]:


y, x = df.loc[
    df["dbd_or_flank_affected"] == "DBD flank affected",
    ["delta_pdi_trunc", "delta_n_K_or_R"],
].values.T

fig, ax = plt.subplots(1, 1)
ax.set_title("r = {:.2f}, p = {:.2f}".format(*stats.pearsonr(x, y)))
ax.scatter(x=x, y=y)
ax.set_xlabel("reference vs alternative (Arg + Lys) count difference")
ax.set_yticks([-1, 0, 1])
ax.set_yticks(np.linspace(-1, 1, 9), minor=True)
ax.set_yticklabels(["-100%", "0", "+≥100%"])
ax.set_ylabel("Change in number of PDI\nin alternative isoform")
fig.savefig(
    "../figures/PDI-change-vs-Arg-Lys-count-diff_DBD-flank-affected-proteins_scatter.pdf",
    bbox_inches="tight",
)


# In[28]:


y, x = (
    df.loc[
        df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform",
        ["delta_pdi_trunc", "f_disorder_delta_aa"],
    ]
    .dropna()
    .values.T
)

fig, ax = plt.subplots(1, 1)
ax.set_title("r = {:.2f}, p = {:.2f}".format(*stats.pearsonr(x, y)))
ax.scatter(x=x, y=y)
ax.set_xlabel("Fraction of alternate regions in predicted disordered regions")
ax.set_yticks([-1, 0, 1])
ax.set_yticks(np.linspace(-1, 1, 9), minor=True)
ax.set_yticklabels(["-100%", "0", "+≥100%"])
ax.set_ylabel("Change in number of PDI\nin alternative isoform")
fig.savefig(
    "../figures/PDI-change-vs-disordered-fraction-of-alternate-regions_full-DBD-and-flanks-proteins_scatter.pdf",
    bbox_inches="tight",
)


# In[29]:


# check for family enrichment of DBD unaffected PDI changes
df.loc[
    (df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform")
    & (df["delta_pdi"] != 0),
    "tf_family",
].value_counts()


# In[30]:


df.loc[
    (df["dbd_or_flank_affected"] == "Full DBD in\nalternative isoform")
    & (df["delta_pdi"] != 0),
    "gene",
].value_counts()


# In[31]:


# PPIs with other TFs as a predictor for non-DBD related PDI changes?


# In[32]:


# 15 aa flanks
" ".join(
    df.loc[
        (df["dbd_fraction"] == 0)
        & (df["dbd_flank_affected"] == False)
        & (df["delta_pdi"] != 0),
        "gene",
    ].unique()
)


# In[33]:


(df["dbd_pct_lost"] > 0).sum()


# In[34]:


fig, ax = plt.subplots(1, 1)
ax.scatter(
    df.loc[
        (df["dbd_pct_lost"] > 0) & (df["dbd_pct_lost"] < 100), "dbd_pct_lost"
    ].values,
    df.loc[
        (df["dbd_pct_lost"] > 0) & (df["dbd_pct_lost"] < 100), "delta_pdi_trunc"
    ].values,
    alpha=1,
)
ax.axhline(0, linestyle="--", color="grey")

ax.set_xlabel("% DBD removed in alternative isoform")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xlim(0, 100)
ax.set_yticks([-1, 0, 1])
ax.set_yticks(np.linspace(-1, 1, 9), minor=True)
ax.set_yticklabels(["-100%", "0", "+≥100%"])
ax.set_ylabel("Change in number of PDI in alternative isoform")
plt.savefig(
    "../figures/DBD_change_vs_PDI_scatter_only_partially_affected.pdf",
    bbox_inches="tight",
)


# In[35]:


df.sort_values("delta_pdi", ascending=False).head()


# In[36]:


# try distance from DBD
# TODO
# check y variable now that we use reference isoform
# horizontal line across whole

fig, axs = plt.subplots(nrows=1, ncols=4, sharey=True, width_ratios=[1, 1, 0.35, 1.5])
fig.set_size_inches(w=8, h=2)
point_size = 7

axs[0].set_title("Full loss of DBD", fontsize=10)
sns.swarmplot(
    data=df,
    y="delta_pdi_trunc",
    x="dbd_or_flank_affected",
    size=point_size,
    order=[
        "Full loss\nof DBD",
    ],
    ax=axs[0],
    color=COLOR_PURPLE,
    alpha=1,
)

axs[1].set_title("Partial loss of DBD", fontsize=10)
axs[1].scatter(
    df.loc[
        (df["dbd_pct_lost"] > 0) & (df["dbd_pct_lost"] < 100), "dbd_pct_lost"
    ].values,
    df.loc[
        (df["dbd_pct_lost"] > 0) & (df["dbd_pct_lost"] < 100), "delta_pdi_trunc"
    ].values,
    alpha=1,
    s=point_size**2
    / 1.5,  # I don't think there should be a divide by anything here....
    color=COLOR_PURPLE,
    clip_on=False,
)
axs[1].set_xlabel("Proportion missing")
axs[1].set_xlim(100, 0)
axs[1].set_xticks([99, 50, 1])
axs[1].set_xticklabels(["{}%".format(x) for x in axs[1].get_xticks()])
axs[1].set_xticks(range(10, 91, 10), minor=True)


axs[2].set_title("Insertion\nwithin DBD", fontsize=10)
axs[2].scatter(
    df.loc[
        (df["dbd_pct_lost"] == 0) & (df["dbd_insertion_n_aa"] > 0), "dbd_insertion_n_aa"
    ].values,
    df.loc[
        (df["dbd_pct_lost"] == 0) & (df["dbd_insertion_n_aa"] > 0), "delta_pdi_trunc"
    ].values,
    alpha=1,
    s=point_size**2
    / 1.5,  # I don't think there should be a divide by anything here....
    color=COLOR_PURPLE,
    clip_on=False,
)
axs[2].set_xlabel("amino acids\ninserted")
axs[2].set_xticks([1, 4])
axs[2].set_xticks(range(1, 6), minor=True)

axs[3].set_title("Full DBD in\nalternative isoform", fontsize=10)
axs[3].scatter(
    df.loc[
        (df["dbd_affected"] == "Full DBD in\nalternative isoform"), "dbd_n_aa_to_change"
    ].values,
    df.loc[
        (df["dbd_affected"] == "Full DBD in\nalternative isoform"), "delta_pdi_trunc"
    ].values,
    alpha=1,
    s=point_size**2
    / 1.5,  # I don't think there should be a divide by anything here....
    color=COLOR_PURPLE,
    clip_on=False,
)
axs[3].set_xlabel("Distance of alternative\nsequence from DBD\n(number of AA)")
# axs[3].set_xlim(100, 0)
# axs[3].set_xticks([99, 50, 1])
# axs[3].set_xticklabels(['{}%'.format(x)for x in axs[1].get_xticks()])
# axs[3].set_xticks(range(10, 91, 10), minor=True)


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
for ax in axs[1:]:
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(which="both", length=0)
for i in [0]:
    axs[i].set_xlabel("")
    axs[i].set_ylabel("")
    axs[i].spines["bottom"].set_visible(False)
    axs[i].xaxis.set_tick_params(length=0)
    axs[i].set_xticks([])
axs[0].set_yticks([-1, 0, 1])
axs[0].set_yticks(np.linspace(-1, 1, 9), minor=True)
axs[0].set_yticklabels(["-100%", "0", "+≥100%"])
axs[0].set_ylabel("Change in number of PDI\nin alternative isoform")

plt.savefig(
    "../figures/DBD_or_flank_change_vs_PDI_composite_alt_with_distance.pdf",
    bbox_inches="tight",
)


# In[37]:


# full DBD in alternative isoform, fraction in disordered
df["f_disorder_difference"] = df.apply(
    lambda x: disordered_fraction_of_different_regions(
        tfs[x["gene"]], x["ref_iso"], x["alt_iso"]
    ),
    axis=1,
)


# In[38]:


# check low values
fig, ax = plt.subplots(1, 1, tight_layout=True)
fig.set_size_inches(h=3, w=3)
sns.swarmplot(
    data=df.loc[
        (df["dbd_affected"] == "Full DBD in\nalternative isoform")
        & (df["delta_pdi_trunc"] != 0),
        :,
    ],
    y="f_disorder_difference",
    color=COLOR_PURPLE,
    ax=ax,
    clip_on=False,
    # alpha=0.5  # TMP - to check no overlap
)
ax.set_xlabel("Full DBD in alternative isoform\n+ difference in DNA binding")
ax.set_ylabel("Fraction of alternative protein\nsequence in disordered regions")
ax.set_xticks([])
for loc in ["top", "bottom", "right"]:
    ax.spines[loc].set_visible(False)
ax.set_yticks(np.linspace(0, 1, 11))
ax.set_yticks(np.linspace(0, 1, 21), minor=True)
ax.set_ylim(0, 1)
ax.set_yticklabels(["{:.0%}".format(y) for y in ax.get_yticks()])
fig.savefig(
    "../figures/disordered-pct-alt-sequence_alt-isoforms-full-DBD-diff-PDI_dotplot.pdf",
    bbox_inches="tight",
)


# In[39]:


# DEBUG
df.sort_values("f_disorder_difference", ascending=True).head(10)


# In[40]:


# could also do number of aa different or aa similarity
df.loc[df["gene"] == "TBX5"]


# In[41]:


df.sort_values("dbd_n_aa_to_change", ascending=False).head(10)


# In[42]:


iso = tfs["HEY1"]["HEY1-2"]
iso.aa_seq[iso.aa_seq_features[0].start : iso.aa_seq_features[0].end]


# In[43]:


iso = tfs["ZIC3"]["ZIC3-2"]
iso.aa_seq[iso.aa_seq_features[0].start : iso.aa_seq_features[0].end]


# In[44]:


len(iso.aa_seq)


# In[45]:


(1060 + 164) / 3


# In[46]:


iso = tfs["ZIC3"]["ZIC3-2"]
(iso.aa_seq_features[0].start, iso.aa_seq_features[0].end)


# In[47]:


iso = tfs["HEY1"]["HEY1-2"]
(iso.aa_seq_features[0].start, iso.aa_seq_features[0].end)


# In[48]:


tfs["HEY1"].pairwise_changes_relative_to_reference("HEY1-2", "HEY1-1").find("I")


# In[49]:


tfs["CREB1"].pairwise_changes_relative_to_reference("CREB1-2", "CREB1-1").find("I")


# In[50]:


iso = tfs["CREB1"]["CREB1-1"]
(iso.aa_seq_features[1].start, iso.aa_seq_features[1].end)


# In[51]:


iso = tfs["KLF7"]["KLF7-1"]
print(iso.aa_seq_features[0].start, iso.aa_seq_features[0].end)
iso.aa_seq[iso.aa_seq_features[0].start : iso.aa_seq_features[0].end]


# In[52]:


iso = tfs["KLF7"]["KLF7-4"]
print(iso.aa_seq_features[0].start, iso.aa_seq_features[0].end)


# In[53]:


print(tfs["KLF7"].pairwise_changes_relative_to_reference("KLF7-1", "KLF7-4").rfind("I"))
tfs["KLF7"].pairwise_changes_relative_to_reference("KLF7-1", "KLF7-4")


# In[54]:


# closest PDB with DNA for HEY1 is 4H10 chain A (ARNTL) sequence ID 59%


# In[55]:


df.loc[
    (df["dbd_pct_lost"] > 0)
    & (df["dbd_pct_lost"] < 100)
    & (df["delta_pdi_trunc"] > -1),
    :,
]


# In[56]:


tfs["ZIC3"]["ZIC3-2"].aa_seq_features


# In[57]:


df.head()


# In[58]:


df.loc[
    (df["dbd_pct_lost"] == 0)
    & (df["dbd_flank_affected"] == False)
    & (df["delta_pdi_trunc"] > -1),
    :,
].sort_values("delta_pdi").tail()


# In[59]:


df.head()


# In[60]:


df.loc[df["dbd_insertion_n_aa"] == 1]


# In[61]:


tfs["EBF3"].exon_diagram()


# In[62]:


m1h = load_m1h_activation_data()


# In[63]:


y1h.head()


# In[64]:


m1h.head()


# In[65]:


# TODO: move to per TF gene notebook
for gene_name, tf in tfs.items():
    if gene_name in y1h["gene_symbol"].values:
        y1h_pdi_per_tf_gene_plot(gene_name, y1h)
        plt.savefig(
            "../figures/per_gene/y1h_profile/{}_y1h_profile.pdf".format(gene_name),
            bbox_inches="tight",
        )
        plt.close(plt.gcf())
    if gene_name in m1h["gene"].values:
        m1h_activation_per_tf_gene_plot(gene_name, m1h)
        plt.savefig(
            "../figures/per_gene/m1h_profile/{}_m1h_profile.pdf".format(gene_name),
            bbox_inches="tight",
        )
        plt.close(plt.gcf())
    tf.exon_diagram()
    plt.savefig(
        "../figures/per_gene/exon_diagram/{}_exon_diagram.pdf".format(gene_name),
        bbox_inches="tight",
    )
    plt.close(plt.gcf())


# In[66]:


tfs["HEY1"].exon_diagram()


# In[67]:


y1h_pdi_per_tf_gene_plot("HEY1", y1h)


# In[68]:


tfs["ZIC3"].exon_diagram()


# In[69]:


y1h_pdi_per_tf_gene_plot("ZIC3", y1h)


# In[70]:


tfs["EBF3"].isoforms[0].dna_binding_domains


# In[71]:


df.loc[df["dbd_affected"].str.startswith("Full"), :].sort_values("delta_pdi")


# In[72]:


tfs["PPARG"].genomic_alignment_of_aa_seqs()


# In[73]:


# load DNA bait sequnces
# look for E-box within HS385 / HS416 / MUT_116
baits = load_Y1H_DNA_bait_sequences()


# In[74]:


y1h = load_y1h_pdi_data()


# In[75]:


# DEBUG
y1h.tail().dtypes


# In[76]:


def baits_set(row):
    return set(row.columns[2:].values[row.iloc[0, 2:].fillna(False).values])


df = (
    y1h.groupby(["gene_symbol", "clone_acc"])
    .apply(baits_set)
    .rename("baits")
    .reset_index()
)


def is_differential(gene):
    non_zero = (frozenset(x) for x in gene["baits"].values if len(x) > 0)
    return len(set(non_zero)) >= 2


df = df.loc[df["gene_symbol"].map(df.groupby("gene_symbol").apply(is_differential))]
df["baits"] = df["baits"].apply(lambda x: " ".join(x))
df.to_csv("../output/Y1H_differential-subset.tsv", index=False, sep="\t")


# In[77]:


df["gene_symbol"].unique()


# In[78]:


# look for E-box within hs385 / hs416 / MUT_116

# Ebox is CANNTG

# re.search('CA..TG', baits['MUT_116'])
re.findall("CA..TG", baits["hs416"])


# In[79]:


re.findall("CA..TG", baits["hs385"])


# In[80]:


baits["hs385"]


# In[81]:


baits["hs416"]


# In[82]:


baits["MUT_116"]


# In[83]:


re.findall("CG..TG", baits["MUT_116"])


# In[84]:


re.findall("CG..TG", baits["hs385"])


# In[85]:


re.findall("CG..TG", baits["hs416"])
