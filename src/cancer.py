# coding: utf-8

# In[1]:


import pandas as pd

from data_loading import load_valid_isoform_clones


# In[2]:


df = pd.read_excel(
    "../data/internal/from_kaia/prioritized_BRCA_isos_to_clone.CDS_seqs.xlsx"
)


# In[3]:


df["seq_cds"]


# In[4]:


iso = load_valid_isoform_clones()
iso["cds_len"] = iso["cds"].str.len()
iso.head()


# In[7]:


# is the matching on CDS working?
df.loc[~df["seq_cds"].str.startswith("ATG"), :]


# In[10]:


iso.loc[~iso["cds"].str.startswith("ATG"), :]


# In[34]:


from data_loading import (
    load_isoform_and_paralog_y2h_data,
    load_y1h_pdi_data,
    load_m1h_activation_data,
)

y2h = load_isoform_and_paralog_y2h_data()
y1h = load_y1h_pdi_data()
m1h = load_m1h_activation_data()


# In[69]:


ppi_per_clone = (
    y2h.loc[(y2h["category"] == "tf_isoform_ppis") & y2h["score"].isin(["1", "0"])]
    .groupby(["ad_gene_symbol", "ad_clone_acc"])
    .size()
)
clones_gte_2ppi = ppi_per_clone[ppi_per_clone >= 2].groupby("ad_gene_symbol").size()
tf_with_ppi_data = set(clones_gte_2ppi[clones_gte_2ppi >= 2].index)
df["at_least_2_ppis_for_two_isoforms"] = df["gene"].isin(tf_with_ppi_data)

y1h["any_pdi"] = y1h.iloc[:, 2:].any(axis=1)
tf_with_pdi_data = (y1h.groupby("gene_symbol").size() >= 2) & y1h.groupby(
    "gene_symbol"
)["any_pdi"].any()
tf_with_pdi_data = set(tf_with_pdi_data[tf_with_pdi_data].index)
df["at_least_1_pdi"] = df["gene"].isin(tf_with_pdi_data)


df["at_least_2_m1h_data"] = df["gene"].isin(
    set(m1h.groupby("gene").size()[m1h.groupby("gene").size() >= 2].index)
)


# In[70]:


# use this. Add tf gene level PPI + PDI + M1H data
df.loc[
    df["index"].apply(lambda x: x.split("_")[-1].startswith("TF"))
    & df["at_least_2_ppis_for_two_isoforms"]
    & df["at_least_1_pdi"]
    & df["at_least_2_m1h_data"],
    :,
]


# In[25]:


# WT1:::GENCPID22245__NA__TF1P0PID361 not matching CDS length in Kaia's table: 867, in isoform table: 864
a = iso.loc[iso["clone_acc"] == "WT1|6/6|10G06"]["cds"].values[0]
b = df.loc[df["index"] == "WT1:::GENCPID22245__NA__TF1P0PID361", "seq_cds"].values[0]
a


# In[77]:


iso.loc[iso["gene"] == "NFIX"]


# ## candidate genes for cancer vinyette
#
# - MAX (lots of isoforms)
# - MLX (lots of isoforms)
# - NFIX?
# - PPARG-1
#     - data looks good...
# - TEAD2-1
#     - data looks good
# - WT1-6
#     - we have a lot of data (already known splicing effects in cancer)
# - ZBTB18
#     - data looks good

# In[26]:


b


# In[27]:


df["seq_cds"].apply(lambda x: x[-3:] in {"TAG", "TAA", "TGA"}).value_counts()


# In[28]:


iso["cds"].apply(lambda x: x[-3:] in {"TAG", "TAA", "TGA"}).value_counts()


# In[30]:


def remove_stop_codon(s):
    if s[-3:] in {"TAG", "TAA", "TGA"}:
        return s[:-3]
    else:
        return s


iso["cds"] = iso["cds"].apply(remove_stop_codon)
df["seq_cds"] = df["seq_cds"].apply(remove_stop_codon)


# In[31]:


pd.merge(df, iso, left_on="seq_cds", right_on="cds")


# In[6]:


iso["gene"].nunique()


# In[7]:


df = pd.read_csv("/Users/lukelambourne/Downloads/cancer_gene_census.csv")


# In[8]:


cgc = set(df["Gene Symbol"].unique())
cgc.intersection(set(iso["gene"].unique()))


# In[9]:


df.head()


# In[10]:


df = df.loc[
    df["Gene Symbol"].isin(iso["gene"].unique()),
    [
        "Gene Symbol",
        "Tumour Types(Somatic)",
        "Tissue Type",
        "Molecular Genetics",
        "Role in Cancer",
        "Mutation Types",
    ],
]


# In[11]:


df.shape


# In[12]:


df["n_cloned_isoforms"] = df["Gene Symbol"].map(iso.groupby("gene").size())


# In[13]:


df["has_novel_isoform"] = df["Gene Symbol"].map(
    iso.groupby("gene")["is_novel_isoform"].any()
)


# In[16]:


# add has useful
df.loc[df["n_cloned_isoforms"] >= 2, :].to_csv(
    "/Users/lukelambourne/Desktop/tmp.csv", index=False
)


# In[17]:


(df["n_cloned_isoforms"] >= 2).sum()
