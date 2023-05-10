
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')


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

from Bio.Seq import Seq
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr

import plotting
from plotting import PAPER_PRESET, PAPER_FONTSIZE, nice_boxplot, nice_violinplot, mimic_r_boxplot


get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'svg'")
mpl.rcParams['figure.autolayout'] = False


# In[3]:


from data_loading import (load_annotated_6k_collection,
                          load_valid_isoform_clones,
                          load_developmental_tissue_expression_remapped,
                          load_gtex_remapped)


# In[4]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[5]:


np.random.seed(2023)


# ## functions

# In[6]:


def calculate_tau(df):
    array = df.values
    
    ## will return NaN as tau for every row that has any NaNs
    array_max = np.max(array, axis=1)
    tmp = array.T / array_max
    tmp = 1 - tmp.T
    nonan_taus = np.sum(tmp, axis=1) / (array.shape[1])
    
    ## will ignore NaNs and compute on the rest of the values
    array_max = np.nanmax(array, axis=1)
    tmp = array.T / array_max
    tmp = 1 - tmp.T
    nan_taus = np.nansum(tmp, axis=1) / np.count_nonzero(~np.isnan(array), axis=1)
    
    
    return nonan_taus, nan_taus, array_max


# In[7]:


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


# ## variables

# In[8]:


dn_cats_f = "../data/processed/DN_cats_Joung.tsv"


# In[9]:


pal = {"ref": sns.color_palette("Set2")[0],
       "ref-v-ref": sns.color_palette("Set2")[0],
       "rewire": sns.color_palette("Set2")[2],
       "DN": sns.color_palette("Set2")[1],
       "NA": "lightgray",
       "likely": "darkgray",
       "both": sns.color_palette("Set2")[3]}


# ## 1. import data

# In[10]:


dn_cats = pd.read_table(dn_cats_f)
dn_cats["dn_cat"].fillna("NA", inplace=True)
dn_cats.dn_cat.value_counts()


# In[11]:


tfs = load_annotated_6k_collection()


# In[12]:


df_gtex, metadata_gtex, genes_gtex = load_gtex_remapped()

exclusion_list_gtex = {'Cells - Leukemia cell line (CML)',
                       'Cells - EBV-transformed lymphocytes',
                       'Cells - Cultured fibroblasts'}

df_gtex = df_gtex.loc[:, ~df_gtex.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex = metadata_gtex.loc[~metadata_gtex['body_site'].isin(exclusion_list_gtex), :]

means_gtex = df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean()


# In[13]:


metadata_gtex_dummy = pd.read_table("../data/processed/metadata_gtex_dummy.csv", sep=",", index_col=0)


# In[14]:


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


# ## 2. down-sample gtex using the same dummy metadata as fig 1

# In[15]:


means_gtex_downsample = df_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean()


# ## 3. calculate expression ratios (copy-pasted from luke's code in expr_data.ipynb)

# In[16]:


per_gene_dev = ((2 ** df_dev - 1)
                .groupby(genes_dev)
                .transform('sum'))
f_dev = (((2 ** df_dev - 1) / per_gene_dev)
        .groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']),
         axis=1)
        .mean())
f_dev = f_dev * ((per_gene_dev.groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']),
                                             axis=1)
                                             .mean() >= 1)
                                         .applymap(lambda x: {False: np.nan, True: 1}[x]))  # only count fractions if gene TPM is >= 1

f_dev = f_dev * 100


# In[17]:


per_gene_gtex = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .transform('sum'))
f_gtex = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
        .mean())
f_gtex = f_gtex * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex = f_gtex * 100


# In[18]:


df_gtex.loc[:,metadata_gtex_dummy.index].shape


# In[19]:


df_gtex.shape


# In[20]:


per_gene_gtex_ds = ((2 ** df_gtex.loc[:,metadata_gtex_dummy.index] - 1)
                   .groupby(genes_gtex)
                   .transform('sum'))

f_gtex_downsample = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1)
        .mean())
f_gtex_downsample = f_gtex_downsample * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex_downsample = f_gtex_downsample * 100


# ## 3. calculate tissue-specificity
# 
# this is currently across all individual samples for genes but not for isos

# In[21]:


gene_dev_nonan_taus, gene_dev_nan_taus, gene_dev_array_max = calculate_tau(per_gene_dev.drop_duplicates())
gene_dev_nan_taus[0:5]


# In[22]:


gene_gtex_nonan_taus, gene_gtex_nan_taus, gene_gtex_array_max = calculate_tau(per_gene_gtex.drop_duplicates())
gene_gtex_nan_taus[0:5]


# In[23]:


gene_gtex_ds_nonan_taus, gene_gtex_ds_nan_taus, gene_gtex_ds_array_max = calculate_tau(per_gene_gtex_ds.drop_duplicates())
gene_gtex_ds_nan_taus[0:5]


# In[24]:


gene_taus = pd.DataFrame()
gene_taus["UID"] = per_gene_dev.drop_duplicates().index
gene_taus["dev_tau"] = gene_dev_nan_taus
gene_taus["gtex_tau"] = gene_gtex_nan_taus
gene_taus["gtex_ds_tau"] = gene_gtex_ds_nan_taus
gene_taus["gene_name"] = gene_taus["UID"].str.split("|", expand=True)[0]
gene_taus.sample(5)


# ## 4. join with DN categories

# In[25]:


indiv_cols = f_dev.columns
dev_ratios = f_dev.reset_index()

dev_ratios["clone_acc"] = dev_ratios["UID"].str.split(" ", expand=True)[0]
dev_ratios.head()


# In[26]:


indiv_cols = f_gtex.columns
gtex_ratios = f_gtex.reset_index()

gtex_ratios["clone_acc"] = gtex_ratios["UID"].str.split(" ", expand=True)[0]


# In[27]:


indiv_cols = f_gtex_downsample.columns
gtex_ds_ratios = f_gtex_downsample.reset_index()

gtex_ds_ratios["clone_acc"] = gtex_ds_ratios["UID"].str.split(" ", expand=True)[0]


# In[28]:


dev_ratios = dev_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id").drop_duplicates()
gtex_ratios = gtex_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id").drop_duplicates()
gtex_ds_ratios = gtex_ds_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id").drop_duplicates()
print(len(dev_ratios))
print(len(gtex_ratios))
print(len(gtex_ds_ratios))


# ## 5. calculate co-expression

# In[29]:


# def count_coex(row, thresh):
#     tot = 0
#     for col in indiv_cols:
#         ref_val = row["%s_ref" % col]
#         alt_val = row["%s_alt" % col]
#         if ref_val > thresh and alt_val > thresh:
#             tot += 1
#     return tot


# In[30]:


# ratios_ref = ratios[ratios["dn_cat"] == "ref"]
# ratios_alt = ratios[ratios["dn_cat"] != "ref"]

# ratios_v = ratios_ref.merge(ratios_alt, on="gene_name", suffixes=("_ref", "_alt"))
# ratios_v["num_coex"] = ratios_v.apply(count_coex, thresh=0.1, axis=1)
# ratios_v["num_coex"].fillna(0, inplace=True)
# ratios_v.sample(5)


# In[31]:


# # do the null test: randomly sample 2 refs 100 times and do the same thing
# null = pd.DataFrame()
# for i in range(100):
#     ref1 = ratios_ref.sample()
#     ref2 = ratios_ref.sample()
    
#     ref1["tmp_gene"] = "tmp_gene"
#     ref2["tmp_gene"] = "tmp_gene"
    
#     mrg = ref1.merge(ref2, on="tmp_gene", suffixes=("_ref", "_alt"))
#     mrg["gene_name"] = mrg["gene_name_ref"] + "-v-" + mrg["gene_name_alt"]
#     mrg["dn_cat_alt"] = "ref-v-ref"
#     mrg["num_coex"] = mrg.apply(count_coex, thresh=0.1, axis=1)
#     null = null.append(mrg)

# null.sample(5)


# In[32]:


# ratios_coex = ratios_v[["gene_name", "tf1p0_id_ref", 
#                         "tf1p0_id_alt", "dn_cat_alt", 
#                         "num_coex", "neglog_diff_pval_ref", 
#                         "neglog_diff_pval_alt"]].append(null[["gene_name", "tf1p0_id_ref", 
#                                                               "tf1p0_id_alt", "dn_cat_alt", "num_coex",
#                                                               "neglog_diff_pval_ref", "neglog_diff_pval_alt"]])
# ratios_coex.dn_cat_alt.value_counts()


# ## 6. make some plots

# In[33]:


dn_cats = dn_cats.merge(gene_taus, on="gene_name").drop_duplicates()
print(len(dn_cats))
dn_cats.head()


# In[34]:


ref_expr = dn_cats.groupby(["gene_name", "family", "dn_cat", "dev_tau",
                            "gtex_tau", "gtex_ds_tau"])["tf1p0_id"].agg("count").reset_index()
ref_expr = ref_expr.pivot(index="gene_name",
                          columns="dn_cat", values="tf1p0_id")
ref_expr.fillna(0, inplace=True)


# In[35]:


def categorize_gene(row):
    if row.DN > 0 and row.rewire > 0:
        return "both"
    elif row.DN > 0:
        return "DN"
    elif row.rewire > 0:
        return "rewire"
    elif row.NA > 0:
        return "NA"
    
ref_expr["gene_cat"] = ref_expr.apply(categorize_gene, axis=1)
ref_expr.reset_index(inplace=True)
ref_expr = ref_expr.merge(dn_cats[["gene_name", "family", "dev_tau", "gtex_tau", "gtex_ds_tau"]],
                          on="gene_name").drop_duplicates()
print(len(ref_expr))
ref_expr.sample(5)


# In[36]:


nice_boxplot(ref_expr, "dev_tau", "gene_cat", pal, ["rewire", "DN", "both", "NA"], 
            [1.04, 1.11, 1.17, 1.02], 0.35, "", ["rewirer", "negative regulator", "both", "NA"], 
            "gene-level tissue specificity (tau)", False, (0.3, 1.23), 
            "developmental gene expression\nclassified TF genes", 
            "../figures/DN_DevTau_Gene_Boxplot.pdf")


# In[39]:


nice_boxplot(ref_expr, "gtex_ds_tau", "gene_cat", pal, ["rewire", "DN", "both", "NA"], 
            [1.05, 1.11, 1.17, 1.02], 0.45, "", ["rewirer", "negative regulator", "both", "NA"], 
            "gene-level tissue specificity (tau)", False, (0.4, 1.23), 
            "GTEx gene expression\nclassified TF genes", 
            "../figures/DN_GTExDsTau_Gene_Boxplot.pdf")


# In[40]:


nice_boxplot(ref_expr, "gtex_tau", "gene_cat", pal, ["rewire", "DN", "both", "NA"], 
            [1.05, 1.11, 1.17, 1.02], 0.45, "", ["rewirer", "negative regulator", "both", "NA"], 
            "gene-level tissue specificity (tau)", False, (0.4, 1.23), 
            "GTEx gene expression\nclassified TF genes", 
            "../figures/DN_GTExTau_Gene_Boxplot.pdf")


# In[41]:


def developmental_tissue_expression_plot(gene_name, figsize, ylim, means, cols, fig_suffix):
    locs = [x for x in list(means.index) if x.split("|")[0] == gene_name]
    n_isos = len(means.loc[locs])
    palette = sns.color_palette("Spectral", as_cmap=False, n_colors=n_isos)
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
    axes[1].legend(loc='lower left', bbox_to_anchor=(1, 0))
    axes[0].axhline(y=1, color='black', linewidth=0.5, linestyle="dashed")
    plt.subplots_adjust(hspace=0.25)
    plt.savefig('../figures/expression_' + gene_name + '_' + fig_suffix + '.pdf',
                bbox_inches='tight')


# In[42]:


notestis_cols = [x for x in means_dev.columns if "testis" not in x]
notestis_cols = [x for x in notestis_cols if "_dev" not in x]
notestis_cols = [x for x in notestis_cols if "max_" not in x]
notestis_cols = [x for x in notestis_cols if "ovary" not in x]
notestis_cols = [x for x in notestis_cols if "brain" not in x]
developmental_tissue_expression_plot("PKNOX1", (5, 1.75), (0, 6), means_dev, notestis_cols, 
                                     "means_dev_notestis")


# In[45]:


notestis_cols = [x for x in means_dev.columns if "testis" not in x]
notestis_cols = [x for x in notestis_cols if "_dev" not in x]
notestis_cols = [x for x in notestis_cols if "max_" not in x]
notestis_cols = [x for x in notestis_cols if "ovary" not in x]
notestis_cols = [x for x in notestis_cols if "brain" not in x]
developmental_tissue_expression_plot("PKNOX1", (7, 1.75), (0, 6), means_dev, notestis_cols, 
                                     "means_dev_notestis_large")


# In[44]:


notestis_cols = [x for x in means_gtex.columns if "Testis" not in x]
notestis_cols = [x for x in notestis_cols if "_gtex" not in x]
notestis_cols = [x for x in notestis_cols if "max_" not in x]
notestis_cols = [x for x in notestis_cols if "Ovary" not in x]
notestis_cols = [x for x in notestis_cols if "Brain" not in x]
developmental_tissue_expression_plot("PKNOX1", (5, 1.75), (0, 6), means_gtex, notestis_cols, 
                                     "means_gtex_notestis")


# In[47]:


notestis_cols = [x for x in means_dev.columns if "testis" not in x]
notestis_cols = [x for x in notestis_cols if "_dev" not in x]
notestis_cols = [x for x in notestis_cols if "max_" not in x]
# notestis_cols = [x for x in notestis_cols if "ovary" not in x]
# notestis_cols = [x for x in notestis_cols if "brain" not in x]
developmental_tissue_expression_plot("KLF7", (8, 1.75), (0, 6), means_dev, notestis_cols, 
                                     "means_dev_notestis_large")


# In[49]:


notestis_cols = [x for x in means_gtex.columns if "Testis" not in x]
notestis_cols = [x for x in notestis_cols if "_dev" not in x]
notestis_cols = [x for x in notestis_cols if "max_" not in x]
# notestis_cols = [x for x in notestis_cols if "ovary" not in x]
# notestis_cols = [x for x in notestis_cols if "brain" not in x]
developmental_tissue_expression_plot("KLF7", (8, 1.75), (0, 6), means_gtex, notestis_cols, 
                                     "means_gtex_notestis_large")


# In[50]:


fig, ax = plt.subplots(figsize=(4.5, 0.75))

tfs["PKNOX1"].exon_diagram(ax=ax)

fig.savefig("../figures/PKNOX1_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[51]:


fig, ax = plt.subplots(figsize=(3, 1))

tfs["PKNOX1"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)

fig.savefig("../figures/PKNOX1_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[56]:


fig, ax = plt.subplots(figsize=(4.5, 1.5))

tfs["KLF7"].exon_diagram(ax=ax)

fig.savefig("../figures/KLF7_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[57]:


fig, ax = plt.subplots(figsize=(4.5, 1.5))

tfs["KLF7"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)

fig.savefig("../figures/KLF7_protein_diagram.pdf", bbox_inches="tight", dpi="figure")

