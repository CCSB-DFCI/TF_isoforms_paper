#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy import stats
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

from plotting import violinplot_reflected
from data_loading import load_ref_vs_alt_isoforms_table


# In[2]:


df = pd.read_excel('../data/internal/TFiso_LLPS_scores_HEK_20230909_V6.xlsx')
df['gene_symbol'] = df['isoform_acc'].map(lambda x: x.split('|')[0])
df['condensates_observed'] = df['Cond_Score'].map({1: True, 0: False})
df['is_cloned_reference'] = df['Ref_isoform'].map({'Reference': True,
                                                   'Alternative': False})
df['HEK_Condensate'] = df['HEK_Condensate'].str.upper().str.strip()
df['HEK_Condensate'] = df['HEK_Condensate'].map(lambda x: {'BOTH(MOST NC)': 'BOTH'}.get(x, x))
if df['isoform_acc'].duplicated().any():
    raise UserWarning('unexpected duplicates')


# In[3]:


df.head()


# In[4]:


df.tail()


# In[5]:


df.isnull().sum()


# In[6]:


df['Mutation_Class'].value_counts()


# In[7]:


print('tested {} isoforms of {} TF genes'.format(df['isoform_acc'].nunique(),
                df['gene_symbol'].nunique())
)
print('{} ({:.0%}) reference isoforms show condensates'.format(
    df.loc[df['is_cloned_reference'], 'condensates_observed'].sum(),
    df.loc[df['is_cloned_reference'], 'condensates_observed'].mean()))
print('{} ({:.0%}) alternative isoforms show condensates'.format(
    df.loc[~df['is_cloned_reference'], 'condensates_observed'].sum(),
    df.loc[~df['is_cloned_reference'], 'condensates_observed'].mean()))
print('{} ({:.0%}) alternative isoforms change condensate formation compared to reference'.format(
    (df.loc[~df['is_cloned_reference'], 'Mutation_Class'] != 'Unchanged').sum(),
    (df.loc[~df['is_cloned_reference'], 'Mutation_Class'] != 'Unchanged').mean()))

print()
print(df.loc[~df['is_cloned_reference'], 'Mutation_Class'].value_counts())
print()


# In[8]:


fig, axs = plt.subplots(nrows=1, ncols=2)
n = df['is_cloned_reference'].sum()
(
    df.loc[df['is_cloned_reference'], 'HEK_Condensate']
    .value_counts()
    .plot.pie(autopct=lambda x: '{:.0f} ({:.0f}%)'.format(x / 100 * n, x),
 ax=axs[0])
)
n = (~df['is_cloned_reference']).sum()
(
    df.loc[~df['is_cloned_reference'], 'HEK_Condensate']
    .value_counts()
    .plot.pie(autopct=lambda x: '{:.0f} ({:.0f}%)'.format(x / 100 * n, x),
 ax=axs[1])
)
axs[0].set_title('Reference isoform condensates')
axs[1].set_title('Alternative isoform condensates')
for ax in axs:
    ax.set_ylabel('')
fig.savefig('../figures/condensate-localisation_ref-vs-alt_pie.pdf',
            bbox_inches='tight')


# In[9]:


fig, ax = plt.subplots(1, 1)
n = (~df['is_cloned_reference']).sum()
(df.loc[~df['is_cloned_reference'], 'Mutation_Class']
 .value_counts()
 .plot.pie(autopct=lambda x: '{:.0f} ({:.0f}%)'.format(x / 100 * n, x),
 ax=ax))
ax.set_ylabel('')
ax.set_title('Differences in alternative isoforms')
fig.savefig('../figures/condensate-change-categories_pie.pdf',
            bbox_inches='tight')


# In[10]:


print('all:')
print(df['HEK_Condensate'].value_counts())
print('\nreference:')
print(df.loc[df['is_cloned_reference'], 'HEK_Condensate'].value_counts())
print('\nalternative')
print(df.loc[~df['is_cloned_reference'], 'HEK_Condensate'].value_counts())


# In[11]:


pairs = load_ref_vs_alt_isoforms_table()
pairs.head()


# In[12]:


df = df.set_index('isoform_acc')


# In[13]:


df.head()


# In[14]:


for x in ['ref', 'alt']:
    for var in ['condensates_observed', 'HEK_Condensate']:
        pairs[var + '_' + x] = pairs['clone_acc_' + x].map(df[var])


# In[15]:


pairs['condensate_cat'] = pairs['clone_acc_alt'].map(df['Mutation_Class'])


# In[16]:


pairs = pairs.loc[pairs['condensates_observed_ref'].notnull(), :]


# In[17]:


pairs.head()


# In[18]:


pairs['condensate_cat'].value_counts()


# In[19]:


pairs['condensate_cat_merged'] = pairs['condensate_cat'].map({
    'Unchanged': 'No difference',
    'LOC': 'Difference',
    'GOC': 'Difference',
    'Changed localization': 'Difference',
    })


# In[20]:


def permutation_test(x, y):
    """
    two-sided
    """
    nx = x.shape[0]
    ny = y.shape[0]
    obs = x.mean() - y.mean()
    merged = np.concatenate([x, y])
    rnd = []
    for _i in range(10000):
        np.random.shuffle(merged)
        rnd.append(merged[:nx].mean() - merged[nx:].mean())
    return (min([sum(r >= obs for r in rnd), sum(r <= obs for r in rnd)]) / len(rnd)) * 2


# In[21]:


pairs.loc[(pairs['n_positive_PPI_ref'] == 0) | (pairs['n_positive_PPI_alt'] == 0),
          'PPI_jaccard'] = np.nan


# In[22]:


var = 'PPI_jaccard'
x = pairs.loc[(pairs['condensate_cat'] == 'Unchanged')
              & pairs[var].notnull(), 
              var].values
y = pairs.loc[(pairs['condensate_cat'] != 'Unchanged')
              & pairs[var].notnull(), 
              var].values

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=4, h=3)
sns.swarmplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              order=['No difference', 'Difference'],
              color='black',
              ax=ax,
              clip_on=False,)
violinplot_reflected(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
                    inner=None,
              #cut=0,
              color='lightgreen',
              order=['No difference', 'Difference'],
              ax=ax,
              )
ax.set_ylim(0, 1)

pval = permutation_test(x, y)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])

ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('PPI Jaccard index')
ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../figures/PPI-Jaccard-vs-condensate-change_violinplot.pdf',
            bbox_inches='tight')


# In[23]:


var = 'PDI_jaccard'
x = pairs.loc[(pairs['condensate_cat'] == 'Unchanged')
              & pairs[var].notnull(), 
              var].values
y = pairs.loc[(pairs['condensate_cat'] != 'Unchanged')
              & pairs[var].notnull(), 
              var].values

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=4, h=3)
sns.swarmplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              order=['No difference', 'Difference'],
              color='black',
              ax=ax,
              clip_on=False,)
violinplot_reflected(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
                    inner=None,
              #cut=0,
              color='lightgreen',
              order=['No difference', 'Difference'],
              ax=ax,
              )
ax.set_ylim(0, 1)

pval = permutation_test(x, y)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])

ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('PDI Jaccard index')
ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../figures/PDI-Jaccard-vs-condensate-change_violinplot.pdf',
            bbox_inches='tight')


# In[24]:


var = 'activation_fold_change_log2'
x = pairs.loc[(pairs['condensate_cat'] == 'Unchanged')
              & pairs[var].notnull(), 
              var].values
y = pairs.loc[(pairs['condensate_cat'] != 'Unchanged')
              & pairs[var].notnull(), 
              var].values

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=4, h=3)
sns.swarmplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              order=['No difference', 'Difference'],
              color='black',
              ax=ax,
              clip_on=False,)
sns.boxplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              #cut=0,
              color='lightgreen',
              order=['No difference', 'Difference'],
              ax=ax,
              )

pval = permutation_test(x, y)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])

ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('Activation log2 fold change')
ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../figures/activation-vs-condensate-change_boxplot.pdf',
            bbox_inches='tight')


# In[25]:


pairs.head()


# In[26]:


# loss vs gain?

# loss should be due to loss of PPIs?
# loss shoudl lead to loss of activation?
var = 'PPI_jaccard'

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=4, h=3)
sns.swarmplot(data=pairs,
                     x='condensate_cat',
                     y=var,
              #order=['Unchanged', 'Difference'],
              color='black',
              ax=ax,
              clip_on=False,)
violinplot_reflected(data=pairs,
                     x='condensate_cat',
                     y=var,
                    inner=None,
              #cut=0,
              color='lightgreen',
              #order=['No difference', 'Difference'],
              ax=ax,
              )
ax.set_ylim(0, 1)
"""
pval = permutation_test(x, y)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])
"""
ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('PPI Jaccard index')
#ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
#                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../figures/PPI-Jaccard-vs-condensate-change_categories_boxplot.pdf',
            bbox_inches='tight')


# In[27]:


# loss vs gain?

# loss should be due to loss of PPIs?
# loss shoudl lead to loss of activation?
var = 'activation_fold_change_log2'

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=4, h=3)
sns.swarmplot(data=pairs,
                     x='condensate_cat',
                     y=var,
              #order=['Unchanged', 'Difference'],
              color='black',
              ax=ax,
              clip_on=False,)
sns.boxplot(data=pairs,
                     x='condensate_cat',
                     y=var,
              color='lightgreen',
              #order=['No difference', 'Difference'],
              ax=ax,
              )
"""
pval = permutation_test(x, y)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])
"""
ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('Activation log2 fold change')
#ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
#                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../figures/activation-vs-condensate-change_categories_boxplot.pdf',
            bbox_inches='tight')


# In[28]:


def detailed_condensate_cat(row):
    a = row['HEK_Condensate_ref']
    if pd.isnull(a):
        a = 'None'
    b = row['HEK_Condensate_alt']
    if pd.isnull(b):
        b = 'None'
    return '{} -> {}'.format(a, b)

pairs['condensate_cat_detailed'] = pairs.apply(detailed_condensate_cat, axis=1)


# In[29]:


pairs.sort_values('activation_fold_change_log2').head()


# In[30]:


pairs.sort_values('activation_fold_change_log2', ascending=False).head()


# In[31]:


# loss vs gain?

# loss should be due to loss of PPIs?
# loss shoudl lead to loss of activation?
var = 'activation_fold_change_log2'

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=9, h=3)
sns.swarmplot(data=pairs,
                     x='condensate_cat_detailed',
                     y=var,
              #order=['Unchanged', 'Difference'],
              color='black',
              ax=ax,
              clip_on=False,)
sns.boxplot(data=pairs,
                     x='condensate_cat_detailed',
                     y=var,
              color='lightgreen',
              #order=['No difference', 'Difference'],
              ax=ax,
              )
"""
pval = permutation_test(x, y)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])
"""
ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('Activation log2 fold change')
#ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
#                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0, rotation=90)
fig.savefig('../figures/activation-vs-condensate-change_categories-detailed_boxplot.pdf',
            bbox_inches='tight')


# - CREB1 - the alternative isoform, with a small insertion, forms nuclear condensates (and doesn't bind DNA...)
# - TBX5 - alternative isoforms form cytoplasmic condensates. TBX5-3, which doesn't activate, doesn't seem to be in the nucleus...
# - ZIC3 - novel isoforms form condensates
# - **PBX1** - ref forms condensates in both nucleus and cytoplasm. alt looses them

# In[32]:


# look for PPIs with known phase separating proteins
from data_loading import load_y2h_isoform_data
y2h = load_y2h_isoform_data()
y2h.head()


# In[33]:


llps_proteins = {'FUS', 'EWS', 'TAF15', 'DDX4',
                 'BRD4', 'MED1',
                 'TFEB',
                 'YAP',}
y2h.loc[y2h['db_gene_symbol'].isin(llps_proteins), :]


# In[34]:


# could also look for low complexity disordered regions


# In[35]:


from data_loading import (load_annotated_TFiso1_collection,
                          load_developmental_tissue_expression_remapped,
                          load_gtex_remapped)


# In[36]:


df_gtex, metadata_gtex, genes_gtex = load_gtex_remapped()

exclusion_list_gtex = {'Cells - Leukemia cell line (CML)',
                       'Cells - EBV-transformed lymphocytes',
                       'Cells - Cultured fibroblasts'}

df_gtex = df_gtex.loc[:, ~df_gtex.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex = metadata_gtex.loc[~metadata_gtex['body_site'].isin(exclusion_list_gtex), :]

means_gtex = df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean()


# In[37]:


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


# In[38]:


means_gtex["max_gtex"] = means_gtex.max(axis=1)
max_gtex = means_gtex[["max_gtex"]] 
max_gtex.head()


# In[39]:


means_dev["max_dev"] = means_dev.max(axis=1)
max_dev = means_dev[["max_dev"]] 
max_dev.head()


# In[40]:


max_tpm = max_gtex.join(max_dev)
max_tpm = np.log2(max_tpm+1)
max_tpm = max_tpm.reset_index()
max_tpm["clone_acc"] = max_tpm.UID.str.split(" ", expand=True)[0]
max_tpm.head()


# In[41]:


pairs = pairs.merge(max_tpm[["clone_acc", "max_gtex", 
                             "max_dev"]],
                    left_on="clone_acc_ref", right_on="clone_acc").merge(max_tpm[["clone_acc", "max_gtex", 
                                                                                  "max_dev"]],
                                                                         left_on="clone_acc_alt", 
                                                                         right_on="clone_acc",
                                                                        suffixes=("_ref", "_alt"))
pairs.head()


# In[48]:


ax = sns.boxplot(data=pairs, x="condensate_cat", y="max_gtex_ref")
sns.swarmplot(data=pairs, x="condensate_cat", y="max_gtex_ref", edgecolor="black", linewidth=1, ax=ax)


# In[50]:


ax = sns.boxplot(data=pairs, x="condensate_cat", y="max_dev_ref")
sns.swarmplot(data=pairs, x="condensate_cat", y="max_dev_ref", edgecolor="black", linewidth=1, ax=ax)


# In[49]:


ax = sns.boxplot(data=pairs, x="condensate_cat", y="max_gtex_alt")
sns.swarmplot(data=pairs, x="condensate_cat", y="max_gtex_alt", edgecolor="black", linewidth=1, ax=ax)


# In[51]:


ax = sns.boxplot(data=pairs, x="condensate_cat", y="max_dev_alt")
sns.swarmplot(data=pairs, x="condensate_cat", y="max_dev_alt", edgecolor="black", linewidth=1, ax=ax)


# In[57]:


max_tpm["gene_symbol"] = max_tpm["clone_acc"].str.split("|", expand=True)[0]
max_tpm_gene = max_tpm[["gene_symbol", "max_gtex", "max_dev"]]
max_tpm_gene = max_tpm_gene.groupby("gene_symbol")[["max_gtex", "max_dev"]].agg("sum").reset_index()
max_tpm_gene.columns = ["gene_symbol", "max_gtex_gene", "max_dev_gene"]
max_tpm_gene.head()


# In[58]:


pairs = pairs.merge(max_tpm_gene, on="gene_symbol")
pairs.head()


# In[59]:


pairs.columns


# In[65]:


ax = sns.boxplot(data=pairs, x="condensates_observed_ref", y="max_dev_gene")
sns.swarmplot(data=pairs, x="condensates_observed_ref", y="max_dev_gene", edgecolor="black", linewidth=1, ax=ax)


# In[66]:


ax = sns.boxplot(data=pairs, x="condensates_observed_alt", y="max_dev_gene")
sns.swarmplot(data=pairs, x="condensates_observed_alt", y="max_dev_gene", edgecolor="black", linewidth=1, ax=ax)


# In[68]:


ax = sns.boxplot(data=pairs, x="condensates_observed_ref", y="max_dev_ref")
sns.swarmplot(data=pairs, x="condensates_observed_ref", y="max_dev_ref", edgecolor="black", linewidth=1, ax=ax)


# In[67]:


ax = sns.boxplot(data=pairs, x="condensates_observed_alt", y="max_dev_alt")
sns.swarmplot(data=pairs, x="condensates_observed_alt", y="max_dev_alt", edgecolor="black", linewidth=1, ax=ax)


# In[69]:


pairs.sort_values(by="max_dev_alt", ascending=True)[["gene_symbol", "clone_acc_ref", "clone_acc_alt", 
                                                     "condensates_observed_ref", "HEK_Condensate_ref",
                                                     "condensates_observed_alt", "HEK_Condensate_alt",
                                                     "max_dev_ref", "max_dev_alt"]]


# In[71]:


df = pd.read_excel('../data/external/Geiger-et-al_MCP_2012_Supplementary-Table-2.xlsx',
                   skiprows=1)
hek_avrg = df[['iBAQ HEK293_1', 'iBAQ HEK293_2', 'iBAQ HEK293_3']].mean(axis=1)
print((hek_avrg > 0).sum(), 'proteins expressed in HEK293 proteome')
hek_expressed_genes = set(df.loc[(hek_avrg > 0) & df['Gene Names'].notnull(),
       'Gene Names'].str.split(';').explode().values)
all_partners = set(y2h['db_gene_symbol'].unique())
print('of {} PPI partners, {} are expressed in HEK293 cells'.format(len(all_partners), 
      len(all_partners.intersection(hek_expressed_genes))))


# In[84]:


hek_prot = df[["Gene Names", "iBAQ HEK293_1", 'iBAQ HEK293_2', 'iBAQ HEK293_3']]
hek_prot["HEK_avrg"] = hek_prot[['iBAQ HEK293_1', 'iBAQ HEK293_2', 'iBAQ HEK293_3']].mean(axis=1)
hek_prot = hek_prot[~pd.isnull(hek_prot["Gene Names"])]
print(len(hek_prot))
hek_prot.sample(5)


# In[85]:


hek_prot["gene_names_list"] = hek_prot["Gene Names"].str.split(';')
hek_prot = hek_prot.explode("gene_names_list")
hek_prot = hek_prot[["gene_names_list", "HEK_avrg"]]
hek_prot.sample(5)


# In[87]:


from data_loading import load_annotated_gencode_tfs
tfs = load_annotated_gencode_tfs()


# In[91]:


hek_prot["is_tf"] = hek_prot["gene_names_list"].isin(list(tfs.keys()))
hek_prot.is_tf.value_counts()


# In[92]:


sns.boxplot(data=hek_prot, x="is_tf", y="HEK_avrg")


# In[94]:


pairs = pairs.merge(hek_prot[["gene_names_list", "HEK_avrg"]], left_on="gene_symbol",
                    right_on="gene_names_list", how="left")
pairs.sample(5)


# In[95]:


ax = sns.boxplot(data=pairs, x="condensates_observed_ref", y="HEK_avrg")
sns.swarmplot(data=pairs, x="condensates_observed_ref", y="HEK_avrg", edgecolor="black", linewidth=1, ax=ax)


# In[96]:


pairs["HEK_avrg_incl0"] = pairs["HEK_avrg"].fillna(0)


# In[97]:


ax = sns.boxplot(data=pairs, x="condensates_observed_ref", y="HEK_avrg_incl0")
sns.swarmplot(data=pairs, x="condensates_observed_ref", y="HEK_avrg_incl0", edgecolor="black", linewidth=1, ax=ax)


# In[ ]:




