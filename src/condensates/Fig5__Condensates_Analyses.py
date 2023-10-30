#!/usr/bin/env python
# coding: utf-8

# # Figure 5: Condensates Analyses

# In[1]:


import numpy as np
import pandas as pd
import seaborn as sns
import sys

from matplotlib import pyplot as plt
from scipy import stats
from statannotations.Annotator import Annotator

sys.path.append("../")

from plotting import violinplot_reflected, mimic_r_boxplot
from data_loading import load_ref_vs_alt_isoforms_table

PAPER_PRESET = {"style": "ticks", "font": "Helvetica", "context": "paper", 
                "rc": {"font.size":7,"axes.titlesize":7,
                       "axes.labelsize":7, 'axes.linewidth':0.5,
                       "legend.fontsize":6, "xtick.labelsize":6,
                       "ytick.labelsize":6, "xtick.major.size": 3.0,
                       "ytick.major.size": 3.0, "axes.edgecolor": "black",
                       "xtick.major.pad": 3.0, "ytick.major.pad": 3.0}}
PAPER_FONTSIZE = 7


# In[2]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[3]:


np.random.seed(2023)


# ## functions

# In[4]:


# def permutation_test(x, y):
#     """
#     two-sided
#     """
#     nx = x.shape[0]
#     ny = y.shape[0]
#     obs = x.mean() - y.mean()
#     merged = np.concatenate([x, y])
#     rnd = []
#     for _i in range(100000):
#         np.random.shuffle(merged)
#         rnd.append(merged[:nx].mean() - merged[nx:].mean())
#     return (min([sum(r >= obs for r in rnd), sum(r <= obs for r in rnd)]) / len(rnd)) * 2


# In[5]:


def permutation_test(sample1, sample2, num_permutations=1000, seed=None, alternative='two-sided'):
    """
    Conduct a permutation test on two samples.

    :param sample1: First sample (array-like)
    :param sample2: Second sample (array-like)
    :param num_permutations: Number of permutations to perform (int)
    :param seed: Seed for random number generator (int)
    :param alternative: Defines the alternative hypothesis. 
                        'two-sided': the distributions are not equal,
                        'less': the distribution of sample1 is less than the distribution of sample2,
                        'greater': the distribution of sample1 is greater than the distribution of sample2
    :return: p-value (float)
    """

    # Ensure reproducibility
    if seed is not None:
        np.random.seed(seed)

    # Combine the samples
    combined = np.concatenate([sample1, sample2])

    # Calculate the observed test statistic
    observed_stat = np.mean(sample1) - np.mean(sample2)

    # Perform the permutations
    count = 0
    for _ in range(num_permutations):
        np.random.shuffle(combined)
        new_sample_1 = combined[:len(sample1)]
        new_sample_2 = combined[len(sample1):]

        # Calculate the new test statistic
        new_stat = np.mean(new_sample_1) - np.mean(new_sample_2)

        # Check if the new test statistic is at least as extreme as the original
        if alternative == 'two-sided':
            count += abs(new_stat) >= abs(observed_stat)
        elif alternative == 'less':
            count += new_stat <= observed_stat
        elif alternative == 'greater':
            count += new_stat >= observed_stat
        else:
            raise ValueError("alternative must be 'two-sided', 'less', or 'greater'")

    # Calculate the p-value
    p_value = (count + 1) / (num_permutations + 1)

    return p_value


# ## 1. import data

# In[6]:


df = pd.read_excel('../../data/internal/TFiso_LLPS_scores_HEK_20230909_V6.xlsx')
df['gene_symbol'] = df['isoform_acc'].map(lambda x: x.split('|')[0])
df['condensates_observed'] = df['Cond_Score'].map({1: True, 0: False})
df['is_cloned_reference'] = df['Ref_isoform'].map({'Reference': True,
                                                   'Alternative': False})
df['HEK_Condensate'] = df['HEK_Condensate'].str.upper().str.strip()
df['HEK_Condensate'] = df['HEK_Condensate'].map(lambda x: {'BOTH(MOST NC)': 'BOTH'}.get(x, x))
if df['isoform_acc'].duplicated().any():
    raise UserWarning('unexpected duplicates')


# In[7]:


df['Mutation_Class'].value_counts()


# In[8]:


df['Loc_Kaia'].value_counts()


# ## 2. summary of data

# In[9]:


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


# In[10]:


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(2.5, 2.5))
n = df['is_cloned_reference'].sum()
(
    df.loc[df['is_cloned_reference'], 'Loc_Kaia']
    .value_counts()
    .plot.pie(autopct=lambda x: '{:.0f}\n({:.0f}%)'.format(x / 100 * n, x),
 ax=axs[0])
)
n = (~df['is_cloned_reference']).sum()
(
    df.loc[~df['is_cloned_reference'], 'Loc_Kaia']
    .value_counts()
    .plot.pie(autopct=lambda x: '{:.0f}\n({:.0f}%)'.format(x / 100 * n, x),
 ax=axs[1])
)
axs[0].set_title('Reference isoform')
axs[1].set_title('Alternative isoform')
for ax in axs:
    ax.set_ylabel('')
fig.savefig('../../figures/fig5/kaia-localisation_ref-vs-alt_pie.pdf',
            bbox_inches='tight')


# In[11]:


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(2.5, 2.5))
n = df['is_cloned_reference'].sum()
(
    df.loc[df['is_cloned_reference'], 'HEK_Condensate']
    .value_counts()
    .plot.pie(autopct=lambda x: '{:.0f}\n({:.0f}%)'.format(x / 100 * n, x),
 ax=axs[0])
)
n = (~df['is_cloned_reference']).sum()
(
    df.loc[~df['is_cloned_reference'], 'HEK_Condensate']
    .value_counts()
    .plot.pie(autopct=lambda x: '{:.0f}\n({:.0f}%)'.format(x / 100 * n, x),
 ax=axs[1])
)
axs[0].set_title('Reference isoform')
axs[1].set_title('Alternative isoform')
for ax in axs:
    ax.set_ylabel('')
fig.savefig('../../figures/fig5/condensate-localisation_ref-vs-alt_pie.pdf',
            bbox_inches='tight')


# In[12]:


fig, ax = plt.subplots(1, 1, figsize=(1.4, 1.4))
n = (~df['is_cloned_reference']).sum()
(df.loc[~df['is_cloned_reference'], 'Mutation_Class']
 .value_counts()
 .plot.pie(autopct=lambda x: '{:.0f}\n({:.0f}%)'.format(x / 100 * n, x),
 ax=ax))
ax.set_ylabel('')
ax.set_title('Differences in alternative isoforms')
fig.savefig('../../figures/fig5/condensate-change-categories_pie.pdf',
            bbox_inches='tight')


# In[13]:


print('all:')
print(df['HEK_Condensate'].value_counts())
print('\nreference:')
print(df.loc[df['is_cloned_reference'], 'HEK_Condensate'].value_counts())
print('\nalternative')
print(df.loc[~df['is_cloned_reference'], 'HEK_Condensate'].value_counts())


# ## 3. merge info with pairs data to see how condensates correlate w PPIs/PDIs/etc

# In[14]:


pairs = load_ref_vs_alt_isoforms_table()
pairs.head()


# In[15]:


df = df.set_index('isoform_acc')


# In[16]:


for x in ['ref', 'alt']:
    for var in ['condensates_observed', 'HEK_Condensate', 'Loc_Kaia']:
        pairs[var + '_' + x] = pairs['clone_acc_' + x].map(df[var])


# In[17]:


pairs['condensate_cat'] = pairs['clone_acc_alt'].map(df['Mutation_Class'])


# In[18]:


pairs = pairs.loc[pairs['condensates_observed_ref'].notnull(), :]


# In[19]:


pairs['condensate_cat'].value_counts()


# In[20]:


pairs['Loc_Kaia_ref'].value_counts()


# In[21]:


pairs['condensate_cat_merged'] = pairs['condensate_cat'].map({
    'Unchanged': 'No difference',
    'LOC': 'Difference',
    'GOC': 'Difference',
    'Changed localization': 'Difference',
    })


# In[22]:


pairs.loc[(pairs['n_positive_PPI_ref'] == 0) | (pairs['n_positive_PPI_alt'] == 0),
          'PPI_jaccard'] = np.nan


# In[23]:


var = 'PPI_jaccard'
x = pairs.loc[(pairs['condensate_cat'] == 'Unchanged')
              & pairs[var].notnull(), 
              var].values
y = pairs.loc[(pairs['condensate_cat'] != 'Unchanged')
              & pairs[var].notnull(), 
              var].values

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2, h=1.75)
sns.swarmplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              order=['No difference', 'Difference'],
              color='white',
              ax=ax,
              edgecolor="black",
              linewidth=0.5,
              clip_on=False,
              size=4)
violinplot_reflected(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
                    inner=None,
              #cut=0,
              color=sns.color_palette("Set2")[0],
              order=['No difference', 'Difference'],
              ax=ax,
              )
ax.set_ylim(0, 1)

pval = permutation_test(x, y, num_permutations=10000, seed=2023)
print(pval)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], 
                      data=pairs, x='condensate_cat_merged', 
                      y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.3f}'.format(pval)])

ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('PPI Jaccard index')
ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../../figures/fig5/PPI-Jaccard-vs-condensate-change_violinplot.pdf',
            bbox_inches='tight')


# In[24]:


var = 'PDI_jaccard'
x = pairs.loc[(pairs['condensate_cat'] == 'Unchanged')
              & pairs[var].notnull(), 
              var].values
y = pairs.loc[(pairs['condensate_cat'] != 'Unchanged')
              & pairs[var].notnull(), 
              var].values

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2, h=1.75)
sns.swarmplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              order=['No difference', 'Difference'],
              color='white',
              edgecolor='black',
              linewidth=0.5,
              ax=ax,
              clip_on=False,
              size=4)
violinplot_reflected(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
                    inner=None,
              #cut=0,
              color=sns.color_palette("Set2")[0],
              order=['No difference', 'Difference'],
              ax=ax,
              )
ax.set_ylim(0, 1)

pval = permutation_test(x, y, num_permutations=10000, seed=2023)
print(pval)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.3f}'.format(pval)])

ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('PDI Jaccard index')
ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../../figures/fig5/PDI-Jaccard-vs-condensate-change_violinplot.pdf',
            bbox_inches='tight')


# In[25]:


var = 'activation_fold_change_log2'
x = pairs.loc[(pairs['condensate_cat'] == 'Unchanged')
              & pairs[var].notnull(), 
              var].values
y = pairs.loc[(pairs['condensate_cat'] != 'Unchanged')
              & pairs[var].notnull(), 
              var].values

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2, h=1.75)
sns.swarmplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              order=['No difference', 'Difference'],
              color='white',
              edgecolor='black',
              alpha=0.5,
              linewidth=0.5,
              ax=ax,
              clip_on=False,
              size=4)
sns.boxplot(data=pairs,
                     x='condensate_cat_merged',
                     y=var,
              #cut=0,
              color=sns.color_palette("Set2")[0],
              order=['No difference', 'Difference'],
              ax=ax,
              )
mimic_r_boxplot(ax)

pval = permutation_test(x, y, num_permutations=10000, seed=2023)
print(pval)
annotator = Annotator(ax=ax, pairs=[('No difference', 'Difference')], data=pairs, x='condensate_cat_merged', y=var, order=['No difference', 'Difference'],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.3f}'.format(pval)])

ax.set_xlabel('Condensate formation between reference and alternative')
for loc in ['right', 'top', 'bottom']:
    ax.spines[loc].set_visible(False)
ax.set_ylabel('Activation log2 fold change')
ax.set_xticklabels(['No difference\n(N = {})'.format(x.shape[0]),
                    'Difference\n(N = {})'.format(y.shape[0])])
ax.xaxis.set_tick_params(length=0)
fig.savefig('../../figures/fig5/activation-vs-condensate-change_boxplot.pdf',
            bbox_inches='tight')


# ## 4. look in more granularity (cytoplasmic v nuclear, etc)

# In[26]:


var = 'PPI_jaccard'

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2.5, h=1.75)
sns.swarmplot(data=pairs,
                     x='condensate_cat',
                     y=var,
              color='white',
              edgecolor='black',
              linewidth=0.5,
              ax=ax,
              clip_on=False,
              size=4)
violinplot_reflected(data=pairs,
                     x='condensate_cat',
                     y=var,
                    inner=None,
              color=sns.color_palette("Set2")[0],
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
fig.savefig('../../figures/fig5/PPI-Jaccard-vs-condensate-change_categories_boxplot.pdf',
            bbox_inches='tight')


# In[27]:


var = 'activation_fold_change_log2'

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=2.5, h=1.75)
sns.swarmplot(data=pairs,
                     x='condensate_cat',
                     y=var,
              color='white',
              edgecolor='black',
              alpha=0.5,
              linewidth=0.5,
              ax=ax,
              clip_on=False,
              size=4)
sns.boxplot(data=pairs,
                     x='condensate_cat',
                     y=var,
              color=sns.color_palette("Set2")[0],
              ax=ax,
              )
mimic_r_boxplot(ax)

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
fig.savefig('../../figures/fig5/activation-vs-condensate-change_categories_boxplot.pdf',
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


var = 'activation_fold_change_log2'

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=5.5, h=1.75)
sns.swarmplot(data=pairs,
                     x='condensate_cat_detailed',
                     y=var,
              color='white',
              edgecolor='black',
              linewidth=0.5,
              ax=ax,
              clip_on=False,
              size=4)
sns.boxplot(data=pairs,
                     x='condensate_cat_detailed',
                     y=var,
              color=sns.color_palette("Set2")[0],
              ax=ax,
              )
mimic_r_boxplot(ax)

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
fig.savefig('../../figures/fig4/activation-vs-condensate-change_categories-detailed_boxplot.pdf',
            bbox_inches='tight')


# - CREB1 - the alternative isoform, with a small insertion, forms nuclear condensates (and doesn't bind DNA...)
# - TBX5 - alternative isoforms form cytoplasmic condensates. TBX5-3, which doesn't activate, doesn't seem to be in the nucleus...
# - ZIC3 - novel isoforms form condensates
# - **PBX1** - ref forms condensates in both nucleus and cytoplasm. alt looses them

# ## 5. are any PPIs in our y2h data well-known LLPS drivers?

# In[32]:


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


# ## 6. examine how expression correlates with condensate formation

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

means_dev["max_dev"] = means_dev.max(axis=1)
max_dev = means_dev[["max_dev"]]

max_tpm = max_gtex.join(max_dev)
max_tpm = np.log2(max_tpm+1)
max_tpm = max_tpm.reset_index()
max_tpm["clone_acc"] = max_tpm.UID.str.split(" ", expand=True)[0]


# In[39]:


max_tpm["gene_symbol"] = max_tpm["clone_acc"].str.split("|", expand=True)[0]
max_tpm_gene = max_tpm[["gene_symbol", "max_gtex", "max_dev"]]
max_tpm_gene = max_tpm_gene.groupby("gene_symbol")[["max_gtex", "max_dev"]].agg("sum").reset_index()
max_tpm_gene.columns = ["gene_symbol", "max_gtex_gene", "max_dev_gene"]


# In[40]:


pairs = pairs.merge(max_tpm[["clone_acc", "max_gtex", 
                             "max_dev"]],
                    left_on="clone_acc_ref", right_on="clone_acc")
pairs.drop("clone_acc", axis=1, inplace=True)
pairs = pairs.merge(max_tpm[["clone_acc", "max_gtex", "max_dev"]],
                    left_on="clone_acc_alt", 
                    right_on="clone_acc",
                    suffixes=("_ref", "_alt"))
pairs.drop("clone_acc", axis=1, inplace=True)
pairs = pairs.merge(max_tpm_gene, on="gene_symbol")
pairs.head()


# In[41]:


fig = plt.figure(figsize=(1.75, 1.75))

x_var = "condensates_observed_ref"
y_var = "max_dev_gene"
data = pairs[["clone_acc_ref", x_var, y_var]].drop_duplicates()
ax = sns.boxplot(data=data, x=x_var, y=y_var,
                 color=sns.color_palette("Set2")[0],
                 fliersize=0)
mimic_r_boxplot(ax)
sns.swarmplot(data=data, x=x_var, y=y_var, 
              edgecolor="black", linewidth=0.5, color=sns.color_palette("Set2")[0], ax=ax,
              size=3, alpha=0.5)

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Condensates observed")
ax.set_ylabel("Maximum log2(tpm)\nDevelopmental RNA-seq")
ax.set_title("Reference isoforms\n\n")

x = data[data[x_var] == False][y_var].values
y = data[data[x_var] == True][y_var].values
pval = permutation_test(x, y, num_permutations=10000, seed=2023)
annotator = Annotator(ax=ax, pairs=[(False, True)], 
                      data=data, x='condensates_observed_ref', 
                      y='max_dev_gene', order=[False, True],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])

fig.savefig("../../figures/fig5/Dev_expr_cond.pdf", dpi="figure", bbox_inches="tight")


# In[42]:


fig = plt.figure(figsize=(1.75, 1.75))

x_var = "condensates_observed_ref"
y_var = "max_gtex_gene"
data = pairs[["clone_acc_ref", x_var, y_var]].drop_duplicates()
ax = sns.boxplot(data=data, x=x_var, y=y_var,
                 color=sns.color_palette("Set2")[0],
                 fliersize=0)
mimic_r_boxplot(ax)
sns.swarmplot(data=data, x=x_var, y=y_var, 
              edgecolor="black", linewidth=0.5, color=sns.color_palette("Set2")[0], ax=ax,
              size=3, alpha=0.5)

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Condensates observed")
ax.set_ylabel("Maximum log2(tpm)\nGTEx RNA-seq")
ax.set_title("Reference isoforms\n\n")

x = data[data[x_var] == False][y_var].values
y = data[data[x_var] == True][y_var].values
pval = permutation_test(x, y, num_permutations=10000, seed=2023)
annotator = Annotator(ax=ax, pairs=[(False, True)], 
                      data=data, x=x_var, 
                      y=y_var, order=[False, True],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])

fig.savefig("../../figures/fig5/GTEx_expr_cond.pdf", dpi="figure", bbox_inches="tight")


# In[43]:


fig = plt.figure(figsize=(1.75, 1.75))

x_var = "condensates_observed_alt"
y_var = "max_dev_alt"
data = pairs[["clone_acc_alt", x_var, y_var]].drop_duplicates()
ax = sns.boxplot(data=data, x=x_var, y=y_var,
                 color=sns.color_palette("Set2")[0],
                 fliersize=0)
mimic_r_boxplot(ax)
sns.swarmplot(data=data, x=x_var, y=y_var, 
              edgecolor="black", linewidth=0.5, color=sns.color_palette("Set2")[0], ax=ax,
              size=3, alpha=0.5)

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Condensates observed")
ax.set_ylabel("Maximum log2(tpm)\nDevelopmental RNA-seq")
ax.set_title("Alternative isoforms\n\n")

x = data[data[x_var] == False][y_var].values
y = data[data[x_var] == True][y_var].values
pval = permutation_test(x, y, num_permutations=10000, seed=2023)
annotator = Annotator(ax=ax, pairs=[(False, True)], 
                      data=data, x=x_var, 
                      y=y_var, order=[False, True],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])

fig.savefig("../../figures/fig5/Dev_expr_cond_alt.pdf", dpi="figure", bbox_inches="tight")


# In[44]:


fig = plt.figure(figsize=(1.75, 1.75))

x_var = "condensates_observed_alt"
y_var = "max_gtex_alt"
data = pairs[["clone_acc_alt", x_var, y_var]].drop_duplicates()
ax = sns.boxplot(data=data, x=x_var, y=y_var,
                 color=sns.color_palette("Set2")[0],
                 fliersize=0)
mimic_r_boxplot(ax)
sns.swarmplot(data=data, x=x_var, y=y_var, 
              edgecolor="black", linewidth=0.5, color=sns.color_palette("Set2")[0], ax=ax,
              size=3, alpha=0.5)

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Condensates observed")
ax.set_ylabel("Maximum log2(tpm)\nGTEx RNA-seq")
ax.set_title("Alternative isoforms\n\n")

x = data[data[x_var] == False][y_var].values
y = data[data[x_var] == True][y_var].values
pval = permutation_test(x, y, num_permutations=10000, seed=2023)
annotator = Annotator(ax=ax, pairs=[(False, True)], 
                      data=data, x=x_var, 
                      y=y_var, order=[False, True],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.2f}'.format(pval)])

fig.savefig("../../figures/fig5/GTEx_expr_cond_alt.pdf", dpi="figure", bbox_inches="tight")


# In[45]:


df = pd.read_excel('../../data/external/Geiger-et-al_MCP_2012_Supplementary-Table-2.xlsx',
                   skiprows=1)
hek_avrg = df[['iBAQ HEK293_1', 'iBAQ HEK293_2', 'iBAQ HEK293_3']].mean(axis=1)
print((hek_avrg > 0).sum(), 'proteins expressed in HEK293 proteome')
hek_expressed_genes = set(df.loc[(hek_avrg > 0) & df['Gene Names'].notnull(),
       'Gene Names'].str.split(';').explode().values)
all_partners = set(y2h['db_gene_symbol'].unique())
print('of {} PPI partners, {} are expressed in HEK293 cells'.format(len(all_partners), 
      len(all_partners.intersection(hek_expressed_genes))))


# In[46]:


hek_prot = df[["Gene Names", "iBAQ HEK293_1", 'iBAQ HEK293_2', 'iBAQ HEK293_3']]
hek_prot["HEK_avrg"] = hek_prot[['iBAQ HEK293_1', 'iBAQ HEK293_2', 'iBAQ HEK293_3']].mean(axis=1)
hek_prot = hek_prot[~pd.isnull(hek_prot["Gene Names"])]
print(len(hek_prot))


# In[47]:


hek_prot["gene_names_list"] = hek_prot["Gene Names"].str.split(';')
hek_prot = hek_prot.explode("gene_names_list")
hek_prot = hek_prot[["gene_names_list", "HEK_avrg"]]


# In[48]:


from data_loading import load_annotated_gencode_tfs
tfs = load_annotated_gencode_tfs()


# In[49]:


hek_prot["is_tf"] = hek_prot["gene_names_list"].isin(list(tfs.keys()))
hek_prot.is_tf.value_counts()


# In[50]:


fig = plt.figure(figsize=(1.5, 1.75))

x_var = "is_tf"
y_var = "HEK_avrg"
ax = sns.boxplot(data=hek_prot, x=x_var, y=y_var,
                 color=sns.color_palette("Set2")[0])
mimic_r_boxplot(ax)

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Is TF?")
ax.set_ylabel("HEK expression\nproteomics data")
ax.set_title("All detected genes\n\n")

x = hek_prot[hek_prot[x_var] == False][y_var].values
y = hek_prot[hek_prot[x_var] == True][y_var].values
pval = permutation_test(x, y, num_permutations=10000, seed=2023)
print(pval)
annotator = Annotator(ax=ax, pairs=[(False, True)], 
                      data=hek_prot, x=x_var, 
                      y=y_var, order=[False, True],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.4f}'.format(pval)])

fig.savefig("../../figures/fig5/HEK_proteomics_expr_TFs.pdf", dpi="figure", bbox_inches="tight")


# In[51]:


pairs = pairs.merge(hek_prot[["gene_names_list", "HEK_avrg"]], left_on="gene_symbol",
                    right_on="gene_names_list", how="left")
pairs.sample(5)


# In[52]:


fig = plt.figure(figsize=(1.75, 1.75))

x_var = "condensates_observed_ref"
y_var = "HEK_avrg"
data = pairs[["clone_acc_ref", x_var, y_var]].drop_duplicates()
ax = sns.boxplot(data=data, x=x_var, y=y_var,
                 color=sns.color_palette("Set2")[0],
                 fliersize=0)
mimic_r_boxplot(ax)
sns.swarmplot(data=data, x=x_var, y=y_var, 
              edgecolor="black", linewidth=0.5, color=sns.color_palette("Set2")[0], ax=ax,
              size=3, alpha=0.5)

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Condensates observed")
ax.set_ylabel("HEK expression\nproteomics")
ax.set_title("Reference isoforms\n\n")

x = data[data[x_var] == False][y_var].values
x = [x for x in x if not pd.isnull(x)]
print(len(x))
y = data[data[x_var] == True][y_var].values
y = [y for y in y if not pd.isnull(y)]
print(len(y))
pval = permutation_test(x, y, num_permutations=10000, seed=2023)
print(pval)
annotator = Annotator(ax=ax, pairs=[(False, True)], 
                      data=pairs, x=x_var, 
                      y=y_var, order=[False, True],)
annotator.configure(loc='outside')
annotator.annotate_custom_annotations(['P = {:.4f}'.format(pval)])

fig.savefig("../../figures/fig5/HEK_proteomics_expr_cond.pdf", dpi="figure", bbox_inches="tight")


# In[53]:


pairs.HEK_Condensate_ref.value_counts()


# In[54]:


pairs[pairs["HEK_Condensate_ref"] == "CC"][["gene_symbol", "clone_acc_ref",
                                            "clone_acc_alt", "condensates_observed_ref",
                                            "condensates_observed_alt", "HEK_Condensate_ref",
                                            "HEK_Condensate_alt", "condensate_cat"]].head(20)


# In[55]:


hpa = pd.read_table("../../data/external/HPA_subcellular_location.tsv", sep="\t")
hpa["all_observed"] = hpa["Approved"].astype(str) + ";" + hpa["Enhanced"].astype(str) + ";" + hpa["Supported"].astype(str) + ";" + hpa["Uncertain"].astype(str)
hpa.head()


# In[56]:


def cytosolic_loc(row, col):
    cytosolic_locs_to_consider = ["Actin filaments", "Cleavage furrow", "Focal adhesion sites",
                                  "Intermediate filaments", "Centriolar satellite", "Centrosome",
                                  "Cytokinetic bridge", "Microtubule ends", "Microtubules",
                                  "Midbody", "Midbody ring", "Mitotic spindle",
                                  "Aggresome", "Cytoplasmic bodies", "Cytosol", "Rods & rings",
                                  "Mitochondria", "Endoplasmic reticulum", "Vesicles",
                                  "Endosomes", "Lipid droplets", "Lysosomes", "Peroxisomes",
                                  "Golgi apparatus", "Cell junctions", "Plasma membrane"]
    
    for loc in cytosolic_locs_to_consider:
        if loc in str(row[col]):
            return True
    return False

hpa["cyto_observed"] = hpa.apply(cytosolic_loc, col="all_observed", axis=1)
len(hpa[hpa["cyto_observed"] == True]["Gene name"].unique())


# In[57]:


hpa["cyto_observed_approved"] = hpa.apply(cytosolic_loc, col="Approved", axis=1)
len(hpa[hpa["cyto_observed_approved"] == True]["Gene name"].unique())


# In[58]:


pairs = pairs.merge(hpa, left_on="gene_symbol", right_on="Gene name", how="left")
pairs.head()


# In[59]:


dd = pairs[["gene_symbol", "clone_acc_ref",
            "condensates_observed_ref",
            "HEK_Condensate_ref",
            "Loc_Kaia_ref",
            "Approved",
            "cyto_observed_approved",
            "all_observed",
            "cyto_observed"]].drop_duplicates()

tot = dd.groupby("Loc_Kaia_ref")["clone_acc_ref"].agg("count").reset_index()
cyto = dd[dd["cyto_observed"] == True].groupby("Loc_Kaia_ref")["clone_acc_ref"].agg("count").reset_index()
cyto_perc = tot.merge(cyto, on="Loc_Kaia_ref")
cyto_perc.columns = ["Loc_Kaia_ref", "tot", "cyto_observed"]
cyto_perc["perc_cyto_observed"] = cyto_perc["cyto_observed"]/cyto_perc["tot"]*100


tot = dd[~pd.isnull(dd["Approved"])].groupby("Loc_Kaia_ref")["clone_acc_ref"].agg("count").reset_index()
cyto_app = dd[dd["cyto_observed_approved"] == True].groupby("Loc_Kaia_ref")["clone_acc_ref"].agg("count").reset_index()
cyto_perc_app = tot.merge(cyto_app, on="Loc_Kaia_ref")
cyto_perc_app.columns = ["Loc_Kaia_ref", "tot", "cyto_observed_approved"]
cyto_perc_app["perc_cyto_observed_approved"] = cyto_perc_app["cyto_observed_approved"]/cyto_perc_app["tot"]*100


# In[60]:


fig = plt.figure(figsize=(1.5, 1.75))

ax = sns.barplot(data=cyto_perc, x="Loc_Kaia_ref", y="perc_cyto_observed",
                 color=sns.color_palette("Set2")[0])

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Subcellular localization\nin our assay")
ax.set_ylabel("% where cytoplasmic loc.\nis observed in HPA")
ax.set_title("Reference isoforms\n")


fig.savefig("../../figures/fig5/hpa.pdf", dpi="figure", bbox_inches="tight")


# In[61]:


fig = plt.figure(figsize=(1.5, 1.75))

ax = sns.barplot(data=cyto_perc_app, x="Loc_Kaia_ref", y="perc_cyto_observed_approved",
                 color=sns.color_palette("Set2")[0])

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Subcellular localization\nin our assay")
ax.set_ylabel("% where cytoplasmic loc.\nis approved in HPA")
ax.set_title("Reference isoforms\n")


fig.savefig("../../figures/fig5/hpa_approved.pdf", dpi="figure", bbox_inches="tight")


# In[62]:


cyto_perc


# In[63]:


cyto_perc_app


# In[64]:


pairs[["clone_acc_ref", "cyto_observed"]].drop_duplicates().cyto_observed.value_counts()


# In[65]:


pairs[pairs["Loc_Kaia_ref"] == "cyto"][["gene_symbol", "clone_acc_ref",
                                            "condensates_observed_ref",
                                            "HEK_Condensate_ref",
                                            "all_observed",
                                            "cyto_observed"]].drop_duplicates().head(20)


# In[66]:


dd[dd["HEK_Condensate_ref"] == "CC"][["gene_symbol", "clone_acc_ref",
                                            "condensates_observed_ref",
                                            "HEK_Condensate_ref",
                                            "all_observed",
                                            "cyto_observed"]].drop_duplicates().head(20)


# In[67]:


pairs["HEK_Condensate_ref_na"] = pairs["HEK_Condensate_ref"].fillna("none")


# In[68]:


fig = plt.figure(figsize=(1.75, 1.75))

x_var = "HEK_Condensate_ref_na"
y_var = "max_dev_gene"
ax = sns.boxplot(data=pairs, x=x_var, y=y_var,
                 color=sns.color_palette("Set2")[0],
                 fliersize=0)
mimic_r_boxplot(ax)
sns.swarmplot(data=pairs, x=x_var, y=y_var, 
              edgecolor="black", linewidth=0.5, color=sns.color_palette("Set2")[0], ax=ax,
              size=3, alpha=0.5)

for loc in ['right', 'top']:
    ax.spines[loc].set_visible(False)
    
ax.set_xlabel("Condensate type")
ax.set_ylabel("Maximum log2(tpm)\nDevelopmental RNA-seq")
ax.set_title("Reference isoforms\n\n")


# In[69]:


pairs[pairs["gene_symbol"] == "PPARG"][["clone_acc_ref", "clone_acc_alt", "HEK_Condensate_ref", 
                                        "HEK_Condensate_alt", "condensate_cat_detailed"]]


# In[70]:


pairs.columns


# In[71]:


pairs[pairs["condensate_cat"] == "LOC"][["gene_symbol", "clone_acc_ref", "clone_acc_alt", "family",
                                         "PPI_jaccard", "PDI_jaccard", "activation_fold_change_log2",
                                         "condensate_cat_detailed"]]


# In[ ]:




