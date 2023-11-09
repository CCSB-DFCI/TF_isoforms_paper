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


# In[51]:


df = pd.read_excel('../data/internal/TFiso_LLPS_scores_HEK_20230909_V6.xlsx')
df['gene_symbol'] = df['isoform_acc'].map(lambda x: x.split('|')[0])
df['condensates_observed'] = df['Cond_Score'].map({1: True, 0: False})
df['is_cloned_reference'] = df['Ref_isoform'].map({'Reference': True,
                                                   'Alternative': False})
df['HEK_Condensate'] = df['HEK_Condensate'].str.upper().str.strip()
df['HEK_Condensate'] = df['HEK_Condensate'].map(lambda x: {'BOTH(MOST NC)': 'BOTH'}.get(x, x))
if df['isoform_acc'].duplicated().any():
    raise UserWarning('unexpected duplicates')

pairs = load_ref_vs_alt_isoforms_table()
df = df.set_index('isoform_acc')
for x in ['ref', 'alt']:
    for var in ['condensates_observed', 'HEK_Condensate', 'Loc_Kaia']:
        pairs[var + '_' + x] = pairs['clone_acc_' + x].map(df[var])
pairs['condensate_cat'] = pairs['clone_acc_alt'].map(df['Mutation_Class'])
pairs = pairs.loc[pairs['condensates_observed_ref'].notnull(), :]
pairs['condensate_cat_merged'] = pairs['condensate_cat'].map({
    'Unchanged': 'No difference',
    'LOC': 'Difference',
    'GOC': 'Difference',
    'Changed localization': 'Difference',
    })
pairs.loc[(pairs['n_positive_PPI_ref'] == 0) | (pairs['n_positive_PPI_alt'] == 0),
          'PPI_jaccard'] = np.nan


# In[3]:


df.isnull().sum()


# In[5]:


df['Mutation_Class'].value_counts()


# In[4]:


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


# In[9]:


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


# In[10]:


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


# In[11]:


print('all:')
print(df['HEK_Condensate'].value_counts())
print('\nreference:')
print(df.loc[df['is_cloned_reference'], 'HEK_Condensate'].value_counts())
print('\nalternative')
print(df.loc[~df['is_cloned_reference'], 'HEK_Condensate'].value_counts())


# In[19]:


pairs['condensate_cat'].value_counts()


# In[21]:


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


# In[23]:


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


# In[24]:


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


# In[25]:


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


# In[27]:


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


# In[28]:


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


# In[16]:


def detailed_condensate_cat(row):
    a = row['HEK_Condensate_ref']
    if pd.isnull(a):
        a = 'None'
    b = row['HEK_Condensate_alt']
    if pd.isnull(b):
        b = 'None'
    return '{} -> {}'.format(a, b)

pairs['condensate_cat_detailed'] = pairs.apply(detailed_condensate_cat, axis=1)


# In[30]:


pairs.sort_values('activation_fold_change_log2').head()


# In[31]:


pairs.sort_values('activation_fold_change_log2', ascending=False).head()


# In[18]:


pairs.sort_values('gene_symbol').to_csv('../output/TF-iso_condensates_ref-vs-alt.tsv', sep='\t', index=False)


# In[32]:


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

# In[5]:


# look at NLS and NES
from data_loading import load_annotated_TFiso1_collection


tfs = load_annotated_TFiso1_collection()
tf = tfs['ATF2']
nls = pd.concat(tf.aa_feature_disruption(tf.cloned_reference_isoform.name) for tf in tfs.values())
nls = nls.loc[nls['category'] == 'UniProt motif', :]
nls['clone_acc_alt'] = nls['alt_iso'].map({iso.name: iso.clone_acc for tf in tfs.values() for iso in tf.cloned_isoforms})
nls['type'] = nls['accession'].apply(lambda x: x.split('_')[0])
nls['affected'] = (nls['deletion'] + nls['insertion'] + nls['frameshift']) > 0 
pairs['NLS_affected'] = pairs['clone_acc_alt'].map(nls.loc[nls['type'] == 'NLS', :].groupby('clone_acc_alt')['affected'].any())
pairs['NES_affected'] = pairs['clone_acc_alt'].map(nls.loc[nls['type'] == 'NES', :].groupby('clone_acc_alt')['affected'].any())


# In[53]:


pairs['NLS_affected'].value_counts()


# In[40]:


pairs['NES_affected'].value_counts()


# In[67]:


pairs.loc[pairs['NLS_affected'].notnull(),
          ['clone_acc_alt',
           'NLS_affected',
           'NES_affected',
           'HEK_Condensate_ref', 'HEK_Condensate_alt',
           'condensate_cat',
           'Loc_Kaia_ref', 'Loc_Kaia_alt']].sort_values('NLS_affected')


# In[65]:


pairs.loc[pairs['NES_affected'].notnull(),
          ['clone_acc_alt',
           'NES_affected',
           'NLS_affected',
           'HEK_Condensate_ref', 'HEK_Condensate_alt',
           'condensate_cat',
           'Loc_Kaia_ref', 'Loc_Kaia_alt']].sort_values('NES_affected')


# In[44]:


stats.fisher_exact([[0, 4], [2, 1]])


# In[33]:


# look for PPIs with known phase separating proteins
from data_loading import load_y2h_isoform_data
y2h = load_y2h_isoform_data()
y2h.head()


# In[34]:


llps_proteins = {'FUS', 'EWS', 'TAF15', 'DDX4',
                 'BRD4', 'MED1',
                 'TFEB',
                 'YAP',}
y2h.loc[y2h['db_gene_symbol'].isin(llps_proteins), :]


# In[35]:


# could also look for low complexity disordered regions

