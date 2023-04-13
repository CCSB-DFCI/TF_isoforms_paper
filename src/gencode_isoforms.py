
# coding: utf-8

# In[33]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import functools

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

from data_loading import load_annotated_gencode_tfs


# In[6]:


tfs = load_annotated_gencode_tfs()


# In[3]:


print(len(tfs), 'TFs')


# In[3]:


n_iso = pd.Series({name: len(tf.orfs) for name, tf in tfs.items()})
n_iso.head()


# In[17]:


print(f'TF with most isoforms is {n_iso.idxmax()} with {n_iso.max()} isoforms')


# In[16]:


# add fraction axis on rhs
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=6, h=2)
xs = range(1, n_iso.max() + 1)
ax.bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])
ax.set_xticks(xs)
for pos in ['top', 'right', 'bottom']:
    ax.spines[pos].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_xlabel('Unique protein isoforms per gene')

def num2pct(y):
    return (y / n_iso.shape[0]) * 100

def pct2num(y):
    return (y / 100) * n_iso.shape[0]


ax.set_ylim(0, 800)
ax.set_yticks(range(0, 800, 100), minor=True)
ax.set_ylabel('TF genes')
for pos in ['top', 'right', 'bottom']:
    ax.spines[pos].set_visible(False)
pctax = ax.secondary_yaxis('right', functions=(num2pct, pct2num))
pctax.set_ylabel('')
pctax.set_yticks(range(0, 46, 5), minor=True)
pctax.set_yticklabels([f'{y:.0f}%' for y in pctax.get_yticks()])
fig.savefig('../figures/n-isoforms-per-gene_GENCODE-v30-TFs_hist.pdf',
            bbox_inches='tight')


# In[58]:


tfdb.loc[tfdb['DBD'] == 'Myb/SANT', :]


# In[56]:


tfdb['DBD'].value_counts().tail(60)


# In[36]:


from data_loading import load_human_tf_db

tfdb = load_human_tf_db()
non_nr_ligand_tfs = {'AHR'} # Hand-curated

fig, ax = plt.subplots(1, 1, sharex=True)
fig.set_size_inches(h=2, w=6)
xs = range(1, df['n_isoforms'].max() + 1)

def num2pct(y, total):
    return (y / total) * 100

def pct2num(y, total):
    return (y / 100) * total

n_iso = df.loc[df['gene'].isin(non_nr_ligand_tfs), 'n_isoforms']
ax.bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])
ax.set_title('Non-nuclear receptor, receptor TFs')
for pos in ['top', 'right']:
    ax.spines[pos].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_ylabel('Genes')

n = df['gene'].isin(non_nr_ligand_tfs).sum()
pctax = ax.secondary_yaxis('right', 
                            functions=(functools.partial(num2pct, total=n),
                                        functools.partial(pct2num, total=n)))
ax.set_ylim(0, 0.75 * n)
p = min(
    (stats.mannwhitneyu(df.loc[df['gene'].isin(non_nr_ligand_tfs), 'n_isoforms'].values,
                        df.loc[~df['gene'].isin(non_nr_ligand_tfs), 'n_isoforms'].values).pvalue
                        ),
                        1
                        )
ax.text(s=f'mean = {n_iso.mean():.2f} isoforms\np = {p:.2E}',
        x=10,
        y=0.5 * n)
pctax.set_yticks(range(0, 80, 20))
pctax.set_yticks(range(0, 80, 5), minor=True)
pctax.set_yticklabels([f'{y:.0f}%' for y in pctax.get_yticks()])

ax.set_xticks(xs)
ax.set_xlabel('Unique protein isoforms per gene')


fig.savefig('../figures/n-isoforms-per-gene_GENCODE-v30-non-NR-receptor-TFs_hist.pdf',
            bbox_inches='tight')


# In[42]:


print(' '.join(non_nr_ligand_tfs))


# In[41]:


tfdb.loc[tfdb['HGNC symbol'] == 'SREBF1']


# In[7]:


# number of isoforms per family
df = pd.DataFrame([(tf.name, len(tf.orfs), tf.tf_family) for name, tf in tfs.items()],
                  columns=['gene', 'n_isoforms', 'family'])
df.head()


# In[8]:


fam_size = df['family'].value_counts()
small_families = set(fam_size[fam_size < 20].index.values)
df['family_merged'] = df['family']
df.loc[df['family'].isin(small_families) | (df['family'] == 'Unknown'), 'family_merged'] = 'other'
df.head()


# In[80]:


# p-value
# make look nice
families = list(df['family_merged'].value_counts().index.values)
families.remove('other')
families += ['other']

fig, axs = plt.subplots(len(families), 1, sharex=True)
fig.set_size_inches(h=2 * len(families), w=6)
xs = range(1, df['n_isoforms'].max() + 1)

def num2pct(y, total):
    return (y / total) * 100

def pct2num(y, total):
    return (y / 100) * total

for fam, ax in zip(families, axs):
    n_iso = df.loc[df['family_merged'] == fam, 'n_isoforms']
    ax.bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])
    ax.set_title(fam)
    for pos in ['top', 'right']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    ax.set_ylabel('Genes')

    n = (df['family_merged'] == fam).sum()
    pctax = ax.secondary_yaxis('right', 
                               functions=(functools.partial(num2pct, total=n),
                                          functools.partial(pct2num, total=n)))
    ax.set_ylim(0, 0.75 * n)
    p = min(
        (stats.mannwhitneyu(df.loc[df['family_merged'] == fam, 'n_isoforms'].values,
                           df.loc[df['family_merged'] != fam, 'n_isoforms'].values).pvalue
                           * len(families)
                           ),
                           1
                           )
    ax.text(s=f'mean = {n_iso.mean():.2f} isoforms\np = {p:.2E}',
            x=10,
            y=0.5 * n)
    pctax.set_yticks(range(0, 80, 20))
    pctax.set_yticks(range(0, 80, 5), minor=True)
    pctax.set_yticklabels([f'{y:.0f}%' for y in pctax.get_yticks()])

axs[-1].set_xticks(xs)
axs[-1].set_xlabel('Unique protein isoforms per gene')


fig.savefig('../figures/n-isoforms-per-gene_GENCODE-v30-TFs_by-family_hist.pdf',
            bbox_inches='tight')


# In[74]:


df['n_isoforms'].mean()


# In[60]:


num2pct['C2H2 ZF'](5)


# In[61]:


num2pct['Ets'](5)


# In[63]:


(df['family_merged'] == 'Ets').sum()


# In[3]:


# types of splicing
# including matrix or venn or something...
tfs['ATF2'].exon_diagram()


# In[6]:


tfs['ATF2'].orfs[0]


# In[22]:


cats = pd.DataFrame([tf.splicing_categories(tf.MANE_select_isoform.name, alt_iso.name) 
              for tf in tfs.values()
              for alt_iso in tf.orfs[1:]
              if tf.has_MANE_select_isoform
              and alt_iso != tf.MANE_select_isoform])


# In[ ]:


# check at least one of something


# In[23]:


cats.mean()


# In[16]:


cats.head()


# In[25]:


cats.loc[cats['intron retention']].head()


# In[30]:


# DEBUG
tfs['FOXN2'].exon_diagram()


# In[28]:


cat_cols = cats.columns[3:]
n_cats = cats[cat_cols].sum(axis=1).value_counts().sort_index()
n_cats


# In[37]:


n_cats / n_cats.sum() * 100


# In[ ]:


cats.groupby(cat_cols)


# In[39]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(h=2.5, w=3)
(cats.mean() * 100).plot.bar(ax=ax)
ax.set_ylabel('Proportion of alternative\nprotein isoforms')
for pos in ['top', 'right', 'bottom']:
    ax.spines[pos].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_yticks(range(0, 61, 10))
ax.set_yticklabels(f'{y:d}%' for y in ax.get_yticks())
ax.set_yticks(range(0, int(ax.get_ylim()[1]), 5),
              minor=True)
fig.savefig('../figures/splice-categories_TFs-in-GENCODE-v30_bar.pdf',
            bbox_inches='tight')


# In[39]:


import venn


# In[29]:


import upsetplot


upsetplot.plot(cats.groupby(list(cat_cols)).size())
plt.savefig('../figures/isoform-categories_gencode_UpSet-plot.pdf',
            bbox_inches='tight')


# In[2]:


from data_loading import load_annotated_6k_collection

tfs_cloned = load_annotated_6k_collection()


# In[14]:


df['n_isoforms_cloned'] = (df['gene'].map({tf.name: len(tf.orfs) for tf in tfs_cloned.values()})
                            .fillna(0)
                            .astype(int))


# In[18]:


# p-value
# make look nice
families = list(df['family_merged'].value_counts().index.values)
families.remove('other')
families += ['other']

fig, axs = plt.subplots(len(families), 1, sharex=True)
fig.set_size_inches(h=2 * len(families), w=6)
xs = range(1, df['n_isoforms_cloned'].max() + 1)

def num2pct(y, total):
    return (y / total) * 100

def pct2num(y, total):
    return (y / 100) * total

for fam, ax in zip(families, axs):
    n_iso = df.loc[df['family_merged'] == fam, 'n_isoforms_cloned']
    ax.bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])
    ax.set_title(fam)
    for pos in ['top', 'right']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    ax.set_ylabel('Genes')

    n = (df.loc[df['n_isoforms_cloned'] >= 1, 'family_merged'] == fam).sum()
    pctax = ax.secondary_yaxis('right', 
                               functions=(functools.partial(num2pct, total=n),
                                          functools.partial(pct2num, total=n)))
    ax.set_ylim(0, 0.75 * n)
    p = min(
        (stats.mannwhitneyu(df.loc[df['family_merged'] == fam, 'n_isoforms_cloned'].values,
                           df.loc[df['family_merged'] != fam, 'n_isoforms_cloned'].values).pvalue
                           * len(families)
                           ),
                           1
                           )
    ax.text(s=f'mean = {n_iso.mean():.2f} isoforms\np = {p:.2E}',
            x=10,
            y=0.5 * n)
    pctax.set_yticks(range(0, 80, 20))
    pctax.set_yticks(range(0, 80, 5), minor=True)
    pctax.set_yticklabels([f'{y:.0f}%' for y in pctax.get_yticks()])

axs[-1].set_xticks(xs)
axs[-1].set_xlabel('Unique protein isoforms per gene')


fig.savefig('../figures/n-isoforms-per-gene_TFiso1-TFs_by-family_hist.pdf',
            bbox_inches='tight')


# In[22]:


# add fraction axis on rhs
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=6, h=2)
xs = range(1, df['n_isoforms_cloned'].max() + 1)
ax.bar(x=xs, height=[df['n_isoforms_cloned'].value_counts().to_dict().get(x, 0) for x in xs])
ax.set_xticks(xs)
for pos in ['top', 'right', 'bottom']:
    ax.spines[pos].set_visible(False)
ax.xaxis.set_tick_params(length=0)
ax.set_xlabel('Unique protein isoforms per gene')

def num2pct(y):
    return (y / df['n_isoforms_cloned'].shape[0]) * 100

def pct2num(y):
    return (y / 100) * df['n_isoforms_cloned'].shape[0]


ax.text(s=f'mean = {df["n_isoforms_cloned"].mean():.2f} isoforms',
        x=10,
        y=0.5 * n)

ax.set_ylim(0, 100)
ax.set_yticks(range(0, 100, 10), minor=True)
ax.set_ylabel('TF genes')
for pos in ['top', 'right', 'bottom']:
    ax.spines[pos].set_visible(False)
pctax = ax.secondary_yaxis('right', functions=(num2pct, pct2num))
pctax.set_ylabel('')
pctax.set_yticks(range(0, 46, 5), minor=True)
pctax.set_yticklabels([f'{y:.0f}%' for y in pctax.get_yticks()])
fig.savefig('../figures/n-isoforms-per-gene_TFiso1-TFs_hist.pdf',
            bbox_inches='tight')

