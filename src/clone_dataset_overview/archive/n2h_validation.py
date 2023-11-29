#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

from plotting import validation_plot, validation_titration_plot
from data_loading import (load_n2h_ppi_validation_data,
                          load_full_y2h_data_including_controls)


# In[2]:


df = load_n2h_ppi_validation_data()


# In[3]:


# TODO: remove this once everything finalized 
# (they should already be removed in the table)

# sequence confirmation was done after the expeirment design
# so need to remove sequence failures
y2h = load_full_y2h_data_including_controls()
y2h = y2h.loc[y2h['category'] == 'tf_isoform_ppis', :]
y2h_positives = y2h.loc[y2h['Y2H_result'] == True, ['ad_orf_id', 'db_orf_id']].values
y2h_positives = set(map(tuple, y2h_positives))
y2h_negatives = y2h.loc[y2h['Y2H_result'] == False, ['ad_orf_id', 'db_orf_id']].values
y2h_negatives = set(map(tuple, y2h_negatives))

# check positives / negatives are in Y2H dataset
def in_y2h_positves(row):
    pair = (row['test_orf_idb'],
            row['test_orf_ida'])
    return pair in y2h_positives


def in_y2h_negatives(row):
    pair = (row['test_orf_idb'],
            row['test_orf_ida'])
    return pair in y2h_negatives


df = df.loc[~((df['source'] == 'isoform positives') & 
            ~df.apply(in_y2h_positves, axis=1)), :]
df = df.loc[~((df['source'] == 'isoform negatives') & 
            ~df.apply(in_y2h_negatives, axis=1)), :]


# In[4]:


COLOR_LIT = (60 / 255, 134 / 255, 184 / 255)
COLOR_HURI = (155 / 255, 97 / 255, 153 / 255)
colors = {    
    'vignettes': 'yellow', 
    'isoform positives': COLOR_HURI,
     'RRS - TF space specific': 'tab:red',
      'Lit-BM - TF space specific': COLOR_LIT,
       'isoform negatives': 'grey',
        'RRS - from HuRI': 'tab:red',
        'Lit-BM-13': COLOR_LIT,
        'PRS - hPRS-v2': COLOR_LIT,
        'RRS - hRRS-v2': 'tab:red'}


# In[5]:


# NOTE: the plate numbers here for the controls for the PRS/RRS are incorrect
fig, axs = plt.subplots(3, 1)
fig.set_size_inches(h=12, w=12)
sns.stripplot(data=df, y='score_pair_log10', x='test_pla', hue='source', ax=axs[0])
sns.stripplot(data=df, y='score_empty-N1_log10', x='test_pla', hue='source', ax=axs[1])
sns.stripplot(data=df, y='score_empty-N2_log10', x='test_pla', hue='source', ax=axs[2])
for ax in axs:
    ax.set_ylim(1.5, 7)
axs[0].legend(bbox_to_anchor=(1., 1.))
axs[1].get_legend().remove()
axs[2].get_legend().remove()
axs[0].set_xlabel('')
axs[1].set_xlabel('')
axs[2].set_xlabel('Plate number')
fig.savefig('../figures/TFv02_scores-by-plate-and-category.pdf',
            bbox_inches='tight')


# In[6]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(h=4, w=12)
sns.stripplot(data=df, y='log2 NLR', x='test_pla', hue='source', ax=ax)
ax.legend(bbox_to_anchor=(1., 1.))
ax.set_xlabel('Plate number')
fig.savefig('../figures/TFv02_NLR-by-plate-and-category.pdf',
            bbox_inches='tight')


# In[7]:


df['source'].value_counts()


# In[8]:


df['source'].unique()


# In[9]:


sources = [
           'PRS - hPRS-v2', 
            'Lit-BM-13', 
            'Lit-BM - TF space specific',
            'RRS - hRRS-v2',
            'RRS - from HuRI', 
            'RRS - TF space specific',
            'isoform positives', 
            'isoform negatives',
            'vignettes',
           ]
fig, ax = plt.subplots(1, 1)
sns.stripplot(data=df, 
              x='source',
              y='log2 NLR', 
              order=sources,
              ax=ax,
              alpha=0.5,
              palette=colors)
ax.tick_params(axis="x", rotation=90)
ax.axhline(y=df.loc[df['source'] == 'RRS - hRRS-v2', 'log2 NLR'].max(), color='black')
ax.set_xlabel('')
fig.savefig('../figures/TFv02_categories_points.pdf',
            bbox_inches='tight')


# In[10]:


sources = [
           'PRS - hPRS-v2', 
            'Lit-BM-13', 
            'Lit-BM - TF space specific',
            'RRS - hRRS-v2',
            'RRS - from HuRI', 
            'RRS - TF space specific',
            'isoform positives', 
            'isoform negatives',
           ]
line_styles = ['-', '--', ':', '-', '--', ':', '-', '-']
fig, ax = plt.subplots(1, 1)
validation_titration_plot(data=df, 
                          selections=[df['source'] == x for x in sources],
                          labels=sources,
                          colors=[colors[x] for x in sources],
                          line_styles=line_styles,
                          score_column='log2 NLR',
                          threshold=df.loc[df['source'] == 'RRS - hRRS-v2', 'log2 NLR'].max(),
                          xmin=3,
                          ax=ax)
ax.set_xlabel('Log2 NLR threshold')
ax.set_yticklabels([f'{x:.0%}' for x in ax.get_yticks()])
fig.savefig('../figures/TFv02_titration.pdf',
            bbox_inches='tight')


# In[11]:


df['source'].unique()


# In[12]:


# bar chart
df['result'] = df['NLR'] > df.loc[df['source'] == 'RRS - hRRS-v2', 'NLR'].max()

fig, ax = plt.subplots(1, 1)
validation_plot(data=df,
                selections=[df['source'] == x for x in sources],
                labels=[str(x) for x in sources],
                colors=[colors[x] for x in sources],
                result_column='result',
                errorbar_capsize=0.1,
                y_max=0.41,
                xlabel_rotation=90)
ax.set_yticklabels([f'{x:.0%}' for x in ax.get_yticks()])
fig.savefig('../figures/TFv02_bar.pdf',
            bbox_inches='tight')


# In[13]:


stats.fisher_exact([[27, 133-27], [5, 133-5]])


# In[14]:


fig, axs = plt.subplots(3, 3)
fig.set_size_inches(12, 12)
s_min = 1.5
s_max = 7
cutoff = df.loc[df['source'] == 'RRS - hRRS-v2', 'log2 NLR'].max()
for source, ax in zip(sources, axs.flatten()):
    ax.scatter(df.loc[df['source'] == source, ['score_empty-N1_log10', 'score_empty-N2_log10']].max(axis=1).values,
               df.loc[df['source'] == source, 'score_pair_log10'].values,
               color=colors[source],
               alpha=0.5)
    ax.set_title(source)
    ax.set_xlim(s_min, s_max)
    ax.set_ylim(s_min, s_max)
    ax.plot([s_min, s_max],
            [s_min + np.log10(2 ** cutoff), s_max + np.log10(2 ** cutoff)],
             color='black')
for ax in axs[-1, :]:
    ax.set_xlabel('max empty control – Log 10 luminescence')
for ax in axs[:, 0]:
    ax.set_ylabel('pair – Log 10 luminescence')
axs[-1, -1].axis('off')
axs[-2, -1].set_xlabel('max empty control – Log 10 luminescence')
fig.savefig('../figures/TFv02_control-vs-pair_scatter.pdf',
            bbox_inches='tight')


# In[15]:


# pair plot
pairs = []
for _i, row_p in df.loc[df['source'] == 'isoform positives', :].iterrows():
    tf_gene = row_p['gene_symbol_tf']
    for _j, row_n in df.loc[(df['source'] == 'isoform negatives') &
                            (df['gene_symbol_tf'] == tf_gene) &
                            (df['test_orf_ida'] == row_p['test_orf_ida']), :].iterrows():
        pairs.append((tf_gene, 
                      row_p['test_orf_idb'],
                      row_n['test_orf_idb'],
                      row_n['test_orf_ida'],
                      row_p['clone_acc'],
                       row_n['clone_acc'],
                       row_n['gene_symbol_partner'],
                      row_p['log2 NLR'],
                      row_n['log2 NLR']))
pairs = pd.DataFrame(pairs, columns=['TF gene', 
                                     'Y2H_positive_orf_id',
                                     'Y2H_negative_orf_id',
                                     'partner_orf_id',
                                     'Y2H_positive_iso_acc',
                                     'Y2H_negative_iso_acc',
                                     'partner gene',
                                     'log2_NLR_Y2H_positive',
                                     'log2_NLR_Y2H_negative'])
fig, ax = plt.subplots(1, 1)
for _i, pair in pairs.iterrows():
    ax.plot([0, 1], pair[['log2_NLR_Y2H_positive', 'log2_NLR_Y2H_negative']].values, '-o', color='grey', alpha=0.5)
ax.axhline(y=cutoff, color='black', linewidth=0.75)


# In[16]:


pairs.sort_values('log2_NLR_Y2H_negative', ascending=False).head()


# In[17]:


vignettes = df.loc[df['source'] == 'vignettes', :].copy()


# In[18]:


vignettes.sort_values(['gene_symbol_tf', 'gene_symbol_partner', 'clone_acc'])


# In[ ]:




