#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from ccsblib import paros_connection
from ccsblib import ccsbplotlib as cplt

from data_loading import (load_isoform_and_paralog_y2h_data, 
                          load_m1h_activation_data,
                          load_y1h_pdi_data,
                          load_valid_isoform_clones,
                          load_ref_vs_alt_isoforms_table)

pd.set_option('display.max_columns', 50)


# In[2]:


# G12 and H12 are positive controls 


# In[3]:


y2h = load_isoform_and_paralog_y2h_data()
iso_gte_1_pos_ppi_iso_data_only = set(y2h.loc[(y2h['category'] == 'tf_isoform_ppis') &
                          (y2h['Y2H_result'] == True), 'ad_clone_acc'].unique()) 
iso_gte_1_pos_ppi_all_data = set(y2h.loc[(y2h['Y2H_result'] == True), 'ad_clone_acc'].unique())
# restict to TF isoform data (i.e. not paralogs etc.)
ppi = y2h.loc[(y2h['category'] == 'tf_isoform_ppis'), 
              ['category',
               'ad_clone_acc',
               'ad_gene_symbol',
               'db_gene_symbol',
               'Y2H_result']].copy()
# at least one positive per PPI partner
ppi = ppi.loc[ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])
                 ['Y2H_result']
                 .transform(lambda row: (row == True).any()),
              :]
# at least one successfully tested PPI per isoform
ppi = ppi.loc[ppi.groupby('ad_clone_acc')
                  ['Y2H_result']
                  .transform(lambda x: (x.notnull().any())),
              :]
# at least two partners per isoform
ppi = ppi.loc[ppi.groupby('ad_gene_symbol')
                 ['ad_clone_acc']
                 .transform(lambda x: x.nunique() >= 2),
              :]
y1h = load_y1h_pdi_data()
m1h = load_m1h_activation_data()
# rna = load_rna_expression_data()
isoforms = load_valid_isoform_clones()
iso_pairs = load_ref_vs_alt_isoforms_table()
iso_pairs['both_iso_y2h_pos'] = (iso_pairs['clone_acc_ref'].isin(iso_gte_1_pos_ppi_iso_data_only) &
                                 iso_pairs['clone_acc_alt'].isin(iso_gte_1_pos_ppi_iso_data_only))
iso_pairs['both_iso_y2h_pos_all_data'] = (iso_pairs['clone_acc_ref'].isin(iso_gte_1_pos_ppi_all_data) &
                                          iso_pairs['clone_acc_alt'].isin(iso_gte_1_pos_ppi_all_data))


# In[4]:


iso_pairs.head()


# In[4]:


clone_to_orf_id = y2h[['ad_clone_acc', 'ad_orf_id']].drop_duplicates().set_index('ad_clone_acc')['ad_orf_id']


# In[5]:


non_zero_iso = set(iso_pairs.loc[(iso_pairs['n_PPI_successfully_tested_in_ref_and_alt'] >= 2)
                             & iso_pairs['both_iso_y2h_pos'], ['clone_acc_ref', 'clone_acc_alt']].values.flatten())
non_zero_iso = {clone_to_orf_id[s] for s in non_zero_iso}


# In[6]:


#TEMP DEBUG
y2h = load_isoform_and_paralog_y2h_data()


# In[7]:


qry = """SELECT a.test_orf_ida,
	   a.test_orf_idb,
	   b.iso_orf_id,
	   b.huri_orf_id,
	   b.source,
	   a.call_1_percent_RRS AS result,
	   a.final_score,
	   c.score,
	   c.empty_n1, c.empty_n2
        FROM tf_validation.validation AS a
        LEFT JOIN tf_validation.validation_source AS b
        ON (a.test_orf_ida = b.orf_id1
        	AND a.test_orf_idb = b.orf_id2)
          OR (a.test_orf_idb = b.orf_id1
        	AND a.test_orf_ida = b.orf_id2)
        LEFT JOIN (select j.score_id, j.score, 
	   						k.empty_n1, k.empty_n2 
							from tf_validation.mn2h_scoring AS j
							left join tf_validation.mn2h_control AS k
							using (plate, well)) AS c
          on a.final_score_id = c.score_id;"""
n2h = pd.read_sql(qry, paros_connection())
n2h['max_control'] = n2h[['empty_n1', 'empty_n2']].max(axis=1)
n2h['non_zero_iso'] = n2h['iso_orf_id'].isin(non_zero_iso)
n_rows_b4 = n2h.shape[0]
n2h = pd.merge(n2h,
         y2h.loc[y2h['category'].isin({'tf_isoform_ppis', 
                                        'tf_paralog_ppis',
                                        'paralog_with_PDI',
                                        'non_paralog_control'}),
				 ['ad_orf_id',
				 'db_orf_id',
				 'ad_gene_symbol',
				 'ad_clone_acc',
				 'db_gene_symbol',
				 'category',
				 'Y2H_result']],
         how='left',
         left_on=['iso_orf_id', 'huri_orf_id'],
         right_on=['ad_orf_id', 'db_orf_id'],
         suffixes=('_n2h', '_y2h'))
if n2h.shape[0] != n_rows_b4:
    raise UserWarning('Problem with table join')


# In[8]:


tcf4_orf_ids = y2h.loc[y2h['ad_clone_acc'].str.startswith('TCF4'), 'ad_orf_id'].unique()
n2h.loc[(n2h['test_orf_ida'].isin(tcf4_orf_ids) | 
        n2h['test_orf_idb'].isin(tcf4_orf_ids)) &
        (n2h['source'] == 'pos'), 'result'].value_counts()


# In[9]:


np.sqrt((0.8*0.2) / 20)


# In[10]:


mismatch = (n2h['source'].isin({'pos', 'neg', 'pos-matched-neg'}) & 
               n2h['Y2H_result'].isnull())
mismatch = mismatch | (n2h['source'].isin({'pos'}) & (n2h['Y2H_result'] != True))
mismatch = mismatch | (n2h['source'].isin({'neg', 'pos-matched-neg'}) & (n2h['Y2H_result'] != False))
n2h['mismatch_current_data'] = mismatch
n2h = n2h.loc[~n2h['mismatch_current_data'], :]


# In[11]:


n2h = pd.merge(n2h.drop(columns=['category']),
         y2h.loc[:, ['ad_orf_id', 'db_orf_id', 'category', 'Y2H_result']],
         how='left',
         left_on=['iso_orf_id', 'huri_orf_id'],
         right_on=['ad_orf_id', 'db_orf_id'])


# In[12]:


n2h.groupby('source')['category'].apply(lambda x: x.isnull().sum())


# In[13]:


for source in n2h['source'].unique():
    print(source)
    print(n2h.loc[n2h['source'] == source, 'category'].value_counts())
    print()


# In[14]:


for source in n2h['source'].unique():
    print(source)
    print(n2h.loc[n2h['source'] == source, 'score'].value_counts())
    print()


# In[15]:


n2h.groupby('source')['non_zero_iso'].value_counts()


# In[16]:


n2h['source'].unique()


# In[17]:


sources = ['huri_rrs', 'TF_RRS', 'lit_bm_2013_rand250', 'litbm',
           'pos', 'pos-matched-neg', 'neg']
cplt.validation_plot(data=n2h,
                     selections=[n2h['source'] == s for s in sources],
                     labels=sources,
                     y_max=0.4,
                    xlabel_rotation=90,
                     errorbar_capsize=0.1)
plt.show()


# In[18]:


samples = {'huri_rrs': n2h['source'] == 'huri_rrs',
           'TF_RRS': n2h['source'] == 'TF_RRS',
           'lit_bm_2013_rand250': n2h['source'] == 'lit_bm_2013_rand250',
           'litbm': n2h['source'] == 'litbm',
           '': pd.Series(index=n2h.index, data=False),
           'pos_non_zero_iso': (n2h['source'] == 'pos') &
                               n2h['non_zero_iso'],
           'pos-matched-neg_non_zero_iso': (n2h['source'] == 'pos-matched-neg') & n2h['non_zero_iso'],
            'neg_non_zero_iso': (n2h['source'] == 'neg') & n2h['non_zero_iso'],
            ' ': pd.Series(index=n2h.index, data=False),
           'pos_other': (n2h['source'] == 'pos') & ~n2h['non_zero_iso'],
           'pos-matched-neg_other': (n2h['source'] == 'pos-matched-neg') & ~n2h['non_zero_iso'],
           'neg_other': (n2h['source'] == 'neg') & ~n2h['non_zero_iso']}

cplt.validation_plot(data=n2h,
                     selections=list(samples.values()),
                     labels=list(samples.keys()),
                     y_max=0.25,
                    xlabel_rotation=90,
                     errorbar_capsize=0.1)
plt.savefig('../figures/N2H_non_zero_isoforms.pdf', bbox_inches='tight')
plt.show()


# In[19]:


n2h.columns


# In[20]:


(n2h.loc[n2h['source'].isin({'pos', 'neg', 'pos-matched-neg'}) &
        n2h['non_zero_iso'],
        ['ad_gene_symbol', 'ad_clone_acc', 'db_gene_symbol',
         'source', 'category', 'result']]
    .rename(columns={'source': 'N2H_category',
                     'category': 'Y2H_category',
                     'result': 'N2H_result'})
    .sort_values(['N2H_category', 'Y2H_category', 'ad_gene_symbol'])
    .to_csv('../output/N2H_non_zero_iso_only.tsv',
            index=False,
            sep='\t'))


# In[21]:


fig, axes = plt.subplots(n2h['source'].nunique())
fig.set_size_inches(4, 4 * axes.shape[0])
for cat, ax in zip(n2h['source'].unique(), axes):
    ax.set_title(cat)
    is_cat = n2h['source'] == cat
    is_pos = n2h['result'] == 1
    is_neg = n2h['result'] == 0
    ax.scatter(x=np.log2(n2h.loc[is_cat & is_neg, 'max_control']),
               y=np.log2(n2h.loc[is_cat & is_neg, 'score']))
    ax.scatter(x=np.log2(n2h.loc[is_cat & is_pos, 'max_control']),
               y=np.log2(n2h.loc[is_cat & is_pos, 'score']))
    ax.set_xlim(2, 18)
    ax.set_ylim(2, 18)
    ax.set_xlabel('Log2 max empty N1/N2 score')
    ax.set_ylabel('Log2 pair score')
plt.subplots_adjust(hspace=0.35)
plt.savefig('../figures/n2h_scatter.pdf', bbox_inches='tight')


# In[22]:


sources = n2h['source'].unique()
cplt.validation_plot(data=n2h,
                     selections=[n2h['source'] == s for s in sources],
                     labels=sources,
                     y_max=0.4,
                    xlabel_rotation=90,
                     errorbar_capsize=0.1)
plt.show()


# In[23]:


# did the TF RRS / Lit-BM really take the longest isoforms?
n2h.loc[n2h['source'].isin(['TF_RRS', 'litbm']), 'iso_orf_id'].map(clone_to_orf_id.reset_index().set_index('ad_orf_id')['ad_clone_acc']).unique()


# In[24]:


from data_loading import load_valid_isoform_clones
clones = load_valid_isoform_clones()


# In[25]:


clones.loc[clones['gene_symbol'] == 'PAX8', :]

