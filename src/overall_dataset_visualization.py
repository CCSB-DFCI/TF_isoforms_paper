
# coding: utf-8

# In[13]:


import pandas as pd

from data_loading import load_isoform_and_paralog_y2h_data, load_y1h_pdi_data, load_valid_isoform_clones


# In[47]:


# table of edges
#    - clone to (edge + clone_id) + to duplicate
# table of nodes
#    - clone to gene
#    - dna vs isoform vs 
ppi = load_isoform_and_paralog_y2h_data()
ppi = ppi.loc[(ppi['category'] == 'tf_isoform_ppis') &
              (ppi['score'] == '1'),
              ['ad_clone_acc', 'ad_gene_symbol', 'db_gene_symbol']]
ppi = ppi.rename(columns={'ad_clone_acc': 'isoform',
                          'db_gene_symbol': 'partner'})
ppi['partner'] = ppi['partner'] + '-' + ppi['ad_gene_symbol']
pdi = pd.read_csv('../data/a2_juan_pdi_w_unique_isoacc.tsv', sep='\t')
clones = load_valid_isoform_clones()
pdi = pdi.loc[pdi['unique_acc'].isin(clones['clone_acc']), :]
pdi['partner'] = pdi['bait'] + '-' + pdi['tf']
pdi['isoform'] = pdi['unique_acc']
edges = pd.concat([ppi.loc[:, ['isoform', 'partner']],
                   pdi.loc[:, ['isoform', 'partner']]])
edges.to_csv('../output/edges.tsv', sep='\t', index=False)

clones = clones.rename(columns={'clone_acc': 'node_id'})
clones['type'] = 'isoform'
dna = pd.DataFrame(data=pdi['partner'].unique(), columns=['node_id'])
dna['type'] = 'DNA'
proteins = pd.DataFrame(data=ppi['partner'].unique(), columns=['node_id'])
proteins['type'] = 'Protein'
nodes = pd.concat([clones, proteins, dna], sort=True)
nodes.to_csv('../output/node_table.tsv', sep='\t', index=False)

