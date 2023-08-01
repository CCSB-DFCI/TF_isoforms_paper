
# coding: utf-8

# # notebook to make nicer PPI/PDI/M1H panels for figures

# In[1]:


import numpy as np
import seaborn as sns
import shutil
from pathlib import Path

from matplotlib import pyplot as plt
import pandas as pd

from data_loading import (load_isoform_and_paralog_y2h_data,
                          load_y1h_pdi_data,
                          load_m1h_activation_data,
                          load_valid_isoform_clones,
                          load_annotated_6k_collection)
from plotting import (y2h_ppi_per_tf_gene_plot,
                      y1h_pdi_per_tf_gene_plot,
                      m1h_activation_per_tf_gene_plot, PAPER_PRESET, PAPER_FONTSIZE)


# In[2]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[3]:


y2h = load_isoform_and_paralog_y2h_data(add_missing_data=True)
y1h = load_y1h_pdi_data(add_missing_data=True)
m1h = load_m1h_activation_data(add_missing_data=True)
isoforms = load_valid_isoform_clones()
y2h = y2h.loc[y2h['ad_clone_acc'].isin(isoforms['clone_acc']).values, :]
y1h = y1h.loc[y1h['unique_acc'].isin(isoforms['clone_acc']).values, :]
m1h = m1h.loc[m1h['clone_acc'].isin(isoforms['clone_acc'].values), :]

tfs = load_annotated_6k_collection()


# In[4]:


len([orf for tf in tfs.values() for orf in tf.orfs])


# ## PDIs

# In[20]:


gene_name = "PKNOX1"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(1.3, 1.0))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax, equal_aspect=False)
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[21]:


gene_name = "GRHL3"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(6, 1.0))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax, equal_aspect=False)
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ## PPIs

# In[13]:


gene_name = "PKNOX1"

fig, ax = plt.subplots(1, 1, figsize=(0.87, 0.95))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h, equal_aspect=False)
plt.savefig('../figures/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[27]:


gene_name = "GRHL3"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(0.9, 1.6))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h, equal_aspect=False)
plt.savefig('../figures/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[31]:


tfs["PKNOX1"].splicing_categories("PKNOX1-1", "PKNOX1-2")

