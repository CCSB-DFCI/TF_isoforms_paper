#!/usr/bin/env python
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
                          load_annotated_TFiso1_collection)
from plotting import (y2h_ppi_per_tf_gene_plot,
                      y1h_pdi_per_tf_gene_plot,
                      m1h_activation_per_tf_gene_plot)


# In[2]:


PAPER_PRESET = {"style": "ticks", "font": "Helvetica", "context": "paper", 
                "rc": {"font.size":10,"axes.titlesize":10,
                       "axes.labelsize":10, 'axes.linewidth':0.5,
                       "legend.fontsize":10, "xtick.labelsize":10,
                       "ytick.labelsize":10, "xtick.major.size": 3.0,
                       "ytick.major.size": 3.0, "axes.edgecolor": "black",
                       "xtick.major.pad": 3.0, "ytick.major.pad": 3.0}}
PAPER_FONTSIZE = 10


# In[3]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[4]:


y2h = load_isoform_and_paralog_y2h_data()
y1h = load_y1h_pdi_data(add_missing_data=True)
m1h = load_m1h_activation_data(add_missing_data=True)
isoforms = load_valid_isoform_clones()
y2h = y2h.loc[y2h['ad_clone_acc'].isin(isoforms['clone_acc']).values, :]
y1h = y1h.loc[y1h['clone_acc'].isin(isoforms['clone_acc']).values, :]
m1h = m1h.loc[m1h['clone_acc'].isin(isoforms['clone_acc'].values), :]

tfs = load_annotated_TFiso1_collection()


# In[5]:


len([orf for tf in tfs.values() for orf in tf.isoforms])


# ## PDIs

# In[6]:


gene_name = "PKNOX1"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(2, 1.3))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax)
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[7]:


gene_name = "HEY1"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(1.3, 1.3))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax, iso_order=["HEY1-2", "HEY1-1"])
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[8]:


gene_name = "KLF7"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(1.6, 1.6))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax)
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[9]:


gene_name = "CREB1"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(1.6, 1.6))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax, iso_order=["CREB1-2", "CREB1-1"])
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[10]:


gene_name = "TBX5"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax)
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[11]:


gene_name = "LHX9"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax)
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[12]:


gene_name = "ZIC3"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(3, 3))
y1h_pdi_per_tf_gene_plot(gene_name, data=y1h, ax=ax)
plt.savefig('../figures/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ## PPIs

# In[13]:


gene_name = "PKNOX1"

fig, ax = plt.subplots(1, 1, figsize=(3, 3))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../figures/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[14]:


gene_name = "GRHL3"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../figures/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[15]:


gene_name = "LHX9"

tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../figures/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ## M1H

# In[16]:


gene_name = "RFX3"

fig, ax = plt.subplots(1, 1, figsize=(3, 1.2))

df = m1h_activation_per_tf_gene_plot("RFX3", data=m1h, ax=ax, xlim=(0, 6))
plt.savefig('../figures/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[17]:


gene_name = "LHX9"

fig, ax = plt.subplots(1, 1, figsize=(3, 1.2))

df = m1h_activation_per_tf_gene_plot("LHX9", data=m1h, ax=ax, xlim=(0, 6))
plt.savefig('../figures/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[18]:


gene_name = "PKNOX1"

fig, ax = plt.subplots(1, 1, figsize=(3, 1.5))

df = m1h_activation_per_tf_gene_plot("PKNOX1", data=m1h, ax=ax, xlim=(0, 6))
plt.savefig('../figures/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[19]:


gene_name = "CREB1"

fig, ax = plt.subplots(1, 1, figsize=(2, 1))

df = m1h_activation_per_tf_gene_plot("CREB1", data=m1h, ax=ax)
plt.savefig('../figures/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[33]:


gene_name = "TGIF1"

fig, ax = plt.subplots(1, 1, figsize=(2, 1))

df = m1h_activation_per_tf_gene_plot("TGIF1", data=m1h, ax=ax)
plt.savefig('../figures/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ## exon diagrams

# In[20]:


fig, ax = plt.subplots(figsize=(7, 1))

tfs["HEY1"].exon_diagram(ax=ax)
fig.savefig("../figures/HEY1_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[21]:


fig, ax = plt.subplots(figsize=(7, 1.5))

tfs["LHX9"].exon_diagram(ax=ax)
fig.savefig("../figures/LHX9_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[22]:


fig, ax = plt.subplots(figsize=(7, 0.75))

tfs["CREB1"].exon_diagram(ax=ax)
fig.savefig("../figures/CREB1_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[23]:


fig, ax = plt.subplots(figsize=(7, 0.75))

tfs["ZNF414"].exon_diagram(ax=ax)
fig.savefig("../figures/ZNF414_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[24]:


fig, ax = plt.subplots(figsize=(4, 2))

tfs["HEY1"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)
fig.savefig("../figures/HEY1_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[25]:


fig, ax = plt.subplots(figsize=(4, 2))

tfs["CREB1"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)
fig.savefig("../figures/CREB1_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[26]:


fig, ax = plt.subplots(figsize=(4, 3))

tfs["ZIC3"].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../figures/ZIC3_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[27]:


fig, ax = plt.subplots(figsize=(3, 2.5))

tfs["ZIC3"].exon_diagram(ax=ax)
fig.savefig("../figures/ZIC3_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[28]:


fig, ax = plt.subplots(figsize=(8, 1))

tfs["RFX3"].exon_diagram(ax=ax)
fig.savefig("../figures/RFX3_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[29]:


fig, ax = plt.subplots(figsize=(8, 1))

tfs["RFX3"].protein_diagram(ax=ax, only_cloned_isoforms=True)
fig.savefig("../figures/RFX3_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[30]:


fig, ax = plt.subplots(figsize=(7, 0.75))

tfs["TBPL1"].exon_diagram(ax=ax)
fig.savefig("../figures/TBPL1_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[31]:


fig, ax = plt.subplots(figsize=(7, 0.75))

tfs["TBPL1"].protein_diagram(ax=ax)
fig.savefig("../figures/TBPL1_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[32]:


fig, ax = plt.subplots(figsize=(7, 0.75))

tfs["TGIF1"].protein_diagram(ax=ax)
fig.savefig("../figures/TGIF1_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[ ]:




