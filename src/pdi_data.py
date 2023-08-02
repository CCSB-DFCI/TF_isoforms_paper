
# coding: utf-8

# In[1]:


import seaborn as sns

from data_loading import load_y1h_pdi_data, load_tf_families


# In[2]:


pdi = load_y1h_pdi_data()
pdi[pdi.columns[2:]].sum().sum()
pdi['tf'].nunique()
(pdi[pdi.columns[2:]].sum() > 0).sum()
len(pdi.columns[2:])


# In[3]:


tf_fam = load_tf_families()


# In[4]:


pdi.head()


# In[5]:


pdi['DBD'] = pdi['tf'].map(tf_fam)


# In[6]:


pdi['n_pdi'] = pdi.iloc[:, 2:-1].sum(axis=1)


# In[7]:


pdi['DBD'].value_counts()


# In[8]:


dbd_to_plot = ['C2H2 ZF', 'bHLH', 'Homeodomain', 'Nuclear receptor']
sns.swarmplot(data=pdi, x='DBD', y='n_pdi', order=dbd_to_plot, alpha=1)

