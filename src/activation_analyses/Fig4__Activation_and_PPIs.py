#!/usr/bin/env python
# coding: utf-8

# # Fig 4: Activation and PPIs/Dimerization
# 
# considering effector domains to just be from the functional studies (no pfam domains). results do not change that much, and makes the code/explanation cleaner to focus on effector domains imo.

# In[1]:


import numpy as np
import matplotlib as mpl
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import seaborn as sns
import sys

from statannotations.Annotator import Annotator
from scipy.stats import wilcoxon

# import utils
sys.path.append("../")
sys.path.append("../data_loading")

import plotting
from plotting import (mimic_r_boxplot,
                      violinplot_reflected,
                      y2h_ppi_per_tf_gene_plot,
                      y1h_pdi_per_tf_gene_plot,
                      m1h_activation_per_tf_gene_plot,
                      annotate_pval)

from data_loading import (load_y2h_isoform_data, 
                          load_m1h_activation_data, 
                          load_ppi_partner_categories, 
                          load_annotated_TFiso1_collection,
                          load_human_tf_db,
                          load_y1h_pdi_data,
                          load_tf_families,
                          DIMERIZING_TF_FAMILIES,
                          load_ref_vs_alt_isoforms_table,
                          load_pfam_domains_horfeome,
                          load_full_y2h_data_including_controls)

from data_loading.isoform_pairwise_metrics import _add_PPI_columns


# In[2]:


PAPER_PRESET = {"style": "ticks", "font": "Helvetica", "context": "paper", 
                "rc": {"font.size":7,"axes.titlesize":7,
                       "axes.labelsize":7, 'axes.linewidth':0.5,
                       "legend.fontsize":6, "xtick.labelsize":6,
                       "ytick.labelsize":6, "xtick.major.size": 3.0,
                       "ytick.major.size": 3.0, "axes.edgecolor": "black",
                       "xtick.major.pad": 3.0, "ytick.major.pad": 3.0}}
PAPER_FONTSIZE = 7


# In[3]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[4]:


np.random.seed(2023)


# ## 1. load TFs, assay data (Y2H, M1H), Pfam domains, and PPI partner categories/cofactors

# In[5]:


tfs = load_annotated_TFiso1_collection()
pairs = load_ref_vs_alt_isoforms_table()

# RORC-1 alt iso is causing an error - filter out here - there's no data for it?
pairs = pairs[pairs["clone_acc_alt"] != "RORC|1/6|05F11"]

pairs['ref_iso'] = pairs['clone_acc_ref'].apply(lambda x: x.split('|')[0] + '-' + x.split('|')[1].split('/')[0])
pairs['alt_iso'] = pairs['clone_acc_alt'].apply(lambda x: x.split('|')[0] + '-' + x.split('|')[1].split('/')[0])
pairs['f_disorder_difference'] = pairs.apply(lambda x: tfs[x['gene_symbol']].disordered_fraction_of_different_regions(x['ref_iso'], x['alt_iso']), axis=1)


# In[6]:


y2h = load_y2h_isoform_data()
m1h = load_m1h_activation_data()
m1h['mean'] = m1h[['M1H_rep1', 'M1H_rep2', 'M1H_rep3']].mean(axis=1)
cats = load_ppi_partner_categories()


# In[7]:


df = pd.read_excel('../../data/external/Geiger-et-al_MCP_2012_Supplementary-Table-2.xlsx',
                   skiprows=1)
hek_avrg = df[['iBAQ HEK293_1', 'iBAQ HEK293_2', 'iBAQ HEK293_3']].mean(axis=1)
print((hek_avrg > 0).sum(), 'proteins expressed in HEK293 proteome')
hek_expressed_genes = set(df.loc[(hek_avrg > 0) & df['Gene Names'].notnull(),
       'Gene Names'].str.split(';').explode().values)
all_partners = set(y2h['db_gene_symbol'].unique())
print('of {} PPI partners, {} are expressed in HEK293 cells'.format(len(all_partners), 
      len(all_partners.intersection(hek_expressed_genes))))


# In[8]:


# now add Pfam AD/RDs
pfam = pd.read_csv('../../data/external/Pfam-A.clans.tsv',
                   sep='\t',
                   names=['pfam_accession', 'clan', 'clan_name', 'short_name', 'name'])


# In[9]:


cof = pd.read_csv('../../data/external/AnimalTFDB3_Homo_sapiens_TF_cofactors.txt',
                 sep='\t')
if cof['Symbol'].duplicated().any():
    raise UserWarning('unexpected duplicates')


# ## 2. categorize effector domain changes between alt/ref iso

# In[10]:


dom = pd.concat([g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()])

# add activation or repression annotation from effector domain studies
effector_domain_type = {}
for tf in tfs.values():
    for d in tf.cloned_reference_isoform.aa_seq_features:
        if d.category == 'effector_domain':
            effector_domain_type[d.accession] = d.name
dom['type'] = dom['accession'].map(effector_domain_type)


# In[11]:


def fraction_of_effector_domains_removed(row, effector_type):
    ds = dom.loc[(dom['alt_iso'] == row['alt_iso']) 
                  & (dom['type'].isin(effector_type)), :]
    if ds.shape[0] == 0:
        return np.nan
    return ds[['deletion', 'frameshift']].sum().sum() / ds['length'].sum()


def insertion_in_effector_domains(row, effector_type):
    ds = dom.loc[(dom['alt_iso'] == row['alt_iso']) 
                  & (dom['type'].isin(effector_type)), :]
    if ds.shape[0] == 0:
        return np.nan
    return ds['insertion'].sum()

def domain_length(row, effector_type):
    ds = dom.loc[(dom['alt_iso'] == row['alt_iso']) 
                  & (dom['type'].isin(effector_type)), :]
    if ds.shape[0] == 0:
        return np.nan
    return ds['length'].sum()


for effector_type, effector_type_name in zip([['AD'], ['RD'], ['Bif'], ['AD', 'RD', 'Bif']], ['AD', 'RD', 'Bif', 'all']):
    pairs['fraction_of_{}_domains_removed'.format(effector_type_name)] = pairs.apply(fraction_of_effector_domains_removed, effector_type=effector_type, axis=1)
    pairs['insertion_in_{}_domains'.format(effector_type_name)] = pairs.apply(insertion_in_effector_domains, effector_type=effector_type, axis=1)
    pairs['length_of_{}_domains'.format(effector_type_name)] = pairs.apply(domain_length, effector_type=effector_type, axis=1)


# ## 3. plot number of annotated domains (of various types) across isoforms

# In[12]:


# plot of number of activation domains per ref iso
# fraction of sequnce within effector domains
def count_effector_domains(gene):
    iso = gene.cloned_reference_isoform
    c = 0
    for d in iso.aa_seq_features:
        if d.category == 'effector_domain':
            c += 1
    return c

n_effector = [count_effector_domains(tf) for tf in tfs.values()]

fig, ax = plt.subplots(1, 1, figsize=(2, 1))
ax.hist(n_effector,
        range=(-0.25, max(n_effector) + 0.25),
          bins=(max(n_effector) * 2 + 1))
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel('Number of genes ({} total)'.format(len(tfs)))
ax.set_xlabel('Effector domains in reference isoform')


# In[13]:


# plot of number of activation domains per ref iso
# fraction of sequnce within effector domains
def count_Soto_effector_domains(gene):
    iso = gene.cloned_reference_isoform
    c = 0
    for d in iso.aa_seq_features:
        if d.category == 'effector_domain' and d.accession.startswith('Soto'):
            c += 1
    return c

n_effector_soto = [count_Soto_effector_domains(tf) for tf in tfs.values()]

fig, ax = plt.subplots(1, 1, figsize=(2, 1))
ax.hist(n_effector_soto,
        range=(-0.25, max(n_effector) + 0.25),
          bins=(max(n_effector) * 2 + 1))
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel('Number of genes ({} total)'.format(len(tfs)))
ax.set_xlabel('Effector domains in reference isoform')
ax.set_title('Soto et al. data')


# In[14]:


# plot of number of activation domains per ref iso
# fraction of sequnce within effector domains
def count_Bintu_effector_domains(gene):
    iso = gene.cloned_reference_isoform
    c = 0
    for d in iso.aa_seq_features:
        if d.category == 'effector_domain' and not d.accession.startswith('Soto'):
            c += 1
    return c

n_effector_bintu = [count_Bintu_effector_domains(tf) for tf in tfs.values()]

fig, ax = plt.subplots(1, 1, figsize=(2, 1))
ax.hist(n_effector_bintu,
        range=(-0.25, max(n_effector) + 0.25),
          bins=(max(n_effector) * 2 + 1))
ax.set_xticks(range(max(n_effector) + 1))
ax.set_ylabel('Number of genes ({} total)'.format(len(tfs)))
ax.set_xlabel('Effector domains in reference isoform')
ax.set_title('Data from Bintu lab papers')


# In[15]:


counter = {'Soto': {'AD': 0, 'RD': 0, 'Bif': 0},
           'Tycko': {'AD': 0, 'RD': 0},
           'DelRosso': {'AD': 0, 'RD': 0}}
for tf in tfs.values():
    has_effector = False
    for domain in tf.cloned_reference_isoform.aa_seq_features:
        if domain.category == 'effector_domain':
            counter[domain.accession.split('_')[0]][domain.name] += 1
counter


# In[16]:


fig, ax = plt.subplots(1, 1, figsize=(2, 1.5))
ax.bar(x=[0.1, 1.0, 2.1, 3.0, 4.1, 5.0],
       height=[counter[x][y] for x in counter.keys() for y in ['AD', 'RD']],
       color=[sns.color_palette("Set2")[0], sns.color_palette("Set2")[1]] * 3)
ax.set_xticks([0.5, 2.5, 4.5])
ax.set_xticklabels(['Soto et al.', 'Tycko et al.', 'DelRosso et al.'], rotation=30, ha="right", va="top")
ax.set_ylabel('Number of effector domains')

# annotate
rects = ax.patches
labels = [counter[x][y] for x in counter.keys() for y in ['AD', 'RD']]
colors = [sns.color_palette("Set2")[0], sns.color_palette("Set2")[1]] * 3
for rect, label, color in zip(rects, labels, colors):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height, label, ha="center", va="bottom", color=color)

colors = {"activation domain": sns.color_palette("Set2")[0], 
          "repression domain": sns.color_palette("Set2")[1]}
labels = list(colors.keys())
handles = [plt.Rectangle((0,0), 1, 1, color=colors[label]) for label in labels]
plt.legend(handles, labels, loc=2, bbox_to_anchor=(1.01, 1))

plt.ylim((0, 400))

for pos in ['top', 'right']:
    ax.spines[pos].set_visible(False)

fig.savefig("../../figures/fig4/annotated_effector_domain_count.pdf", dpi="figure", bbox_inches="tight")


# ## 4. summary plot looking at presence of activ/repr domains and activity

# In[17]:


# define activity above baseline as >= 1 (absolute value)
m1h['gte_2_fold'] = (m1h['mean'].abs() >= 1)
m1h['gte_above'] = (m1h['mean'] >= 1)
m1h['gte_below'] = (m1h['mean'] <= -1)
m1h.head()


# In[18]:


ads = []
rds = []
bifs = []

for i, row in m1h.iterrows():
    gene = row.gene_symbol
    iso = row.clone_acc
    try:
        tf_iso = tfs[gene][iso]
    except KeyError:
        print("missing tf %s or iso %s" % (gene, iso))
        ads.append(np.nan)
        rds.append(np.nan)
        bifs.append(np.nan)
        continue
    
    has_ad = False
    has_rd = False
    has_bif = False
    for domain in tf_iso.aa_seq_features:
        if domain.category == 'effector_domain':
            if domain.name == "AD":
                has_ad = True
            if domain.name == "RD":
                has_rd = True
            if domain.name == "Bif":
                has_bif = True
    
    ads.append(has_ad)
    rds.append(has_rd)
    bifs.append(has_bif)


# In[19]:


m1h['has_ad'] = ads
m1h['has_rd'] = rds
m1h['has_bif'] = bifs

# filter out the above isos that got removed from clone collection post-hoc
m1h = m1h[~pd.isnull(m1h['has_ad'])]
len(m1h)


# In[20]:


len(m1h.gene_symbol.unique())


# In[21]:


def cat_dom(row):
    if row.has_ad == True and row.has_rd == False and row.has_bif == False:
        return "activation domain"
    elif row.has_ad == False and row.has_rd == True and row.has_bif == False:
        return "repression domain"
    elif row.has_ad == False and row.has_rd == False and row.has_bif == True:
        return "bifunctional domain"
    elif row.has_ad == False and row.has_rd == False and row.has_bif == False:
        return "no annotated domains"
    else:
        return "combination of domains"

m1h["cat_dom"] = m1h.apply(cat_dom, axis=1)
m1h.cat_dom.value_counts()


# In[22]:


def cat_gte(row):
    if row.gte_above:
        return "above baseline"
    elif row.gte_below:
        return "below baseline"
    else:
        return "NA"
    
m1h["cat_gte"] = m1h.apply(cat_gte, axis=1)
m1h.cat_gte.value_counts()


# In[23]:


m1h_filt = m1h[(m1h["cat_gte"] != "NA") & (m1h["cat_dom"].isin(["activation domain", "repression domain"]))]
m1h_filt = pd.pivot_table(m1h_filt, index="cat_dom", columns="cat_gte", values='clone_acc', aggfunc='count')
print(m1h_filt.sum())
m1h_filt = m1h_filt/m1h_filt.sum()
m1h_filt


# In[24]:


palette = {"activation domain": sns.color_palette("Set2")[0],
           "repression domain": sns.color_palette("Set2")[1]}
sns.palplot(palette.values())


# In[25]:


ax = m1h_filt.T.plot.bar(stacked=True, color=palette.values(), figsize=(1, 1.5))

ax.set_ylabel("% of isoforms")
ax.set_xlabel("")

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))
ax.set_xticklabels(["above M1H baseline", "below M1H baseline"], ha="right", va="top", rotation=30)

for pos in ['top', 'right']:
    ax.spines[pos].set_visible(False)

plt.savefig('../../figures/fig4/m1h_baseline_doms.pdf',
            bbox_inches='tight')


# ### how many TFs have isos that are both strong activators and strong repressors?

# In[26]:


m1h_stack = pd.pivot_table(m1h, values="clone_acc", index="gene_symbol", columns="cat_gte", aggfunc='count').reset_index()
m1h_stack.fillna(0, inplace=True)
m1h_stack[(m1h_stack["above baseline"] > 0) & (m1h_stack["below baseline"] > 0)]


# ## 5. summary plot looking at gain/loss of domains and activity

# In[27]:


pairs['m1h_gte_2_fold_at_least_one_iso_per_gene'] = pairs['gene_symbol'].map(m1h.groupby('gene_symbol')
                                                                             ['gte_2_fold']
                                                                             .any())
pairs['abs_activation_fold_change_log2'] = pairs['activation_fold_change_log2'].abs()


# In[28]:


# create a color map of domain length
# sum up lengths of all domains (plot only includes examples w 1 type of domain)
pairs['length_of_all_domains'].fillna(0, inplace=True)
t_dom_length = pairs.loc[:,'length_of_all_domains'].values
t_dom_length = t_dom_length[t_dom_length > 0]

# using min and max makes colors too hard too read - cut off
cmap = sns.color_palette("flare", as_cmap=True)
norm = plt.Normalize(25, 250)
palette_dom_length = {value: cmap(norm(value)) for value in t_dom_length}

def re_color(row, palette):
    if row['length_of_all_domains'] == 0:
        color = sns.color_palette("flare")[0]
    else:
        color = palette[row['length_of_all_domains']]
    return color

pairs["color_dom_length"] = pairs.apply(re_color, axis=1, palette=palette_dom_length)


# In[29]:


df = pairs.copy()
df = df.loc[df['activation_fold_change_log2'].notnull() & df['m1h_gte_2_fold_at_least_one_iso_per_gene'], :]


# In[30]:


fig = plt.figure(figsize=(2, 1.5))
ax = sns.histplot(data=df, x="activation_fold_change_log2", color="slategrey")
ax.set_xlabel("M1H activation foldchange\n(log2(alt/ref))")
ax.set_ylabel("count")

ax.axvline(x=-1, linestyle="dashed", color="black", linewidth=1)
ax.axvline(x=1, linestyle="dashed", color="black", linewidth=1)

# annotate
n_less_neg1 = len(df[df["activation_fold_change_log2"] <= -1])
n_greater_1 = len(df[df["activation_fold_change_log2"] >= 1])

ax.text(-1.2, 40, "%s pairs\n(%s%%)" % (n_less_neg1, np.round(n_less_neg1/len(df)*100, 1)), ha="right", va="center")
ax.text(1.2, 40, "%s pairs\n(%s%%)" % (n_greater_1, np.round(n_greater_1/len(df)*100, 1)), ha="left", va="center")

for pos in ['top', 'right']:
    ax.spines[pos].set_visible(False)

fig.savefig("../../figures/fig4/m1h_alt_ref_dist.pdf", dpi="figure", bbox_inches="tight")


# In[31]:


print("NUM ALT ISOS WITH AT LEAST 2-FOLD CHANGE COMPARED TO REF: %s" % (n_less_neg1 + n_greater_1))
print("PERCENT ALT ISOS WITH AT LEAST 2-FOLD CHANGE COMPARED TO REF: %s" % ((n_less_neg1 + n_greater_1)/len(df)*100))


# In[32]:


# combining full + partial loss since we only have 1 full loss activ. domain
palette = palette_dom_length
hue = 'length_of_all_domains'
color = 'color_dom_length'
t = t_dom_length

gs_kw = dict(width_ratios=[1, 1, 0.8, 1.5, 1])

fig, axs = plt.subplots(1, 5, sharey=True, gridspec_kw=gs_kw)
fig.set_size_inches(w=8.2, h=2)

point_size = 5

######### activation domain (full and partial) #########
tot_n_part_loss_activ = df.loc[(df['fraction_of_AD_domains_removed'] > 0) 
                             & (df['fraction_of_RD_domains_removed'].isnull() | 
                               (df['fraction_of_RD_domains_removed'] == 0))
                             & (df['fraction_of_Bif_domains_removed'].isnull() | 
                               (df['fraction_of_Bif_domains_removed'] == 0)), :]
axs[0].set_title('activation\ndomain')
axs[0].scatter(tot_n_part_loss_activ.loc[:, 'fraction_of_AD_domains_removed'].values,
               tot_n_part_loss_activ.loc[:, 'activation_fold_change_log2'].values,
               alpha=1,
               s=point_size**2,
               c=tot_n_part_loss_activ.loc[:, color].values,
               linewidth=1,
               edgecolor="black",
               clip_on=False)
axs[0].set_xticks([])
axs[0].set_xlabel('Proportion missing')
axs[0].set_xlim(1, 0)
axs[0].set_xticks([1, 0.5, 0.01])
axs[0].set_xticklabels([f'{x:.0%}' for x in axs[0].get_xticks()])

tbx5_y = df.loc[(df["clone_acc_alt"] == "TBX5|3/3|08H01"), 'activation_fold_change_log2'].values[0]
for point in axs[0].collections:
    for x, y in point.get_offsets():
        if np.isclose(tbx5_y, y):
            print("found: %s, %s" % (x, y))
            axs[0].annotate("TBX5-3", xy=(x, y), xytext=(5, -10), textcoords='offset points',
                            arrowprops = dict(arrowstyle="-", connectionstyle="arc3,rad=0.3",
                                              color='black'), 
                            ha="left", va="top", fontsize=7,
                            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))

######### repression domain (full and partial) #########
tot_n_part_loss_repr = df.loc[(df['fraction_of_RD_domains_removed'] > 0)
                            & (df['fraction_of_AD_domains_removed'].isnull() | 
                              (df['fraction_of_AD_domains_removed'] == 0))
                            & (df['fraction_of_Bif_domains_removed'].isnull() | 
                              (df['fraction_of_Bif_domains_removed'] == 0)), :]
axs[1].set_title('repression\ndomain')
axs[1].scatter(tot_n_part_loss_repr.loc[:, 'fraction_of_RD_domains_removed'].values,
               tot_n_part_loss_repr.loc[:, 'activation_fold_change_log2'].values,
               alpha=1,
               s=point_size**2,
               c=tot_n_part_loss_repr.loc[:, color].values,
               linewidth=1,
               edgecolor="black",
               clip_on=False)
axs[1].set_xticks([])
axs[1].set_xlabel('Proportion missing')
axs[1].set_xlim(1, 0)
axs[1].set_xticks([1, 0.5, 0.01])
axs[1].set_xticklabels([f'{x:.0%}' for x in axs[1].get_xticks()])

dlx1_y = df.loc[(df["clone_acc_alt"] == "DLX1|2/2|07E09"), 'activation_fold_change_log2'].values[0]
print(dlx1_y)
for point in axs[1].collections:
    for x, y in point.get_offsets():
        if np.isclose(dlx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[1].annotate("DLX1-2", xy=(x, y), xytext=(5, 10), textcoords='offset points',
                            arrowprops = dict(arrowstyle="-", connectionstyle="arc3,rad=0.3",
                                              color='black'), 
                            ha="left", va="bottom", fontsize=7,
                            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))

######### both activ/repression domain OR bif (full and partial) #########
tot_n_part_loss_both = df.loc[((df['fraction_of_AD_domains_removed'] > 0) &
                              (df['fraction_of_RD_domains_removed'] > 0)) |
                              (df['fraction_of_Bif_domains_removed'] > 0), :]
axs[2].set_title('both activ. &\nrepr. domains')
axs[2].scatter(tot_n_part_loss_both.loc[:, 'fraction_of_all_domains_removed'].values,
               tot_n_part_loss_both.loc[:, 'activation_fold_change_log2'].values,
               alpha=1,
               s=point_size**2,
               c=tot_n_part_loss_both.loc[:, color].values,
               linewidth=1,
               edgecolor="black",
               clip_on=False)
axs[2].set_xticks([])
axs[2].set_xlabel('Proportion missing')
axs[2].set_xlim(1, 0)
axs[2].set_xticks([1, 0.5, 0.01])
axs[2].set_xticklabels([f'{x:.0%}' for x in axs[2].get_xticks()])

e2f3_y = df.loc[(df["clone_acc_alt"] == "E2F3|3/4|10B08"), 'activation_fold_change_log2'].values[0]
for point in axs[2].collections:
    for x, y in point.get_offsets():
        if np.isclose(e2f3_y, y):
            print("found: %s, %s" % (x, y))
            axs[2].annotate("E2F3-3", xy=(x, y), xytext=(5, 15), textcoords='offset points',
                            arrowprops = dict(arrowstyle="-", connectionstyle="arc3,rad=0.3",
                                              color='black'), 
                            ha="left", va="top", fontsize=7,
                            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))



######### everything retained #########
all_retained = df.loc[((df['fraction_of_AD_domains_removed'] == 0) |
                       (df['fraction_of_RD_domains_removed'] == 0) |
                       (df['fraction_of_Bif_domains_removed'] == 0))
                     & (df['fraction_of_AD_domains_removed'].isnull() | 
                        (df['fraction_of_AD_domains_removed'] == 0)) 
                     & (df['fraction_of_RD_domains_removed'].isnull() | 
                        (df['fraction_of_RD_domains_removed'] == 0)) 
                     & (df['fraction_of_Bif_domains_removed'].isnull() | 
                        (df['fraction_of_Bif_domains_removed'] == 0)), :]
axs[3].set_title('All effector domains\nin alt. iso.')
sns.swarmplot(data=all_retained,
              y='activation_fold_change_log2', 
              x='m1h_gte_2_fold_at_least_one_iso_per_gene',
              size=point_size,
              clip_on=False,
              ax=axs[3],
              linewidth=1,
              edgecolor="black",
              alpha=1,
              hue=hue,
              palette=palette)
axs[3].set_xticks([])
axs[3].set_xlabel('')
axs[3].get_legend().remove()

# annotate pbx1 and rfx3
pbx1_y = df.loc[(df["clone_acc_alt"] == "PBX1|2/2|02C05"), 'activation_fold_change_log2'].values[0]
creb5_y = df.loc[(df["clone_acc_alt"] == "CREB5|2/3|08A12"), 'activation_fold_change_log2'].values[0]
for point in axs[3].collections:
    for x, y in point.get_offsets():
        if np.isclose(pbx1_y, y):
            print("found: %s, %s" % (x, y))
            axs[3].annotate("PBX1-2", xy=(x, y), xytext=(-5, -10), textcoords='offset points',
                            arrowprops = dict(arrowstyle="-", connectionstyle="arc3,rad=0.3",
                                              color='black'), 
                            ha="center", va="top", fontsize=7,
                            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))
        if np.isclose(creb5_y, y):
            print("found: %s, %s" % (x, y))
            axs[3].annotate("CREB5-2", xy=(x, y), xytext=(10, -20), textcoords='offset points',
                            arrowprops = dict(arrowstyle="-", connectionstyle="arc3,rad=0.3",
                                              color='black'), 
                            ha="center", va="top", fontsize=7,
                            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))

# missing stuff
incl = pd.concat([tot_n_part_loss_activ, 
                  tot_n_part_loss_repr, 
                  tot_n_part_loss_both, 
                  all_retained])

######### no annotated domains #########
no_annot = df.loc[(~df.index.isin(incl.index.values)) & (pd.isnull(df["fraction_of_AD_domains_removed"])) &
                  (pd.isnull(df["fraction_of_RD_domains_removed"])) & 
                  (pd.isnull(df["fraction_of_Bif_domains_removed"]))]
axs[4].set_title('No annotated\neffector domains')
sns.swarmplot(data=no_annot,
              y='activation_fold_change_log2', 
              x='m1h_gte_2_fold_at_least_one_iso_per_gene',
              size=point_size,
            clip_on=False,
              ax=axs[4],
              linewidth=1,
               edgecolor="black",
              alpha=1,
              color=sns.color_palette("flare")[0])
axs[4].set_xticks([])
axs[4].set_xlabel('')

rfx3_y = df.loc[(df["clone_acc_alt"] == "RFX3|3/5|08G08"), 'activation_fold_change_log2'].values[0]
rfx4_y = df.loc[(df["clone_acc_alt"] == "RFX3|4/5|11D09"), 'activation_fold_change_log2'].values[0]
for point in axs[4].collections:
    for x, y in point.get_offsets():
        if np.isclose(rfx3_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate("RFX3-3", xy=(x, y), xytext=(0, -12), textcoords='offset points',
                            arrowprops = dict(arrowstyle="-", connectionstyle="arc3,rad=0.3",
                                              color='black'), 
                            ha="left", va="top", fontsize=7,
                            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))
        if np.isclose(rfx4_y, y):
            print("found: %s, %s" % (x, y))
            axs[4].annotate("RFX3-4", xy=(x, y), xytext=(-5, -8), textcoords='offset points',
                            arrowprops = dict(arrowstyle="-", connectionstyle="arc3,rad=-0.3",
                                              color='black'), 
                            ha="center", va="top", fontsize=7,
                            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'))


# add colorbar
# mirror figure
fig2, axs2 = plt.subplots(1, 5, sharey=True, gridspec_kw=gs_kw)
fig2.set_size_inches(w=8.2, h=2)
map1 = axs2[4].imshow(np.stack([t, t]), cmap="flare", vmin=25, vmax=250)
cbar = fig.colorbar(map1, ax=axs[4], aspect=30, pad=0.2)
cbar.set_ticks([25, 75, 150, 250])
cbar.set_ticklabels(["<=25", "75", "150", ">=250"])
cbar.set_label("# AA in annotated domain", labelpad=0)


for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim(-7.5, 7.5)
    ax.axhline(y=0, color='black', linewidth=1, linestyle='dashed', zorder=1)
for ax in axs[1:]:
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_tick_params(which='both', length=0)
    ax.set_ylabel("")
axs[0].set_ylabel("log2(activation fold change)")
fig.savefig('../../figures/fig4/activation_vs_domain_removal_colored_by_dom_length.collapsed.pdf', bbox_inches='tight')


# In[33]:


## calculate p-value using wilcoxon
def paired_pval(x, y):
    
    # make sure x, y are filtered for nan while maintaining pair relationship
    x_filt = []
    y_filt = []
    for x_, y_ in zip(x, y):
        if pd.isnull(x_) or pd.isnull(y_):
            continue
        else:
            x_filt.append(x_)
            y_filt.append(y_)

    try:
        stat, p = wilcoxon(x_filt, y_filt)
        return p
    except:
        return np.nan


# In[34]:


x = tot_n_part_loss_activ.activation_ref.values
y = tot_n_part_loss_activ.activation_alt.values
print("PAIRED WILCOXON P-VALUE COMPARING M1H REF V ALT FOR ALTS THAT LOSE ANNOTATED ACTIV DOMAIN: %s" % paired_pval(x, y))


# ## 6. pie chart of PPI categories

# In[35]:


len(cats)


# In[36]:


len(cats.gene_symbol_partner.unique())


# In[37]:


cats.category.value_counts()


# In[38]:


cats[(~pd.isnull(cats["cofactor_type"])) & (cats["cofactor_type"] != "unknown")]


# In[39]:


y2h_nonan = y2h[~pd.isnull(y2h["Y2H_result"])]
print("# of unique PPI partners found to interact w/ at least 1 TF iso")
len(y2h_nonan.db_gene_symbol.unique())


# In[40]:


y2h_true = y2h_nonan[y2h_nonan["Y2H_result"]][["ad_gene_symbol", "db_gene_symbol"]].drop_duplicates()
print("# of total True gene-gene PPIs profiled")
len(y2h_true)


# In[41]:


print("# of unique TF isoforms found to have at least 1 PPI")
len(y2h_nonan.ad_clone_acc.unique())


# In[42]:


print("# of unique TF genes found to have at least 1 PPI")
len(y2h_nonan.ad_gene_symbol.unique())


# In[43]:


ggi = y2h_nonan[["ad_gene_symbol", "db_gene_symbol"]].drop_duplicates()
ggi


# In[44]:


# limiting df to those that are in the y2h iso data
cats_y2h = cats[cats["gene_symbol_partner"].isin(ggi["db_gene_symbol"])]
len(cats_y2h)


# In[45]:


cats_dupe = cats_y2h.groupby("gene_symbol_partner")["category"].agg("count").reset_index()
cats_dupe[cats_dupe["category"] > 1].head()


# gene partners are now in mutually exclusive categories

# In[46]:


ys = np.array([len(cats_y2h[cats_y2h["category"] == "TF"]), len(cats_y2h[cats_y2h["category"] == "cofactor"]),
               len(cats_y2h[cats_y2h["category"] == "signaling"]),
               len(cats_y2h[cats_y2h["category"] == "other"])])
labels = ["TF", "cofactor", "signaling", "other"]
colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], sns.color_palette("Set2")[5], "darkgray"]

fig, ax = plt.subplots(figsize=(1.2, 1.2), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.0f%%', startangle=-45, 
                    explode=(0.02, 0.2, 0.05, 0.05))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")
    n.set_fontsize(6)

fig.savefig("../../figures/fig4/PPIs-gene-level-manual-categories_simplified.pdf", dpi="figure", bbox_inches="tight")


# In[47]:


cofactor_partners = set(cats_y2h.loc[cats_y2h['category'] == 'cofactor', 'gene_symbol_partner'].unique())
signaling_partners = set(cats_y2h.loc[cats_y2h['category'] == 'signaling', 'gene_symbol_partner'].unique())
other_partners = set(cats_y2h.loc[cats_y2h['category'] == 'other', 'gene_symbol_partner'].unique())
tf_gene_symbols = set(load_human_tf_db()['HGNC symbol'].values)

coactivator_partners = set(cats_y2h.loc[cats_y2h['cofactor_type'] == 'coactivator', 'gene_symbol_partner'].unique())
corepressor_partners = set(cats_y2h.loc[cats_y2h['cofactor_type'] == 'corepressor', 'gene_symbol_partner'].unique())


# In[48]:


cats_y2h.cofactor_type.value_counts()


# In[49]:


list(signaling_partners)[25:35]


# In[50]:


# make a similar pie chart but this time focusing on protein categories that get rewired across isoforms


# In[51]:


ppi = load_y2h_isoform_data(require_at_least_one_ppi_per_isoform=True)
ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])['Y2H_result'].apply(lambda x: (x == False).sum()).reset_index()


# In[52]:


a = ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])['Y2H_result'].apply(lambda x: (x == False).sum()).reset_index()
tot = ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])['Y2H_result'].apply(lambda x: (x.notnull().sum())).reset_index()
rw = tot.merge(a, on=["ad_gene_symbol", "db_gene_symbol"], how="left")
rw["rewiring_score"] = rw["Y2H_result_y"]/rw["Y2H_result_x"]
rw.sample(5)


# In[53]:


# anything w a rewiring score > 0 is rewired
rewired_ppis = rw[rw["rewiring_score"] > 0]
cats_y2h_rw = cats_y2h[cats_y2h["gene_symbol_partner"].isin(rewired_ppis["db_gene_symbol"])]
len(cats_y2h_rw)


# In[54]:


ys = np.array([len(cats_y2h_rw[cats_y2h_rw["category"] == "TF"]), 
               len(cats_y2h_rw[cats_y2h_rw["category"] == "cofactor"]),
               len(cats_y2h_rw[cats_y2h_rw["category"] == "signaling"]),
               len(cats_y2h_rw[cats_y2h_rw["category"] == "other"])])
labels = ["TF", "cofactor", "signaling", "other"]
colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], sns.color_palette("Set2")[5], "darkgray"]

fig, ax = plt.subplots(figsize=(1.2, 1.2), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.0f%%', startangle=-45, 
                    explode=(0.02, 0.2, 0.05, 0.05))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")
    n.set_fontsize(6)

fig.savefig("../../figures/fig4/PPIs-gene-level-manual-categories_simplified.rewired_only.pdf", dpi="figure", bbox_inches="tight")


# ## 7. plot the relationship between gain/loss of PPIs and changes in activity

# In[55]:


def add_restricted_ppi_columns(pairs, rows, label):
    pairs_cf = pairs[['clone_acc_ref', 'clone_acc_alt']].copy()
    _add_PPI_columns(df=pairs_cf, y2h=y2h.loc[rows, :])
    return pd.merge(pairs, 
                    pairs_cf,
                    how='left',
                    on=['clone_acc_ref', 'clone_acc_alt'],
                    suffixes=('', '_' + label))


# In[56]:


pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(cofactor_partners),
                           label='cofactors'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(cofactor_partners) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='cofactors_HEK'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(coactivator_partners),
                           label='coactivators'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(coactivator_partners) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='coactivators_HEK'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(corepressor_partners),
                           label='corepressors'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(corepressor_partners) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='corepressors_HEK'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(signaling_partners),
                           label='signaling'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(signaling_partners) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='signaling_HEK'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(tf_gene_symbols),
                           label='tfs'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(tf_gene_symbols) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='tfs_HEK'
)

pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(other_partners) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='other_HEK'
)


# In[57]:


def bar_activation_vs_ppi(x, y, pairs=pairs, x_label=None, y_label=None, color=None):
    """
    TODO:
        - calculate p-value properly
            - this requires permuting in some smart way
            - one question is whether the genes are the number of independent data points or the isoforms are
            - I think the answer is the isoforms are
    
    """
    df = pairs.copy()
    if x_label is None:
        x_label = x
    if y_label is None:
        y_label = y
    if color is None:
        color = sns.color_palette("Set2")[1]
    fig, ax = plt.subplots(1, 1, figsize=(1.15, 1.5))

    def bin_delta_ppi(delta_ppi):
        if pd.isnull(delta_ppi):
            return np.nan
        if delta_ppi < 0:
            return 'change'
        elif delta_ppi > 0:
            return 'change'
        elif delta_ppi == 0:
            return 'equal'
        else:
            raise ValueError(delta_ppi)


    df[x + '_binned'] = df[x].apply(bin_delta_ppi)
    sns.stripplot(data=df,
                  x=x + '_binned',
                  y=y,
                  order=['equal', 'change'],
                  alpha=0.75,
                  color=color,
                  linewidth=1,
                  edgecolor="black",
                  ax=ax)
    if False:
        sns.pointplot(data=df,
                    x=x + '_binned',
                    y=y,
                    order=['equal', 'change'],
                    alpha=0.5,
                    color='black',
                    ax=ax)
    if True:
        sns.boxplot(data=df,
                    x=x + '_binned',
                    y=y,
                    order=['equal', 'change'],
                    fliersize=0,
                    color=color,
                    ax=ax)
        mimic_r_boxplot(ax)
    else:
        sns.violinplot(data=df,
                    x=x + '_binned',
                    y=y,
                    order=['equal', 'change'],
                    color='lightgrey',
                    ax=ax)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    with_data = (df[x].notnull() & df[y].notnull())
    n_pair = with_data.sum()
    n_iso = len(set(df.loc[with_data, ['clone_acc_ref', 'clone_acc_alt']].values.flatten()))
    n_gene = df.loc[with_data, 'gene_symbol'].nunique()
    
    dist_a = df[df[x + '_binned'] == 'equal'][y].values
    dist_b = df[df[x + '_binned'] == 'change'][y].values
    u, p = stats.mannwhitneyu(dist_a, dist_b, alternative="less")
    plotting.annotate_pval(ax, 0.2, 0.8, 2.4, 0, 2.4, p, 7)
    
    ax.text(s=f'{n_pair:d} pairs\n{n_iso:d} isoforms\n{n_gene:d} genes\np = {p:.2f}',
            x=1.03,
            y=0.95,
            ha='left',
            va='top',
            transform=ax.transAxes)
    #ax.set_ylim(-4, 4) # NOTE cuts outlier TODO add broken axis

    for pos in ['top', 'bottom', 'right']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    fig.savefig(f'../../figures/fig4/{x}-vs-{y}_scatter.pdf',
                bbox_inches='tight')


# In[58]:


pairs['activation_abs_fold_change'] = pairs['activation_fold_change_log2'].abs()


# In[59]:


# limit to pairs w signal in m1h
df = pairs.copy()
df = df.loc[df['activation_fold_change_log2'].notnull() & df['m1h_gte_2_fold_at_least_one_iso_per_gene'], :]
len(df)


# In[60]:


bar_activation_vs_ppi(
    pairs=df.loc[df['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_other_HEK',
    y='activation_abs_fold_change',
    x_label='Net difference in number of PPIs',
    y_label='|log2(activation foldchange)|',
    color="darkgrey")


# In[61]:


bar_activation_vs_ppi(
    pairs=df.loc[df['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_cofactors_HEK',
    y='activation_abs_fold_change',
    x_label='Difference in number of cofactor PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='|log2(activation foldchange)|',
    color=sns.color_palette("Set2")[1])


# In[62]:


bar_activation_vs_ppi(
    pairs=df.loc[df['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_signaling_HEK',
    y='activation_abs_fold_change',
    x_label='Difference in number of signaling PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='|log2(activation foldchange)|',
    color=sns.color_palette("Set2")[5])


# In[63]:


bar_activation_vs_ppi(
    pairs=df.loc[df['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_tfs_HEK',
    y='activation_abs_fold_change',
    x_label='Difference in number of TF PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='|log2(activation foldchange)|',
    color=sns.color_palette("Set2")[2])


# ## 7. dimerizing TFs: plot conservation of interactions across isoforms

# In[64]:


ppi = load_y2h_isoform_data(require_at_least_one_ppi_per_isoform=True)
ppi.head()


# In[65]:


n_iso_per_ppi = (ppi.groupby(['ad_gene_symbol', 'db_gene_symbol']))['ad_clone_acc'].agg("count").reset_index()
n_true_iso_per_ppi = (ppi[ppi['Y2H_result'] == True].groupby(['ad_gene_symbol', 
                                                              'db_gene_symbol']))['ad_clone_acc'].agg("count").reset_index()
n_iso_per_ppi = n_iso_per_ppi.merge(n_true_iso_per_ppi, on=['ad_gene_symbol', 'db_gene_symbol'], how='left')
n_iso_per_ppi['f_iso_positive'] = n_iso_per_ppi['ad_clone_acc_y']/n_iso_per_ppi['ad_clone_acc_x']
n_iso_per_ppi = n_iso_per_ppi[['ad_gene_symbol', 'db_gene_symbol', 'f_iso_positive']]
n_iso_per_ppi


# In[66]:


# isoform-variable PPIs are ones that the fraction positive is < 1
iso_variable = len(n_iso_per_ppi[n_iso_per_ppi['f_iso_positive'] < 1])
print("NUMBER OF PPIs THAT VARY ACROSS ISOFORMS: %s" % iso_variable)
print("PERCENT OF PPIS THAT VARY ACROSS ISOFORMS: %s" % (iso_variable/len(n_iso_per_ppi)*100))


# In[67]:


n_iso_per_ppi = pd.merge(n_iso_per_ppi,
                         ppi.groupby(['ad_gene_symbol',
                                      'db_gene_symbol'])
                                ['Y2H_result']
                                .apply(lambda x: x.notnull().sum())
                                .rename('n_iso_successfully_tested')
                                .reset_index(),
                            how='left',
                            on=['ad_gene_symbol', 'db_gene_symbol'],
                            )
if n_iso_per_ppi['n_iso_successfully_tested'].isnull().any():
    raise UserWarning('unexpected missing values')
n_iso_per_ppi = n_iso_per_ppi.loc[n_iso_per_ppi['n_iso_successfully_tested'] >= 2, :]


# In[68]:


tf_fam = load_tf_families()


# In[69]:


n_iso_per_ppi['db_is_tf'] = n_iso_per_ppi['db_gene_symbol'].isin(tf_fam)
n_iso_per_ppi['ad_tf_family'] = n_iso_per_ppi['ad_gene_symbol'].map(tf_fam)
n_iso_per_ppi['db_tf_family'] = n_iso_per_ppi['db_gene_symbol'].map(tf_fam)


# In[70]:


def tf_tf_dimer_ppi_catagories(row):
    is_dimer_ad = row['ad_tf_family'] in DIMERIZING_TF_FAMILIES
    if pd.isnull(row['db_tf_family']):
        if is_dimer_ad:
            return 'Obligate dimer TF / other'
        else:
            return 'Non obligate dimer TF / other'
    else:  # TF-TF PPI
        if is_dimer_ad:
            if row['db_tf_family'] == row['ad_tf_family']:
                return 'Obligate dimer TF / within-family TF'
            else:
                return 'Obligate dimer TF / other'
        else:
            if row['db_tf_family'] == row['ad_tf_family']:
                return 'Non obligate dimer TF / within-family TF'
            else:
                return 'Non obligate dimer TF / other'


n_iso_per_ppi['dimer_cat'] = n_iso_per_ppi.apply(tf_tf_dimer_ppi_catagories,
                                                 axis=1)

n_iso_per_ppi.head()


# In[71]:


cats = [
 'Obligate dimer TF / within-family TF',
 'Obligate dimer TF / other',
 'Non obligate dimer TF / within-family TF',
 'Non obligate dimer TF / other',
 ]


# In[72]:


n_iso_per_ppi.dimer_cat.value_counts()


# In[73]:


for cat in cats:
    print("%s | # unique TF genes: %s" % (cat, 
                                          len(n_iso_per_ppi[n_iso_per_ppi["dimer_cat"] == cat].ad_gene_symbol.unique())))


# In[74]:


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


# In[75]:


def bootstrap_99_ci(series, n_bootstraps=1000):
    """
    Calculate the 95% confidence interval using bootstrapping.
    
    :param series: Data points as a pandas Series.
    :param n_bootstraps: Number of bootstrap iterations.
    :return: A tuple containing the lower and upper confidence interval bounds.
    """
    bootstrapped_means = []
    for _ in range(n_bootstraps):
        sample = series.sample(n=len(series), replace=True)  # Sampling with replacement
        sample_mean = sample.mean()
        bootstrapped_means.append(sample_mean)

    lower, upper = np.percentile(bootstrapped_means, [0.5, 99.5])  # 99% confidence interval
    return lower, upper


# In[76]:


# n for each bar
# axis labels etc.
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=1.7, h=1.7)
violinplot_reflected(data=n_iso_per_ppi,
              x='dimer_cat',
              y='f_iso_positive',
              inner=None,
              cut=0,
              color=sns.color_palette("Set2")[0],
              order=cats,
              alpha=0.75,
              linewidth=0.5,
              ax=ax)

# add markers for mean + 95% ci [avoids using pointplot which is not that customizable, aesthetically]
means = n_iso_per_ppi.groupby("dimer_cat")["f_iso_positive"].agg("mean").reset_index()
cis = n_iso_per_ppi.groupby("dimer_cat")["f_iso_positive"].apply(bootstrap_99_ci).reset_index()
cis[['lower', 'upper']] = pd.DataFrame(cis['f_iso_positive'].tolist(), index=cis.index)

# reorder so they're in the same order
cis.set_index("dimer_cat", inplace=True)
means_tmp = means.set_index("dimer_cat")
cis = cis.loc[cats, :]
means_tmp = means_tmp.loc[cats, :]

# add points
sns.swarmplot(data=means, x='dimer_cat', y='f_iso_positive', order=cats, ax=ax,
              color="white", edgecolor="black", linewidth=0.25, size=4)
ax.errorbar(x=range(len(cats)), y=means_tmp["f_iso_positive"], 
            yerr=[means_tmp["f_iso_positive"] - cis["lower"], cis["upper"] - means_tmp["f_iso_positive"]], 
            fmt='o', color="black")

# now add p-vals
pairs_to_test = [(cats[0], cat) for cat in cats[1:]]
for i, pair in enumerate(pairs_to_test):
    x = n_iso_per_ppi[n_iso_per_ppi["dimer_cat"] == pair[0]]["f_iso_positive"].values
    y = n_iso_per_ppi[n_iso_per_ppi["dimer_cat"] == pair[1]]["f_iso_positive"].values
    p = permutation_test(x, y, num_permutations=10000, alternative="two-sided", seed=2023)
    print(p)
    annotate_pval(ax, 0, 1+i, 1.02+(i*0.08), 0, 1.02+(i*0.08), p, fontsize-1, text="< 1e-4")
    print("pair: %s | p-val: %s" % (pair, p))

ax.set_ylabel('Fraction isoforms interacting', position=(0, 0.425))
ax.set_xlabel('')

ax.set_xticklabels(cats, rotation=30, ha='right', va='top')

ax.set_ylim(0, 1.2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# manually set left axis so it stops at 1.0
ax.spines['left'].set_visible(False)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
axes_to_data = ax.transAxes + ax.transData.inverted()
left_spine_in_data_coords = axes_to_data.transform((0, 0))
ax.plot([left_spine_in_data_coords[0], left_spine_in_data_coords[0]], [left_spine_in_data_coords[0], 1],
        color=ax.spines['bottom'].get_edgecolor(), linewidth=ax.spines['bottom'].get_linewidth())



fig.savefig('../../figures/fig4/n_iso_ppi_by_dimer_cat.pdf',
            bbox_inches='tight')


# In[77]:


n_iso_per_ppi[n_iso_per_ppi["ad_gene_symbol"] == "ATF2"]


# ## 8. dimerizing TFs: plot M1H change

# In[78]:


len(pairs)


# In[79]:


len(pairs[~pd.isnull(pairs['activation_abs_fold_change'])])


# In[80]:


# limit to pairs w signal in m1h
df = pairs.copy()
df = df.loc[df['activation_fold_change_log2'].notnull() & df['m1h_gte_2_fold_at_least_one_iso_per_gene'], :]
len(df)


# In[81]:


df[~pd.isnull(df['activation_abs_fold_change'])].groupby("is_dimerizing_TF_family")["activation_abs_fold_change"].agg("mean")


# In[82]:


df[~pd.isnull(df['activation_abs_fold_change'])].groupby("is_dimerizing_TF_family")["activation_abs_fold_change"].agg("median")


# In[83]:


df_nonan = df[~pd.isnull(df["activation_abs_fold_change"])]
x = df_nonan[df_nonan["is_dimerizing_TF_family"] == True]["activation_abs_fold_change"].values
y = df_nonan[df_nonan["is_dimerizing_TF_family"] == False]["activation_abs_fold_change"].values

permutation_test(x, y, num_permutations=10000, alternative="two-sided")


# In[84]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(w=1, h=2)

ax = sns.boxplot(data=df, x="is_dimerizing_TF_family", y="activation_abs_fold_change", 
                 flierprops={"marker": "o"}, color="slategrey")
mimic_r_boxplot(ax)

# p-value
df_nonan = df[~pd.isnull(df["activation_abs_fold_change"])]
x = df_nonan[df_nonan["is_dimerizing_TF_family"] == True]["activation_abs_fold_change"].values
y = df_nonan[df_nonan["is_dimerizing_TF_family"] == False]["activation_abs_fold_change"].values

pval = permutation_test(x, y, num_permutations=10000, alternative="two-sided")
plotting.annotate_pval(ax, 0.2, 0.8, 4, 0, 4, pval, 7)

ax.set_xlabel("")
ax.set_xticklabels(["non-dimerizing TFs", "dimerizing TFs"], rotation=30, ha="right", va="top")
ax.set_ylabel("|log2 fold-change in activation|")

fig.savefig("../../figures/fig4/M1H_activ_v_dimerization.pdf", dpi="figure", bbox_inches="tight")


# ## 9. dimerizing TFs: retained interactions heatmap

# In[85]:


# TF isoforms are on the 'ad' side; ORFeome on 'db' side
ppi.head()


# In[86]:


# TF-TF binding
tftf = ppi.loc[ppi['Y2H_result'] == True, ['ad_gene_symbol', 'db_gene_symbol']].drop_duplicates().copy()
tftf['ad_dbd'] = tftf['ad_gene_symbol'].map(tf_fam)
tftf['db_dbd'] = tftf['db_gene_symbol'].map(tf_fam)
tftf = tftf.dropna()
tftf.head()


# In[87]:


# TF-TF rewiring
tftf = pd.merge(tftf, 
        (ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])['Y2H_result'].apply(lambda x: (x == False).sum()) / 
 ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])['Y2H_result'].apply(lambda x: (x.notnull().sum()))).reset_index(),
        how='left',
        on=['ad_gene_symbol', 'db_gene_symbol'])


# In[88]:


tftf.loc[tftf['db_dbd'] == 'bHLH', 'ad_gene_symbol'].value_counts()


# In[89]:


DIMERIZING_TF_FAMILIES


# In[90]:


fams = tf_fam.value_counts().index
fams = list(filter(lambda x: x in tftf['ad_dbd'].unique() or x in tftf['db_dbd'].unique(), fams))
len(fams)


# In[91]:


num_pairs = [((tftf['ad_dbd'] == x) & (tftf['db_dbd'] == y)).sum() for x in fams for y in fams]
num_pairs = np.reshape(num_pairs, (-1, len(fams)))
num_pairs = pd.DataFrame(num_pairs, index=fams, columns=fams)
num_pairs.head()


# In[92]:


rewiring = [tftf.loc[(tftf['ad_dbd'] == x) & (tftf['db_dbd'] == y), 'Y2H_result'].mean() for x in fams for y in fams]
rewiring = np.reshape(rewiring, (-1, len(fams)))
rewiring = pd.DataFrame(rewiring, index=fams, columns=fams)
rewiring.head()


# In[93]:


# since rewiring score is 1-fraction isos interacting, we decided to just go with that
rewiring = 1-rewiring
rewiring.head()


# In[94]:


tftf[(tftf["ad_dbd"] == "Homeodomain") & (tftf["db_dbd"] == "bHLH")].mean()


# In[95]:


len(tftf[(tftf["ad_dbd"] == "Homeodomain") & (tftf["db_dbd"] == "bHLH")])


# In[96]:


tftf[(tftf["ad_dbd"] == "bHLH") & (tftf["db_dbd"] == "Homeodomain")].mean()


# In[97]:


len(tftf[(tftf["ad_dbd"] == "bHLH") & (tftf["db_dbd"] == "Homeodomain")])


# confirming that AD is on the rows and DB is on the columns as expected

# In[98]:


# limit to fams w 5 DB and remove 'unknown' DBDs
num_pairs_sum = num_pairs.sum(axis=0)
filt_fams = num_pairs_sum[num_pairs_sum >= 5]
filt_fams = filt_fams[filt_fams.index != "Unknown"]
filt_fams


# In[99]:


rewiring_filt = rewiring.loc[list(filt_fams.index), list(filt_fams.index)]
rewiring_filt


# In[100]:


fig = plt.figure(figsize=(2.1, 1.8))

g = sns.heatmap(rewiring_filt, cmap="viridis", cbar_kws={"label": "mean fraction isoforms interacting"})

# highlight the squares corresponding to dimerizing pairs
for i, fam in enumerate(list(filt_fams.index)):
    if fam in DIMERIZING_TF_FAMILIES:
        g.add_patch(Rectangle((i, i), 1, 1, fill=False, edgecolor='black', lw=1))
        
g.set_ylabel("TF isoform family")
g.set_xlabel("Y2H partner TF family")
g.xaxis.tick_top()
g.xaxis.set_label_position('top')
g.set_xticklabels(g.get_xticklabels(), rotation=90, ha="left", va="bottom")
fig.savefig("../../figures/fig4/Dimerizing_PPI_heatmap.small.pdf", dpi="figure", bbox_inches="tight")


# In[101]:


filt_fams = list(rewiring_filt.columns)
filt_fams


# In[102]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(2, 2)
scaling = 5
sc = ax.scatter(x=[x for x in filt_fams for y in filt_fams],
                y=[y for x in filt_fams for y in filt_fams],
                s=[((tftf['db_dbd'] == x) & (tftf['ad_dbd'] == y)).sum() * scaling for x in filt_fams for y in filt_fams],
                c=[tftf.loc[(tftf['db_dbd'] == x) & (tftf['ad_dbd'] == y), 'Y2H_result'].mean() for x in filt_fams for y in filt_fams],
                cmap='viridis_r',
                vmin=0,
                vmax=7/8,
                clip_on=False)

dm = ax.scatter(x=[x for x in filt_fams for y in filt_fams if x in DIMERIZING_TF_FAMILIES and y in DIMERIZING_TF_FAMILIES and x == y],
                y=[y for x in filt_fams for y in filt_fams if x in DIMERIZING_TF_FAMILIES and y in DIMERIZING_TF_FAMILIES and x == y],
                s=[((tftf['db_dbd'] == x) & (tftf['ad_dbd'] == y)).sum() * scaling for x in filt_fams for y in filt_fams if x in DIMERIZING_TF_FAMILIES and y in DIMERIZING_TF_FAMILIES and x == y],
                c=[tftf.loc[(tftf['db_dbd'] == x) & (tftf['ad_dbd'] == y), 'Y2H_result'].mean() for x in filt_fams for y in filt_fams if x in DIMERIZING_TF_FAMILIES and y in DIMERIZING_TF_FAMILIES and x == y],
                cmap='viridis_r',
                vmin=0,
                vmax=7/8,
                edgecolor="black",
                linewidth=1,
                clip_on=False)

ax.xaxis.set_tick_params(rotation=90)
for s in ax.spines.values():
    s.set_visible(False)
    
ax.set_xlim(-0.5, len(filt_fams) - 0.5)
ax.set_ylim(len(filt_fams) - 0.5, -0.5)

ax.set_xlabel('Y2H partner TF family',
              fontsize=fontsize,
              labelpad=5)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_ylabel('TF isoform family',
              fontsize=fontsize,
              labelpad=5)

ax.legend(*sc.legend_elements('sizes',
                              num=[1, 10, 50],
                              func=lambda x: x / scaling),
          bbox_to_anchor=[1.01, 1],
          labelspacing=1.1,
          title='# PPIs',
          frameon=False,
          title_fontsize=fontsize)

cbax = ax.inset_axes([1.07, 0, 0.05, 0.42], transform=ax.transAxes)
fig.colorbar(mpl.cm.ScalarMappable(cmap='viridis', norm=mpl.colors.Normalize(vmin=0, vmax=7/8)),
             ax=ax,
             cax=cbax,
             orientation='vertical',
             fraction=0.02,
             label='fraction of isoforms\ninteracting (mean)')
fig.savefig('../../figures/fig4/Dimerizing_PPI_correlogram.small.pdf', bbox_inches='tight')


# In[103]:


fig = plt.figure(figsize=(4.5, 4.5))
g = sns.heatmap(rewiring.T, cmap="viridis", cbar_kws={"label": "mean fraction of isoforms interacting", 
                                                      "aspect": 40},
                annot=num_pairs.T, annot_kws={"fontsize": 6})

# highlight the squares corresponding to dimerizing pairs
for i, fam in enumerate(list(rewiring.index)):
    if fam in DIMERIZING_TF_FAMILIES:
        g.add_patch(Rectangle((i, i), 1, 1, fill=False, edgecolor='black', lw=1))
        
g.set_ylabel("TF isoform family")
g.set_xlabel("Y2H partner TF family")
g.xaxis.tick_top()
g.xaxis.set_label_position('top')
g.set_xticklabels(g.get_xticklabels(), rotation=90, ha="left", va="bottom")
fig.savefig("../../figures/fig4/Dimerizing_PPI_heatmap.all.pdf", dpi="figure", bbox_inches="tight")


# ## 10. domain-domain interactions plot

# In[104]:


def load_3did_DDIs():
    fpath = '../../data/external/3did_flat_2022-05.txt'
    domain_pairs = []
    for line in open(fpath, 'r'):
        if line.startswith('#=ID'):
            domain_pairs.append((line.split()[3][1:8], line.split()[4][:7]))
    df = pd.DataFrame(data=domain_pairs, columns=['pfam_a', 'pfam_b'])
    if df.duplicated().any():
        raise UserWarning('unexpected duplicates')
    return df


# In[105]:


pfam = load_pfam_domains_horfeome()
pfam['orf_id'] = pfam['orf_id'].astype(int)
pfam.head()


# In[106]:


# TODO: add in pfam AC for ZF for ZF array
tf_pfam_domains = {tf.name: {dom.accession for iso in tf.cloned_isoforms
                             for dom in iso.aa_seq_features 
                             if dom.category == 'Pfam_domain'}
                   for tf in tfs.values()}


# In[107]:


print(len(ppi))

ppi['gene_level_pair'] = ppi['ad_gene_symbol'] + '_' + ppi['db_gene_symbol']
ppi['ad_iso_id'] = ppi['ad_clone_acc'].apply(lambda x: x.split('|')[0] + '-' + x.split('|')[1].split('/')[0])
ppi.head()


# In[108]:


# code requires the orf id so map this using the full df
ppi_full = load_full_y2h_data_including_controls()
ppi = ppi.merge(ppi_full[['db_gene_symbol', 'db_orf_id']], on='db_gene_symbol', how='left')
ppi.sample(5)


# In[109]:


# check no-overlap with domain motif
ddi = load_3did_DDIs()

def matching_DDIs(row):
    if (
        pd.isnull(row['partner_domains']) or
        pd.isnull(row['tf_domains'])
        ):
        return np.nan
    if (
        len(row['partner_domains']) == 0 or
        len(row['tf_domains']) == 0
    ):
        return np.nan
    matches = ddi.loc[(ddi['pfam_a'].isin(row['partner_domains']) &
                    ddi['pfam_b'].isin(row['tf_domains'])) |
                   (ddi['pfam_a'].isin(row['tf_domains']) &
                    ddi['pfam_b'].isin(row['partner_domains']))].values
    if len(matches) == 0:
        return np.nan
    return frozenset((a, b) for a, b in matches)

partner_domains = pfam.groupby('orf_id')['pfam_accession'].apply(set)
ppi['partner_domains'] = ppi['db_orf_id'].map(partner_domains)
ppi['tf_domains'] = ppi['ad_gene_symbol'].map(tf_pfam_domains)
ppi['matching_DDI'] = ppi.apply(matching_DDIs, axis=1)
ppi.head()


# In[110]:


ggi = ppi.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'db_orf_id', 'matching_DDI']].drop_duplicates()
print(len(ggi))
ggi.head()


# In[111]:


# domain removal
# for each alt iso, for each DDI, calc fraction of domain
# removed and fraction of PPIs retained
dom = pd.concat([g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()])


# In[112]:


# filter and consolidate DDIs

ddi_annot = pd.read_csv('../../data/internal/DDI_manual_annotation.tsv', sep='\t')

valid_ddi_pairs = {frozenset(p) for p in ddi_annot.loc[ddi_annot['to_use'], ['pfam_a', 'pfam_b']].values}
def filter_ddi(pairs):
    if pd.isnull(pairs):
        return np.nan
    passed = frozenset(p for p in pairs if frozenset(p) in valid_ddi_pairs)
    if len(passed) == 0:
        return np.nan
    return passed

ddi_to_merge = {
    ('PF00170', 'PF07716'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF00170', 'PF03131'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF07716', 'PF07716'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF03131', 'PF03131'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF00046', 'PF05920'): ('PF00046',	'PF00046'),  # homeobox-homeobox
    ('PF05920', 'PF05920'): ('PF00046',	'PF00046'),  # homeobox-homeobox
    ('PF00989', 'PF14598'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF00989', 'PF08447'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF00989', 'PF13426'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF14598', 'PF14598'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF08447', 'PF14598'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF08447', 'PF08447'): ('PF00989', 'PF00989'),  # PAS-PAS
}


# In[113]:


df = ppi.loc[ppi['matching_DDI'].notnull(), :].copy()
ref_isos = {tf.cloned_reference_isoform.clone_acc for tf in tfs.values()}
positive_in_ref = df.loc[df['ad_clone_acc'].isin(ref_isos) & 
                         (df['Y2H_result'] == True), 'gene_level_pair'].unique()
df = df.loc[df['gene_level_pair'].isin(positive_in_ref) &
            ~df['ad_clone_acc'].isin(ref_isos), :]
df['matching_DDI'] = df['matching_DDI'].apply(filter_ddi)
df['matching_DDI'] = df['matching_DDI'].apply(lambda x: frozenset(ddi_to_merge.get(tuple(sorted(p)), tuple(sorted(p))) for p in x) if pd.notnull(x) else x)
df = df.loc[df['matching_DDI'].notnull(), :]


# In[114]:


merged_domains = {'PF07716': 'PF00170',  # bZIP
 'PF03131': 'PF00170',  # bZIP
 'PF05920': 'PF00046',  # homeobox
'PF14598': 'PF00989',  # PAS
'PF08447': 'PF00989',  # PAS
'PF13426': 'PF00989',  # PAS
}
df['tf_domains'] = df['tf_domains'].apply(lambda x: {merged_domains.get(d, d) for d in x})
dom['accession'] = dom['accession'].apply(lambda x: merged_domains.get(x, x))


def pick_the_one_domain(row):
    ddi_domains = {x for a, b in row['matching_DDI'] for x in [a, b]}
    ds = {d for d in row['tf_domains'] if d in ddi_domains}
    if len(ds) == 0:
        print(row)
        raise UserWarning('something wrong')
    return ds

df['tf_domains'] = df.apply(pick_the_one_domain, axis=1)

def filter_for_domain_in_cloned_reference_isoform(row):
    in_ref = {d.accession for d in tfs[row['ad_gene_symbol']].cloned_reference_isoform.aa_seq_features}
    return {d for d in row['tf_domains'] if d in in_ref}


df['tf_domains'] = df.apply(filter_for_domain_in_cloned_reference_isoform, axis=1)
df = df.loc[df['tf_domains'].map(lambda x: len(x) > 0), :]

def fraction_of_DDI_domains_removed(row):
    ds = dom.loc[(dom['alt_iso'] == row['ad_iso_id']) 
                  & dom['accession'].isin(row['tf_domains']), :]
    if ds.shape[0] == 0:
        print(row)
        raise UserWarning('something wrong')
    return ds[['deletion', 'frameshift']].sum().sum() / ds['length'].sum()

def insertion_in_DDI_domains(row):
    ds = dom.loc[(dom['alt_iso'] == row['ad_iso_id']) 
                  & dom['accession'].isin(row['tf_domains']), :]
    return ds['insertion'].sum()

df['fraction_of_DDI_domains_removed'] = df.apply(fraction_of_DDI_domains_removed, axis=1)
df['insertion_in_DDI_domains'] = df.apply(insertion_in_DDI_domains, axis=1)

df['tf_domains'] = df['tf_domains'].apply(frozenset)


# In[115]:


data = (df.loc[df['Y2H_result'].notnull(), :]
          .groupby(['ad_iso_id', 'tf_domains'])
          ['fraction_of_DDI_domains_removed']
          .mean()
          .to_frame())

nonan = df.loc[df['Y2H_result'].notnull(), :]
nonan['Y2H_result'] = nonan['Y2H_result'].astype(int)

data['Y2H_result_mean'] = nonan.groupby(['ad_iso_id', 'tf_domains'])['Y2H_result'].mean()
data['insertion_in_DDI_domains'] = nonan.groupby(['ad_iso_id', 'tf_domains'])['insertion_in_DDI_domains'].mean()


# In[116]:


# add distance from domain
# TODO move to isolib
def n_aa_change_from_feature(gene, ref_iso_name, alt_iso_name, domain_start, domain_end):
    algn = gene.pairwise_changes_relative_to_reference(ref_iso_name, alt_iso_name)

    def _coords_transform_aa_seq_to_alignment(i, alignment):
        if i > len(alignment.replace("I", "")):
            raise ValueError("position is not in isoform AA sequence")
        aa_seq_indices = [
            "" if c == "I" else len(alignment[:j].replace("I", ""))
            for j, c in enumerate(alignment)
        ]
        return aa_seq_indices.index(i)
    
    start = _coords_transform_aa_seq_to_alignment(domain_start, algn)
    end = _coords_transform_aa_seq_to_alignment(domain_end - 1, algn) + 1

    if not all(x == 'M' for x in algn[start:end]):
        return 0  # change is within the domain
    
    big_number = 9999999999999999999999999
    c_dist = big_number
    n_dist = big_number
    for i, l in enumerate(reversed(algn[:start])):
        if l != 'M':
            c_dist = i + 1
            break
    for i, l in enumerate(algn[end:]):
        if l != 'M':
            n_dist = i + 1
            break
    if c_dist == big_number and n_dist == big_number:
        raise UserWarning('problem calculating distance')
    return min([c_dist, n_dist])


def n_aa_to_all_features(self, ref_iso_name):
    results = []
    ref_iso = self._iso_dict[ref_iso_name]
    row = {"gene": self.name, "ref_iso": ref_iso_name}
    for aa_feature in ref_iso.aa_seq_features:
        for alt_iso_name, alt_iso in self._iso_dict.items():
            if alt_iso_name == ref_iso_name:
                continue
            row.update(
                {
                    "alt_iso": alt_iso_name,
                    "accession": aa_feature.accession,
                    "category": aa_feature.category,
                    "start_in_ref_iso": aa_feature.start,
                    "end_in_ref_iso": aa_feature.end,
                    "length": aa_feature.end - aa_feature.start,
                }
            )
            row.update({"n_aa_change_to_domain": n_aa_change_from_feature(self, ref_iso_name, alt_iso_name, aa_feature.start, aa_feature.end)})
            results.append(row.copy())
    results = pd.DataFrame(results)
    return results

dist = pd.concat([n_aa_to_all_features(g, g.cloned_reference_isoform.name) for g in tfs.values()])

def get_dist(row):
    alt_iso, dom_accessions = row.name
    return dist.loc[(dist['alt_iso'] == alt_iso) & (dist['accession'].isin(dom_accessions)),
                    'n_aa_change_to_domain'].min()

data['domain_n_aa_to_change'] = data.apply(get_dist, axis=1)


# In[117]:


COLOR_PURPLE = (155 / 255, 97 / 255, 153 / 255)


gs_kw = dict(width_ratios=[0.7, 1, 1.9])

fig, axs = plt.subplots(1, 3, sharey=True, gridspec_kw=gs_kw)
fig.set_size_inches(w=4, h=1.5)

point_size = 6

axs[0].set_title('Full loss\nof domain',
fontsize=PAPER_FONTSIZE)
sns.swarmplot(data=data.loc[data['fraction_of_DDI_domains_removed'] == 1, :],
              y='Y2H_result_mean', 
              x='fraction_of_DDI_domains_removed',
              size=point_size,
         #     order=[
         #            'Full loss\nof DBD',
         #            ],
            clip_on=False,
              ax=axs[0],
              color=COLOR_PURPLE,
              edgecolor='black',
              linewidth=1,
              alpha=1)

axs[1].set_title('Partial loss\nof domain',
fontsize=PAPER_FONTSIZE)
partial_loss = (data['fraction_of_DDI_domains_removed'] > 0) & (data['fraction_of_DDI_domains_removed'] < 1)
axs[1].scatter(data.loc[partial_loss, 'fraction_of_DDI_domains_removed'].values,
               data.loc[partial_loss, 'Y2H_result_mean'].values,
           alpha=1,
           s=point_size**2,  
            color=COLOR_PURPLE,
               edgecolor='black',
               linewidth=1,
           clip_on=False)
axs[1].set_xlabel('Proportion missing')
axs[1].set_xlim(1, 0)
axs[1].set_xticks([0.99, 0.5, 0.01])
axs[1].set_xticklabels([f'{x:.0%}' for x in axs[1].get_xticks()])
#axs[1].set_xticks(range(10, 91, 10), minor=True)



axs[2].set_title('Full domain in\nalternative isoform', fontsize=PAPER_FONTSIZE)
axs[2].scatter(data.loc[(data['fraction_of_DDI_domains_removed'] == 0), 'domain_n_aa_to_change'].values,
               data.loc[(data['fraction_of_DDI_domains_removed'] == 0), 'Y2H_result_mean'].values,
           alpha=1,
           s=point_size**2,
            color=COLOR_PURPLE,
               edgecolor='black',
               linewidth=1,
           clip_on=False)
axs[2].set_xlabel('Distance of alternative\nsequence from domain\n(number of AA)')

for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim(0, 1)
for ax in axs[1:]:
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_tick_params(which='both', length=0)
for i in [0]:
    axs[i].set_xlabel('')
    axs[i].set_ylabel('')
    axs[i].spines['bottom'].set_visible(False)
    axs[i].xaxis.set_tick_params(length=0)
    axs[i].set_xticks([])
axs[0].set_yticks([0, 0.25, 0.5, 0.75, 1])
axs[0].set_yticks(np.linspace(0, 1, 21), minor=True)
axs[0].set_yticklabels(['{:.0%}'.format(y) for y in axs[0].get_yticks()])
axs[0].set_ylabel('Fraction of domain-domain mediated\nPPIs with alternative isoform')
fig.savefig('../../figures/fig4/PPI_vs_domain_removal.pdf', bbox_inches='tight')


# In[118]:


print("NUMBER OF DOMAIN-DOMAIN PPIS: %s" % len(ggi[~pd.isnull(ggi["matching_DDI"])]))
print("PERCENT OF DOMAIN-DOMAIN PPIS: %s" % (len(ggi[~pd.isnull(ggi["matching_DDI"])])/len(ggi)*100))
len(ggi)


# ## 11. isoform example vignettes

# In[119]:


# reload data since we edited dfs above
y2h = load_y2h_isoform_data()
y1h = load_y1h_pdi_data(add_missing_data=True)
m1h = load_m1h_activation_data(add_missing_data=True)
tfs = load_annotated_TFiso1_collection()


# ### RFX3

# In[120]:


gene_name = "RFX3"


# In[121]:


fig, ax = plt.subplots(1, 1, figsize=(1, 0.8))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(0, 6))
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[122]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[123]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[124]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# ### PBX1

# In[125]:


gene_name = "PBX1"


# In[126]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../../figures/fig4/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[127]:


cats_y2h[cats_y2h["gene_symbol_partner"].isin(["PIN1", "TMF1"])]


# In[128]:


pairs[pairs["gene_symbol"] == "PBX1"]


# In[129]:


fig, ax = plt.subplots(1, 1, figsize=(2, 0.5))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(-0.1, 3))
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[130]:


fig, ax = plt.subplots(figsize=(5, 2))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[131]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# ### CREB5

# In[132]:


gene_name = "CREB5"


# In[133]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(0.8, 0.8))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../../figures/fig4/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[134]:


fig, ax = plt.subplots(1, 1, figsize=(2, 0.6))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(-2.2, 2.2))
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[135]:


fig, ax = plt.subplots(figsize=(5, 1))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[136]:


fig, ax = plt.subplots(figsize=(5, 0.7))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[137]:


tfs["CREB5"]["CREB5-204"].exons


# In[138]:


tfs["CREB5"]["CREB5-202"]


# ### DLX1

# In[139]:


gene_name = "DLX1"
fig, ax = plt.subplots(figsize=(3, 0.6))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[140]:


fig, ax = plt.subplots(figsize=(5, 0.7))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[141]:


fig, ax = plt.subplots(1, 1, figsize=(1, 0.5))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(0, 4.2))
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ### ATF2

# In[142]:


gene_name = "ATF2"


# In[143]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(1.25, 1.25))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../../figures/fig4/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[144]:


fig, ax = plt.subplots(figsize=(6, 2))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[145]:


fig, ax = plt.subplots(figsize=(6, 2))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# ### TBX5

# In[146]:


gene_name = "TBX5"


# In[147]:


fig, ax = plt.subplots(figsize=(4, 1))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[148]:


fig, ax = plt.subplots(1, 1, figsize=(1, 0.6))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(-0.5, 4.1))
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ### TGIF1

# In[149]:


gene_name = "TGIF1"


# In[150]:


fig, ax = plt.subplots(figsize=(4, 1.5))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[151]:


fig, ax = plt.subplots(figsize=(4, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[152]:


fig, ax = plt.subplots(1, 1, figsize=(2, 0.6))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(0, 4.1))
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ### E2F3

# In[153]:


gene_name = "E2F3"


# In[154]:


fig, ax = plt.subplots(figsize=(4, 1.5))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[155]:


fig, ax = plt.subplots(figsize=(4, 1))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[156]:


fig, ax = plt.subplots(1, 1, figsize=(1, 0.6))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(-1.5, 7))
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# ### FOXP3

# In[157]:


gene_name = "FOXP3"


# In[158]:


fig, ax = plt.subplots(figsize=(4, 2))

tfs[gene_name].exon_diagram(ax=ax)
fig.savefig("../../figures/fig4/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[159]:


fig, ax = plt.subplots(figsize=(4, 2))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig4/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[160]:


fig, ax = plt.subplots(1, 1, figsize=(1, 1.2))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(-2, 2), iso_order=["FOXP3-3", "FOXP3-1",
                                                                                          "FOXP3-4", "FOXP3-5",
                                                                                          "FOXP3-6"])
plt.savefig('../../figures/fig4/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[ ]:




