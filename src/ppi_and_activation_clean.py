
# coding: utf-8

# In[1]:


# TODO
# p-value / statistic
# finalize plot


import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

from plotting import mimic_r_boxplot

from data_loading import (load_y2h_isoform_data, 
    load_m1h_activation_data, 
    load_ppi_partner_categories, 
    load_valid_isoform_clones,
    load_human_tf_db,
    load_y1h_pdi_data)
from isoform_pairwise_metrics import pairs_of_isoforms_comparison_table


# In[2]:


pd.set_option('display.max_columns', 100)


# In[3]:


y2h = load_y2h_isoform_data()
m1h = load_m1h_activation_data()
cats = load_ppi_partner_categories()


# In[4]:


cof = pd.read_csv('../data/external/AnimalTFDB3_Homo_sapiens_TF_cofactors.txt',
                 sep='\t')
if cof['Symbol'].duplicated().any():
    raise UserWarning('unexpected duplicates')


# In[5]:


cof.tail()


# In[6]:


y2h.head()


# In[7]:


# number of cofactors
# number of PPI with cofactors
# split by family
print(cof.shape[0], 'cofactors')
print(y2h.loc[y2h['db_gene_symbol'].isin(cof['Symbol']), 'db_gene_symbol'].nunique(),
      ' of which we have PPIs for')
print(y2h['db_gene_symbol'].isin(cof['Symbol']).sum(),
      ' cofactor PPIs')
cof_with_ppi = set(y2h.loc[y2h['db_gene_symbol'].isin(cof['Symbol']), 'db_gene_symbol'].unique())
print(cof.loc[cof['Symbol'].isin(cof_with_ppi), 'Family'].value_counts())


# In[8]:


df = pd.read_excel('../data/external/Geiger-et-al_MCP_2012_Supplementary-Table-2.xlsx',
                   skiprows=1)
hek_avrg = df[['iBAQ HEK293_1', 'iBAQ HEK293_2', 'iBAQ HEK293_3']].mean(axis=1)
print((hek_avrg > 0).sum(), 'proteins expressed in HEK293 proteome')
hek_expressed_genes = set(df.loc[(hek_avrg > 0) & df['Gene Names'].notnull(),
       'Gene Names'].str.split(';').explode().values)
all_partners = set(y2h['db_gene_symbol'].unique())
print('of {} PPI partners, {} are expressed in HEK293 cells'.format(len(all_partners), 
      len(all_partners.intersection(hek_expressed_genes))))


# ## pie charts

# In[17]:


len(cats)


# In[16]:


len(cats.partner.unique())


# In[22]:


cats.category.value_counts()


# In[35]:


y2h_nonan = y2h[~pd.isnull(y2h["Y2H_result"])]
len(y2h_nonan.db_gene_symbol.unique())


# In[37]:


ggi = y2h_nonan[["ad_gene_symbol", "db_gene_symbol"]].drop_duplicates()
ggi


# In[38]:


# limiting df to those that are in the y2h iso data
cats_y2h = cats[cats["partner"].isin(ggi["db_gene_symbol"])]
len(cats_y2h)


# In[40]:


cats_dupe = cats_y2h.groupby("partner")["category"].agg("count").reset_index()
cats_dupe[cats_dupe["category"] > 1].head()


# In[41]:


# keeping partners that are assigned more than 1 category for now...
cats_y2h[cats_y2h["partner"] == "BARD1"]


# In[42]:


def categorize_PPI_partner(row):
    if row['category'] == "TF":
        return 'TF'
    elif row['category'] == "cofactor":
        return 'cofactor'
    elif row['category'] == "signaling":
        return 'signaling'
    else:
        return 'other'
    
cats_y2h["partner_category"] = cats_y2h.apply(categorize_PPI_partner, axis=1)
cats_y2h.partner_category.value_counts()


# In[83]:


ys = np.array([len(cats[cats["partner_category"] == "TF"]), len(cats[cats["partner_category"] == "cofactor"]),
               len(cats[cats["partner_category"] == "signaling"]),
              len(cats[cats["partner_category"] == "other"])])
labels = ["TF", "cofactor", "signaling", "other"]
colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], sns.color_palette("Set2")[5], "darkgray"]

fig, ax = plt.subplots(figsize=(1.6, 1.6), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.0f%%', startangle=-45, 
                    explode=(0.02, 0.2, 0.05, 0.05))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")

fig.savefig("../figures/PPIs-gene-level-manual-categories_simplified.pdf", dpi="figure", bbox_inches="tight")


# In[53]:


cofactor_partners = set(cats.loc[cats['category'] == 'cofactor', 'partner'].unique())
coactivator_partners = set(cats.loc[cats['cofactor_type'] == 'coactivator', 'partner'].unique())
corepressor_partners = set(cats.loc[cats['cofactor_type'] == 'corepressor', 'partner'].unique())
both_partners = set(cats.loc[cats['cofactor_type'] == 'both', 'partner'].unique())
signaling_partners = set(cats.loc[cats['category'] == 'signaling', 'partner'].unique())
cofactor_animal_db = set(cof['Symbol'].unique())
tf_gene_symbols = set(load_human_tf_db()['HGNC symbol'].values)


# In[50]:


def categorize_cofactors(row):
    if row["partner_category"] == "cofactor":
        if row["partner"] in both_partners:
            return "both"
        elif row["partner"] in coactivator_partners:
            return "coactivator"
        elif row["partner"] in corepressor_partners:
            return "corepressor"
        else:
            return "unknown"
    else:
        return "NA"
    
cats_y2h["cofactor_type"] = cats_y2h.apply(categorize_cofactors, axis=1)
cats_y2h.cofactor_type.value_counts()


# In[72]:


cofacs = cats_y2h[cats_y2h["partner_category"] == "cofactor"]

ys = np.array([len(cofacs[cofacs["cofactor_type"] == "coactivator"]), 
               len(cofacs[cofacs["cofactor_type"] == "corepressor"]),
               len(cofacs[cofacs["cofactor_type"] == "both"]),
               len(cofacs[cofacs["cofactor_type"] == "unknown"])])
labels = ["coactivator", "corepressor", "both", "unknown"]
colors = [sns.color_palette("Set2")[0], sns.color_palette("Set2")[3], sns.color_palette("Set2")[4], "darkgray"]

fig, ax = plt.subplots(figsize=(1.9, 1.9), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.0f%%', startangle=60, 
                    explode=(0.05, 0.05, 0.05, 0.1))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")

fig.savefig("../figures/PPIs-gene-level-manual-categories_cofactors.pdf", dpi="figure", bbox_inches="tight")


# In[54]:


pairs = pairs_of_isoforms_comparison_table(isoforms=load_valid_isoform_clones(),
                                           y2h=y2h,
                                           m1h=m1h,
                                           y1h=load_y1h_pdi_data())


def add_restricted_ppi_columns(pairs, rows, label):
    pairs_cf = pairs_of_isoforms_comparison_table(isoforms=load_valid_isoform_clones(),
                                                  y2h=y2h.loc[rows, :])
    return pd.merge(pairs, 
                    pairs_cf.loc[:, [c for c in pairs_cf.columns if c.startswith('ppi')]], 
                    how='left', left_index=True, right_index=True,
                    suffixes=('', '_' + label))

pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(cofactor_partners),
                           label='cofactors'
)
pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(cofactor_animal_db),
                           label='cofactors_animal_db'
)
pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(coactivator_partners) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='coactivators_HEK'
)
pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(coactivator_partners) &
                                 ~y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='coactivators_not_HEK'
)
pairs = add_restricted_ppi_columns(pairs, 
                           rows=~y2h['db_gene_symbol'].isin(cofactor_partners),
                           label='not_cofactors'
)
pairs = add_restricted_ppi_columns(pairs, 
                           rows=y2h['db_gene_symbol'].isin(coactivator_partners),
                           label='coactivators'
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
                           rows=y2h['db_gene_symbol'].isin(tf_gene_symbols),
                           label='tfs'
)
pairs = add_restricted_ppi_columns(pairs, 
                           rows=(y2h['db_gene_symbol'].isin(signaling_partners) &
                                 y2h['db_gene_symbol'].isin(hek_expressed_genes)),
                           label='signaling_HEK'
)


# In[77]:


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
    fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.25))

    def bin_delta_ppi(delta_ppi):
        if pd.isnull(delta_ppi):
            return np.nan
        if delta_ppi < 0:
            return 'loss'
        elif delta_ppi > 0:
            return 'gain'
        elif delta_ppi == 0:
            return 'equal'
        else:
            raise ValueError(delta_ppi)


    df[x + '_binned'] = df[x].apply(bin_delta_ppi)
    sns.stripplot(data=df,
                  x=x + '_binned',
                  y=y,
                  order=['loss', 'equal', 'gain'],
                  alpha=0.75,
                  color=color,
                  linewidth=1,
                  edgecolor="black",
                  ax=ax)
    if False:
        sns.pointplot(data=df,
                    x=x + '_binned',
                    y=y,
                    order=['loss', 'equal', 'gain'],
                    alpha=0.5,
                    color='black',
                    ax=ax)
    if True:
        sns.boxplot(data=df,
                    x=x + '_binned',
                    y=y,
                    order=['loss', 'equal', 'gain'],
                    fliersize=0,
                    color=color,
                    ax=ax)
        mimic_r_boxplot(ax)
    else:
        sns.violinplot(data=df,
                    x=x + '_binned',
                    y=y,
                    order=['loss', 'equal', 'gain'],
                    color='lightgrey',
                    ax=ax)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    with_data = (df[x].notnull() & df[y].notnull())
    n_pair = with_data.sum()
    n_iso = len(set(df.loc[with_data, ['clone_acc_a', 'clone_acc_b']].values.flatten()))
    n_gene = df.loc[with_data, 'tf_gene_symbol'].nunique()
    scc, p_scc = stats.spearmanr(df.loc[df[x].notnull() & df[y].notnull(), x].values,
                    df.loc[df[x].notnull() & df[y].notnull(), y].values)
    ax.text(s=f'{n_pair:d} pairs\n{n_iso:d} isoforms\n{n_gene:d} genes\nSpearman r = {scc:.2f}\np = {p_scc:.2f}',
            x=1.03,
            y=0.95,
            ha='left',
            va='top',
            transform=ax.transAxes)
    #ax.set_ylim(-4, 4) # NOTE cuts outlier TODO add broken axis
    ax.axhline(y=0, color='black', linestyle="dashed", linewidth=1)
    for pos in ['top', 'bottom', 'right']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    fig.savefig(f'../figures/{x}-vs-{y}_scatter.pdf',
                bbox_inches='tight')


# In[78]:


bar_activation_vs_ppi(
    pairs=pairs.loc[(pairs['m1h_min'] < -1) | (pairs['m1h_max'] > 1), :],
    x='ppi_delta_n_coactivators_HEK',
    y='activation_fold_change',
    x_label='Difference in number of coactivator PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='log-2 fold-change in activation',
    color=sns.color_palette("Set2")[0])


# In[84]:


# TEST what is the outlier on the left
(pairs.loc[((pairs['m1h_min'] < -1) | 
           (pairs['m1h_max'] > 1))
           & (pairs['ppi_delta_n_coactivators_HEK'].notnull()), 
              :]
           .sort_values('activation_fold_change',
                        ascending=False))[["tf_gene_symbol", "clone_acc_a", "clone_acc_b",
                                           "ppi_delta_n_coactivators_HEK", "activation_fold_change"]]


# In[80]:


bar_activation_vs_ppi(
    pairs=pairs.loc[(pairs['m1h_min'] < -1) | (pairs['m1h_max'] > 1), :],
    x='ppi_delta_n_corepressors_HEK',
    y='activation_fold_change',
    x_label='Difference in number of corepressor PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='log-2 fold-change in activation',
    color=sns.color_palette("Set2")[3])


# In[81]:


bar_activation_vs_ppi(
    pairs=pairs.loc[(pairs['m1h_min'] < -1) | (pairs['m1h_max'] > 1), :],
    x='ppi_delta_n_signaling_HEK',
    y='activation_fold_change',
    x_label='Difference in number of signaling PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='log-2 fold-change in activation',
    color=sns.color_palette("Set2")[5])


# In[89]:


# TEST what is the outlier on the left
(pairs.loc[((pairs['m1h_min'] < -1) | 
           (pairs['m1h_max'] > 1))
           & (pairs['ppi_delta_n_signaling_HEK'].notnull()), 
              :]
           .sort_values('activation_fold_change',
                        ascending=True))[["tf_gene_symbol", "clone_acc_a", "clone_acc_b",
                                           "ppi_delta_n_coactivators_HEK", "activation_fold_change"]].head(30)


# In[90]:


pairs[pairs["tf_gene_symbol"] == "PBX1"][["tf_gene_symbol", "clone_acc_a", "clone_acc_b",
                                           "ppi_delta_n_signaling_HEK", "activation_fold_change"]]


# In[92]:


cats[cats["partner"] == "PIN1"]

