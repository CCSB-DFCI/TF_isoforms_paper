
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
    load_y1h_pdi_data,
    load_ref_vs_alt_isoforms_table,
    load_annotated_TFiso1_collection)
from data_loading.isoform_pairwise_metrics import _add_PPI_columns


# In[2]:


pd.set_option('display.max_columns', 100)


# In[3]:


y2h = load_y2h_isoform_data()
m1h = load_m1h_activation_data()
cats = load_ppi_partner_categories()


# In[4]:


# TODO: move to data_loading
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


# In[9]:


pairs = load_ref_vs_alt_isoforms_table()
tfs = load_annotated_TFiso1_collection()


def add_restricted_ppi_columns(pairs, rows, label):
    pairs_cf = pairs[['clone_acc_ref', 'clone_acc_alt']].copy()
    _add_PPI_columns(df=pairs_cf, y2h=y2h.loc[rows, :])
    return pd.merge(pairs, 
                    pairs_cf,
                    how='left',
                    on=['clone_acc_ref', 'clone_acc_alt'],
                    suffixes=('', '_' + label))


cofactor_partners = set(cats.loc[cats['category'] == 'cofactor', 'gene_symbol_partner'].unique())
coactivator_partners = set(cats.loc[cats['cofactor_type'] == 'coactivator', 'gene_symbol_partner'].unique())
corepressor_partners = set(cats.loc[cats['cofactor_type'] == 'corepressor', 'gene_symbol_partner'].unique())
cofactor_animal_db = set(cof['Symbol'].unique())
tf_gene_symbols = set(load_human_tf_db()['HGNC symbol'].values)
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


# In[10]:


pairs['activation_abs_fold_change'] = pairs['activation_fold_change_log2'].abs()


# In[11]:


def scatter_with_correlation(x, y, pairs=pairs, x_label=None, y_label=None):
    """
    TODO:
        - calculate p-value properly
            - this requires permuting in some smart way
            - one question is whether the genes are the number of independent data points or the isoforms are
            - I think the answer is the isoforms are
    
    """
    if x_label is None:
        x_label = x
    if y_label is None:
        y_label = y
    fig, ax = plt.subplots(1, 1)
    pairs.plot.scatter(x=x, y=y, ax=ax, alpha=0.3)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    with_data = (pairs[x].notnull() & pairs[y].notnull())
    n_pair = with_data.sum()
    n_iso = len(set(pairs.loc[with_data, ['clone_acc_ref', 'clone_acc_alt']].values.flatten()))
    n_gene = pairs.loc[with_data, 'gene_symbol'].nunique()
    scc, p_scc = stats.spearmanr(pairs.loc[pairs[x].notnull() & pairs[y].notnull(), x].values,
                    pairs.loc[pairs[x].notnull() & pairs[y].notnull(), y].values)
    ax.text(s=f'{n_pair:d} pairs\n{n_iso:d} isoforms\n{n_gene:d} genes\nSpearman r = {scc:.2f}\np = {p_scc:.2f}',
            x=0.95,
            y=0.95,
            ha='right',
            va='top',
            transform=ax.transAxes)
    fig.savefig(f'../figures/{x}-vs-{y}_scatter.pdf',
                bbox_inches='tight')

scatter_with_correlation(
    x='n_PPI_diff',
    y='activation_abs_fold_change',
    x_label='Number of different PPIs',
    y_label='Absolute log-2 fold-change in activation')


# In[12]:


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
    fig, ax = plt.subplots(1, 1, figsize=(3.25, 2.25))

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
    n_iso = len(set(df.loc[with_data, ['clone_acc_ref', 'clone_acc_alt']].values.flatten()))
    n_gene = df.loc[with_data, 'gene_symbol'].nunique()
    scc, p_scc = stats.spearmanr(df.loc[df[x].notnull() & df[y].notnull(), x].values,
                    df.loc[df[x].notnull() & df[y].notnull(), y].values)
    ax.text(s=f'{n_pair:d} pairs\n{n_iso:d} isoforms\n{n_gene:d} genes\nSpearman r = {scc:.2f}\np = {p_scc:.2f}',
            x=1.03,
            y=0.95,
            ha='left',
            va='top',
            transform=ax.transAxes)
    ax.set_ylim(-4, 4) # NOTE cuts outlier TODO add broken axis
    ax.axhline(y=0, color='black', linestyle="dashed", linewidth=1)
    for pos in ['top', 'bottom', 'right']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    fig.savefig(f'../figures/{x}-vs-{y}_scatter.pdf',
                bbox_inches='tight')

bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n',
    y='activation_fold_change_log2',
    x_label='Net difference in number of PPIs',
    y_label='log-2 fold-change in activation')


# In[13]:


bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_coactivators',
    y='activation_fold_change_log2',
    x_label='Difference in number of coactivator PPIs',
    y_label='log-2 fold-change in activation')


# In[14]:


bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_corepressors',
    y='activation_fold_change_log2',
    x_label='Difference in number of corepressor PPIs',
    y_label='log-2 fold-change in activation')


# In[15]:


bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_cofactors_animal_db',
    y='activation_fold_change_log2',
    x_label='Difference in number of cofactor PPIs',
    y_label='log-2 fold-change in activation')


# In[16]:


bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_coactivators_HEK',
    y='activation_fold_change_log2',
    x_label='Difference in number of coactivator PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='log-2 fold-change in activation',
    color=sns.color_palette("Set2")[0])


# In[17]:


# TEST what is the outlier on the left
(pairs.loc[(pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True)
           & (pairs['PPI_delta_n_coactivators_HEK'].notnull()), 
              :]
           .sort_values('activation_fold_change_log2',
                        ascending=False)).head()


# In[18]:


bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_corepressors_HEK',
    y='activation_fold_change_log2',
    x_label='Difference in number of corepressor PPIs\nrestricted to those expressed in HEK293 cells',
    y_label='log-2 fold-change in activation',
    color=sns.color_palette("Set2")[3])


# In[19]:


# this is not correct
# what I want is the pairs where all the lost cofactors are not expressed
bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_coactivators_not_HEK',
    y='activation_fold_change_log2',
    x_label='Difference in number of coactivator PPIs\nrestricted to those NOT expressed in HEK283 cells',
    y_label='log-2 fold-change in activation')


# In[20]:


bar_activation_vs_ppi(
    pairs=pairs.loc[pairs['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] == True, :],
    x='PPI_delta_n_tfs',
    y='activation_fold_change_log2',
    x_label='Difference in number of TF-TF PPIs',
    y_label='log-2 fold-change in activation')


# In[21]:


scatter_with_correlation(
    x='PPI_delta_n',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs',
    y_label='log-2 fold-change in activation')


# In[22]:


scatter_with_correlation(
    x='n_PPI_diff_cofactors',
    y='activation_abs_fold_change',
    x_label='Number of different PPIs\nrestricting to cofactors',
    y_label='Absolute log-2 fold-change in activation')


# In[23]:


scatter_with_correlation(
    x='PPI_delta_n_cofactors',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs\nrestricting to cofactors',
    y_label='log-2 fold-change in activation')


# In[26]:


scatter_with_correlation(
    x='n_PPI_diff_not_cofactors',
    y='activation_abs_fold_change',
    x_label='Number of different PPIs\nrestricting to everything other than cofactors',
    y_label='Absolute log-2 fold-change in activation')


# In[27]:


scatter_with_correlation(
    x='PPI_delta_n_not_cofactors',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs\nrestricting to everything other than cofactors',
    y_label='log-2 fold-change in activation')


# In[30]:


pairs.loc[(pairs['n_PPI_successfully_tested_in_ref_and_alt_coactivators'] >= 1) & (pairs['n_PPI_diff_cofactors'] == 0), :]


# In[34]:


# restrict to no change with cofactors
# BUG: file name could overwrite?
scatter_with_correlation(
    pairs=pairs.loc[(pairs['n_PPI_successfully_tested_in_ref_and_alt_coactivators'] >= 1) & (pairs['n_shared_PPI_cofactors'] == pairs['n_PPI_successfully_tested_in_ref_and_alt_coactivators']), :],
    x='PPI_delta_n_not_cofactors',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs\nrestricting to everything other than cofactors',
    y_label='log-2 fold-change in activation')


# In[38]:


scatter_with_correlation(
    pairs=pairs.loc[(pairs['n_PPI_diff_coactivators'] == 0), :],
    x='PPI_delta_n',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs\nrestricting to everything other than cofactors',
    y_label='log-2 fold-change in activation')


# In[40]:


scatter_with_correlation(
    x='n_PPI_diff_coactivators',
    y='activation_fold_change_log2',
    x_label='Number of different PPIs\nrestricting to coactivators',
    y_label='Absolute log-2 fold-change in activation')


# In[41]:


scatter_with_correlation(
    x='PPI_delta_n_coactivators',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs\nrestricting to coactivators',
    y_label='log-2 fold-change in activation')


# In[43]:


scatter_with_correlation(
    x='PPI_delta_n_tfs',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs\nrestricting to TFs',
    y_label='log-2 fold-change in activation')


# In[45]:


scatter_with_correlation(
    x='n_PPI_diff_corepressors',
    y='activation_fold_change_log2',
    x_label='Number of different PPIs\nrestricting to corepressors',
    y_label='Absolute log-2 fold-change in activation')


# In[46]:


scatter_with_correlation(
    x='PPI_delta_n_corepressors',
    y='activation_fold_change_log2',
    x_label='Difference in number of PPIs\nrestricting to corepressors',
    y_label='log-2 fold-change in activation')


# In[47]:


# co-repressors minus coactivators
pairs['PPI_delta_n_coactivators_minus_corepressors'] = pairs['PPI_delta_n_coactivators'].fillna(0) - pairs['PPI_delta_n_corepressors'].fillna(0)
pairs.loc[pairs['PPI_delta_n_coactivators'].isnull() & pairs['PPI_delta_n_corepressors'].isnull(),
          'PPI_delta_n_coactivators_minus_corepressors'] = np.nan
scatter_with_correlation(
    x='PPI_delta_n_coactivators_minus_corepressors',
    y='activation_fold_change_log2',
    x_label='number of coactivator PPIs\n- number of corepressor PPIs',
    y_label='log-2 fold-change in activation')


# In[51]:


# look at the biggest difference examples
(pairs.loc[pairs['n_PPI_diff'] > 0, 
['gene_symbol', 'clone_acc_ref', 'clone_acc_alt', 'activation_abs_fold_change', 'n_PPI_successfully_tested_in_ref_and_alt', 'n_PPI_diff', 'n_PPI_diff_cofactors']]
.sort_values('activation_abs_fold_change', ascending=False).head(25)
)


# In[54]:


# look for repression differences
(pairs.loc[(pairs['n_PPI_diff'] > 0) & pairs['gene_symbol'].isin(m1h.loc[(m1h['M1H_rep1'] < -1), 'gene_symbol'].unique()), 
['gene_symbol', 'clone_acc_ref', 'clone_acc_alt', 'activation_abs_fold_change', 'n_PPI_successfully_tested_in_ref_and_alt', 'n_PPI_diff', 'n_PPI_diff_cofactors']]
.sort_values('activation_abs_fold_change', ascending=False).head(25)
)


# In[57]:


# do category plot on some vs no difference
pairs['PPI_any_diff'] = (pairs['n_PPI_diff'] != 0)
pairs.loc[pairs['n_PPI_diff'].isnull(), 'PPI_any_diff'] = np.nan
x = 'PPI_any_diff'
y = 'activation_fold_change_log2'
sns.swarmplot(data=pairs, x='PPI_any_diff', y='activation_fold_change_log2')
sns.pointplot(data=pairs, x='PPI_any_diff', y='activation_fold_change_log2', estimator=np.median)
stats.mannwhitneyu(pairs.loc[pairs[x] == True, y], pairs.loc[pairs[x] == False, y])


# In[58]:


pairs['n_PPI_diff'].value_counts().sum()


# In[59]:


pairs['PPI_any_diff'].value_counts()


# In[60]:


pairs['PPI_any_diff_cofactors'] = pairs['n_PPI_diff_cofactors'] != 0
pairs.loc[pairs['n_PPI_diff_cofactors'].isnull(), 'PPI_any_diff_cofactors'] = np.nan
x = 'PPI_any_diff_cofactors'
y = 'activation_fold_change_log2'
sns.swarmplot(data=pairs, x='PPI_any_diff_cofactors', y='activation_fold_change_log2')
sns.pointplot(data=pairs, x='PPI_any_diff_cofactors', y='activation_fold_change_log2', estimator=np.median)
stats.mannwhitneyu(pairs.loc[pairs[x] == True, y], pairs.loc[pairs[x] == False, y])


# In[61]:


x = 'PPI_any_diff_not_cofactors'
y = 'activation_fold_change_log2'
pairs['PPI_any_diff_not_cofactors'] = pairs['n_PPI_diff_not_cofactors'] != 0
pairs.loc[pairs['n_PPI_diff_not_cofactors'].isnull(), 'PPI_any_diff_not_cofactors'] = np.nan
sns.swarmplot(data=pairs, x='PPI_any_diff_not_cofactors', y='activation_fold_change_log2')
sns.pointplot(data=pairs, x='PPI_any_diff_not_cofactors', y='activation_fold_change_log2', estimator=np.median)
stats.mannwhitneyu(pairs.loc[pairs[x] == True, y], pairs.loc[pairs[x] == False, y])


# In[62]:


y2h.head()


# In[63]:


# what is the best metric for PPI partners whose loss is most associciated with change in
# activation?

def ppi_tf_gene(data, gene_name):
    tf = data.loc[(data['ad_gene_symbol'] == gene_name),
                  ['ad_clone_acc', 'db_gene_symbol', 'Y2H_result']].copy()
    tf = tf.pivot(index='ad_clone_acc',
                  columns='db_gene_symbol',
                  values='Y2H_result')
    return tf

# for each tf, for each partner, if differentially bound, absolute change in activation
data = []
for tf_gene in y2h['ad_gene_symbol'].unique():
    tf_ppi = ppi_tf_gene(y2h, tf_gene)
    tf_ppi = tf_ppi.loc[tf_ppi.any(axis=1), :]
    if tf_ppi.shape[0] < 2:
        continue
    tf_ppi = tf_ppi.loc[:, (tf_ppi == True).any(axis=0) & (tf_ppi == False).any(axis=0)]
    for partner in tf_ppi.columns:
        binds = tf_ppi.index[tf_ppi[partner] == True].values
        not_binds = tf_ppi.index[tf_ppi[partner] == False].values
        actv_diff = (m1h.loc[m1h['clone_acc'].isin(binds), ['M1H_rep1', 'M1H_rep2', 'M1H_rep3']].mean().mean()
                     - m1h.loc[m1h['clone_acc'].isin(not_binds), ['M1H_rep1', 'M1H_rep2', 'M1H_rep3']].mean().mean())
        data.append([tf_gene, partner, actv_diff])
df = pd.DataFrame(data=data, columns=['gene_symbol', 'gene_symbol_partner', 'activation_diff']).dropna()
df['num_points'] = df['gene_symbol_partner'].map(df['gene_symbol_partner'].value_counts())
df = df.sort_values(['num_points', 'gene_symbol_partner', 'activation_diff'], ascending=False)


# In[64]:


df.head(12)


# In[65]:


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(16, 4)
sns.stripplot(data=df.loc[df['num_points'] > 2, :],
              x='gene_symbol_partner',
              y='activation_diff',
              ax=ax)
ax.axhline(0, linestyle='--', color='grey')
ax.set_xlabel('')
ax.set_ylabel('log2 activation fold change')
ax.tick_params(axis='x', rotation=90, length=0)
plt.savefig('../figures/activation_change_per_ppi_partner.pdf',
            bbox_inches='tight')

