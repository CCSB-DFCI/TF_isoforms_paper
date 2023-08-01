
# coding: utf-8

# Make a table in relation to this dominant negative idea.
# 
# Variables:
# 
# - DBD intact?
# - Dimerizing PPI intact?
# - M1H - loss of activation
# - Y1H - loss of PDIs

# In[1]:


import numpy as np
import pandas as pd

from data_loading import (load_annotated_6k_collection,
                          load_y2h_isoform_data,
                          load_y1h_pdi_data,
                          load_m1h_activation_data,
                          load_valid_isoform_clones)
from isoform_pairwise_metrics import pairs_of_isoforms_comparison_table


# In[2]:


# TODO: move to data loading
DIMERIZING_TF_FAMILIES = {'bHLH',
              'bZIP',
              'Nuclear receptor',
              'E2F',
              'CENPB',
              'MADS box',
              'Grainyhead',
              'SAND',
              'Rel',
              'EBF1',
              'STAT',
              'IRF',
              'RFX',
              'HSF',
              'p53',
              'AP-2',
              'GCM',
              'BED ZF',
              'MADF',
              'ARID/BRIGHT',
              'Myb/SANT',
              'SMAD'}


# In[3]:


tfs = load_annotated_6k_collection()


# In[4]:


m1h = load_m1h_activation_data()
y1h = load_y1h_pdi_data()
y2h = load_y2h_isoform_data()


# In[5]:


from isoform_pairwise_metrics import _pairs_comparison_table


def pairs_of_ref_vs_alt_isoforms_comparison_table(tfs, y2h=None, y1h=None, m1h=None):
    iso_pairs = []
    for tf in tfs.values():
        ref = tf.cloned_reference_isoform
        for alt in tf.cloned_isoforms:
            if alt.name == ref.name:
                continue
            iso_pairs.append((tf.name,
                              tf.ensembl_gene_id,
                              tf.tf_family,
                              tf.tf_family in DIMERIZING_TF_FAMILIES,
                              ref.clone_acc,
                              alt.clone_acc,
                              '|'.join(ref.ensembl_transcript_ids) if ref.ensembl_transcript_ids is not None else np.nan,
                              '|'.join(alt.ensembl_transcript_ids) if alt.ensembl_transcript_ids is not None else np.nan,
                              ref.is_novel_isoform(),
                              alt.is_novel_isoform(),
                              tf.cloned_MANE_select_isoform,
                              len(ref.aa_seq),
                              len(alt.aa_seq),
                              len(ref.exons),
                              len(alt.exons),
                              tf.splicing_categories(ref.name, alt.name)["alternative N-terminal"],
                              tf.splicing_categories(ref.name, alt.name)["alternative C-terminal"],
                              tf.splicing_categories(ref.name, alt.name)["alternative internal exon"],
                              tf.splicing_categories(ref.name, alt.name)["alternative 5' splice site"],
                              tf.splicing_categories(ref.name, alt.name)["alternative 3' splice site"],
                              tf.splicing_categories(ref.name, alt.name)["exon skipping"],
                              tf.splicing_categories(ref.name, alt.name)["mutually exclusive exons"],
                              tf.splicing_categories(ref.name, alt.name)["intron retention"],

                              ))
    iso_pairs = pd.DataFrame(
        data=iso_pairs,
        columns=["gene_symbol",
                 "Ensembl_gene_ID",
                 "family",
                 "is_dimerizing_TF_family",
                 "clone_acc_ref",
                 "clone_acc_alt",
                 "Ensembl_transcript_IDs_ref",
                 "Ensembl_transcript_IDs_alt",
                 "is_ref_novel_isoform",
                 "is_alt_novel_isoform",
                 "is_MANE_select_isoform_cloned",
                 "n_aa_ref",
                 "n_aa_alt",
                 "n_exons_ref",
                 "n_exons_alt",
                 "is_alternative_N_terminal",
                 "is_alternative_C_terminal",
                 "is_alternative_internal_exon",
                 "is_alternative_5_prime_donor",
                 "is_alternative_3_prime_acceptor",
                 "is_exon_skipping",
                 "is_mutually_exclusive_exons",
                 "is_intron_retention",
                 ]
    )
    return iso_pairs



df = pairs_of_ref_vs_alt_isoforms_comparison_table(tfs, y2h=y2h, y1h=y1h, m1h=m1h)


# In[6]:


# DBD intact
from data_loading import load_dbd_accessions

def load_dbd_affected():
    df = pd.concat([g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()])
    df['is_DBD'] = df['accession'].isin(load_dbd_accessions())
    df_new = (df.loc[df['is_DBD'], :]
        .groupby(['gene', 'ref_iso', 'alt_iso'])
        [['deletion', 'frameshift']].sum()
        .sum(axis=1) / df.loc[df['is_DBD'], :]
        .groupby(['gene', 'ref_iso', 'alt_iso'])
        ['length'].sum()).to_frame(name='dbd_fraction')
    df_new['dbd_insertion_n_aa'] = (df.loc[df['is_DBD'], :]
                                  .groupby(['gene', 'ref_iso', 'alt_iso'])
                                  ['insertion']
                                  .sum())
    df = df_new.reset_index()
    df['dbd_pct_lost'] = df['dbd_fraction'] * 100.
    df = df.drop(columns=['dbd_fraction'])
    return df


dbd = load_dbd_affected()
dbd['clone_acc_ref'] = dbd['ref_iso'].map({iso.name: iso.clone_acc for tf in tfs.values() for iso in tf.cloned_isoforms})
dbd['clone_acc_alt'] = dbd['alt_iso'].map({iso.name: iso.clone_acc for tf in tfs.values() for iso in tf.cloned_isoforms})
dbd = dbd.drop(columns=['gene', 'ref_iso', 'alt_iso'])
df = pd.merge(df, dbd, how='left', on=['clone_acc_ref', 'clone_acc_alt'])
df['dbd_affected'] = df['dbd_pct_lost'] > 0


# In[7]:


from data_loading import load_seq_comparison_data

aa_ident = load_seq_comparison_data()
df["aa_seq_pct_id"] = df.apply(
    lambda x: "_".join(sorted([x["clone_acc_ref"], x["clone_acc_alt"]])), axis=1
).map(aa_ident)
if df['aa_seq_pct_id'].isnull().any():
    raise UserWarning('Unexpected missing sequence similarity values')


# In[8]:


# now do the assays

# y2h n_positive_ref n_positive_ref_filtered n_shared_ref_alt
# y2h n successfully tested in both
# M1H at least one isoform of gene has |activation| >= 2 fold

y2h_complete = load_y2h_isoform_data(require_at_least_one_ppi_per_isoform=False)
n_ppi = y2h_complete.loc[(y2h_complete['Y2H_result'] == True), :].groupby('ad_clone_acc').size()
df['n_positive_PPI_ref'] = df['clone_acc_ref'].map(n_ppi)
df['n_positive_PPI_alt'] = df['clone_acc_alt'].map(n_ppi)
# BUG MISSING 0's here!
df.loc[df['n_positive_PPI_ref'].isnull() &
       df['clone_acc_ref'].isin(y2h_complete.loc[(y2h_complete['Y2H_result'] == False), 
                              'ad_clone_acc'].unique()),
       'n_positive_PPI_ref'] = 0
df.loc[df['n_positive_PPI_alt'].isnull() &
       df['clone_acc_alt'].isin(y2h_complete.loc[(y2h_complete['Y2H_result'] == False), 
                              'ad_clone_acc'].unique()),
       'n_positive_PPI_alt'] = 0


# In[9]:


from isoform_pairwise_metrics import (
 number_tested_partners,
  number_shared_partners,
  jaccard_index)

def ppi_metric(row, data, function, suffixes=('_a', '_b')):
    ad_a = row["clone_acc" + suffixes[0]]
    ad_b = row["clone_acc" + suffixes[1]]
    pair = data.loc[data["ad_clone_acc"].isin([ad_a, ad_b]), :].pivot(
        values="Y2H_result", index="db_gene_symbol", columns="ad_clone_acc"
    )
    if ad_a not in pair.columns or ad_b not in pair.columns:
        return np.nan
    # remove any partner with AA / NC / NS / NaN in either
    pair = pair.loc[pair.notnull().all(axis=1), :].astype(int).astype(bool)
    # remove partners that tested negative in both
    pair = pair.loc[pair.any(axis=1), :]
    if pair.shape[0] > 0:
        return function(
            set(pair.index[pair[ad_a]].values), set(pair.index[pair[ad_b]].values)
        )
    else:
        return np.nan

df['n_PPI_successfully_tested_in_ref_and_alt'] = df.apply(
            ppi_metric, data=y2h_complete, suffixes=("_ref", "_alt"), function=number_tested_partners, axis=1
        )
df['n_positive_PPI_ref_filtered'] = df.apply(
            ppi_metric, data=y2h_complete, suffixes=("_ref", "_alt"), function=lambda a, b: len(a), axis=1
        )
df['n_positive_PPI_alt_filtered'] = df.apply(
            ppi_metric, data=y2h_complete, suffixes=("_ref", "_alt"), function=lambda a, b: len(b), axis=1
        )
df['n_shared_PPI'] = df.apply(
            ppi_metric, data=y2h_complete, suffixes=("_ref", "_alt"), function=number_shared_partners, axis=1
        )
df['PPI_jaccard'] = df.apply(
            ppi_metric, data=y2h_complete, suffixes=("_ref", "_alt"), function=jaccard_index, axis=1
        )


# In[10]:


from data_loading import load_tf_families
from data_loading import load_human_tf_db
from data_loading import load_ppi_partner_categories


ppi_partner_cats = load_ppi_partner_categories()
tfdb = load_human_tf_db()
fam = load_tf_families()
y2h['ad_tf_family'] = y2h['ad_gene_symbol'].map(fam)
y2h['db_tf_family'] = y2h['db_gene_symbol'].map(fam)
y2h['is_dimerizing_ppi'] = (y2h['ad_tf_family'].isin(DIMERIZING_TF_FAMILIES) &
                        (y2h['ad_tf_family'] == y2h['db_tf_family']))
y2h['is_tf_tf_ppi'] = y2h['db_gene_symbol'].isin(tfdb['HGNC symbol'].unique())

# of reference dimer PPI, are all lost, some lost, none lost
def ppi_pertubation(row, ppi):
    ref_clone_acc = row['clone_acc_ref']
    alt_clone_acc = row['clone_acc_alt']
    if ref_clone_acc not in ppi['ad_clone_acc'].unique() or alt_clone_acc not in ppi['ad_clone_acc'].unique():
        return np.nan
    df = (ppi.loc[ppi['ad_clone_acc'].isin([ref_clone_acc, alt_clone_acc]),
                  ['ad_clone_acc', 'db_gene_symbol', 'Y2H_result']]
            .pivot(values='Y2H_result', index='db_gene_symbol', columns='ad_clone_acc')
            .dropna())
    df = df.loc[df.any(axis=1), :]
    if df.shape[0] == 0:
        return np.nan
    if df.all().all():
        return 'retains all'
    elif not df[alt_clone_acc].any():
        return 'loses all'
    elif df[alt_clone_acc].sum() > df[ref_clone_acc].sum():
        return 'gains some'
    else:
        return 'loses some'

df['dimer_ppi'] = df.apply(ppi_pertubation, 
                           ppi=y2h.loc[y2h['is_dimerizing_ppi'], :],
                           axis=1)
df['other_than_dimer_ppi'] = df.apply(ppi_pertubation, 
                           ppi=y2h.loc[~y2h['is_dimerizing_ppi'], :],
                           axis=1)
df['tf_tf_ppi'] = df.apply(ppi_pertubation, 
                           ppi=y2h.loc[y2h['is_tf_tf_ppi'], :],
                           axis=1)

                           
df['tf_cofactor_ppi'] = df.apply(ppi_pertubation, 
                           ppi=y2h.loc[y2h['db_gene_symbol'].isin(ppi_partner_cats.loc[ppi_partner_cats['category'] == 'cofactor', 'partner'].unique()) &
                                       ~y2h['is_tf_tf_ppi'], :],
                           axis=1)
df['tf_signalling_ppi'] = df.apply(ppi_pertubation, 
                           ppi=y2h.loc[y2h['db_gene_symbol'].isin(ppi_partner_cats.loc[ppi_partner_cats['category'].isin(['cell cycle', 'protein traffiking', 'protein turnover', 'signaling', 'cell-cell signaling']), 'partner'].unique()) &
                                       ~y2h['is_tf_tf_ppi'] &
                                       ~y2h['db_gene_symbol'].isin(ppi_partner_cats.loc[ppi_partner_cats['category'] == 'cofactor', 'partner'].unique()), :],
                           axis=1)


# In[11]:


# TODO: move code to data loading
y1h = y1h.drop_duplicates(subset=['unique_acc'])


# In[16]:


y1h.head()


# In[18]:


y1h = y1h.set_index('unique_acc')


# In[20]:


def pdi_metric(row, data, function, suffixes=('_a', '_b')):
    clone_acc_a = row["clone_acc" + suffixes[0]]
    clone_acc_b = row["clone_acc" + suffixes[1]]
    df = data.loc[
        (data.index == clone_acc_a)
        | (data.index == clone_acc_b),
        data.columns[1:],
    ].copy()
    if df.shape[0] < 2:
        return np.nan
    df = df.loc[[clone_acc_a, clone_acc_b], df.any(axis=0)]
    if df.shape[1] == 0:
        return np.nan
    a = set(df.columns[df.iloc[0]])
    b = set(df.columns[df.iloc[1]])
    return function(a, b)

n_pdi = y1h.drop(columns=['tf']).sum(axis=1)
df['n_positive_PDI_ref'] = df['clone_acc_ref'].map(n_pdi)
df['n_positive_PDI_alt'] = df['clone_acc_alt'].map(n_pdi)
df['n_PDI_successfully_tested_in_ref_and_alt'] = df.apply(
            pdi_metric, data=y1h, suffixes=("_ref", "_alt"), function=number_tested_partners, axis=1
        )

df['n_positive_PDI_ref_filtered'] = df.apply(
            pdi_metric, data=y1h, suffixes=("_ref", "_alt"), function=lambda a, b: len(a), axis=1
        )
df['n_positive_PDI_alt_filtered'] = df.apply(
            pdi_metric, data=y1h, suffixes=("_ref", "_alt"), function=lambda a, b: len(b), axis=1
        )


df['n_shared_PDI'] = df.apply(
            pdi_metric, data=y1h, suffixes=("_ref", "_alt"), function=number_shared_partners, axis=1
        )
df['PDI_jaccard'] = df.apply(
            pdi_metric, data=y1h, suffixes=("_ref", "_alt"), function=jaccard_index, axis=1
        )


# In[22]:


m1h['mean'] = m1h[['M1H_rep1', 'M1H_rep2', 'M1H_rep3']].mean(axis=1)
m1h['abs_mean'] = m1h['mean'].abs()
df['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] = df['gene_symbol'].map(m1h.groupby('gene')['abs_mean'].max() >= 1)
df['activation_ref'] = df['clone_acc_ref'].map(m1h.set_index('clone_acc')['mean'])
df['activation_alt'] = df['clone_acc_alt'].map(m1h.set_index('clone_acc')['mean'])
df['activation_fold_change_log2'] = (df['activation_alt'] - df['activation_ref'])


# In[23]:


df.to_csv('../output/TF-iso_ref-vs-alt.tsv', sep='\t', index=False)


# In[4]:


df = pairs_of_isoforms_comparison_table(isoforms=load_valid_isoform_clones(),
                                        m1h=load_m1h_activation_data(),
                                        y1h=load_y1h_pdi_data(),
                                        y2h=load_y2h_isoform_data())


# In[5]:


ref_clones = {tf.cloned_reference_isoform.clone_acc for tf in tfs.values()}
df = df.loc[df['clone_acc_a'].isin(ref_clones) | df['clone_acc_b'].isin(ref_clones)]


# In[6]:


df['clone_acc_ref'] = (df['clone_acc_a'] * df['clone_acc_a'].isin(ref_clones) +
                       df['clone_acc_b'] * df['clone_acc_b'].isin(ref_clones))
df['clone_acc_alt'] = (df['clone_acc_b'] * df['clone_acc_a'].isin(ref_clones) +
                       df['clone_acc_a'] * df['clone_acc_b'].isin(ref_clones))


# In[10]:


from data_loading import load_tf_families

y2h = load_y2h_isoform_data()
fam = load_tf_families()
# NOTE: this is not complete. See issue #69
y2h['ad_tf_family'] = y2h['ad_gene_symbol'].map(fam)
y2h['db_tf_family'] = y2h['db_gene_symbol'].map(fam)
y2h['is_dimerizing_ppi'] = (y2h['ad_tf_family'].isin(DIMERIZING_TF_FAMILIES) &
                        (y2h['ad_tf_family'] == y2h['db_tf_family']))


# In[11]:


# of reference dimer PPI, are all lost, some lost, none lost
def ppi_pertubation(row, ppi):
    ref_clone_acc = row['clone_acc_ref']
    alt_clone_acc = row['clone_acc_alt']
    if ref_clone_acc not in ppi['ad_clone_acc'].unique() or alt_clone_acc not in ppi['ad_clone_acc'].unique():
        return np.nan
    df = (ppi.loc[ppi['ad_clone_acc'].isin([ref_clone_acc, alt_clone_acc]),
                  ['ad_clone_acc', 'db_gene_symbol', 'Y2H_result']]
            .pivot(values='Y2H_result', index='db_gene_symbol', columns='ad_clone_acc')
            .dropna())
    df = df.loc[df.any(axis=1), :]
    if df.shape[0] == 0:
        return np.nan
    if df.all().all():
        return 'retains all'
    elif not df[alt_clone_acc].any():
        return 'loses all'
    elif df[alt_clone_acc].sum() > df[ref_clone_acc].sum():
        return 'gains some'
    else:
        return 'loses some'

df['dimer_ppi'] = df.apply(ppi_pertubation, 
                           ppi=y2h.loc[y2h['is_dimerizing_ppi'], :],
                           axis=1)
df['other_ppi'] = df.apply(ppi_pertubation, 
                           ppi=y2h.loc[~y2h['is_dimerizing_ppi'], :],
                           axis=1)


# In[13]:


df['dimer_ppi'].value_counts()


# In[117]:


df.loc[(df['dimer_ppi'] == 'loses some'), :]


# In[149]:


df.loc[(df['dimer_ppi'] == 'loses all'), :]


# In[120]:


df.loc[(df['dimer_ppi'] == 'gains some'), :]


# In[158]:


# retains dimer PPI and loses PDI or M1H or DBD
df.loc[(df['dimer_ppi'] == 'retains all') & 
       (df['activation_abs_fold_change'] >= 1)].shape


# In[157]:


df.loc[(df['dimer_ppi'] == 'retains all') & 
       ~(
       (df['activation_abs_fold_change'] >= 1) |
       (df['pdi_n_diff'].abs() > 0) |
       (df['dbd_pct_lost'] > 0) |
       (df['other_ppi'].isin(['loses some', 'loses all']))
       )]


# In[115]:


df.loc[(df['dimer_ppi'] == 'retains all') & 
       (df['pdi_n_diff'].abs() > 0)]


# In[151]:


# retains dimer PPI loses other PPI
df.loc[(df['dimer_ppi'] == 'retains all') & 
       (df['other_ppi'].isin(['loses some', 'loses all']))]


# In[14]:


# DBD intact
from data_loading import load_dbd_accessions

def load_dbd_affected():
    df = pd.concat([g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()])
    df['is_DBD'] = df['accession'].isin(load_dbd_accessions())
    df_new = (df.loc[df['is_DBD'], :]
        .groupby(['gene', 'ref_iso', 'alt_iso'])
        [['deletion', 'frameshift']].sum()
        .sum(axis=1) / df.loc[df['is_DBD'], :]
        .groupby(['gene', 'ref_iso', 'alt_iso'])
        ['length'].sum()).to_frame(name='dbd_fraction')
    df_new['dbd_insertion_n_aa'] = (df.loc[df['is_DBD'], :]
                                  .groupby(['gene', 'ref_iso', 'alt_iso'])
                                  ['insertion']
                                  .sum())
    df = df_new.reset_index()
    df['dbd_pct_lost'] = df['dbd_fraction'] * 100.
    df = df.drop(columns=['dbd_fraction'])
    return df


dbd = load_dbd_affected()
dbd['clone_acc_ref'] = dbd['ref_iso'].map({iso.name: iso.clone_acc for tf in tfs.values() for iso in tf.cloned_isoforms})
dbd['clone_acc_alt'] = dbd['alt_iso'].map({iso.name: iso.clone_acc for tf in tfs.values() for iso in tf.cloned_isoforms})
df = pd.merge(df, dbd, how='left', on=['clone_acc_ref', 'clone_acc_alt'])


# In[15]:


df.loc[(df['dimer_ppi'] == 'retains all') & 
       (df['dbd_pct_lost'] > 0)]


# In[137]:


df['dbd_pct_lost'].isnull().sum()


# In[161]:


df['tf_gene_symbol'].map(fam).isin(DIMERIZING_TF_FAMILIES).value_counts()


# In[16]:


df['tf_family'] = df['tf_gene_symbol'].map(fam)
df['is_dimer_fam'] = df['tf_gene_symbol'].map(fam).isin(DIMERIZING_TF_FAMILIES)
df['dbd_affected'] = df['dbd_pct_lost'] > 0


# In[18]:


# try and understand differences in DBD affected between dimer and non-dimer DBD affected fraction
# confounding factors:
# length of DBD / fraction of protein of DBD
# identity of isoforms
df.groupby('is_dimer_fam')['dbd_affected'].mean()


# In[41]:


df.groupby('tf_family')['dbd_affected'].mean()


# In[165]:


# tidy up table and send
# reorder columns
# add loss of ligand binding domain?
# can we add dimer domains?
df.head()


# In[19]:


# add loses all PPIs
# how is it done in the website figures?
y2h_complete = load_y2h_isoform_data(require_at_least_one_ppi_per_isoform=False)


# In[20]:


def ppi_total_loss(row, ppi):
    ref_clone_acc = row['clone_acc_ref']
    alt_clone_acc = row['clone_acc_alt']
    if ref_clone_acc not in ppi['ad_clone_acc'].unique() or alt_clone_acc not in ppi['ad_clone_acc'].unique():
        return np.nan
    df = (ppi.loc[ppi['ad_clone_acc'].isin([ref_clone_acc, alt_clone_acc]),
                  ['ad_clone_acc', 'db_gene_symbol', 'Y2H_result']]
            .pivot(values='Y2H_result', index='db_gene_symbol', columns='ad_clone_acc')
            .dropna())
    df = df.loc[df.any(axis=1), :]
    if df.shape[0] == 0:
        return np.nan
    return not df[alt_clone_acc].any()


df['ppi_alternative_loses_all'] = df.apply(ppi_total_loss, 
                                           ppi=y2h_complete,
                                           axis=1)


# In[21]:


df['ppi_alternative_loses_all'].sum()


# In[22]:


# check there is at least one assay
# but need to include all zero ppi pairs
no_assay_data = (df['ppi_n_tested'].isnull() &
                 ~(df['ppi_alternative_loses_all'] == True) &
                 df['activation_fold_change'].isnull() &
                 df['pdi_jaccard'].isnull())
df = df.loc[~no_assay_data, :]


# In[28]:


clone_acc_to_enst = {iso.clone_acc: '|'.join(iso.ensembl_transcript_ids)
                     if iso.ensembl_transcript_ids is not None else np.nan
                     for tf in tfs.values()
                     for iso in tf.cloned_isoforms}
df['ENST_ref'] = df['clone_acc_ref'].map(clone_acc_to_enst)
df['ENST_alt'] = df['clone_acc_alt'].map(clone_acc_to_enst)


# In[31]:


# change to ppi_n_ref / alt / shared
# add has PPI/PDI/M1H data
# add M1H absolute highest across the whole gene
# add splicing type?
# add length in aa of ref and alt
# add ref is MANE
# add ref(?) / alt is novel

# add expression data?

cols = ['tf_gene_symbol',
        'tf_family',
        'is_dimer_fam',
        'clone_acc_ref',
        'clone_acc_alt',
        'ENST_ref',
        'ENST_alt',
        'aa_seq_pct_id',
        'dbd_affected',
        'dbd_insertion_n_aa',
        'dbd_pct_lost',
        'ppi_n_tested',
        'ppi_jaccard',
        'dimer_ppi',
        'other_ppi',
        'ppi_alternative_loses_all',
        'pdi_n_union',
        'pdi_n_shared',
        'pdi_jaccard',
        'm1h_min',
        'm1h_max',
        'activation_abs_fold_change']

df = df.rename(columns={'pdi_n_tested': 'pdi_n_union'})
df.loc[:, cols].to_csv('../output/TF-iso_ref-vs-alt.tsv', sep='\t', index=False)


# In[18]:


df.to_csv('../output/TF-iso_ref-vs-alt.tsv', sep='\t', index=False)


# In[40]:


y2h.loc[y2h['db_gene_symbol'] == 'ID3']


# In[138]:


# DEBUG
df.loc[df['dbd_pct_lost'].isnull(), :]

