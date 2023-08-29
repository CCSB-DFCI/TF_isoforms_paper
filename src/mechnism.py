
# coding: utf-8

# Make the `TF-iso_ref-vs-alt.tsv` table

# In[1]:


import numpy as np
import pandas as pd

from data_loading import (load_annotated_6k_collection,
                          load_y2h_isoform_data,
                          load_y1h_pdi_data,
                          load_m1h_activation_data,
                          load_valid_isoform_clones,
                          DIMERIZING_TF_FAMILIES)
from isoform_pairwise_metrics import pairs_of_isoforms_comparison_table


# In[2]:


tfs = load_annotated_6k_collection()

m1h = load_m1h_activation_data()
m1h['mean'] = m1h[['M1H_rep1', 'M1H_rep2', 'M1H_rep3']].mean(axis=1)
m1h['abs_mean'] = m1h['mean'].abs()

y1h = load_y1h_pdi_data()
y1h = y1h.drop_duplicates(subset=['unique_acc'])  # TODO: move code to data loading
y1h = y1h.set_index('unique_acc')  # TODO: move code to data loading

y2h = load_y2h_isoform_data()


# In[3]:


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


# In[4]:


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


# In[5]:


from data_loading import load_seq_comparison_data

aa_ident = load_seq_comparison_data()
df["aa_seq_pct_id"] = df.apply(
    lambda x: "_".join(sorted([x["clone_acc_ref"], x["clone_acc_alt"]])), axis=1
).map(aa_ident)
if df['aa_seq_pct_id'].isnull().any():
    raise UserWarning('Unexpected missing sequence similarity values')


# In[6]:


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


# In[7]:


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


# In[8]:


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


# In[10]:


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
    # kaia edited these 2 lines to drop any baits with NA
    # my version of pandas wasn't auto-filtering these out, think it got fixed in later versions
    a = set(df.columns[df.iloc[0].fillna(False)])
    b = set(df.columns[df.iloc[1].fillna(False)])
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


# In[11]:


df['at_least_one_isoform_in_gene_abs_activation_gte_2fold'] = df['gene_symbol'].map(m1h.groupby('gene')['abs_mean'].max() >= 1)
df['activation_ref'] = df['clone_acc_ref'].map(m1h.set_index('clone_acc')['mean'])
df['activation_alt'] = df['clone_acc_alt'].map(m1h.set_index('clone_acc')['mean'])
df['activation_fold_change_log2'] = (df['activation_alt'] - df['activation_ref'])


# In[12]:


df.to_csv('../output/TF-iso_ref-vs-alt.tsv', sep='\t', index=False)

