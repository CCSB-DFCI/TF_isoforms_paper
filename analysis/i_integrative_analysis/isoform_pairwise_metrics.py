from itertools import combinations

import numpy as np
import pandas as pd

from data_loading import (load_seq_comparison_data,
                          load_isoforms_of_paralogs_pairs)


def paralog_pair_ppi_table(data, tf_gene_a, tf_gene_b):
    gene_a_partners = data.loc[(data['ad_gene_symbol'] == tf_gene_a) &
                               (data['category'] == 'tf_isoform_ppis') &
                               (data['score'] == '1'),
                               'db_gene_symbol'].unique()
    gene_b_partners = data.loc[(data['ad_gene_symbol'] == tf_gene_b) &
                               (data['category'] == 'tf_isoform_ppis') &
                               (data['score'] == '1'),
                               'db_gene_symbol'].unique()
    partners = np.concatenate([gene_a_partners, gene_b_partners])
    tf = data.loc[((data['ad_gene_symbol'] == tf_gene_a) |
                  (data['ad_gene_symbol'] == tf_gene_b))
                  & data['category'].isin(['tf_isoform_ppis',
                                           'tf_paralog_ppis',
                                           'non_paralog_control',
                                           'paralog_with_PDI'])
                  & data['db_gene_symbol'].isin(partners),
                  ['ad_clone_acc', 'db_gene_symbol', 'score']].copy()
    return tf


def pairs_of_paralogs_and_isoforms_comparison_table(isoforms,
                                                    paralog_pairs,
                                                    y2h,
                                                    y1h,
                                                    m1h):
    """[summary]

    Restricted to isoforms pairs in the paralogs dataset.

    Args:
        isoforms (bool): [description]
        paralog_pairs ([type]): [description]
        y2h ([type]): [description]
        y1h ([type]): [description]
        m1h ([type]): [description]

    Returns:
        [type]: [description]

    """
    pairs = load_isoforms_of_paralogs_pairs(paralog_pairs, isoforms)
    pairs['category'] = pairs['is_paralog_pair'].map({True: 'paralogs', False: 'non-paralog-control'})
    pairs = pairs.drop(columns=['aa_seq_a', 'aa_seq_b', 'is_paralog_pair'])
    iso_pairs = []
    for tf_gene in isoforms['gene'].unique():
        tf_iso = isoforms.loc[isoforms['gene'] == tf_gene,
                              'clone_acc'].values
        for iso_a, iso_b in combinations(tf_iso, 2):
            iso_pairs.append((tf_gene, tf_gene, iso_a, iso_b))
    iso_pairs = pd.DataFrame(data=iso_pairs,
                             columns=['tf_gene_a', 'tf_gene_b', 'clone_acc_a', 'clone_acc_b'])
    iso_pairs['category'] = 'isoforms'
    iso_pairs = iso_pairs.loc[iso_pairs['tf_gene_a'].isin(pairs['tf_gene_a']) |
                              iso_pairs['tf_gene_a'].isin(pairs['tf_gene_b']), :]
    pairs = pd.concat([pairs, iso_pairs])
    pairs['pair'] = pairs.apply(lambda x: '_'.join(sorted([x['clone_acc_a'],
                                                           x['clone_acc_b']])),
                                axis=1)
    return _pairs_comparison_table(pairs, y2h, y1h, m1h)


def pairs_of_isoforms_comparison_table(isoforms, y2h, y1h, m1h):
    iso_pairs = []
    for tf_gene in isoforms['gene'].unique():
        tf_iso = isoforms.loc[isoforms['gene'] == tf_gene,
                              'clone_acc'].values
        for iso_a, iso_b in combinations(tf_iso, 2):
            iso_pairs.append((tf_gene, iso_a, iso_b))
    iso_pairs = pd.DataFrame(data=iso_pairs,
                             columns=['tf_gene_symbol', 'clone_acc_a', 'clone_acc_b'])
    iso_pairs['pair'] = iso_pairs.apply(lambda x: '_'.join(sorted([x.clone_acc_a, x.clone_acc_b])), axis=1)
    return _pairs_comparison_table(iso_pairs, y2h, y1h, m1h)


def _pairs_comparison_table(pairs, y2h, y1h, m1h):
    pairs['ppi_n_tested'] = pairs.apply(ppi_metric,
                                        data=y2h,
                                        function=number_tested_partners,
                                        axis=1)
    pairs['ppi_n_shared'] = pairs.apply(ppi_metric,
                                        data=y2h,
                                        function=number_shared_partners,
                                        axis=1)
    pairs['ppi_n_min'] = pairs.apply(ppi_metric,
                                     data=y2h,
                                     function=number_min_partners,
                                     axis=1)
    pairs['ppi_n_min_diff'] = pairs.apply(ppi_metric,
                                          data=y2h,
                                          function=min_difference,
                                          axis=1)
    pairs['ppi_jaccard'] = pairs.apply(ppi_metric,
                                       data=y2h,
                                       function=jaccard_index,
                                       axis=1)
    pairs['ppi_simpson'] = pairs.apply(ppi_metric,
                                       data=y2h,
                                       function=simpsons_index,
                                       axis=1)
    pairs['ppi_n_diff'] = pairs['ppi_n_tested'] - pairs['ppi_n_shared']

    pairs['pdi_n_tested'] = pairs.apply(pdi_metric,
                                        data=y1h,
                                        function=number_tested_partners,
                                        axis=1)
    pairs['pdi_n_shared'] = pairs.apply(pdi_metric,
                                        data=y1h,
                                        function=number_shared_partners,
                                        axis=1)
    pairs['pdi_n_min'] = pairs.apply(pdi_metric,
                                     data=y1h,
                                     function=number_min_partners,
                                     axis=1)
    pairs['pdi_n_min_diff'] = pairs.apply(pdi_metric,
                                          data=y1h,
                                          function=min_difference,
                                          axis=1)
    pairs['pdi_jaccard'] = pairs.apply(pdi_metric,
                                       data=y1h,
                                       function=jaccard_index,
                                       axis=1)
    pairs['pdi_simpson'] = pairs.apply(pdi_metric,
                                       data=y1h,
                                       function=simpsons_index,
                                       axis=1)
    pairs['pdi_n_diff'] = pairs['pdi_n_tested'] - pairs['pdi_n_shared']

    pairs['activation_fold_change'] = pairs.apply(fold_change_m1h, data=m1h, axis=1)

    aa_ident = load_seq_comparison_data()
    pairs['aa_seq_pct_id'] = (pairs.apply(lambda x: '_'.join(sorted([x['clone_acc_a'],
                                                                     x['clone_acc_b']])),
                                          axis=1)
                                   .map(aa_ident))
    if pairs['aa_seq_pct_id'].isnull().any():
        raise UserWarning('Problem with sequence similarity data')
    return pairs


def ppi_metric(row, data, function):
    ad_a = row['clone_acc_a']
    ad_b = row['clone_acc_b']
    pair = (data.loc[data['ad_clone_acc'].isin([ad_a, ad_b]), :]
                .pivot(values='score',
                       index='db_gene_symbol',
                       columns='ad_clone_acc'))
    if ad_a not in pair.columns or ad_b not in pair.columns:
        return np.nan
    # remove any partner with AA / NC / NS / NaN in either
    pair = pair.loc[pair.isin(['0', '1']).all(axis=1), :].astype(int).astype(bool)
    # remove partners that tested negative in both
    pair = pair.loc[pair.any(axis=1), :]
    if pair.shape[0] > 0:
        return function(set(pair.index[pair[ad_a]].values),
                        set(pair.index[pair[ad_b]].values))
    else:
        return np.nan


def pdi_metric(row, data, function):
    df = data.loc[(data['unique_acc'] == row['clone_acc_a']) |
                  (data['unique_acc'] == row['clone_acc_b']),
                  data.columns[2:]].copy()
    if df.shape[0] < 2:
        return np.nan
    df = df.loc[:, df.any(axis=0)]
    if df.shape[1] == 0:
        return np.nan
    a = set(df.columns[df.iloc[0]])
    b = set(df.columns[df.iloc[1]])
    return function(a, b)


def jaccard_index(a, b):
    return len(a.intersection(b)) / float(len(a.union(b)))


def simpsons_index(a, b):
    min_size =  min(len(a), len(b))
    if min_size == 0:
        return np.nan
    else:
        return len(a.intersection(b)) / float(min_size)


def number_tested_partners(a, b):
    """Comes up with nan when it should be 0?"""
    return len(a.union(b))


def number_shared_partners(a, b):
    return len(a.intersection(b))


def number_min_partners(a, b):
    return min(len(a), len(b))


def min_difference(a, b):
    return min(len(a.difference(b)), len(b.difference(a)))


def fold_change_m1h(row, data):
    if (row['clone_acc_a'] not in data['clone_acc'].values or
        row['clone_acc_b'] not in data['clone_acc'].values):
        return np.nan
    a = data.loc[data['clone_acc'] == row['clone_acc_a'],
                [c for c in data.columns if c.startswith('M1H_rep')]].mean(axis=1).values[0]
    b = data.loc[data['clone_acc'] == row['clone_acc_b'],
                [c for c in data.columns if c.startswith('M1H_rep')]].mean(axis=1).values[0]
    return max(a, b) - min(a, b)
