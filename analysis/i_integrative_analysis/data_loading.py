from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


DATA_DIR = Path(__file__).resolve().parents[2] / 'data'


def load_valid_isoform_clones():
    """The subset of TF isoform clones that have passed the stringent
       annotation process (at-length AA match to GENCODE, or approved by
       Gloria/GENCODE team).

    Returns:
        set(str): clone accession IDs

    """
    df = pd.read_csv(DATA_DIR / 'internal/isoform_clones.tsv', sep='\t')
    y2h = load_isoform_and_paralog_y2h_data(filter_for_valid_clones=False)
    y1h = load_y1h_pdi_data()
    m1h = load_m1h_activation_data()
    df['in_m1h'] = df['clone_acc'].isin(m1h['clone_acc'])
    df['in_y1h'] = df['clone_acc'].isin(y1h['unique_acc'])
    df['in_y2h'] = df['clone_acc'].isin(y2h.loc[(y2h['category'] == 'tf_isoform_ppis') &
                                                y2h['score'].isin(['0', '1']),
                                                'ad_clone_acc'])
    # dropping duplicates with identical AA seqs, keeping those with M1H data
    df = df.sort_values(['gene', 'in_m1h', 'in_y2h', 'in_y1h'], ascending=[True] + [False] * 3)
    df = (df.loc[~df['dup_idx'].duplicated() | df['dup_idx'].isnull(), ['gene', 'clone_acc']]
            .sort_values(['gene', 'clone_acc']))
    aa_seq_file = DATA_DIR / 'internal/j2_6k_unique_isoacc_and_prot_seqs.fa'
    aa = {r.id.split('xxx')[1]: str(r.seq) for r in
          SeqIO.parse(aa_seq_file, format='fasta')}
    df['aa_seq'] = df['clone_acc'].map(aa)
    df['num_aa'] = df['aa_seq'].str.len()
    iso_annot = pd.read_csv(DATA_DIR / 'internal/c_conso_annot_table_man_annot.tsv',
                            sep='\t')
    df['is_novel_isoform'] = df['clone_acc'].map(iso_annot.set_index('unique_acc')['gc_match'] == 0)
    return df


def load_tf_isoform_y2h_screen_results():
    """There were two screens performed:

    The cloned TF isoforms as AD-fusions against DB-fusions of:
    (1) ORFeome 9.1
    (2) Subset of TFs and co-factors

    Returns:
        pandas.DataFrame: for each pair, was it found in the first and second screens

    """
    df = pd.read_csv(DATA_DIR / 'internal/tf_isoform_y2h_screen.tsv', sep='\t')
    if not (df['in_orfeome_screen'] | df['in_focussed_screen']).all():
        raise UserWarning('Something went wrong...')
    return df


def load_isoform_and_paralog_y2h_data(add_missing_data=False, filter_for_valid_clones=True):
    """
    - NS: sequencing failed
    - NC: no call (e.g., mis-spotting)
    - AA: autoactivator

    """
    if add_missing_data and not filter_for_valid_clones:
        raise ValueError('This combination of arguments will not work')
    df = pd.read_csv(DATA_DIR / 'internal/y2h_data.tsv',
                     sep='\t',
                     na_values=[''],
                     keep_default_na=False)
    if filter_for_valid_clones:
        valid_clones = load_valid_isoform_clones()
        df = df.loc[df['ad_clone_acc'].isin(valid_clones['clone_acc']), :]
    if add_missing_data:
        all_possible_ints = (pd.merge(df.loc[df['category'] == 'tf_isoform_ppis',
                                             ['ad_gene_symbol', 'db_gene_symbol', 'category']]
                                        .drop_duplicates(),
                                      valid_clones,
                                      left_on='ad_gene_symbol',
                                      right_on='gene')
                             .drop(columns='gene')
                             .rename(columns={'clone_acc': 'ad_clone_acc'}))
        df = pd.merge(df,
                      all_possible_ints,
                      on=['ad_gene_symbol', 'db_gene_symbol', 'ad_clone_acc', 'category'],
                      how='outer')
        df['score'] = df['score'].fillna('NA')
    cat_info = load_ppi_partner_categories()
    cats = cat_info.groupby('category')['partner'].apply(set).to_dict()
    for cat, members in cats.items():
        df['is_partner_category_' + '_'.join(cat.split())] = df['db_gene_symbol'].isin(members)
    cofac_type = cat_info.groupby('cofactor_type')['partner'].apply(set).to_dict()
    for subtype, members in cofac_type.items():
        df['is_cofactor_subtype_' + subtype] = df['db_gene_symbol'].isin(members)
    return df


def load_y1h_pdi_data(add_missing_data=False):
    df = pd.read_csv(DATA_DIR / 'internal/a2_juan_pdi_w_unique_isoacc.tsv', sep='\t')
    df = (pd.concat([df.loc[:, ['tf', 'unique_acc']],
                     pd.get_dummies(df['bait'])],
                    axis=1)
            .groupby(['tf', 'unique_acc']).sum() > 0).reset_index()
    zeros = pd.read_csv(DATA_DIR / 'internal/a2_juan_isoforms_wo_pdi.tsv', sep='\t')
    df = pd.concat([df, zeros], axis=0, sort=False).reset_index(drop=True).fillna(False)
    if add_missing_data:
        isoforms = load_valid_isoform_clones()
        df = (pd.merge(df,
                       isoforms.rename(columns={'gene': 'tf', 'clone_acc': 'unique_acc'}),
                       on=['tf', 'unique_acc'],
                       how='outer'))
    df = df.sort_values(['tf', 'unique_acc'])
    return df


def load_m1h_activation_data(add_missing_data=False):
    """


    """
    df = pd.read_csv(DATA_DIR / 'internal/a_m1h_final_table.tsv', sep='\t')
    df = df.rename(columns={'pos_acc': 'clone_acc'})
    for column in df.columns:
        if column.startswith('M1H_rep'):
            df[column] = np.log2(df[column])
    if add_missing_data:
        isoforms = load_valid_isoform_clones()
        df = pd.merge(df, isoforms, on=['gene', 'clone_acc'], how='outer')
    df = df.sort_values(['gene', 'clone_acc'])
    return df


def load_rna_expression_data():
    """
    Isoform clone expression across HPA tissues.
    """
    df = pd.read_table(DATA_DIR / 'internal/a_kallisto_hpa_tpms_prot_seq_grouped_w_lambert.tsv')
    df = df.loc[df.target_id.str.contains('/'), :]
    idxs = [x for x in df.columns if not x.startswith('ERR')]
    df = df.set_index(idxs)
    df2 = df.stack().reset_index()
    df2[['err', 'ers']] = df2.level_5.str.split('|', expand=True)

    # sample manifest
    dfm = pd.read_table(DATA_DIR / 'internal/b_sample_manifest_E-MTAB-2836.sdrf.txt')
    dfm = dfm[['Comment[ENA_RUN]', 'Source Name']]
    dfm.columns = ['err', 'tiss']
    dfm['tiss'] = dfm['tiss'].apply(lambda x: x.split('_')[0])

    # add tissue type to the expression matrix
    df3 = df2.merge(dfm, how='left', on='err')

    df3 = df3[['gene', 'target_id', 0, 'tiss']]
    df3.columns = ['gene', 'isoform', 'tpm', 'tiss']
    df4 = df3.groupby(['gene', 'isoform', 'tiss']).agg({'tpm': ['mean', 'std']}).reset_index()
    df4.columns = df4.columns.get_level_values(0)
    df4.columns = ['gene', 'isoform', 'tiss', 'tpm', 'tpm_stdev']
    df4 = df4[['gene', 'isoform', 'tiss', 'tpm', 'tpm_stdev']]
    df4['isoacc'] = df4['isoform'].str.split('_').str.get(0)
    df4 = df4[['gene', 'isoacc', 'tiss', 'tpm', 'tpm_stdev']]
    # write out table
    # df4.to_csv('expression_table_tfisoclones.tsv', sep='\t', index=False)
    return df4


def load_seq_comparison_data():
    """
    Pairwise sequence comparisons of AA.
    Needleman algorithm, global alignment.
    Note - I checked and there are no duplicate rows.
    """
    df = pd.read_table(DATA_DIR / 'processed/a_2019-09-10_AA_seq_identity_Isoform_series_for_all_6Kpairs_unique_acc.txt')
    df['pair'] = df.apply(lambda x: '_'.join(sorted([x.iso1, x.iso2])), axis=1)
    df = df[['pair', 'AAseq_identity%']]
    df.columns = ['pair', 'aa_seq_pct_id']

    df_b = pd.read_csv(DATA_DIR / 'processed/paralog_non_paralog_seq_id.tsv',
                       sep='\t').drop_duplicates()
    df_b['pair'] = df_b.apply(lambda x: '_'.join(sorted([x['clone_acc_a'], x['clone_acc_b']])), axis=1)
    df_b = df_b.loc[:, ['pair', 'aa_seq_pct_id']]
    df = pd.concat([df, df_b])

    if df['pair'].duplicated().any():
        raise UserWarning('Unexpected duplicates')
    df = df.set_index('pair')
    if (df['aa_seq_pct_id'] < 0).any() or (df['aa_seq_pct_id'] > 100).any():
        raise UserWarning('Percent values outside 0-100')
    return df['aa_seq_pct_id']


def load_paralog_pairs():
    """Pairs of TF gene paralogs and non-paralogs that were tested

    Returns:
        pandas.DataFrame: one row for each pair

    """
    df = pd.read_csv(DATA_DIR / 'internal/a_tf_iso_paralog_nonparalogs_tested.tsv', sep='\t')
    df['is_paralog_pair'] = (df['cat2'] == 'paralog')
    aa = pd.read_csv(DATA_DIR / 'processed/b_2018-11-30_AA_seq_identity_Paralog_comparisons_unique_acc.txt',
                     sep='\t')
    if not (df['tf1'] < df['tf2']).all():
        raise UserWarning('Expected genes to be ordered')
    (aa[['gene1', 'gene2']].min(axis=1) + '_' + aa[['gene1', 'gene2']].max(axis=1)).duplicated().any()
    aa['tf1'] = aa[['gene1', 'gene2']].min(axis=1)
    aa['tf2'] = aa[['gene1', 'gene2']].max(axis=1)
    df = (pd.merge(df, aa, how='left', on=['tf1', 'tf2'])
            .loc[:, ['tf1', 'tf2', 'is_paralog_pair', 'AAseq_identity%']]
            .rename(columns={'tf1': 'tf_gene_a',
                             'tf2': 'tf_gene_b',
                             'AAseq_identity%': 'pct_aa_seq_identity'}))
    df = df.drop_duplicates()
    return df


def load_isoforms_of_paralogs_pairs(pairs, isoforms):
    pairs = pd.merge(pairs.loc[:, ['tf_gene_a', 'tf_gene_b', 'is_paralog_pair']],
                     isoforms.loc[:, ['gene', 'clone_acc', 'aa_seq']],
                     how='inner',
                     left_on='tf_gene_a',
                     right_on='gene').rename(columns={'aa_seq': 'aa_seq_a', 'clone_acc': 'clone_acc_a'})
    pairs = pd.merge(pairs,
                     isoforms.loc[:, ['gene', 'clone_acc', 'aa_seq']],
                     how='inner',
                     left_on='tf_gene_b',
                     right_on='gene').rename(columns={'aa_seq': 'aa_seq_b', 'clone_acc': 'clone_acc_b'})
    pairs = pairs.drop(columns=['gene_x', 'gene_y'])
    return pairs


def load_ppi_partner_categories():
    """Juan's manual classification of the PPI interaction partners.

     Note that a gene can be in multiple categories.

    Returns:
        pandas.DataFrame: gene and category

    """
    df = pd.read_excel(DATA_DIR / 'internal/20191028- Uniprot functions for interactors.xlsx',
                       sheet_name='Final')
    if df['Function class'].isnull().any():
        raise UserWarning('Unexpected missing values')
    if df['partner'].duplicated().any():
        raise UserWarning('Unexpected duplicate entries')
    df = df.set_index('partner')
    cofac_type = df['Cofactor type?'].copy()
    df = df['Function class'].str.split(', ', expand=True).stack().reset_index().loc[:, ['partner', 0]].rename(columns={0: 'category'})
    df['category'] = df['category'].str.strip()
    cf_rows = df['category'] == 'cofactor'
    df.loc[cf_rows, 'cofactor_type'] = (df.loc[cf_rows, 'partner']
                                          .map(cofac_type))
    if not (df['partner'] == df['partner'].str.strip()).all():
        raise UserWarning('Possibly something wrong with gene names column')
    return df


def load_tf_families():
    """From the Lambert et al. review in Cell 2018

    Returns:
        pandas.Series: HGNC gene symbol to TF family

    """
    tf_fam = pd.read_csv(DATA_DIR / 'external/Human_TF_DB_v_1.01.csv')
    tf_fam = tf_fam.loc[tf_fam['Is TF?'] == 'Yes', ['HGNC symbol', 'DBD']]
    if tf_fam['HGNC symbol'].duplicated().any():
        raise UserWarning('Unexpected duplicates')
    tf_fam = tf_fam.set_index('HGNC symbol')['DBD']
    return tf_fam
