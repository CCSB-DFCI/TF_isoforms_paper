import numpy as np
import pandas as pd
import ccsblib



def load_valid_isoform_clones():
    """The subset of TF isoform clones that have passed the stringent
       annotation process (at-length AA match to GENCODE, or approved by
       Gloria/GENCODE team).

    Returns:
        set(str): clone accession IDs

    """
    df = pd.read_sql("""SELECT symbol as gene,
                               unique_acc as clone_acc,
                               dup_idx
                          FROM tf_screen.iso6k_annotation
                         ORDER BY gene, clone_acc;""",
                       ccsblib.paros_connection())
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
    return df


def load_tf_isoform_screen_results():
    """There were two screens performed:

    The cloned TF isoforms as AD-fusions against DB-fusions of:
    (1) ORFeome 9.1
    (2) Subset of TFs and co-factors

    Returns:
        pandas.DataFrame: for each pair, was it found in the first and second screens

    """
    qry = """SELECT ad_orf_id,
                    db_orf_id,
                    pool_name
               FROM swimseq_run.INGS_IST a,
                    swimseq_run.NGS_POOL b
              WHERE a.pool_id = b.pool_id
                AND a.pool_id in (787, 788)
                AND ist_score>=0.2;"""
    df = pd.read_sql(qry, ccsblib.paros_connection())
    df = (pd.get_dummies(df, columns=['pool_name']).groupby(['ad_orf_id', 'db_orf_id']).sum() > 0).reset_index()
    df = df.rename(columns={'pool_name_ds20180213_TFisoS04': 'in_orfeome_screen',
                            'pool_name_ds20180213_TFisoS05': 'in_focussed_screen'})
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
    qry_a = """select a.category,
                      a.ad_orf_id,
                      b.unique_acc AS ad_clone_acc,
                      a.ad_symbol AS ad_gene_symbol,
                      a.db_orf_id,
                      c.symbol AS db_gene_symbol,
                      a.final_score AS score
                 FROM tf_screen.tf_isoform_final AS a
                 LEFT JOIN tf_screen.iso6k_sequences AS b
                   ON a.ad_orf_id = b.orf_id
                 LEFT JOIN horfeome_annotation_gencode27.orf_class_map_ensg AS c
                   ON a.db_orf_id = c.orf_id;"""
    df_a = pd.read_sql(qry_a, ccsblib.paros_connection())
    if filter_for_valid_clones:
        valid_clones = load_valid_isoform_clones()
        df_a = df_a.loc[df_a['ad_clone_acc'].isin(valid_clones['clone_acc']), :]

    # remove duplicate ORF for gene DDX39B, where sequencing mostly failed
    df_a = df_a.loc[df_a['db_orf_id'] != 3677, :]

    df_a['category'] = df_a['category'].map({'ppi': 'tf_isoform_ppis',
                                             'ng_stem_cell_factor': 'tf_isoform_ppis',
                                             'rrs': 'rrs_isoforms',
                                             'litbm': 'lit_bm_isoforms'})
    qry_b = """select a.simple_category AS category,
                      a.ad_orf_id,
                      b.unique_acc AS ad_clone_acc,
                      a.ad_symbol AS ad_gene_symbol,
                      a.db_orf_id,
                      c.symbol AS db_gene_symbol,
                      a.final_score AS score
                 FROM tf_screen.paralog_final AS a
                 LEFT JOIN tf_screen.iso6k_sequences AS b
                   ON a.ad_orf_id = b.orf_id
                 LEFT JOIN horfeome_annotation_gencode27.orf_class_map_ensg AS c
                   ON a.db_orf_id = c.orf_id;"""
    df_b = pd.read_sql(qry_b, ccsblib.paros_connection())
    if filter_for_valid_clones:
        df_b = df_b.loc[df_b['ad_clone_acc'].isin(valid_clones), :]
    df_b['category'] = df_b['category'].map({'paralog': 'tf_paralog_ppis',
                                             'PDI_PPI': 'paralog_with_PDI',
                                             'nonparalog': 'non_paralog_control',
                                             'rrs': 'rrs_paralogs',
                                             'litbm': 'lit_bm_paralogs'})
    df = pd.concat([df_a, df_b])
    # drop cases where single orf ID mapped to multiple gene symbols
    df = df.drop_duplicates(['category', 'ad_orf_id', 'db_orf_id'])
    # drop interaction partner ORFs whose sequence does not map to an ensembl gene
    df = df.dropna(subset=['db_gene_symbol'])

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

    """
    # Need to map the screen data to the gene level first
    screen = load_tf_isoform_screen_results()
    pd.merge(df,
            screen,
            how='left',
            on=['ad_orf_id', 'db_orf_id']).sort_values(['ad_gene_symbol', 'db_gene_symbol'])
    """
    return df


def load_y1h_pdi_data(add_missing_data=False):
    df = pd.read_csv('../../data/a2_juan_pdi_w_unique_isoacc.tsv', sep='\t')
    df = (pd.concat([df.loc[:, ['tf', 'unique_acc']],
                     pd.get_dummies(df['bait'])],
                    axis=1)
            .groupby(['tf', 'unique_acc']).sum() > 0).reset_index()
    zeros = pd.read_csv('../../data/a2_juan_isoforms_wo_pdi.tsv', sep='\t')
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
    df = pd.read_csv('../../data/a_m1h_final_table.tsv', sep='\t')
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
    df = pd.read_table('../../data/tf_iso_expression/a_kallisto_hpa_tpms_prot_seq_grouped_w_lambert.tsv')
    df = df.loc[df.target_id.str.contains('/'), :]
    idxs = [x for x in df.columns if not x.startswith('ERR')]
    df = df.set_index(idxs)
    df2 = df.stack().reset_index()
    df2[['err', 'ers']] = df2.level_5.str.split('|', expand=True)

    # sample manifest
    dfm = pd.read_table('../../data/tf_iso_expression/b_sample_manifest_E-MTAB-2836.sdrf.txt')
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
    df = pd.read_table('../../data/tf_AA_seq_identities/a_2019-09-10_AA_seq_identity_Isoform_series_for_all_6Kpairs_unique_acc.txt')
    df['pair'] = df.apply(lambda x: '_'.join(sorted([x.iso1, x.iso2])), axis=1)
    df = df[['pair', 'AAseq_identity%']]
    df.columns = ['pair', 'aa_seq_pct_id']
    if df['pair'].duplicated().any():
        raise UserWarning('Unexpected duplicates')
    df = df.set_index('pair')
    if (df['aa_seq_pct_id'] < 0).any() or (df['aa_seq_pct_id'] > 100).any():
        raise UserWarning('Percent values outside 0-100')
    return df['aa_seq_pct_id']
