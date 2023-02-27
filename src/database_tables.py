from pathlib import Path

import pandas as pd

from ccsblib import paros_connection
from ccsblib import huri


def isoform_clones():
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
                     paros_connection())
    return df


def tf_isoform_y2h_screen_results():
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
    df = pd.read_sql(qry, paros_connection())
    df = (pd.get_dummies(df, columns=['pool_name']).groupby(['ad_orf_id', 'db_orf_id']).sum() > 0).reset_index()
    df = df.rename(columns={'pool_name_ds20180213_TFisoS04': 'in_orfeome_screen',
                            'pool_name_ds20180213_TFisoS05': 'in_focussed_screen'})
    if not (df['in_orfeome_screen'] | df['in_focussed_screen']).all():
        raise UserWarning('Something went wrong...')
    return df


def _annotate_ppi_source_info(df):

    def load_orf_id_to_tf_gene():
        qry = """SELECT orf_id,
                        symbol AS tf_gene_symbol
                FROM tf_screen.iso6k_sequences;"""
        df = pd.read_sql(qry, paros_connection())
        df = df.drop_duplicates()
        if df['orf_id'].duplicated().any():
            raise UserWarning('Unexpected duplicate ORF IDs')
        return df.set_index('orf_id')['tf_gene_symbol']

    screen = tf_isoform_y2h_screen_results()
    orf_id_to_tf_gene = load_orf_id_to_tf_gene()
    screen['ad_gene_symbol'] = screen['ad_orf_id'].map(orf_id_to_tf_gene)
    screen = (screen.dropna(subset=['ad_gene_symbol'])
                    .drop(columns=['ad_orf_id']))
    screen = screen.groupby(['ad_gene_symbol', 'db_orf_id']).any().reset_index()
    df = pd.merge(df,
                  screen,
                  how='left',
                  on=['ad_gene_symbol', 'db_orf_id'])
    df['in_orfeome_screen'] = df['in_orfeome_screen'].fillna(False)
    df['in_focussed_screen'] = df['in_focussed_screen'].fillna(False)
    hiu = huri.load_nw_hi_union(id_type='orf_id')
    id_map = huri.load_id_map('orf_id', 'hgnc_symbol', via='ensembl_gene_id')
    hiu_pairs = set((pd.merge(hiu, id_map, how='inner',
                    left_on='orf_id_a', right_on='orf_id')
                    [['hgnc_symbol', 'orf_id_b']]
                    .drop_duplicates())
                    .apply(lambda x: x['hgnc_symbol'] + '-' + str(x['orf_id_b']),
                           axis=1))
    hiu_pairs = hiu_pairs.union((pd.merge(hiu, id_map, how='inner',
                                          left_on='orf_id_b', right_on='orf_id')
                                 [['hgnc_symbol', 'orf_id_a']]
                                 .drop_duplicates())
                                .apply(lambda x: x['hgnc_symbol'] + '-' + str(x['orf_id_a']),
                                axis=1))
    df['in_hi_union'] = (df['ad_gene_symbol'] + '-' + df['db_orf_id'].astype(str)).isin(hiu_pairs)
    litbm = huri.load_nw_lit_bm_17(id_type='orf_id')
    litbm_pairs = set((pd.merge(litbm, id_map, how='inner',
                                left_on='orf_id_a', right_on='orf_id')
                      [['hgnc_symbol', 'orf_id_b']]
                      .drop_duplicates())
                      .apply(lambda x: x['hgnc_symbol'] + '-' + str(x['orf_id_b']),
                             axis=1))
    litbm_pairs = litbm_pairs.union((pd.merge(litbm, id_map, how='inner',
                                              left_on='orf_id_b', right_on='orf_id')
                                    [['hgnc_symbol', 'orf_id_a']]
                                    .drop_duplicates())
                                    .apply(lambda x: x['hgnc_symbol'] + '-' + str(x['orf_id_a']),
                                    axis=1))
    df['in_lit_bm'] = (df['ad_gene_symbol'] + '-' + df['db_orf_id'].astype(str)).isin(litbm_pairs)
    return df


def isoform_and_paralog_y2h_data():
    """
    - NS: sequencing failed
    - NC: no call (e.g., mis-spotting)
    - AA: autoactivator

    """
    qry_a = """SELECT a.category,
                      a.ad_orf_id,
                      b.unique_acc AS ad_clone_acc,
                      a.ad_symbol AS ad_gene_symbol,
                      a.db_orf_id,
                      c.symbol AS db_gene_symbol,
                      a.final_score AS score,
                      d.standard_batch,
                      d.retest_pla,
                      d.retest_pos
                 FROM tf_screen.tf_isoform_final AS a
                 LEFT JOIN tf_screen.iso6k_sequences AS b
                   ON a.ad_orf_id = b.orf_id
                 LEFT JOIN horfeome_annotation_gencode27.orf_class_map_ensg AS c
                   ON a.db_orf_id = c.orf_id
                 LEFT JOIN tf_screen.retest AS d
                   ON a.retest_id = d.retest_id;"""
    df_a = pd.read_sql(qry_a, paros_connection())

    # remove duplicate ORF for gene DDX39B, where sequencing mostly failed
    df_a = df_a.loc[df_a['db_orf_id'] != 3677, :]

    df_a['category'] = df_a['category'].map({'ppi': 'tf_isoform_ppis',
                                             'ng_stem_cell_factor': 'tf_isoform_ppis',
                                             'rrs': 'rrs_isoforms',
                                             'litbm': 'lit_bm_isoforms'})
    qry_b = """SELECT a.simple_category AS category,
                      a.ad_orf_id,
                      b.unique_acc AS ad_clone_acc,
                      a.ad_symbol AS ad_gene_symbol,
                      a.db_orf_id,
                      c.symbol AS db_gene_symbol,
                      a.final_score AS score,
                      d.standard_batch,
                      d.retest_pla,
                      d.retest_pos
                 FROM tf_screen.paralog_final AS a
                 LEFT JOIN tf_screen.iso6k_sequences AS b
                   ON a.ad_orf_id = b.orf_id
                 LEFT JOIN horfeome_annotation_gencode27.orf_class_map_ensg AS c
                   ON a.db_orf_id = c.orf_id
                 LEFT JOIN tf_screen.retest AS d
                   ON a.retest_id = d.retest_id;"""
    df_b = pd.read_sql(qry_b, paros_connection())
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

    non_protein_coding_genes = """RBMY2FP
                                  REXO1L6P
                                  SUMO1P1
                                  SLCO4A1-AS1
                                  ZNF213-AS1
                                  AC023509.3
                                  KRT8P41
                                  LINC01658""".split()
    df = df.loc[~df['db_gene_symbol'].isin(non_protein_coding_genes), :]

    # remove unexpected duplicates, some of which have inconsistent scores
    # see Issue #28
    df = df.loc[~(df.duplicated(subset=['ad_clone_acc', 'db_gene_symbol'],
                                keep=False)
                  & (df['category'] == 'paralog_with_PDI')), :]

    df = _annotate_ppi_source_info(df)
    return df


func_to_file_name = [(isoform_clones, 'isoform_clones'),
                     (tf_isoform_y2h_screen_results, 'tf_isoform_y2h_screen'),
                     (isoform_and_paralog_y2h_data, 'y2h_data')]
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'data/internal'


def write_all_tables():
    for func, filename in func_to_file_name:
        func().to_csv(DATA_DIR / (filename + '.tsv'),
                      sep='\t',
                      index=False)


def check_consistency():
    all_good = True
    for func, filename in func_to_file_name:
        a = func()
        b = pd.read_csv(DATA_DIR / (filename + '.tsv'),
                        sep='\t',
                        na_values=[''],
                        keep_default_na=False)
        if not a.equals(b):
            all_good = False
            print(filename + ': inconsitent with database')
    if all_good:
        print('All files consistent with database')
