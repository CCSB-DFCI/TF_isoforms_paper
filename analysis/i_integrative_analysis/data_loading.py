import os
import sys
from pathlib import Path
import itertools
import warnings
from collections import defaultdict
import functools

import numpy as np
import pandas as pd
from Bio import SeqIO, Align
import tqdm

import isolib

DATA_DIR = Path(__file__).resolve().parents[2] / 'data'
CACHE_DIR = Path(__file__).resolve().parents[2] / 'cache'

sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)), '../..'))
from isomodules import isocreate, isofunc


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


    nt_seq_file = DATA_DIR / 'internal/j2_6k_unique_isoacc_and_nt_seqs.fa'
    nt = {r.id: str(r.seq) for r in SeqIO.parse(nt_seq_file, format='fasta')}
    df['cds'] = df['clone_acc'].map(nt)

    aa_seq_file = DATA_DIR / 'internal/j2_6k_unique_isoacc_and_prot_seqs.fa'
    aa = {r.id.split('xxx')[1]: str(r.seq) for r in
          SeqIO.parse(aa_seq_file, format='fasta')}
    df['aa_seq'] = df['clone_acc'].map(aa)
    df['num_aa'] = df['aa_seq'].str.len()
    iso_annot = pd.read_csv(DATA_DIR / 'internal/c_conso_annot_table_man_annot.tsv',
                            sep='\t')
    df['is_novel_isoform'] = df['clone_acc'].map(iso_annot.set_index('unique_acc')['gc_match'] == 0)
    exclude_both_strands_genes = ['FOXD4L3', 'NANOG']
    df = df.loc[~df['gene'].isin(exclude_both_strands_genes), :]
    return df


def load_aligned_aa_seqs(gene_name):
    """
    TODO: this should be deleted, in favour of the genomic_alignment_of_aa_seqs method. Check for uses first

    """
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '../../data')
    path_6k_gtf = os.path.join(data_dir,
                               'hTFIso6K_valid_isoforms/c_6k_unique_acc_aligns.gtf')
    path_6k_fa = os.path.join(data_dir,
                              'hTFIso6K_valid_isoforms/j2_6k_unique_isoacc_and_nt_seqs.fa')
    orf_seqs_6k = isofunc.oc_fasta_to_orf_seq_dict(path_6k_fa)
    gd = isocreate.init_gen_obj(path_6k_gtf, [gene_name])
    gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs_6k)

    gene = gd[gene_name]

    gene_coords = sorted(list(set([pos.coord for pos in gene.poss])))
    tracks = {}
    for orf in sorted(gene.orfs):
        orf_aa = {pos.coord: pos.res.aa for pos in orf.chain if pos.res is not None}
        tracks[orf.name] = ''.join([orf_aa.get(i, '-') for i in gene_coords[::3]])
    if gene.strand == '-':
        tracks = {k: v[::-1] for k, v in tracks.items()}
    return tracks


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


def load_y2h_isoform_data(require_at_least_one_ppi_per_isoform=True,
                          add_missing_data=False,
                          filter_for_valid_clones=True):
    """

    filtered dataset:
        remove partners with no positive tests per tf gene
        remove isoforms with no successful tests
        remove genes with < 2 isoforms

    Args:
        require_at_least_one_ppi_per_isoform (bool, optional): Defaults to True.
        add_missing_data (bool, optional): Defaults to False.
        filter_for_valid_clones (bool, optional): Defaults to True.

    Returns:
        pandas.DataFrame

    """
    y2h = load_isoform_and_paralog_y2h_data(add_missing_data, filter_for_valid_clones)
    ppi = y2h.loc[(y2h['category'] == 'tf_isoform_ppis'),
                  ['ad_clone_acc',
                   'ad_orf_id',
                   'ad_gene_symbol',
                   'db_gene_symbol',
                   'db_orf_id',
                   'score']].copy()
    # at least one positive with an isoform per partner, for each TF gene
    ppi = ppi.loc[ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])
                  ['score']
                  .transform(lambda row: (row == '1').any()),
                  :]
    # at least one succussful test per isoform
    ppi = ppi.loc[ppi.groupby('ad_clone_acc')
                  ['score']
                  .transform(lambda x: (x.isin(['0', '1']).any())),
                  :]
    if require_at_least_one_ppi_per_isoform:
        ppi = ppi.loc[ppi.groupby('ad_clone_acc')
                      ['score']
                      .transform(lambda x: (x == '1').any()),
                      :]
    # at least two isoforms per TF gene
    ppi = ppi.loc[ppi.groupby('ad_gene_symbol')
                  ['ad_clone_acc']
                  .transform(lambda x: x.nunique() >= 2),
                  :]
    # successful tests with at least two isoforms of a TF gene, per partner
    ppi = ppi.loc[ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])['score']
                  .transform(lambda x: x.isin(['0', '1']).sum() >= 2),
                  :]
    # require at least two successfully tested partners per TF gene
    ppi = ppi.loc[ppi['ad_gene_symbol'].map((ppi.loc[ppi['score'].isin(['0', '1'])]
                  .groupby(['ad_gene_symbol', 'db_gene_symbol']).size() >= 2)
                  .groupby('ad_gene_symbol').sum() >= 2),
                  :]
    return ppi


def load_y2h_paralogs_additional_data():
    """Pairs tested in Y2H for the paralogs data, in addition to isoform pairs.

    """
    y2h = load_isoform_and_paralog_y2h_data()
    pairs = load_paralog_pairs()
    y2h_paralog = y2h.loc[y2h['category'].isin(['tf_paralog_ppis',
                                                'paralog_with_PDI',
                                                'non_paralog_control']), :].copy()
    pair_map = defaultdict(set)
    for _i, row in pairs.iterrows():
        a, b = row['tf_gene_a'], row['tf_gene_b']
        pair_map[a].add(b)
        pair_map[b].add(a)

    def find_matching_gene(row):
        matches = pair_map[row['ad_gene_symbol']]
        matches = set(y2h.loc[(y2h['category'] == 'tf_isoform_ppis') &
                              y2h['ad_gene_symbol'].isin(matches) &
                              (y2h['db_gene_symbol'] == row['db_gene_symbol']), 
                              'ad_gene_symbol'].unique())
        if len(matches) == 0:
            return np.nan
        else:
            return '|'.join(matches)

    y2h_paralog['paired_tf_gene'] = y2h_paralog.apply(find_matching_gene, axis=1)

    gte2iso = (y2h.loc[y2h['category'] == 'tf_isoform_ppis', :]
                  .groupby('ad_gene_symbol')['ad_clone_acc']
                  .nunique() >= 2)
    gte2iso = set(gte2iso.index[gte2iso])
    y2h_paralog['at_least_2_isoforms'] = (y2h_paralog['ad_gene_symbol'].isin(gte2iso) &
                                          y2h_paralog['paired_tf_gene'].apply(lambda x: any(g in gte2iso for g in x.split('|')) if pd.notnull(x) else False))

    gte2partner = (y2h.loc[y2h['category'] == 'tf_isoform_ppis', :]
                    .groupby('ad_gene_symbol')['db_gene_symbol']
                    .nunique() >= 2)
    gte2partner = set(gte2partner.index[gte2partner])
    y2h_paralog['at_least_2_partners'] = (y2h_paralog['ad_gene_symbol'].isin(gte2partner) &
                                          y2h_paralog['paired_tf_gene'].apply(lambda x: any(g in gte2partner for g in x.split('|')) if pd.notnull(x) else False))

    return y2h_paralog


def load_isoform_and_paralog_y2h_data(add_missing_data=False, filter_for_valid_clones=True, add_partner_cateogories=False):
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
    if add_partner_cateogories:
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
                       isoforms.loc[:, ['gene', 'clone_acc']]
                               .rename(columns={'gene': 'tf',
                                                'clone_acc': 'unique_acc'}),
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
    print('writing out')
    df4.to_csv('expression_table_tfisoclones.tsv', sep='\t', index=False)
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


def load_paralog_pairs(filter_for_valid_clones=True):
    """Pairs of TF gene paralogs and non-paralogs that were tested in Y2H pairwise tests.

    WARNING: the aa sequence identity is just whatever is first in the file.
    Need to settle on which one to use.

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
    if filter_for_valid_clones:
        valid_clones = load_valid_isoform_clones()
        df = df.loc[df['tf_gene_a'].isin(valid_clones['gene']) &
                    df['tf_gene_b'].isin(valid_clones['gene']), :]
    if (df['tf_gene_a'] == df['tf_gene_b']).any():
        raise ValueError('Same gene twice, should be two different genes')
    return df


def load_isoforms_of_paralogs_pairs(pairs, isoforms):
    """[summary]

    Args:
        pairs ([type]): [description]
        isoforms (bool): [description]

    Returns:
        [type]: [description]
    """
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


def convert_old_id_to_new_id(old_id):
    if len(old_id.split('xxx')) != 3:
        raise UserWarning('Unrecognized old isoform clone ID')
    return old_id.split('xxx')[1]


def read_hmmer3_domtab(filepath):
    """Parser for HMMER 3 hmmscan with the --domtblout option.
    Args:
        filepath (str): location of file
    Returns:
        pandas.DataFrame: see the HMMER User's Guide for a description of the
                          columns.
    """
    columns = ['target name',
               'target accession',
               'tlen',
               'query name',
               'query accession',
               'qlen',
               'E-value',
               'score',
               'bias',
               '#',
               'of',
               'c-Evalue',
               'i-Evalue',
               'full_sequence_score',
               'full_sequence_bias',
               'hmm_coord_from',
               'hmm_coord_to',
               'ali_coord_from',
               'ali_coord_to',
               'env_coord_from',
               'env_coord_to',
               'acc',
               'description of target']
    ncols = len(columns)
    # can't just use pandas.read_csv because spaces are used both as the
    # delimeter and are present in the values of the last columns
    with open(filepath, 'r') as f:
        data = []
        for line in f:
            if line.startswith('#'):
                continue
            s = line.split()
            data.append(s[:ncols - 1] + [' '.join(s[ncols - 1:])])
    df = pd.DataFrame(data=data, columns=columns)
    int_cols = ['tlen', 'qlen', '#', 'of',
                'hmm_coord_from', 'hmm_coord_to',
                'ali_coord_from', 'ali_coord_to',
                'env_coord_from', 'env_coord_to']
    float_cols = ['E-value', 'score', 'bias', 'c-Evalue', 'i-Evalue',
                  'full_sequence_score', 'full_sequence_bias', 'acc']
    df[int_cols] = df[int_cols].astype(int)
    df[float_cols] = df[float_cols].astype(float)
    return df


def _load_pfam_domains(fpath, remove_overlapping=True, cutoff=1E-5):
    pfam = read_hmmer3_domtab(fpath)
    pfam['pfam_ac'] = pfam['target accession'].str.replace(r'\..*', '')
    pfam = pfam.loc[pfam['E-value'] < cutoff, :]
    pfam = _remove_overlapping_domains_from_same_clan(pfam)
    return pfam


def load_pfam_domains_6k():
    filtered_pfam_path = CACHE_DIR / 'pfam_6K.tsv'
    if filtered_pfam_path.exists():
        return pd.read_csv(filtered_pfam_path, sep='\t')
    pfam = _load_pfam_domains(DATA_DIR / 'processed/6K_hmmer_2020-04-16_domtabl.txt')
    pfam['query name'] = pfam['query name'].map(convert_old_id_to_new_id)
    pfam.to_csv(filtered_pfam_path, index=False, sep='\t')
    return pfam


def load_pfam_domains_gencode():
    filtered_pfam_path = CACHE_DIR / 'pfam_gencode.v30.tsv'
    if filtered_pfam_path.exists():
        pfam = pd.read_csv(filtered_pfam_path, sep='\t')
        pfam['query name'] = pfam['query name'].str.slice(8)
        pfam['query name'] = pfam['query name'].apply(lambda x: '|'.join(sorted(x.split('|'))))
        return pfam
    pfam = _load_pfam_domains(DATA_DIR / 'processed/2020-04-23_GC30_6K_domtabl.txt')
    pfam.to_csv(filtered_pfam_path, index=False, sep='\t')
    pfam['query name'] = pfam['query name'].str.slice(8)  # remove 'GC_grp:_'
    pfam['query name'] = pfam['query name'].apply(lambda x: '|'.join(sorted(x.split('|'))))
    return pfam


def _is_overlapping(dom_a, dom_b):
    if dom_b['env_coord_from'] > dom_a['env_coord_to'] or dom_a['env_coord_from'] > dom_b['env_coord_to']:
        return False
    else:
        return True


def load_pfam_clans():
    clans = pd.read_csv(DATA_DIR / 'external/Pfam-A.clans.tsv',
                        sep='\t',
                        header=None)
    if clans.loc[:, 0].duplicated().any():
        raise UserWarning()
    clans = clans.iloc[:, [0, 1]].dropna().set_index(0)[1].to_dict()
    return clans


def _remove_overlapping_domains_from_same_clan(pfam_in):
    """
    Trying to follow the method in the pfam_scan pearl script from EBI.
    """
    pfam = pfam_in.copy()
    clans = load_pfam_clans()
    to_remove = set()
    pfam = pfam.sort_values(['query name', 'E-value'])
    # This is a bit slow
    for iso_id in tqdm.tqdm(pfam['query name'].unique()):
        doms = pfam.loc[(pfam['query name'] == iso_id) &
                        (pfam['pfam_ac'].isin(clans)), :]
        for i, j in itertools.combinations(doms.index, 2):
            if i in to_remove or j in to_remove:
                continue
            dom_a = doms.loc[i, :]
            dom_b = doms.loc[j, :]
            if clans[dom_a['pfam_ac']] == clans[dom_b['pfam_ac']]:
                if _is_overlapping(dom_a, dom_b):
                    to_remove.add(j)
    pfam = pfam.drop(to_remove)
    return pfam


def load_DNA_binding_domains(add_additional_domains=True):
    dbd = pd.read_csv(DATA_DIR / 'internal/a2_final_list_of_dbd_pfam_and_names_ZF_marked.txt',
                      sep='\t')
    clans = load_pfam_clans()
    dbd['clan'] = dbd['pfam'].map(clans)
    if not add_additional_domains:
        return dbd

    # hand adding missing domain see issue #61
    dbd = dbd.append({'dbd': 'BTD',
                      'pfam': 'PF09270',
                      'clan': clans.get('PF09270', np.nan)},
                     ignore_index=True)

    dbd_clans = {'CL0361',  # C2H2-ZF
                'CL0012',  # Histone (mostly DNA binding...)
                'CL0274',  # WRKY-GCM1
                'CL0114',  # HMG-box
                'CL0081',  # MBD-like
                'CL0073',  # P53-like
                'CL0407',  # TATA-Binding Protein like
                'CL0018'}  # bZIP
    for pfam_id in (pfam_id for pfam_id, clan_id in clans.items() if clan_id in dbd_clans):
        if pfam_id not in dbd['pfam'].values:
            dbd = dbd.append({'dbd': pfam_id,
                              'pfam': pfam_id,
                              'clan': clans.get(pfam_id, np.nan)},
                             ignore_index=True)
    return dbd


@functools.lru_cache()
def load_dbd_accessions():
    dbd = load_DNA_binding_domains()
    return set(dbd['pfam'].values).union({'C2H2_ZF_array_' + str(i) for i in range(2, 30)})


def load_annotated_6k_collection():
    path_6k_gtf = DATA_DIR / 'internal/c_6k_unique_acc_aligns.gtf'
    path_6k_fa = DATA_DIR / 'internal/j2_6k_unique_isoacc_and_nt_seqs.fa'
    path_gencode_aa_seq = DATA_DIR / 'external/gencode.v30.pc_translations.fa'
    import pyranges
    # note that pyranges switches the indexing to python 0-indexed, half-open interval
    algn = pyranges.read_gtf(path_6k_gtf).df
    algn = algn.loc[algn['Start'] < algn['End'], :]  # filter out problematic rows
    nt_seq = {r.name: r for r in SeqIO.parse(path_6k_fa, format='fasta')}
    pfam = load_pfam_domains_6k()
    clones = load_valid_isoform_clones()
    algn = algn.loc[algn['transcript_id'].isin(clones['clone_acc'].unique()), :]
    nt_seq = {k: v for k, v in nt_seq.items() if k in clones['clone_acc'].unique()}
    genes = {}
    for gene_name in algn['gene_id'].unique():
        if gene_name == 'PCGF6':  # has a 6nt insertion that doesn't map to reference genome
            continue
        orf_ids = algn.loc[algn['gene_id'] == gene_name, 'transcript_id'].unique()
        missing = [orf for orf in orf_ids if orf not in nt_seq]
        if missing:
            raise ValueError(', '.join(missing) + ' not in ' + path_6k_fa)
        extra = [orf for orf in nt_seq.keys() if orf.split('|')[0] == gene_name and orf not in orf_ids]
        if extra:
            warnings.warn(', '.join(extra) + ' in fasta but not in gtf')
        isoforms = []
        for orf_id in orf_ids:
            exons = []
            columns = ['gene_id', 'transcript_id', 'Chromosome', 'Strand', 'Start', 'End']
            for _i, row in algn.loc[algn['transcript_id'] == orf_id, columns].iterrows():
                exons.append(isolib.Exon(*row.values))
            orf_name = orf_id.split("|")[0] + '-' + orf_id.split("|")[1].split("/")[0]
            isoforms.append(isolib.Cloned_Isoform(clone_name=orf_name,
                                                  exons=exons,
                                                  clone_nt_seq=str(nt_seq[orf_id].seq),
                                                  clone_acc=orf_id))
        genes[gene_name] = isolib.Gene(gene_name, isoforms)
    pfam['gene_name'] = pfam['query name'].apply(lambda x: x.split('|')[0])
    for _i, row in pfam.iterrows():
        gene_name = row['gene_name']
        iso_name = (row['query name'].split('|')[0] + '-' +
                    row['query name'].split('|')[1].split('/')[0])
        if gene_name not in genes or iso_name not in genes[gene_name]:
            continue
        genes[gene_name][iso_name].add_aa_seq_feature(category='Pfam_domain',
                                                      name=row['target name'],
                                                      accession=row['pfam_ac'],
                                                      start=row['env_coord_from'] - 1,
                                                      end=row['env_coord_to'])
    _make_c2h2_zf_arrays(genes)
    _add_dbd_flanks(genes)

    ensembl_proteins = defaultdict(lambda: defaultdict(list))
    for record in SeqIO.parse(path_gencode_aa_seq, 'fasta'):
        ids = record.id.split('|')
        ensembl_protein_id = ids[0].split('.')[0]
        ensembl_transcript_id = ids[1].split('.')[0]
        ensembl_trancript_name = ids[5]
        ensembl_gene_name = ids[6]
        ensembl_proteins[ensembl_gene_name][str(record.seq)].append((ensembl_protein_id,
                                                                     ensembl_transcript_id,
                                                                     ensembl_trancript_name))
    ensembl_proteins['ZNF223'] = ensembl_proteins['AC092072.1']  # HACK fix for gene that has been renamed
    for tf in genes.values():
        if tf.name not in ensembl_proteins:
            raise UserWarning('{} not in gencode AA sequence file'.format(tf.name))
        mapping = _match_clones_to_ensembl_transcripts({iso.name: iso.aa_seq for iso in tf.orfs},
                                                       {tuple(v): k for k, v in ensembl_proteins[tf.name].items()})
        for isoform_name, identifiers in mapping.items():
            tf[isoform_name].ensembl_protein_ids = [ids[0] for ids in identifiers]
            tf[isoform_name].ensembl_transcript_ids = [ids[1] for ids in identifiers]
            tf[isoform_name].ensembl_transcript_names = [ids[2] for ids in identifiers]
    return genes


def _seq_identity(seq_a, seq_b):
    if seq_a == seq_b:
        return 1.
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.gap_score = -1
    alignment = aligner.align(seq_a, seq_b)[0].__str__().split()[1]
    return alignment.count('|') / len(alignment)


MIN_AA_SEQ_ID_TO_MATCH_ENSEMBL_TO_CLONE = 0.98
def _match_clones_to_ensembl_transcripts(clone_seqs, ensembl_seqs, cutoff=MIN_AA_SEQ_ID_TO_MATCH_ENSEMBL_TO_CLONE):
    """
    return mapping of clone to ensembl seqs

    Arguments:
        2 x dict('id': aa_seq)

    Returns:
        dict(in_id -> out_id)
    
    """
    if len(clone_seqs) == 0:
        return ValueError('Empty clone sequences')
    if len(ensembl_seqs) == 0:
        return ValueError('Empty ensembl sequences')

    ids = np.array([[_seq_identity(a, b) for a in clone_seqs.values()] for b in ensembl_seqs.values()])
    best_match_for_clone = ids == ids.max(axis=0)
    best_match_for_ensembl = ids == ids.max(axis=1).reshape(ids.shape[0], 1)
    matches = best_match_for_clone & best_match_for_ensembl & (ids >= cutoff)
    mapping = {a: b for i, a in enumerate(clone_seqs.keys()) 
                    for j, b in enumerate(ensembl_seqs.keys()) 
                    if matches[j, i]}
    return mapping


def _filter_gencode_gtf(out_file_path, genes_subset):
    print('Filtering GENCODE .gtf file. Takes a few minutes but only needs to be done once.')
    import pyranges
    path_gencode_gtf = DATA_DIR / 'external/gencode.v30.annotation.gtf'
    algn = pyranges.read_gtf(path_gencode_gtf, duplicate_attr=True)
    algn = algn[(algn.Feature == 'CDS') &
                (algn.gene_type == 'protein_coding') &
                (algn.transcript_type == 'protein_coding')]
    algn = algn[algn.tag.str.contains('basic')]
    algn = algn[algn.gene_id.str.replace(r'\..*', '', regex=True).isin(genes_subset)]
    algn.to_gtf(out_file_path)


def _filter_gencode_fasta(out_file_path, genes_subset):
    path_gencode_fa = DATA_DIR / 'external/gencode.v30.pc_transcripts.fa'
    out_records = filter(lambda x: x.name.split('|')[1].split('.')[0] in genes_subset,
                         SeqIO.parse(path_gencode_fa, format='fasta'))
    SeqIO.write(out_records, out_file_path, 'fasta')


def load_annotated_gencode_tfs():
    import pyranges  # this import is hidden as it's causing installation problems
    path_filtered_gencode_gtf = CACHE_DIR / 'filtered_CDS_PC_basic_TFs.gencode.v30.annotation.gtf'
    path_filtered_gencode_fa = CACHE_DIR / 'filtered_only_TFs.gencode.v30.pc_transcripts.fa'
    path_gencode_aa_seq = DATA_DIR / 'external/gencode.v30.pc_translations.fa'
    tf_ensembl_gene_ids = set(pd.read_csv(DATA_DIR / 'external/Human_TF_DB_v_1.01.csv')['Ensembl ID'].values)
    if not path_filtered_gencode_gtf.exists():
        _filter_gencode_gtf(path_filtered_gencode_gtf, tf_ensembl_gene_ids)
    if not path_filtered_gencode_fa.exists():
        _filter_gencode_fasta(path_filtered_gencode_fa, tf_ensembl_gene_ids)
    # note that pyranges switches the indexing to python 0-indexed, half-open interval
    algn = pyranges.read_gtf(path_filtered_gencode_gtf, duplicate_attr=True).df
    algn['gene_id'] = algn['gene_id'].str.replace(r'\..*', '', regex=True)
    algn['transcript_id'] = algn['transcript_id'].str.replace(r'\..*', '', regex=True)
    # remove Y-chromosome copies of PAR region, leaving the X copies
    algn = algn.loc[~(algn['tag'].str.contains('PAR')), :]
    nt_seq = {r.name.split('|')[0].split('.')[0]: r for r in SeqIO.parse(path_filtered_gencode_fa, format='fasta')}
    if not algn['transcript_id'].isin(nt_seq).all():
        missing = algn.loc[~algn['transcript_id'].isin(nt_seq), 'transcript_id'].values
        raise ValueError(', '.join(missing) + ' not in ' + path_filtered_gencode_fa)

    tf_gene_names = set(algn['gene_name'].unique())
    valid_transcipts = set(algn['transcript_name'].unique())
    aa_seqs = defaultdict(dict)
    duplicates = {}
    for record in SeqIO.parse(path_gencode_aa_seq, 'fasta'):
        transcript_name, gene_name = record.name.split('|')[5:7]
        if gene_name not in tf_gene_names or transcript_name not in valid_transcipts:
            continue
        aa_seqs[gene_name][transcript_name] = record.seq
    for gene_transcripts in aa_seqs.values():
        for transcript_name, seq in gene_transcripts.items():
            duplicates[transcript_name] = '|'.join(sorted([k for k, v in gene_transcripts.items()
                                                           if v == seq]))
    unique_pc_transcripts = set([v.split('|')[0] for v in duplicates.values()])

    pfam = load_pfam_domains_gencode()
    genes = {}
    for gene_id in tqdm.tqdm(algn['gene_id'].unique()):
        if gene_id == 'ENSG00000163602':
            continue  # RYBP has a sequencing error so don't have full CDS coordinates
        transcript_ids = algn.loc[algn['gene_id'] == gene_id, 'transcript_id'].unique()
        isoforms = []
        for transcript_id in transcript_ids:
            transcript_name, gene_name = nt_seq[transcript_id].name.split('|')[4:6]
            if transcript_name not in unique_pc_transcripts:
                continue
            cds = _extract_region(nt_seq[transcript_id], 'CDS', raise_error=True)
            utr5 = _extract_region(nt_seq[transcript_id], 'UTR5', raise_error=False)
            utr3 = _extract_region(nt_seq[transcript_id], 'UTR3', raise_error=False)
            exons = []
            columns = ['gene_id', 'transcript_id', 'Chromosome', 'Strand', 'Start', 'End']
            for _i, row in algn.loc[algn['transcript_id'] == transcript_id, columns].iterrows():
                exons.append(isolib.Exon(*row.values))
            isoforms.append(isolib.Isoform(transcript_name,
                                           exons,
                                           CDS_nt_seq=cds,
                                           aa_seq=str(aa_seqs[gene_name][transcript_name]),
                                           UTR_5prime_nt_seq=utr5,
                                           UTR_3prime_nt_seq=utr3,
                                           ensembl_transcript_id=transcript_id))
        genes[gene_name] = isolib.Gene(gene_name, isoforms)
    pfam['gene_name'] = pfam['query name'].apply(lambda x: '-'.join(x.split('|')[0].split('-')[:-1]))
    for _i, row in pfam.iterrows():
        gene_name = row['gene_name']
        transcript_name = row['query name'].split('|')[0]
        if gene_name not in genes or transcript_name not in genes[gene_name]:
            continue
        genes[gene_name][transcript_name].add_aa_seq_feature(category='Pfam_domain',
                                                             name=row['target name'],
                                                             accession=row['pfam_ac'],
                                                             start=row['env_coord_from'] - 1,
                                                             end=row['env_coord_to'])
    _make_c2h2_zf_arrays(genes)
    _add_dbd_flanks(genes)
    return genes


def _extract_region(rec, region, raise_error=True):
    """
    Get subset of nucleotide sequence that correpsonds to different
    regions (e.g. CDS, 3'UTR) when bounds of regions are specified
    in desciption in fasta file, like this:
    >HES4-204|HES4|667|UTR5:1-8|CDS:9-578|UTR3:579-667|
    """
    bounds = [x for x in rec.name.split('|') if x.startswith(region + ':')]
    if len(bounds) == 0:
        if raise_error:
            raise UserWarning('Missing {} coordinates\n'.format(region) + rec.name)
        else:
            return None
    if len(bounds) > 1:
        raise UserWarning('More than one set of {} coordinates\n'.format(region) + rec.name)
    start, stop = bounds[0][len(region) + 1:].split('-')
    return str(rec.seq)[int(start) - 1:int(stop)]


def _make_c2h2_zf_arrays(tfs, MAX_NUM_AA_C2H2_ZF_SEPERATION=10):
    clans = load_pfam_clans()
    C2H2_ZF_PFAM_CLAN_AC = 'CL0361'
    c2h2_zf_pfam_ids = {k for k, v in clans.items() if v == C2H2_ZF_PFAM_CLAN_AC}
    for gene in tfs.values():
        for orf in gene.orfs:
            zfs = [d for d in orf.aa_seq_features if d.accession in c2h2_zf_pfam_ids]
            if len(zfs) < 2:
                continue
            zfs = sorted(zfs, key=lambda x: x.start)
            array = [zfs[0]]
            for zf in zfs[1:]:
                if zf.start <= (array[-1].end - 1) + MAX_NUM_AA_C2H2_ZF_SEPERATION:
                    array.append(zf)
                else:
                    if len(array) > 1:
                        orf.add_aa_seq_feature('ZF_array',
                                                'C2H2_ZF_array_' + str(len(array)),
                                                'C2H2_ZF_array_' + str(len(array)),
                                                array[0].start,
                                                array[-1].end)
                        for dom in array:
                            orf.remove_aa_seq_feature(dom.accession, dom.start, dom.end)
                    array = [zf]
            if len(array) > 1:
                orf.add_aa_seq_feature('ZF_array',
                                        'C2H2_ZF_array_' + str(len(array)),
                                        'C2H2_ZF_array_' + str(len(array)),
                                        array[0].start,
                                        array[-1].end)
                for dom in array:
                    orf.remove_aa_seq_feature(dom.accession, dom.start, dom.end)


def _add_dbd_flanks(genes):
    for gene in genes.values():
        for isoform in gene.orfs:
            for dbd in isoform.dna_binding_domains:
                start_n = max(0, dbd.start - 15)
                end_n = dbd.start
                start_c = dbd.end
                end_c = min(len(isoform.aa_seq), dbd.end + 15)
                if dbd.start > 0:
                    isoform.add_aa_seq_feature('DBD_flank',
                                               dbd.accession + '_flank_N',
                                               'N_DBD_flank',
                                               start_n,
                                               end_n)
                if dbd.end < len(isoform.aa_seq):
                    isoform.add_aa_seq_feature('DBD_flank',
                                               dbd.accession + '_flank_C',
                                               'C_DBD_flank',
                                               start_c,
                                               end_c)