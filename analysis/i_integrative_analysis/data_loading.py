import os
import sys
from pathlib import Path
import itertools
import warnings

import numpy as np
import pandas as pd
from Bio import SeqIO
import pyranges
import tqdm

sys.path.append('../..')
import isolib

DATA_DIR = Path(__file__).resolve().parents[2] / 'data'

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
                'ad_gene_symbol',
                'db_gene_symbol',
                'score']].copy()
    ppi = ppi.loc[ppi.groupby(['ad_gene_symbol', 'db_gene_symbol'])
                    ['score']
                    .transform(lambda row: (row == '1').any()),
                :]
    ppi = ppi.loc[ppi.groupby('ad_clone_acc')
                    ['score']
                    .transform(lambda x: (x.isin(['0', '1']).any())),
                :]
    if require_at_least_one_ppi_per_isoform:
        ppi = ppi.loc[ppi.groupby('ad_clone_acc')['score'].transform(lambda x: (x == '1').any()),
                    :]
    ppi = ppi.loc[ppi.groupby('ad_gene_symbol')
                    ['ad_clone_acc']
                    .transform(lambda x: x.nunique() >= 2),
                :]
    return ppi


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


def load_pfam_domains_for_6k(remove_overlapping=True, cutoff=1E-5):
    fpath = DATA_DIR / 'processed/6K_hmmer_2020-04-16_domtabl.txt'
    pfam = read_hmmer3_domtab(fpath)
    pfam['pfam_ac'] = pfam['target accession'].str.replace(r'\..*', '')
    pfam['query name'] = pfam['query name'].map(convert_old_id_to_new_id)
    pfam = pfam.loc[pfam['E-value'] < cutoff, :]
    pfam = _remove_overlapping_domains_from_same_clan(pfam)
    return pfam


def _is_overlapping(dom_a, dom_b):
    if dom_b['env_coord_from'] > dom_a['env_coord_to'] or dom_a['env_coord_from'] > dom_b['env_coord_to']:
        return False
    else:
        return True


def _remove_overlapping_domains_from_same_clan(pfam_in):
    """
    Trying to follow the method in the pfam_scan pearl script from EBI.
    """
    pfam = pfam_in.copy()
    clans = pd.read_csv(DATA_DIR / 'external/Pfam-A.clans.tsv',
                        sep='\t',
                        header=None)
    if clans.loc[:, 0].duplicated().any():
        raise UserWarning()
    clans = clans.iloc[:, [0, 1]].dropna().set_index(0)[1].to_dict()
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


def load_DNA_binding_domains():
    return pd.read_csv(DATA_DIR / 'internal/a2_final_list_of_dbd_pfam_and_names_ZF_marked.txt',
                       sep='\t')


def load_annotated_6k_collection():
    path_6k_gtf = DATA_DIR / 'internal/c_6k_unique_acc_aligns.gtf'
    path_6k_fa = DATA_DIR / 'internal/j2_6k_unique_isoacc_and_nt_seqs.fa'
    # note that pyranges switches the indexing to python 0-indexed, half-open interval
    algn = pyranges.read_gtf(path_6k_gtf).df
    algn = algn.loc[algn['Start'] < algn['End'], :]  # filter out problematic rows
    nt_seq = {r.name: r for r in SeqIO.parse(path_6k_fa, format='fasta')}
    pfam = load_pfam_domains_for_6k()
    dbd = load_DNA_binding_domains()
    clones = load_valid_isoform_clones()
    algn = algn.loc[algn['transcript_id'].isin(clones['clone_acc'].unique()), :]
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
            isoforms.append(isolib.ORF(orf_id, exons, str(nt_seq[orf_id].seq)))
        genes[gene_name] = isolib.Gene(gene_name, isoforms)
    pfam['gene_name'] = pfam['query name'].apply(lambda x: x.split('|')[0])
    for _i, row in pfam.iterrows():
        gene_name = row['gene_name']
        iso_name = row['query name']
        if gene_name not in genes or iso_name not in genes[gene_name]:
            continue
        genes[gene_name][iso_name].add_aa_seq_feature(category='Pfam_domain', 
                                                                    name=row['target name'],
                                                                    accession=row['pfam_ac'],
                                                                    start=row['env_coord_from'] - 1,
                                                                    end=row['env_coord_to'])
    return genes
