import os
import sys

import numpy as np
from matplotlib import pyplot as plt

from ccsblib import ccsbplotlib as cplt

sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)), '../..'))
from isomodules import isocreate, isoimage, isofunc


def isoform_display_name(s):
       """Convert clone accession ID to display friendly format"""
       return s.split('|')[0] + '-' + s.split('|')[1].split('/')[0]


def strikethrough(s):
    return ''.join([c + '\u0336' for c in s])


def isoform_box_and_line_drawing(gene_name, clone_accs, ax=None):
    if ax is None:
        ax = plt.gca()
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '../../data')
    path_6k_gtf = os.path.join(data_dir,
                               'hTFIso6K_valid_isoforms/c_6k_unique_acc_aligns.gtf')
    path_6k_fa = os.path.join(data_dir,
                              'hTFIso6K_valid_isoforms/j2_6k_unique_isoacc_and_nt_seqs.fa')
    orf_seqs_6k = isofunc.oc_fasta_to_orf_seq_dict(path_6k_fa)
    gd = isocreate.init_gen_obj(path_6k_gtf, [gene_name])
    gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs_6k)
    isoimage.render_iso_image([orf for orf in list(gd[gene_name].orfs) if orf.name in clone_accs],
                              ax=ax)


def y2h_ppi_per_tf_gene_plot(gene_name, 
                             data,
                             ax=None, 
                             min_n_isoforms=1,
                             min_n_partners=1):
    tf = data.loc[(data['category'] == 'tf_isoform_ppis') &
                (data['ad_gene_symbol'] == gene_name),
                ['ad_clone_acc', 'db_gene_symbol', 'score']].copy()
    tf['score'] = tf['score'].map({'1': True,
                               '0': False,
                               'AA': np.nan,
                               'NC': np.nan})
    tf['ad_clone_acc'] = tf['ad_clone_acc'].apply(isoform_display_name)
    tf = tf.pivot(index='ad_clone_acc',
                  columns='db_gene_symbol',
                  values='score')
    if ax is None:
        ax = plt.gca()
    if tf.shape[0] < min_n_isoforms or tf.shape[1] < min_n_partners:
        ax.set_axis_off()
        ax.text(0.5, 0.5,
                'No PPI data available',
                ha='center', va='center',
                fontsize=30,
                fontweight='bold',
                color='grey')
        return
    cplt.binary_profile_matrix(tf, ax=ax, column_label_rotation=90)
    ax.set_yticklabels([strikethrough(name) if all_na else name
                        for name, all_na in tf.isnull().all(axis=1).iteritems()])


def y1h_pdi_per_tf_gene_plot(gene_name, 
                             data,
                             ax=None,
                             min_n_isoforms=1,
                             min_n_partners=1):
    tf = data.loc[data['tf'] == gene_name,
                  data.columns[1:]].copy().set_index('unique_acc')
    tf.index = tf.index.map(isoform_display_name)
    tf = tf.loc[:, tf.any(axis=0)]
    if ax is None:
        ax = plt.gca()
    if tf.shape[0] < min_n_isoforms or tf.shape[1] < min_n_partners:
        ax.set_axis_off()
        ax.text(0.5, 0.5,
                'No PDI data available',
                ha='center', va='center',
                fontsize=30,
                fontweight='bold',
                color='grey')
        return
    cplt.binary_profile_matrix(tf, ax=ax, column_label_rotation=90)
    ax.set_yticklabels([strikethrough(name) if all_na else name
                        for name, all_na in tf.isnull().all(axis=1).iteritems()])

def m1h_activation_per_tf_gene_plot(tf_gene_name, data, ax=None):
    if ax is None:
        ax = plt.gca()
    rep_columns = [c for c in data.columns if c.startswith('M1H_rep')]
    is_all_na = data[rep_columns].isnull().groupby(data['gene']).all().all(axis=1)[tf_gene_name]
    if tf_gene_name not in data['gene'].values or is_all_na:
        ax.set_axis_off()
        ax.text(0.5, 0.5,
                'No activation data available',
                ha='center', va='center',
                fontsize=30,
                fontweight='bold',
                color='grey')
        return
    clones = [isoform_display_name(acc) for acc in data.loc[data['gene'] == tf_gene_name, 'clone_acc'].values
              for __ in range(len(rep_columns))]
    values = data.loc[data['gene'] == tf_gene_name, rep_columns].values.flatten()
    ax.scatter(
        y=clones,
        x=values,
        alpha=0.5
        )
    n_reps = len(rep_columns)
    ax.set_yticks(clones[::n_reps])  # needed to avoid truncating clones with missing data
    ax.set_yticklabels([c if not np.isnan(v) else strikethrough(c)
                        for c, v in zip(clones[::n_reps], values[::n_reps])])
    ax.set_xlim(-3, 12)
    ax.set_xlabel('Log2 M1H readout')
    ax.axvline(0, linestyle='-', color='grey')
    ax.axvline(-1, linestyle='--', color='grey')
    ax.axvline(1, linestyle='--', color='grey')
    ax.invert_yaxis()