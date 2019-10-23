#!/usr/bin/env python
#
# title           :driver.py
# description     :Instantiate gen_obj for 1600 TFs. Write-out same-prot. seq.
#                  info to a table.
# author          :Gloria Sheynkman
# date            :May 27th, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================

import os
from isomodules import isocreate
from isomodules import isoclass
from isomodules import isofunc
from isomodules import isogroup
from isomodules import isowrangle

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master_take_3/data_input/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa')
path_tf_list_gloria = os.path.join(data_dir, 'human_tfs_annotation_full_and_stringent.tsv')
path_tf_list_sachi = os.path.join(data_dir, '1332_highconfidence_TFs_with_ensembl_updated_2017-08-08.txt')
path_tf_list_lambert = os.path.join(data_dir, 'Lambert_TFs_v_1.01.txt')


# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
tf_sachi, tf_gloria, tf_lambert = isofunc.load_tf_list(path_tf_list_sachi,
                                     path_tf_list_gloria, path_tf_list_lambert)
# 1332 in Sachi list, 1579 in Gloria list, 1639 in lambert list

# make gene objects
genes = list(set(tf_sachi + tf_gloria + tf_lambert))

d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes, verbose=True)

with open('a_orig_genename_in_gencode30.tsv', 'w') as ofile:
    ofile.write('gname\tpresence\tsachi\tgloria\tlambert\n')
    for gene in genes:
        if gene not in d_gc:
            status = 'absent'
        else:
            status = 'present'
        sachi = '1' if gene in tf_sachi else '0'
        gloria = '1' if gene in tf_gloria else '0'
        lambert = '1' if gene in tf_lambert else '0'
        ofile.write('\t'.join([gene, status, sachi, gloria, lambert]) + '\n')
