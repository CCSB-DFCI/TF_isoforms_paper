#!/usr/bin/env python
#
# title           :driver.py
# description     :Work out protein/genome-based alignment of protein isoforms.
# author          :Gloria Sheynkman
# date            :June 3rd, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================

import os
from isomodules import isocreate
from isomodules import isocreatealign
from isomodules import isoclass
from isomodules import isofunc
from isomodules import isofeature
from isomodules import isomap
from isomodules import isoimage
from isomodules import isogroup
from isomodules import isoalign
from isomodules import isowrangle
from isomodules import isocreatefeat

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master/data/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf.toy')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa.toy')
path_tf_names = os.path.join(data_dir, 'tf_genenames/b_orig_genename_in_gencode30_man_annot.txt')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
genes = isofunc.load_tf_genenames(path_tf_names)  # ~1700 genenames

genes = ['NFYA', 'PAX5']

gd = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs)

for symbol, gene in gd.items():
    grps = isowrangle.get_same_prot_seq_grps_of_orfs_from_gene_(gene)
    repr_orfs = set(grp.repr_orf for grp in grps)  # uniq-prot-seq orfs
    gene.redundant_seq_orfs = gene.orfs - repr_orfs
    gene.orfs = repr_orfs

all_grps = []
for symbol, gene in gd.items():
    grps = isocreatealign.create_and_map_splice_based_align_obj(gene.ref_alt_pairs)
    for grp in grps:
        isocreatefeat.create_and_map_frame_objects(grp)
    all_grps.append(grps)

# %%

for grps in all_grps:
    for grp in grps:
        orf2 = grp.other_orf
        orf2.current_grp = grp
        orf2.current_feat = grp.frmf
