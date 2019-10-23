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

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master_take_3/data_input/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf.toy')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa.toy')
path_tf_names = os.path.join(data_dir, 'tf_genenames/b_orig_genename_in_gencode30_man_annot.txt')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
genes = isofunc.load_tf_genenames(path_tf_names)  # ~1700 genenames

# %%

reload(isocreate)
reload(isocreatealign)
reload(isoclass)
reload(isofunc)
reload(isofeature)
reload(isomap)
reload(isoimage)
reload(isogroup)
reload(isoalign)
reload(isowrangle)

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
    all_grps.append(grps)

# %%

import itertools

def infer_frame_of_non_overlapping_seq(frame_blocks):
    """Infer the frame of sequence regions of alt. iso. that are not in ref.
       Algorithm - look to the upstream-most overlapping segment. If doesn't
       exist (b/c is at N-terminus), look to downstream.
    """
    inferred_frame_blocks = []
    infer_status_chain = []
    for i, block in enumerate(frame_blocks):
        blen = len(block)
        if i == 0 and '*' in block:
            inferred_frm = frame_blocks[1][0]
            infer_status = True
        elif '*' in block:
            inferred_frm = frame_blocks[i - 1][0]
            infer_status = True
        else:
            # don't need to infer frame b/c propogated from ref.
            inferred_frm = frame_blocks[i][0]
            infer_status = False
        infer_status = infer_status * 1  # convert to no. for writeout
        inferred_frame_blocks.append(inferred_frm * blen)
        infer_status_chain.extend([infer_status] * blen)
    frm_list = list(''.join(inferred_frame_blocks))
    infer_list = list(''.join(map(str, infer_status_chain)))
    return infer_list, frm_list

def create_list_of_frms_and_ifrms(inferred_blocks):
    ifrm_str = ''.join(inferred_blocks)
    return ifrm_str

reload(isofeature)

for grps in all_grps:
    for grp in grps:
        orf1 = grp.repr_orf
        orf2 = grp.other_orf
        print grp.alnf.full
        print orf1, orf2
        fstr = ''.join([res.rfrm for res in orf2.res_chain])
        frame_blocks = [''.join(g) for _, g in itertools.groupby(fstr)]
        infer_list, frm_list = infer_frame_of_non_overlapping_seq(frame_blocks)
        featf = isofeature.FeatureFull(grp, orf2, 'frms')
        print orf2.frms
        for i, frm in enumerate(frm_list):
            infer = infer_list[i]
            res = orf2.res_chain[i]
            if not infer:
                frmr = isofeature.FeatureResidue(featf, frm, res, 'direct', 'frms')
            else:
                frmr = isofeature.FeatureResidue(featf, frm, res, 'inferred', 'frms')
            featf.chain.append(frmr)




    break
