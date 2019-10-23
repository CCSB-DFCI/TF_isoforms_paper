#!/usr/bin/env python
#
# title           :driver.py
# description     :Starting template to iso-objs and do analysis.
# author          :Gloria Sheynkman
# date            :May 6th, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================

import os
from isomodules import isocreate
from isomodules import isoclass
from isomodules import isofunc
from isomodules import isofeature
from isomodules import isomap
from isomodules import isoimage
from isomodules import isogroup
from isomodules import isoalign
from isomodules import isowrangle
from isomodules import isocompare

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master_take_3/data_input/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf.toy')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa.toy')
path_gc_aln = os.path.join(data_dir, 'clustal_align/b_align_results_from_sachi/2019-05-20_alignment_table_for_Gloria.txt')
path_gc_aln_toy = os.path.join(data_dir, 'clustal_align/c_make_toy_align_NFYA_PAX5/clustal_pairwise_alns_gc_only.tsv.toy')
path_gc_aln_blocks = os.path.join(data_dir, 'clustal_align/b_align_results_from_sachi/2019-05-20_pairwise_orf_alignment_6KandGencode.txt')
path_gc_aln_blocks_toy = os.path.join(data_dir, 'clustal_align/c_make_toy_align_NFYA_PAX5/clustal_pairwise_block_cats_gc_only.tsv.toy')


# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
aln_data = isofunc.load_clustal_alignments(path_gc_aln)
aln_blocks = isofunc.load_alignment_block_calls(aln_data, path_gc_aln_blocks)


reload(isocreate)
reload(isoclass)
reload(isofunc)
reload(isofeature)
reload(isomap)
reload(isoimage)
reload(isogroup)
reload(isoalign)
reload(isowrangle)

genes = ['NFYA', 'PAX5']

d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
d_gc_w_seqs = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs)
gd = d_gc_w_seqs  # gen_obj dict


# clustal alignments into OOP objects
grps = isocreate.create_and_map_alignments(gd, aln_data, aln_blocks)
for grp in grps:
    print grp.alnf
    grp.enter()
    print grp.anchor_orf.current_grp
    alnf = grp.anchor_orf.aln
    print ''.join([alnr.alnb.cat[1] if alnr not in [alnr.alnb.first, alnr.alnb.last] else '|' for alnr in alnf.chain])
    print ''.join([alnr.alnsb.cat[1] if alnr not in [alnr.alnsb.first, alnr.alnsb.last] else '|' for alnr in alnf.chain])
    print ''.join([str(alnr.alnsb.cds1.ord)[0] for alnr in alnf.chain])
    print ''.join([alnr.res1.aa for alnr in alnf.chain])
    print ''.join([alnr.cat for alnr in alnf.chain])
    print ''.join([alnr.res2.aa for alnr in alnf.chain])
    print ''.join([str(alnr.alnsb.cds2.ord)[0] for alnr in alnf.chain])
    print '\n'
    grp.exit()
