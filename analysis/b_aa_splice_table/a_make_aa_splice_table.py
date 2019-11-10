#!/usr/bin/env python

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
from isomodules import isocreatealign

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master/data/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf.toy')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa.toy')
path_gc_domain_toy = os.path.join(data_dir, 'gencode.v30.domains.toy')
path_gc_aln = os.path.join(data_dir, 'clustal_align/b_align_results_from_sachi/2019-05-20_alignment_table_for_Gloria.txt')
path_gc_aln_toy = os.path.join(data_dir, 'clustal_align/c_make_toy_align_NFYA_PAX5/clustal_pairwise_alns_gc_only.tsv.toy')
path_gc_aln_blocks = os.path.join(data_dir, 'clustal_align/b_align_results_from_sachi/2019-05-20_pairwise_orf_alignment_6KandGencode.txt')
path_gc_aln_blocks_toy = os.path.join(data_dir, 'clustal_align/c_make_toy_align_NFYA_PAX5/clustal_pairwise_block_cats_gc_only.tsv.toy')
path_6k_gtf_toy = os.path.join(data_dir, 'hTFIso6K_gtfhg38.star.gtf.toy')
path_6k_fa_toy = os.path.join(data_dir, 'hTFIso6K_sequences_w_orfid.fasta.toy')
path_hg38_fa = os.path.join(data_dir, 'GRCh38.primary_assembly.genome.fa')
path_tf_list_gloria = os.path.join(data_dir, 'human_tfs_annotation_full_and_stringent.tsv')
path_tf_list_sachi = os.path.join(data_dir, '1332_highconfidence_TFs_with_ensembl_updated_2017-08-08.txt')
path_tf_list_lambert = os.path.join(data_dir, 'tf_annotations/Lambert_TFs_v_1.01.txt')




# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
tf_sachi, tf_gloria, tf_lambert = isofunc.load_tf_list(path_tf_list_sachi,
                                     path_tf_list_gloria, path_tf_list_lambert)
hg38_dict = isofunc.load_hg38(path_hg38_fa)

# %%


reload(isocreate)
reload(isoclass)
reload(isofunc)
reload(isofeature)
reload(isomap)
reload(isoimage)
reload(isogroup)
reload(isoalign)
reload(isowrangle)
reload(isoalign)

genes = ['NFYA', 'PAX5']

d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
d_gc_w_seqs = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs)
d_gc_w_seqs_w_juncs = isocreate.create_and_link_junction_related_obj(d_gc_w_seqs, hg38_dict)
gd = d_gc_w_seqs_w_juncs  # gen_obj dict

# Now, create an aln object, which can bring out a "full" display.
with open('a_full_display_toy.txt', 'w') as ofile:
    for symbol, gene in gd.items():
        for pair in gene.ref_alt_pairs:
            aln_grps = isocreatealign.create_and_map_splice_based_align_obj([pair])
            alnf = aln_grps[0].alnf
            ofile.write(alnf.full)
