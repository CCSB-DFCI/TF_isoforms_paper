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
from isomodules import isofunc
from isomodules import isowrangle

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master_take_3/data_input/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)

# load gene list, with updated manual corrections
genes = []
for line in open('b_orig_genename_in_gencode30_man_annot.txt').readlines()[1:]:
    wds = line.rstrip().split('\t')
    gene = wds[5]
    if gene == '-':
        # skip over DUX genes, which are achromosomal
        continue
    genes.append(gene)


# make gene objects for gloria + sachi + lambert tfs
d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes, verbose=True)
d_gc_w_seqs = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs, verbose=True)
gd = d_gc_w_seqs

# %%

# merge same-protein-seq. orfs, write out info. to table
with open('d_gencode30_uniq_prot_and_annot_sachi_gloria_lambert.tsv', 'w') as ofile:
    ofile.write('orig_gene\tnew_gene\tin_sachi1332\tin_gloria1579\tin_lambert1639\tensg\t'
                'isoname\tenst\tensp\taa_seq\tsame_prot_seq_name\t'
                'is_repr\tappris\ttranscript_level\ttranscript_support_level\t'
                'isoname_derived_rank\tenst_derived_idx\tis_basic\t'
                'is_cds_start_NF\tis_cds_end_NF\tappris_rank\n')
    for line in open('b_orig_genename_in_gencode30_man_annot.txt').readlines()[1:]:
        wds = line.rstrip().split('\t')
        orig_gene, presence, sachi, gloria, lambert, gname = wds
        prefix = [orig_gene, gname, sachi, gloria, lambert]
        if gname not in gd:
            print gname
            continue
        gen_obj = gd[gname]
        grps = isowrangle.get_same_prot_seq_grps_of_orfs_from_gene_(gen_obj)
        for grp in grps:
            orf = grp.repr_orf
            odata = (prefix
                    + [orf.gene.ensg, orf.name, orf.enst, orf.ensp, orf.tr, grp, '1']
                    + orf.flags
                    + [orf.is_basic * 1, orf.cds_start_nf * 1, orf.cds_end_nf * 1, 1])
            ofile.write('\t'.join(map(str, odata)) + '\n')
            tmp = []  # store the out data for all other orfs, sort by appris/evidence to get ranking of ensts
            for orf in grp.other_orfs:
                odata = (prefix
                        + [orf.gene.ensg, orf.name, orf.enst, orf.ensp, orf.tr, grp, '0']
                        + orf.flags
                        + [orf.is_basic * 1, orf.cds_start_nf * 1, orf.cds_end_nf * 1])
                tmp.append(odata)
            tmp = sorted(tmp, key = lambda x: (x[-8], x[-7], x[-6], x[-5], x[-4]))
            for i, odata in enumerate(tmp):
                # i + 2 is the ranking, starting at 2 b/c # 1 is appris (already written out above)
                ofile.write('\t'.join(map(str, odata + [i+2])) + '\n')
