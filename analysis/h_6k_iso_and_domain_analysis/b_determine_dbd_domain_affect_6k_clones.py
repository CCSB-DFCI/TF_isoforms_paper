#!/usr/bin/env python
# Map regulatory domains to 6K isoforms. Compare domains and write-out table.
# ==============================================================================

import os
import re
from isomodules import isocreate
from isomodules import isocreatealign
from isomodules import isocreatefeat
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
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa')
path_6k_gtf = os.path.join(data_dir, 'hTFIso6K_gtfhg38.star.gtf')
path_6k_fa = os.path.join(data_dir, 'hTFIso6K_sequences_w_orfid.fasta')
path_pfam_gc = os.path.join(data_dir, 'pfam/2019-06-14_HMMER_mappings_to_all_GENCODEv30_LambertTFs_isoforms.txt')
path_pfam_6k = os.path.join(data_dir, 'pfam/2019-06-14_HMMER_mappings_to_all_6Kisoforms.txt')
path_pfam_names = os.path.join(data_dir, 'pfam/pfam_a_names.tsv')
path_appris = os.path.join(data_dir, 'principle_isoforms/determined_principle_isoform_for_1700tfs.tsv')
path_isoacc = os.path.join(data_dir, 'hTFIso6K_assigned_isoacc/1_assigned_human_readable_isoacc.tsv')


# load data - isoform sequences and domains
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
orf_seqs_6k = isofunc.oc_fasta_to_orf_seq_dict(path_6k_fa)
appris_orfs = isofunc.load_appris_principle_isonames(path_appris)
domains = isofunc.load_domain_mappings(path_pfam_gc, path_pfam_names)
isoacc_map = isofunc.load_6k_isoacc_map(path_isoacc)


# genes = ['NFYA', 'PAX5']  # toy example
genes = []  # when running full
cl_gd = isocreate.init_gen_obj(path_6k_gtf, gene_list=genes, assign_isoacc=True, isoacc_map=isoacc_map)
cl_gd = isocreate.create_and_link_seq_related_obj(cl_gd, orf_seqs_6k)

genes = cl_gd.keys()  # ~400 tfs

# make gene obj., but only for appris principle orfs and only for 6k tfs
gc_gd = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes, iso_list=appris_orfs)
gc_gd = isocreate.create_and_link_seq_related_obj(gc_gd, orf_seqs)

# %%

def domain_info_is_a_DBD(domain_info):
    if domain_info[2] == 'DBD':
        return True

def all_6k_orfs_of_this_gene(gd, gene):
    return gd[gene].orfs

def get_stats_for_domain_that_is_entirely_subst(dom_alnrs):
    alnpb = dom_alnrs[0].alnpb
    aa1 = alnpb.aa1.replace('-', '')
    aa2 = alnpb.aa2.replace('-', '')
    len_aa1 = len(aa1)
    len_aa2 = len(aa2)
    len_diff = len_aa1 - len_aa2
    perc_diff = float(len_diff) / len_aa1
    return [aa1, aa2, len_aa1, len_aa2, len_diff, perc_diff]

def compute_the_perc_aa_affected(prot_cat_str, dom_len):
    """From the protein category (from alnpb) string, calc. % AA preserved.
       Calculating two ways - using num. AA in orig. denominator
                              using num. splice-aligned AA (alnr) in denom.
    """
    num_matches = prot_cat_str.count('M')
    perc_denom_dom = num_matches / float(dom_len) * 100
    return num_matches, perc_denom_dom


with open('b_6k_isoforms_affect_of_mapped_dbd_domains.tsv', 'w') as ofile:
    ofile.write('gene\tisoacc\tisoacc2\tdom_name\tdom_cat\tdom_start_idx\tdom_end_idx\tprot_cmp_block_string\tprot_cmp_string\tdom_len\tnum_match\t%AA_preserved_denom_dom_len\n')
    for symbol, gene in gc_gd.items():
        appris_orf = gene.appris_orf
        for dom_info in domains[appris_orf.ensp]:
            if domain_info_is_a_DBD(dom_info):
                anchor_domain = isocreatefeat.create_and_map_domain(appris_orf, dom_info)
                for clone_orf in all_6k_orfs_of_this_gene(cl_gd, symbol):
                    # align 6k orf to gencode appris orf (create an aln grp)
                    grp = isocreatealign.create_and_map_splice_based_align_obj([[appris_orf, clone_orf]])[0]
                    grp.enter()
                    alnf = grp.alnf
                    chain = alnf.chain
                    # chain of alnrs corresponding to domain
                    dom_alnrs = alnf.get_subset_aln_chain(anchor_domain.first.res.aln, anchor_domain.last.res.aln)
                    cat_str = ''.join([alnr.alnpb.cat for alnr in dom_alnrs])
                    cat_block = ''.join([x[0] for x in re.split('(M+)', cat_str) if x])
                    dom_len = len(anchor_domain.chain)
                    # for now, deletions and subst. treated as 'removals' or 'R'
                    cat_str = cat_str.replace('D', 'R').replace('S', 'R')
                    cat_block = cat_block.replace('D', 'R').replace('S', 'R')
                    num_matches, perc_aa_den_dom = compute_the_perc_aa_affected(cat_str, dom_len)
                    start_idx, end_idx = dom_alnrs[0].res2.idx, dom_alnrs[-1].res2.idx
                    start_idx = 'NA' if not start_idx else start_idx
                    end_idx = 'NA' if not end_idx else end_idx
                    odata = [symbol, clone_orf.name, clone_orf.name2, anchor_domain.acc, anchor_domain.cat, start_idx, end_idx, cat_block, cat_str, dom_len, num_matches, perc_aa_den_dom]
                    ofile.write('\t'.join(map(str, odata)) + '\n')
