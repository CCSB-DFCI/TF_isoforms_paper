#!/usr/bin/env python
# Map regulatory domains to 6K isoforms. Compare domains and write-out table.
# ==============================================================================

import os
import re
import pandas as pd
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
path_gc_gtf = os.path.join(data_dir, 'gencode_corr_genename_for_uniprot_dbd_mapping/gencode.v30.annotation_corr_ZNF223.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode_corr_genename_for_uniprot_dbd_mapping/gencode.v30.pc_transcripts_corr_ZNF223.fa')
path_6k_gtf = os.path.join(data_dir, 'hTFIso6K_valid_isoforms/c_6k_unique_acc_aligns.gtf')
path_6k_fa = os.path.join(data_dir, 'hTFIso6K_valid_isoforms/j2_6k_unique_isoacc_and_nt_seqs.fa')
path_uniprot_dbd_to_best_iso = '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/project_tf_isoforms/analysis/8_map_domains/4_map_uniprot_domain_to_gc_6k/d_best_iso_per_domain_mapping.tsv'
path_dbd_names = '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/project_tf_isoforms/analysis/8_map_domains/4_map_uniprot_domain_to_gc_6k/0_DBD_curation_table_from_sachi/a_final_list_of_dbd_pfam_and_names.tsv'


reload(isocreate)
reload(isofunc)

## load data

# only DBD domains
domains = isofunc.load_uniprot_dbd_mappings(path_uniprot_dbd_to_best_iso, path_dbd_names)

# %%

# %%
non_zf_domains = isofunc.load_uniprot_dbd_mappings(path_uniprot_dbd_to_best_iso, path_dbd_names, non_zf_only=True)

orf_seqs_gc = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
orf_seqs_6k = isofunc.oc_fasta_to_orf_seq_dict(path_6k_fa)

# list of genes to analyze - those with a DBD mapping
#  1732 genes
genes = list(set(pd.read_table(path_uniprot_dbd_to_best_iso)['gene']))

cl_gd = isocreate.init_gen_obj(path_6k_gtf, gene_list=genes)
cl_gd = isocreate.create_and_link_seq_related_obj(cl_gd, orf_seqs_6k, verbose=True)

genes = cl_gd.keys()

gc_gd = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
gc_gd = isocreate.create_and_link_seq_related_obj(gc_gd, orf_seqs_gc, verbose=True)

# %%

def return_orf_obj_from_enst(gene_gc, enst):
    """Given a gencode gene objection, return the orf_obj corr. to the enst."""
    for orf in gene_gc.orfs:
        if orf.enst == enst:
            return orf

def all_other_gc_6k_orfs_of_gene(gc_gd, cl_gd, symbol, anchor_orf):
    """Grab all other orfs from a gene_obj."""
    other_orfs = []
    for orf in gc_gd[symbol]:
        if orf.enst != anchor_orf.enst:
            other_orfs.append(orf)
    if symbol in cl_gd:
        for orf in cl_gd[symbol]:
            other_orfs.append(orf)
    return other_orfs

def all_other_6k_orfs(cl_gd, symbol):
    return cl_gd[symbol].orfs

def incorporate_flanks_in_domain_range(dom_info, anchor_orf, flank_status):
    """Extend the domain ranges by 15 AA down and upstream."""
    if flank_status == 'with_flank':
        pfam, name, cat, start, end = dom_info
        len_prot = len(anchor_orf.tr)
        extended_start = int(start) - 15
        extended_end = int(end) + 15
        if extended_start < 1:
            extended_start = 1
        if extended_end > len_prot:
            extended_end = len_prot
        return pfam, name, cat, extended_start, extended_end
    else:
        return dom_info

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


# list of symbols already processed (part 1)
df = pd.read_table('e_GC_6K_affect_of_mapped_uniprot_dbds_with_flank.tsv')
already_processed = set(df.gene)


# only process results for 6K TFs
for flank_status in ['with_flank']:

    with open('e_GC_6K_affect_of_mapped_uniprot_dbds_{}_part2.tsv'.format(flank_status), 'w') as ofile:
        ofile.write('gene\tisoacc\tdbd_pfam_acc\tdbd_name\tdom_start_idx\tdom_end_idx\tprot_cmp_block_string\tprot_cmp_string\tdom_len\tnum_match\t%AA_preserved_denom_dom_len\n')

        for symbol in non_zf_domains:
            if symbol in already_processed: continue
            if symbol not in gc_gd:
                print 'symbol not in gencode gene dict :' + symbol
                continue
            gene_gc = gc_gd[symbol]

            # for now, juan needs data, only process results for 6K TFs
            if symbol not in cl_gd: continue

            for enst, dom_infos in non_zf_domains[symbol].items():
                print symbol, enst
                # for now, give Juan non-ZF DBD mappings
                anchor_orf = return_orf_obj_from_enst(gene_gc, enst)

                for other_orf in all_other_6k_orfs(cl_gd, symbol):
                    # align anchor orf to all other orfs (GC and 6K orfs) (create an aln grp)
                    grp = isocreatealign.create_and_map_splice_based_align_obj([[anchor_orf, other_orf]])[0]
                    grp.enter()
                    alnf = grp.alnf
                    chain = alnf.chain
                    for dom_info in dom_infos:
                        dom_info = incorporate_flanks_in_domain_range(dom_info, anchor_orf, flank_status)
                        anchor_domain = isocreatefeat.create_and_map_domain(anchor_orf, dom_info, with_evals=False)
                        dom_len = len(anchor_domain.chain)
                        # # write out info for anchor domain # not writing out the anchor information for now
                        # start_idx, end_idx = anchor_domain.first.res.idx, anchor_domain.last.res.idx
                        # odata = [symbol, anchor_orf.name, anchor_domain.acc, anchor_domain.desc, start_idx, end_idx, 'M', 'M'*dom_len, dom_len, dom_len, 100]
                        # ofile.write('\t'.join(map(str, odata)) + '\n')
                        # chain of alnrs corresponding to domain
                        dom_alnrs = alnf.get_subset_aln_chain(anchor_domain.first.res.aln, anchor_domain.last.res.aln)
                        cat_str = ''.join([alnr.alnpb.cat for alnr in dom_alnrs])
                        cat_block = ''.join([x[0] for x in re.split('(M+)', cat_str) if x])
                        # for now, deletions and subst. treated as 'removals' or 'R'
                        cat_str = cat_str.replace('D', 'R').replace('S', 'R')
                        cat_block = cat_block.replace('D', 'R').replace('S', 'R')
                        num_matches, perc_aa_den_dom = compute_the_perc_aa_affected(cat_str, dom_len)
                        start_idx, end_idx = dom_alnrs[0].res2.idx, dom_alnrs[-1].res2.idx
                        start_idx = 'NA' if not start_idx else start_idx
                        end_idx = 'NA' if not end_idx else end_idx
                        odata = [symbol, other_orf.name, anchor_domain.acc, anchor_domain.desc, start_idx, end_idx, cat_block, cat_str, dom_len, num_matches, perc_aa_den_dom]
                        ofile.write('\t'.join(map(str, odata)) + '\n')
