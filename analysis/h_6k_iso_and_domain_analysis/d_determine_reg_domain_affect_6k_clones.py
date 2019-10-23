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
path_reg_dom_6k = os.path.join(data_dir, 'regulatory_domains/e_regulatory_domains_mapped_to_6k_exact_and_clustal_updated.tsv')
# new filepaths for gtf and fa that uses the new unique_acc
path_6k_gtf = os.path.join(data_dir, 'hTFIso6K_valid_isoforms/c_6k_unique_acc_aligns.gtf')
path_6k_fa = os.path.join(data_dir, 'hTFIso6K_valid_isoforms/j2_6k_unique_isoacc_and_nt_seqs.fa')


# load data - isoform sequences and domains
orf_seqs_6k = isofunc.oc_fasta_to_orf_seq_dict(path_6k_fa)
domains = isofunc.load_regulatory_domain_mappings(path_reg_dom_6k)


reload(isocreate)
reload(isoclass)
reload(isofunc)
reload(isofeature)
reload(isomap)
reload(isoimage)
reload(isogroup)
reload(isoalign)
reload(isowrangle)
reload(isocreatealign)
reload(isocreatefeat)

# genes = ['NFYA', 'PAX5']  # toy example
# genes = ['AEBP2']  # when running full
genes = []
gd = isocreate.init_gen_obj(path_6k_gtf, genes)
gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs_6k)

# %%


def get_anchor_domain_mapping_info(gene):
    """From domain mappings of a gene, define anchor domains.
        Returns list:
            [[anchor_dom, anchor_orf, [other_orfs]]]
            Note - one uniq. anchor_dom acc represented per orf
    """
    anchor_domains = {}
    for orf in gene.orfs_len_ordered_desc:
        other_orfs = gene.get_other_orfs(orf)
        for dom in orf.doms:
            if dom.acc not in anchor_domains:
                anchor_domains[dom.acc] = [dom, orf, other_orfs]
    return anchor_domains.values()

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

with open('d_6k_isoforms_affect_of_mapped_reg_domains.tsv', 'w') as ofile:
    ofile.write('gene\tisoacc\tisoacc2\tdom_name\tdom_cat\tdom_start_idx\tdom_end_idx\tprot_cmp_block_string\tprot_cmp_string\tdom_len\tnum_match\t%AA_preserved_denom_dom_len\n')
    for symbol, gene in gd.items():
        for orf in gene:
            dom_infos = domains[orf.name]
            for dom_info in dom_infos:
                isocreatefeat.create_and_map_domain(orf, dom_info, with_evals=False)
        anchor_domain_infos = get_anchor_domain_mapping_info(gene)
        for anchor_domain, anchor_orf, other_orfs in anchor_domain_infos:
            first_time_writing_out_anchor_orf_mapping = True
            for other_orf in other_orfs:
                try:
                    # create an aln grp
                    grp = isocreatealign.create_and_map_splice_based_align_obj([[anchor_orf, other_orf]])[0]
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
                    if first_time_writing_out_anchor_orf_mapping:
                        start_idx, end_idx = anchor_domain.first.res.idx, anchor_domain.last.res.idx
                        odata = [symbol, anchor_orf.name, '', anchor_domain.acc, anchor_domain.cat, start_idx, end_idx, 'M', 'M'*dom_len, dom_len, dom_len, 100]
                        ofile.write('\t'.join(map(str, odata)) + '\n')
                        first_time_writing_out_anchor_orf_mapping = False
                    num_matches, perc_aa_den_dom = compute_the_perc_aa_affected(cat_str, dom_len)
                    start_idx, end_idx = dom_alnrs[0].res2.idx, dom_alnrs[-1].res2.idx
                    start_idx = 'NA' if not start_idx else start_idx
                    end_idx = 'NA' if not end_idx else end_idx
                    odata = [symbol, other_orf.name, '', anchor_domain.acc, anchor_domain.cat, start_idx, end_idx, cat_block, cat_str, dom_len, num_matches, perc_aa_den_dom]
                    ofile.write('\t'.join(map(str, odata)) + '\n')
                except:
                    print other_orf, anchor_domain
