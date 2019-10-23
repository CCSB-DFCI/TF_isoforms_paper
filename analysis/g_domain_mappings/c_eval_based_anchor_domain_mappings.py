#!/usr/bin/env python
#
# title           :driver.py
# description     :Work out domain mappings to protein isoforms.
# author          :Gloria Sheynkman
# date            :June 17th, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================


import os
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
from bx.intervals.intersection import Interval, IntervalTree
import math

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master_take_3/data_input/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf.toy')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa.toy')
path_tf_names = os.path.join(data_dir, 'tf_genenames/b_orig_genename_in_gencode30_man_annot.txt')
# path_pfam_gc = os.path.join(data_dir, 'pfam/2019-06-14_HMMER_mappings_to_all_GENCODEv30_LambertTFs_isoforms.txt')
# path_pfam_6k = os.path.join(data_dir, 'pfam/2019-06-14_HMMER_mappings_to_all_6Kisoforms.txt')
path_pfam_eval = os.path.join(data_dir, 'pfam/2019-07-04_HMMER_domain_mappings_to_GS_fasta_file.txt')
path_pfam_names = os.path.join(data_dir, 'pfam/pfam_a_names.tsv')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
# genes_all = isofunc.load_tf_genenames(path_tf_names)  # ~1700 genenames

reload(isocreate)
reload(isoclass)
reload(isocreatefeat)
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
# genes = genes_all

gd = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs)

for symbol, gene in gd.items():
    grps = isowrangle.get_same_prot_seq_grps_of_orfs_from_gene_(gene)
    repr_orfs = set(grp.repr_orf for grp in grps)  # uniq-prot-seq orfs
    gene.redundant_seq_orfs = gene.orfs - repr_orfs
    gene.orfs = repr_orfs


domains = isofunc.load_eval_domain_mappings(path_pfam_eval, path_pfam_names)



for symbol, gene in gd.items():
    for orf in gene:
        domain_infos = domains[orf.name]
        for info in domain_infos:
            isocreatefeat.create_and_map_domain(orf, info)



# %%

def is_eval_for_anchor_best_or_equal(eval_anchor, eval_other):
    """Determine if the evalue for domain mapped to anchor isoform is better
       or at least relatively equal to the domain mapped to other isoform.
    """
    if eval_anchor < eval_other:
        return True
    log10_perc_diff = math.log10(eval_anchor)/math.log10(eval_other)
    if 0.9 < log10_perc_diff < 1.1:
        return True
    return False

with open('c_gc_ref_vs_other_domain_eval_compare.tsv', 'w') as ofile:
    ofile.write('gene\tanchor_orf\tother_orf\tdomain_cat\tpfam\tdomain_name\tref_eval\tother_eval\tis_best\n')
    for symbol, gene in gd.items():
        # store ranges of domains mapped to 'other' orfs
        tree = IntervalTree()  # holds domain intervals, for 'other' domains
        for orf in gene.other_orfs:
            for feat in orf.doms:
                start, end = sorted([feat.first.res.p1.coord, feat.last.res.p3.coord])
                tree.insert_interval(Interval(start, end, value=feat))
        # for each domain mapped to 'anchor' orf (here, appris-best), find overlps
        for anchor_feat in gene.repr_orf.doms:
            start, end = sorted([feat.first.res.p1.coord, feat.last.res.p3.coord])
            for hit in tree.find(start, end):
                other_feat = hit.value
                if anchor_feat.acc == other_feat.acc:
                    if is_eval_for_anchor_best_or_equal(anchor_feat.eval, other_feat.eval):
                        is_best_eval = '1'
                    else:
                        is_best_eval = '0'
                    odata = [symbol, gene.repr_orf.name, other_feat.orf.name, anchor_feat.cat, anchor_feat.acc, anchor_feat.name, anchor_feat.eval, other_feat.eval, is_best_eval]
                    ofile.write('\t'.join(map(str, odata)) + '\n')





# %%

for grps in all_grps:
    first_grp = True
    for grp in grps:
        if first_grp:
            orf1 = grp.repr_orf
            domain_infos = domains[orf1.ensp]
            for info in domain_infos:
                isocreatefeat.create_and_map_domain(orf1, info)
            first_grp = False
        orf2 = grp.other_orf
        domain_infos = domains[orf2.ensp]
        for info in domain_infos:
            isocreatefeat.create_and_map_domain(orf2, info)


# %%

def return_paired_domain(feat, other_orf):
    """Given a feat mapping to the anchor orf, see if there is an overlapping
       feature that should be paired with the anchor feat.
    """
    alt_ress = [featr.res.aln.res2 for featr in feat.chain]
    for alt_res in alt_ress:
        for alt_feat in alt_res.doms:
            if feat.desc == alt_feat.featf.desc:
                return alt_feat.featf
    return None

def get_paired_info(grp):
    """Based on protein alignment, determine overlapping domains.
       Return 'paired domains' between orf1/orf2, as well as cases where domain
       only maps to ref (orf1) or alt (orf2).
    """
    orf1 = grp.repr_orf
    orf2 = grp.other_orf
    paired_feats = []  # [[ref_feat, alt_feat], etc...]
    paired_off_alt_feats = []  # for bookkeeping
    unpaired_ref_feats = []
    unpaired_alt_feats = []
    for feat in orf1.doms:
        paired_alt_feat = return_paired_domain(feat, orf2);
        if paired_alt_feat:
            paired_feats.append([feat, paired_alt_feat])
            paired_off_alt_feats.append(paired_alt_feat)
        else:
            unpaired_ref_feats.append(feat)
    # process alt. feats
    for afeat in orf2.doms:
        if afeat not in paired_off_alt_feats:
            unpaired_alt_feats.append(afeat)
    # print grp
    # print 'paired domains', paired_feats
    # print 'paired alts   ', paired_off_alt_feats
    # print 'unpaired refs ', unpaired_ref_feats
    # print 'unpaired alts ', unpaired_alt_feats
    return paired_feats, unpaired_ref_feats, unpaired_alt_feats



# find all paired domains
# paired domain - domain of same-pfam-name, overlapping by at least 1+ AAs
for grps in all_grps:
    for grp in grps:
        grp.enter()
        paired_feats, unpaired_ref_feats, unpaired_alt_feats = get_paired_info(grp)
        # TODO - code to decide anchor domain to transfer from, based on best
        #        mapping to ref or alt, for now move fwd with ref-mapped feat

        grp.exit()
