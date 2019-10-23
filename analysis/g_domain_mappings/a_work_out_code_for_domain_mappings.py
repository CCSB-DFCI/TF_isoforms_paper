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

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master_take_3/data_input/'
          )

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa')
path_tf_names = os.path.join(data_dir, 'tf_genenames/b_orig_genename_in_gencode30_man_annot.txt')
path_pfam_gc = os.path.join(data_dir, 'pfam/2019-06-14_HMMER_mappings_to_all_GENCODEv30_LambertTFs_isoforms.txt')
path_pfam_6k = os.path.join(data_dir, 'pfam/2019-06-14_HMMER_mappings_to_all_6Kisoforms.txt')
path_pfam_names = os.path.join(data_dir, 'pfam/pfam_a_names.tsv')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
genes_all = isofunc.load_tf_genenames(path_tf_names)  # ~1700 genenames

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

# genes = ['NFYA', 'PAX5']
genes = genes_all

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



reload(isofunc)
reload(isofeature)
reload(isocreatefeat)

domains = isofunc.load_domain_mappings(path_pfam_gc, path_pfam_names)

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

# is now a method of alnf, but 1700TFs processed w version without method
def get_subset_aln_chain(self, first_alnr, last_alnr):
    """Return a list of alnr between first/last_alnr.
       Note - assumes that both alnr are in self.chain and first < last.
    """
    subchain = []
    in_subset = False
    for alnr in self.chain:
        if alnr == first_alnr:
            subchain.append(alnr)
            in_subset = True
        elif alnr == last_alnr:
            subchain.append(alnr)
            in_subset = False
        elif in_subset:
            subchain.append(alnr)
        else:
            pass
    return subchain

import re

# write out the alignment cateogry change for domain-domain comparisons
with open('a_domain_per_res_prot_block_change.txt', 'w') as ofile:
    ofile.write('ref_orf\talt_orf\tprot_cmp_block_string\tdom_name\tprot_cmp_string\t'
                'sub_aa1\tsub_aa2\tsub_len_aa1\tsub_len_aa2\tsub_len_diff\t'
                'sub_perc_diff\n')
    for grps in all_grps:
        for grp in grps:
            grp.enter()
            ref_orf = grp.repr_orf
            alt_orf = grp.other_orf
            chain = grp.alnf.chain
            for ref_featf in ref_orf.doms:
                alnf = grp.alnf
                # chain of alnrs corresponding to domain
                dom_alnrs = get_subset_aln_chain(alnf, ref_featf.first.res.aln, ref_featf.last.res.aln)
                cat_str = ''.join([alnr.alnpb.cat for alnr in dom_alnrs])
                cat_block = ''.join([x[0] for x in re.split('(M+)', cat_str) if x])
                extra = ['-'] * 5
                if cat_block == 'S':
                    alnpb = dom_alnrs[0].alnpb
                    aa1 = alnpb.aa1.replace('-', '')
                    aa2 = alnpb.aa2.replace('-', '')
                    len_aa1 = len(aa1)
                    len_aa2 = len(aa2)
                    len_diff = len_aa1 - len_aa2
                    perc_diff = float(len_diff) / len_aa1
                    extra = [aa1, aa2, len_aa1, len_aa2, len_diff, perc_diff]
                odata = [ref_orf.name, alt_orf.name, ref_featf.desc, cat_block, cat_str] + extra
                ofile.write('\t'.join(map(str, odata)) + '\n')
            grp.exit()

# %%


# compare domain len in ref/alt, for overlapping and same-name domain
# print out the lengths

def get_the_ress_in_alt_orf_aligned_to_ref_domain(ref_featf):
    return [featr.res.aln.res2 for featr in ref_featf.chain]

def get_the_domains_mapped_to_alt_ress(alt_ress):
    alt_feats = set()
    for alt_res in alt_ress:
        for featr in alt_res.doms:
            if featr.featf:
                alt_feats.add(featr.featf)
    return alt_feats

def if_exists_grab_alt_domain_mapped(ref_feat, alt_feats):
    """Matches are based on pfam accession."""
    ref_pfam = ref_feat.acc
    for alt_feat in alt_feats:
        if alt_feat.acc == ref_pfam:
            return alt_feat
    return None

with open('a_ref_alt_domain_len_compare.txt', 'w') as ofile:
    ofile.write('ref_orf\talt_orf\tpfam\tref_len\talt_pfam\talt_len\tlen_diff\talt_is_longer\n')
    grp_to_print_out = []
    for grps in all_grps:
        for grp in grps:
            grp.enter()
            ref_orf = grp.repr_orf
            alt_orf = grp.other_orf
            chain = grp.alnf.chain
            for ref_featf in ref_orf.doms:
                alt_ress = get_the_ress_in_alt_orf_aligned_to_ref_domain(ref_featf)
                alt_feats = get_the_domains_mapped_to_alt_ress(alt_ress)
                alt_feat_match = if_exists_grab_alt_domain_mapped(ref_featf, alt_feats)
                if alt_feat_match:
                    ref_len = len(ref_featf.chain)
                    alt_len = len(alt_feat_match.chain)
                    alt_is_longer = True if alt_len > ref_len else False
                    if alt_is_longer:
                        grp_to_print_out.append(grp)
                    len_diff = ref_len - alt_len
                    odata = (ref_orf.name, alt_orf.name, ref_featf.name,
                             ref_len, alt_feat_match.name,
                             alt_len, len_diff, alt_is_longer)
                    ofile.write('\t'.join(map(str, odata)) + '\n')
            grp.exit()

# %%

# find cases in the 1700TFs where a domain maps in the alt. but not ref.

with open('a_domain_in_alt_but_not_ref.txt', 'w') as ofile:
    ofile.write('ref_orf\talt_orf\tpfam_in_alt_only\n')
    grp_to_print_out_later = []
    for grps in all_grps:
        for grp in grps:
            grp.enter()
            ref_orf = grp.repr_orf
            alt_orf = grp.other_orf
            pfam_in_ref = set([ref_featf.desc for ref_featf in ref_orf.doms])
            pfam_in_alt = set([alt_featf.desc for alt_featf in alt_orf.doms])
            pfams_in_alt_only = pfam_in_alt - pfam_in_ref
            if pfams_in_alt_only:
                grp_to_print_out_later.append(grp)
                for alt_only_pfam in pfams_in_alt_only:
                    odata = (ref_orf.name, alt_orf.name, alt_only_pfam)
                ofile.write('\t'.join(map(str, odata)) + '\n')
            grp.exit()

# %%

def print_full_tracks_for_alignment_group(grp):
    """print all domain mappings for anchor/other isoform"""
    ostr = ''
    chain = grp.alnf.chain
    for orf, reslookup in [(grp.repr_orf, 'res1'), (grp.other_orf, 'res2')]:
        feats = sorted([(feat.desc, feat) for feat in orf.doms])
        for name, feat in feats:
            orf.current_feat = feat
            feat_str = ''.join(['X' if getattr(aln, reslookup).dom else '-' for aln in chain])
            ostr += '{:16s}{}\n'.format(name, feat_str)
            orf.current_feat = None
        aa_str = ''.join([getattr(aln, reslookup).aa for aln in chain])
        ostr += '{:16s}{}\n'.format(orf.name, aa_str)
        # cds_str = ''.join(['|' if getattr(aln, reslookup).is_at_cds_edge else str(getattr(aln, reslookup).cds.ord)[0] for aln in chain])
        # ostr += '{:16s}{}\n'.format(orf.name, cds_str)
        ostr += '\n'
    return ostr

with open('a_domain_in_alt_but_not_ref_align_tracks.txt', 'w') as ofile:
    for grp in grp_to_print_out_later:
        ofile.write(print_full_tracks_for_alignment_group(grp) + '\n')

# %%




# print all domain mappings for anchor/other isoform
with open('a_ref_alt_next_150_tf_domain_mappings_DBDs.txt', 'w') as ofile:
    for grps in all_grps:
        for grp in grps:
            ostr = ''
            chain = grp.alnf.chain
            all_feats = []
            for orf, reslookup in [(grp.repr_orf, 'res1'), (grp.other_orf, 'res2')]:
                feats = sorted([(feat.desc, feat) for feat in orf.doms if feat.cat == 'DBD' and 'ZF' not in feat.desc.upper()])
                all_feats.extend(feats)
            if len(all_feats) == 0:
                continue
            for orf, reslookup in [(grp.repr_orf, 'res1'), (grp.other_orf, 'res2')]:
                feats = sorted([(feat.desc, feat) for feat in orf.doms if feat.cat == 'DBD' and 'ZF' not in feat.desc.upper()])
                for name, feat in feats:
                    orf.current_feat = feat
                    feat_str = ''.join(['X' if getattr(aln, reslookup).dom else '-' for aln in chain])
                    ostr += '{:16s}{}\n'.format(name, feat_str)
                    orf.current_feat = None
                aa_str = ''.join([getattr(aln, reslookup).aa for aln in chain])
                ostr += '{:16s}{}\n'.format(orf.name, aa_str)
                # cds_str = ''.join(['|' if getattr(aln, reslookup).is_at_cds_edge else str(getattr(aln, reslookup).cds.ord)[0] for aln in chain])
                # ostr += '{:16s}{}\n'.format(orf.name, cds_str)
                ostr += '\n'
            ostr += '\n'
            ofile.write(ostr)
        ofile.write('*'*100 + '\n')


# %%




for dom in orf1.doms:
    orf1.current_feat = dom
    print dom.full
    # print [res.dom for res in orf1.res_chain]
    # print [cds for cds in orf1.cdss]
    # print [cds.dom for cds in orf1.cdss]
    orf1.current_feat = None

print '\n'

for dom in orf2.doms:
    orf2.current_feat = dom
    print dom.full
    orf2.current_feat = None































# comment
