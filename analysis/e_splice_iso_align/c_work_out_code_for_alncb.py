#!/usr/bin/env python
# test

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
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa')
path_tf_names = os.path.join(data_dir, 'tf_genenames/b_orig_genename_in_gencode30_man_annot.txt')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
genes = isofunc.load_tf_genenames(path_tf_names)  # ~1700 genenames

# %%

# this took a long time to run, few hours
# once loaded into memory, wrote-out tables below

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

genes = ['NFYA', 'PAX5', 'CXXC1']

gd = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs)

for symbol, gene in gd.items():
    grps = isowrangle.get_same_prot_seq_grps_of_orfs_from_gene_(gene)
    repr_orfs = set(grp.repr_orf for grp in grps)  # uniq-prot-seq orfs
    gene.redundant_seq_orfs = gene.orfs - repr_orfs
    gene.orfs = repr_orfs

all_grps = []
i = 1
for symbol, gene in gd.items():
    print symbol, i
    i += 1
    grps = isocreatealign.create_and_map_splice_based_align_obj(gene.ref_alt_pairs)
    all_grps.append(grps)

# %%

# write out table with subst. and the occurance of frameshifts



# %%

# write out table with lengths of insertions/deletions
with open('d_subst_underlying_indel_vs_frmshift_info.tsv', 'w') as ofile:
    ofile.write('gene\torf1\torf2\tseg1\tseg2\tis_ragged\tseg_len\n')
    for grps in all_grps:
        for grp in grps:
            aln = grp.alnf
            for alnpb in aln.protblocks:
                if alnpb.cat == cat:
                    is_ragged = '0'
                    if cat == 'I':
                        if len(alnpb.aa1.replace('-','')) > 0:
                            is_ragged = '1'
                        slen = len(alnpb.aa2)
                    elif cat == 'D':
                        if len(alnpb.aa2.replace('-','')) > 0:
                            is_ragged = '1'
                        slen = len(alnpb.aa1)
                    odata = [aln.orf1.gene.name, aln.orf2.name,
                             aln.other_orf.name, alnpb.aa1, alnpb.aa2,
                             is_ragged, slen]
                    ofile.write('\t'.join(map(str, odata)) + '\n')

# %%

# write out table of tf isoform numbers stats (w filter for evidence codes, NF)
def are_both_ends_found(orf):
    if not orf.cds_start_nf and not orf.cds_end_nf:
        return '1'
    else:
        return '0'

with open('f_lambert_tfs_ref_alt_w_gencode_stats.tsv', 'w') as ofile:
    ofile.write('symbol\tisoname\tstatus\tis_basic\tappris\tlevel\ttsl\tends_found\tstart_nf\tend_nf\n')
    for grps in all_grps:
        for grp in grps:
            orf = grp.repr_orf
            odata = [orf.gene.name, orf.name, 'ref', orf.is_basic * 1,
                     orf.appris, orf.level, orf.tsl, are_both_ends_found(orf),
                     orf.cds_start_nf * 1, orf.cds_end_nf * 1]
            ofile.write('\t'.join(map(str, odata)) + '\n')
            break
        for grp in grps:
            alts = grp.other_orfs
            for orf in alts:
                odata = [orf.gene.name, orf.name, 'alt', orf.is_basic * 1,
                         orf.appris, orf.level, orf.tsl, are_both_ends_found(orf),
                         orf.cds_start_nf * 1, orf.cds_end_nf * 1]
                ofile.write('\t'.join(map(str, odata)) + '\n')



# %%

# write out table with the protein block and splice block calls
with open('c_align_block_calls_splice_and_protein.txt', 'w') as ofile:
    for grps in all_grps:
        for grp in grps:
            self = grp.alnf
            pbl = ''.join([alnr.alnpb.cat for alnr in self.chain])
            bl = ''.join([alnr.alnb.cat for alnr in self.chain])
            aa1 = ''.join([alnr.res1.aa for alnr in self.chain])
            aa2 = ''.join([alnr.res2.aa for alnr in self.chain])
            frm2 = ''.join([alnr.res2.rfrm for alnr in self.chain])
            frm2 = frm2.replace('1', ' ').replace('-', ' ').replace('*', ' ')
            ostr = ('{gene} {strand}\n'
                    '{alnpb:16s}PBL:{pbl}\n'
                    '{alnb:16s} BL:{bl}\n'
                    '{orf1:16s} AA:{aa1}\n'
                    '{orf2:16s} AA:{aa2}\n'
                    '                 FM:{frm2}')
            ostr = ostr.format(gene=self.orf1.gene.name, strand=self.orf1.strand,
                               alnpb='prot. blocks', pbl=pbl, alnb='splice blocks', bl=bl,
                               orf1=self.orf1.name, aa1=aa1, orf2=self.orf2.name,
                               aa2=aa2, frm2=frm2)
            ofile.write(ostr + '\n\n')

# %%

# output a table listing all protein segment calls
with open('c_align_block_calls_individual_blocks.txt', 'w') as ofile:
    for grps in all_grps:
        for grp in grps:
            ofile.write(str(grp)+ '\n')
            alnf = grp.alnf
            for alnpb in alnf.protblocks:
                orow = [grp.anchor_orf.gene, grp.anchor_orf, grp.other_orf,
                        alnpb.cat, alnpb.aa1, alnpb.aa2]
                ostr = '{:12s}\t{:12s}\t{:12s}\t{:4s}\t{}\t{}\n'.format(*orow)
                ofile.write(ostr)

# %%

# output a table listing all protein segment calls
with open('d_substitution_stats.tsv', 'w') as ofile:
    ofile.write('tf\torf1\torf2\tseg1_len\tseg2_len\tlen_diff\tperc_change\n')
    for grps in all_grps:
        for grp in grps:
            alnf = grp.alnf
            for alnpb in alnf.protblocks:
                if alnpb.cat == 'S':
                    anchor_len = len(alnpb.aa1)
                    other_len = len(alnpb.aa2)
                    len_diff = anchor_len - other_len
                    perc_change = float(other_len)/anchor_len * 100
                    orow = [grp.anchor_orf.gene, grp.anchor_orf, grp.other_orf,
                            anchor_len, other_len, len_diff, perc_change]
                    ostr = '\t'.join(map(str, orow)) + '\n'
                    ofile.write(ostr)

# %%

import subprocess

def get_clustal_scores(seq1, seq2):
    """Run clustalw and get alignment scores."""
    tmp_file = './c_prot_segments_clustal_temp_folder/tmp.fa'
    with open(tmp_file, 'w') as ofile:
        ofile.write('>orf1\n{}\n>orf2\n{}\n'.format(seq1, seq2))
    ostr = subprocess.check_output(['clustalw', '-GAPOPEN=80', '-INFILE={}'.format(tmp_file)])
    score1 = ostr.split('\n')[12].split('Score:')[1].strip()
    score2 = ostr.split('\n')[20].split('Score')[1].strip()
    return score1, score2


def get_clustal_alignment():
    tmp_file = './c_prot_segments_clustal_temp_folder/tmp.aln'
    lines = open(tmp_file).readlines()
    orf1 = ''.join([x.strip()[16:] for x in lines[3::4]])
    orf2 = ''.join([x.strip()[16:] for x in lines[4::4]])
    aln = ''.join([x[16:].rstrip('\n') for x in lines[5::4]])
    return orf1, orf2, aln


def get_perc_frameshift_of_the_substitution(alnpb):
    len_alnpb = len(alnpb.chain)
    blocks = alnpb.blocks
    len_frm = 0
    for block in blocks:
        if block.cat == 'F':
            len_frm += len(block.chain)
    perc_frm = float(len_frm) / len_alnpb * 100
    return len_frm, len_alnpb, perc_frm


def get_perc_identity_for_alignment(alnpb, aln):
    smallest_segment_length = min(len(alnpb.aa1), len(alnpb.aa2))
    num_matches = aln.count('*')
    perc_ident = float(num_matches)/smallest_segment_length
    return perc_ident



# output a table listing substitution protein segment and clustal align
# perc. frame is perc. of AA with frameshift compared to longest segment
with open('d_substitution_stats_w_clustal_seq_compare.txt', 'w') as ofile:
    ofile.write('tf\torf1\torf2\tseg1_len\tseg2_len\tlen_diff\tperc_change\t'
                'len_frameshift\tlen_pblock\tperc_frameshift\tperc_ident\t'
                'clustal_score1\tclustal_score2\tclustal_aln\n')
    for grps in all_grps:
        for grp in grps:
            alnf = grp.alnf
            for alnpb in alnf.protblocks:
                if alnpb.cat == 'S':
                    anchor_len = len(alnpb.aa1)
                    other_len = len(alnpb.aa2)
                    len_diff = anchor_len - other_len
                    perc_change = float(other_len) / anchor_len * 100
                    len_frm, len_pblock, perc_frm = get_perc_frameshift_of_the_substitution(alnpb)
                    score1, score2 = get_clustal_scores(alnpb.aa1, alnpb.aa2)
                    orf1, orf2, aln = get_clustal_alignment()
                    perc_ident = get_perc_identity_for_alignment(alnpb, aln)
                    orow = [grp.anchor_orf.gene, grp.anchor_orf, grp.other_orf,
                            anchor_len, other_len,
                            len_diff, perc_change,
                            len_frm, len_pblock, perc_frm, perc_ident,
                            score1, score2, orf1, orf2, aln]
                    ostr = '\t'.join(map(str, orow)) + '\n'
                    ofile.write(ostr)


# %%

def get_perc_frameshift_of_the_substitution(alnpb):
    len_alnpb = len(alnpb.chain)
    blocks = alnpb.blocks
    len_frm = 0
    for block in blocks:
        if block.cat == 'F':
            len_frm += len(block.chain)
    perc_frm = float(len_frm) / len_alnpb * 100
    return len_frm, len_alnpb, perc_frm



# as above, but only clustal scores
with open('d_substitution_stats_w_clustal_scores.txt', 'w') as ofile:
    ofile.write('tf\torf1\torf2\tseg1_len\tseg2_len\tlen_diff\tperc_change\t'
                'len_frameshift\tlen_pblock\tperc_frameshift\t'
                'clustal_score1\tclustal_score2\n')
    for grps in all_grps:
        for grp in grps:
            alnf = grp.alnf
            for alnpb in alnf.protblocks:
                if alnpb.cat == 'S':
                    anchor_len = len(alnpb.aa1)
                    other_len = len(alnpb.aa2)
                    len_diff = anchor_len - other_len
                    perc_change = float(other_len) / anchor_len * 100
                    score1, score2 = get_clustal_scores(alnpb.aa1, alnpb.aa2)
                    len_frm, len_pblock, perc_frm = get_perc_frameshift_of_the_substitution(alnpb)
                    orow = [grp.anchor_orf.gene, grp.anchor_orf, grp.other_orf,
                            anchor_len, other_len,
                            len_diff, perc_change,
                            len_frm, len_pblock, perc_frm, score1, score2]
                    ostr = '\t'.join(map(str, orow)) + '\n'
                    ofile.write(ostr)







# for grps in all_grps:
#     for grp in grps:
#         print grp
#         print grp.alnf.full
#         print ''.join([alnr.alnb.cat for alnr in grp.alnf.chain])
#         print ''.join([alnr.alnpb.cat for alnr in grp.alnf.chain])
#         for alnpb in grp.alnf.protblocks:
#             print alnpb.cat
#             print 'anchor sub:' + ''.join([res.aa for res in alnpb.res_chain1])
#             print 'other  sub:' + ''.join([res.aa for res in alnpb.res_chain2])
#             print '\n'
#         break
