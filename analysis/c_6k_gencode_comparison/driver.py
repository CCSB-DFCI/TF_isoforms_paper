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
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v30.pc_transcripts.fa')
path_6k_gtf = os.path.join(data_dir, 'hTFIso6K_gtfhg38.star.gtf')
path_6k_fa = os.path.join(data_dir, 'hTFIso6K_sequences_w_orfid.fasta')
# genenames changed in GC30, manually updated genenames
path_gene_map = '../a_gencode_repr_orf_of_same_prot_seq/b_orig_genename_in_gencode30_man_annot.txt'


# load sequence data
orf_seqs_gc = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
orf_seqs_6k = isofunc.oc_fasta_to_orf_seq_dict(path_6k_fa)

# load corrected genenames (orig_symbol, and in_GC30), {orig_symb: gc30_symb}
symbol_map = isofunc.get_updated_genenames(path_gene_map)

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
reload(isocompare)

# # toy
# genes = ['NFYA']

# old genenames, in 6k
genes = isofunc.get_list_of_6k_genes(path_6k_fa)
# same as 6k, but with gc30 updated genes
genes_gc30_converted = [symbol_map[gene] for gene in genes]

d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes_gc30_converted)
d_gc = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs_gc)

d_6k = isocreate.init_gen_obj(path_6k_gtf, gene_list=genes)
d_6k = isocreate.create_and_link_seq_related_obj(d_6k, orf_seqs_6k)
# updated the old 6k genenames to gc30, to allow for comparison below
d_6k_updated = {}  # {symbol: gen_obj}
for symbol, gen_obj in d_6k.items():
    new_symbol = symbol_map[symbol]
    gen_obj.name = new_symbol
    d_6k_updated[new_symbol] = gen_obj

# %%

reload(isocompare)

# compare 6KOC sequences to GC
# generate novelty table (from 6K comparison to GC)
def main():
    ofile = open('6K_vs_GC30_cds_coords_matches.tsv', 'w')
    ofile.write('gene\tisoacc\tfrac_tr\tgc_match\tgc_isoname\tgc_enst'
                '\tappris\tlevel\ttsl\n')
    for name, gen_obj in d_6k_updated.items():
        for orf in gen_obj:
            matches_to_gc_orfs = isocompare.get_all_matches_to_orfs(orf, d_gc,
                                                match_mode='cds_intron_chain')
            if there_are_no_matches(matches_to_gc_orfs):
                odata = [orf.gene.name, orf.name, orf.frac_tr, '0', '-', '-',
                         '-', '-', '-']
            else:
                gorf = get_repr_gc_orf_that_this_orf_matched(matches_to_gc_orfs)
                odata = [orf.gene.name, orf.name, orf.frac_tr, '1', gorf.name,
                         gorf.enst, gorf.appris, gorf.level, gorf.tsl]
            odata = map(str, odata)
            ofile.write('\t'.join(odata) + '\n')
    ofile.close()

def there_are_no_matches(matches):
    if len(matches) == 0:
        return True

def get_repr_gc_orf_that_this_orf_matched(matches):
    """Sometimes 6K ORF matches to multiple, identical GC orfs with same prot.
       seq. Return the one representative match.
    """
    if len(matches) > 1:
        # get the best representative GC orf
        repr_orf = isogroup.SameSequenceGencodeGroup(matches).repr_orf
        return repr_orf
    else:
        # return the one GC orf
        return next(iter(matches))

if __name__ == '__main__':
    main()
