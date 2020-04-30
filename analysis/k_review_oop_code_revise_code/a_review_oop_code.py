#!/usr/bin/env python
# Reviewing my OOP code while Luke goes through it.
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

reload(isoclass)
reload(isocreate)

data_dir = (
          '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/'
          'project_tf_isoforms/iso_master/data/'
          )

# filepaths
path_6k_gtf = os.path.join(data_dir, 'hTFIso6K_valid_isoforms/c_6k_unique_acc_aligns.gtf')
path_6k_fa = os.path.join(data_dir, 'hTFIso6K_valid_isoforms/j2_6k_unique_isoacc_and_nt_seqs.fa')
path_hg38 = os.path.join(data_dir, 'GRCh38.primary_assembly.genome.fa')

## load data
orf_seqs_6k = isofunc.oc_fasta_to_orf_seq_dict(path_6k_fa)
genes_6k = list(set([x.split('|')[0] for x in orf_seqs_6k.keys()]))

# luke - changed this to all genes if you want to make iso-matrices for all 366 TFs
# here I'm only taking the first two
genes_6k = genes_6k[0:2]

# make gene dictionary (gd)
gd = isocreate.init_gen_obj(path_6k_gtf, gene_list=genes_6k)
gd = isocreate.create_and_link_junct_and_ss_objs(gd, hg38)
gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs_6k, verbose=True)



for gname, gene in gd.items():
    for orf in gene:
        for junc in orf.juncs:
            print junc.full
            break
        break
# %%

def get_ordered_list_of_genome_coords(poss):
    """From a set/list of position objects, return an ordered list of all used
       genomic coordinates.
       Input - poss, a set of pos_obj for a gene
       Return - an ordered list of absolute genomic coordinates
    """
    coords = sorted(list(set([pos.coord for pos in poss])))
    return coords

def get_genome_coord_alignment_track(orf, coords):
    """For an isoform (orf), and the full set of genomic coords, including
       the spacers, return a string of the amino acids.
       For example: MMMAAATTT---AAA
    """
    track = ''  # a track of the orf AAs, aligned to genome
    orf_coords = get_genome_coords_for_orf(orf)
    for coord in coords:
        if coord in orf_coords:
            pos = orf_coords[coord]
            track += pos.res.aa
        else:
            track += '-'
    return track

def get_genome_coords_for_orf(orf):
    orf_coords = {pos.coord:pos for pos in orf.chain}
    return orf_coords


# for each gene, make an isoform alignment matrix, based on genome-coords
with open('a_isoform_align_matrices.tsv', 'w') as ofile:
    ofile.write('orf_name\tgeneome_aligned_track\n')
    for gname, gene in gd.items():
        coords = get_ordered_list_of_genome_coords(gene.poss)
        for orf in gene.orfs:
            track = get_genome_coord_alignment_track(orf, coords)
            ofile.write(orf.name + '\t' + track + '\n')
