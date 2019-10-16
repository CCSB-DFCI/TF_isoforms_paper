# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 11:59:09 2016

@author: wentingbian
"""

import MySQLdb
import getpass
import re
from Bio import SeqIO

def rev_comp(seq):
    seq = seq.upper()
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([seq_dict[base] if base in seq_dict else "N" for base in reversed(seq)])

bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

orf_seq = {}
for record in SeqIO.parse("all_seq.fa", "fasta"):
    orf_id = int(re.findall("orf_id=(\d+)", record.id)[0])
    try:
        cds_start = int(re.findall("cds_start=(\d+)", record.id)[0])
    except:
        cds_start = 0
    try:
        cds_end = int(re.findall("cds_end=(\d+)", record.id)[0])
    except:
        cds_end = 0
    	
    if cds_end == 0:
        orf_seq[int(orf_id)] = str(record.seq[cds_start:].upper())
    else:
        orf_seq[int(orf_id)] = str(record.seq[cds_start:-cds_end].upper())



con = MySQLdb.connect(host = "paros", user = "user", passwd = getpass.getpass())
cursor = con.cursor()
cursor.execute("select orf_id, gi, chrom, strand, ref_pos, ref_base, cdna_pos, cdna_base from horfeome_annotation.horfeome_cdna_snp_hg38 where mut_type is null")
records = cursor.fetchall()

a = 0
b = 0
c = 0
for orf_id, gi, chrom, strand, ref_pos, ref_base, cdna_pos, cdna_base in records:
    if strand == "-":
        ref_base = rev_comp(ref_base)
        cdna_base = rev_comp(cdna_base)
    if cdna_base == "N" and orf_seq[orf_id][cdna_pos - 1] == "Y":
        cursor.execute("update horfeome_annotation.horfeome_cdna_snp_hg38 set mut_type = 'Syn' where orf_id = %s and gi = %s and ref_pos = %s", (orf_id, gi, ref_pos))
        continue
    assert orf_seq[orf_id][cdna_pos - 1] == cdna_base
        #print strand, orf_seq[orf_id][cdna_pos - 2: cdna_pos + 1], cdna_base
    if cdna_pos % 3 == 0:
        cdna_codon = orf_seq[orf_id][cdna_pos - 3: cdna_pos]
        ref_codon = orf_seq[orf_id][cdna_pos - 3: cdna_pos - 1] + ref_base
    elif cdna_pos % 3 == 1:
        cdna_codon = orf_seq[orf_id][cdna_pos - 1: cdna_pos + 2]
        ref_codon = ref_base + orf_seq[orf_id][cdna_pos: cdna_pos + 2]
    else:
        cdna_codon = orf_seq[orf_id][cdna_pos - 2: cdna_pos + 1]
        #print cdna_pos, len(orf_seq[orf_id])
        ref_codon = orf_seq[orf_id][cdna_pos - 2] + ref_base + orf_seq[orf_id][cdna_pos]
    if len(cdna_codon) != 3:
        continue
    if codon_table[cdna_codon] == codon_table[ref_codon]:
        cursor.execute("update horfeome_annotation.horfeome_cdna_snp_hg38 set mut_type = 'Syn' where orf_id = %s and gi = %s and ref_pos = %s", (orf_id, gi, ref_pos))
    elif codon_table[cdna_codon] == "*":
        cursor.execute("update horfeome_annotation.horfeome_cdna_snp_hg38 set mut_type = 'Nonsense' where orf_id = %s and gi = %s and ref_pos = %s", (orf_id, gi, ref_pos))
    elif codon_table[ref_codon] == "*":
        cursor.execute("update horfeome_annotation.horfeome_cdna_snp_hg38 set mut_type = 'Other' where orf_id = %s and gi = %s and ref_pos = %s", (orf_id, gi, ref_pos))
    else:
        cursor.execute("update horfeome_annotation.horfeome_cdna_snp_hg38 set mut_type = 'Nonsyn' where orf_id = %s and gi = %s and ref_pos = %s", (orf_id, gi, ref_pos))
