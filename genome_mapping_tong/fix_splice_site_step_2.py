# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 14:57:09 2016

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

def hamming_dist(str1, str2):
    assert len(str1) == len(str2)
    distance = 0
    for i in xrange(len(str1)):
        if str1[i] != str2[i]:
            distance += 1
    return distance
   
def clean_splice_site(cds_sub, ref_up, ref_down, strand):
    result = []
    if strand == "+":
        conan1, conan2 = "GT", "AG"
    elif strand == "-":
        conan1, conan2 = "CT", "AC"
    for i in re.finditer(conan1, ref_up):
        for j in re.finditer(conan2, ref_down):
            up_seq = ref_up[:i.start()]
            down_seq = ref_down[j.end():]
            ref_seq = up_seq + down_seq
            if len(ref_seq) == len(cds_sub) and hamming_dist(ref_seq, cds_sub) <= 3:
                #print up_seq, down_seq, ref_seq
                #print ref_up, ref_down, cds_sub
                #print i.start(), j.end()
                #assert 0
                result.append([i.start(), j.end()])
    ## TODO what if there are several options
    #if len(result) > 1:
    #    print cds_sub, ref_up, ref_down, strand
    return result


ref_seq = {}
for record in SeqIO.parse("/data/bioinfo/ngs/helper/hg38/hg38.fa", "fasta"):
    ref_seq[record.id] = record.seq

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
        orf_seq[int(orf_id)] = record.seq[cds_start:]
    else:
        orf_seq[int(orf_id)] = record.seq[cds_start:-cds_end]

con = MySQLdb.connect(host = "paros", user = "user", passwd = getpass.getpass())
cursor = con.cursor()
cursor.execute("select orf_id, gi, chrom, strand, ref_start, ref_end, cdna_start, cdna_end, id from horfeome_annotation.horfeome_cdna_mapping_hg38 ")
mapping_info = {}
for record in cursor:
    if (record[0], record[1]) not in mapping_info:
        mapping_info[(record[0], record[1])] = []
    mapping_info[(record[0], record[1])].append(record)

for orf_id, entrez in mapping_info:
    records = mapping_info[(orf_id, entrez)]
    ## TODO remove some special cases
    splice_site_pos = {}    
    for i in xrange(1, len(records)):
            ## there maybe some bug near the splice sites
        if records[i][6] - records[i - 1][7] > 1 and records[i][4] - records[i - 1][5] > 1:
            ref_up = str(ref_seq[records[i][2]][records[i - 1][5] - 10: records[i - 1][5] + 10])
            ref_down = str(ref_seq[records[i][2]][records[i][4] - 11: records[i][4] + 9])
            cds_seq = orf_seq[orf_id] if records[i][3] == "+" else rev_comp(orf_seq[orf_id])
            cds_sub = str(cds_seq[records[i - 1][7] - 10: records[i][6] + 9])
            tmp = clean_splice_site(cds_sub, ref_up, ref_down, records[i][3])
            if len(tmp) > 1:
                print orf_id, entrez, i
            elif len(tmp) == 1:
                splice_site_pos[i] = tmp[0]
                up_splice_site, down_splice_site = splice_site_pos[i]
                up_site_move = up_splice_site - 10
                down_site_move = down_splice_site - 10
                print records[i - 1][5], records[i - 1][7], records[i][4], records[i][6], up_splice_site, down_splice_site
                print records[i - 1][5] + up_site_move, records[i - 1][7] + up_site_move, records[i][4] + down_site_move, records[i][6] + down_site_move
                cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_end = %s, cdna_end = %s where id = %s", (records[i - 1][5] + up_site_move, records[i - 1][7] + up_site_move, records[i - 1][8]))
                cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_start = %s, cdna_start = %s where id = %s", (records[i][4] + down_site_move, records[i][6] + down_site_move, records[i][8]))
                
con.close()
