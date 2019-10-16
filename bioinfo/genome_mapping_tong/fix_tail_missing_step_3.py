# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:40:50 2016

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
cursor.execute("select orf_id, gi, chrom, strand, ref_start, ref_end, cdna_start, cdna_end, id from horfeome_annotation.horfeome_cdna_mapping_hg38")
mapping_info = {}
for record in cursor:
    if (record[0],  record[2], record[3]) not in mapping_info:
        mapping_info[(record[0], record[2], record[3])] = []
    mapping_info[(record[0],  record[2] , record[3])].append(record)
  
for orf_id, chrom, strand in mapping_info:
    records = mapping_info[(orf_id, chrom, strand)]  
   # chrom, strand = records[0][2], records[0][3]
    if records[0][6] > 1:
        if strand == "-":
            ref = rev_comp(str(ref_seq[chrom][records[0][4] - records[0][6] : records[0][4] - 1]))
            cds = str(orf_seq[orf_id][-records[0][6] + 1:])
            if ref != cds:
                if len(ref) < 3:
                    ref = rev_comp(str(ref_seq[chrom][records[0][4] - records[0][6]: records[0][4] - records[0][6] + 3]))
                    cds = str(orf_seq[orf_id][-3:])
                    if ref in ["TGA", "TAA", "TAG"] and cds in ["TGA", "TAA", "TAG"]:
                        continue
                        print records[0][4], records[0][6], records[0][4] - records[0][6] + 1, 1, records[0][8]
                        cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_start = %s, cdna_start = %s where id = %s", (records[0][4] - records[0][6] + 1, 1, records[0][8]))
                else:
                    ref_last = ref[:-3]
                    cds_last = cds[:-3]
                    if ref_last == cds_last and ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"]: # there is no such case.
                        continue
                    elif ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"] and len(ref_last) <= 2:
                        print ref, cds, records[0][-1]
                        cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_start = %s, cdna_start = %s where id = %s", (records[0][4] - records[0][6] + 1, 1, records[0][8]))
                    else:
                        if len(cds_last)> 0 and hamming_dist(ref_last, cds_last) / float(len(cds_last)) <= 0.25 and ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"]:
                            #continue
                            print records[0][4], records[0][6], records[0][4] - records[0][6] + 1, 1, records[0][8]
                            cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_start = %s, cdna_start = %s where id = %s", (records[0][4] - records[0][6] + 1, 1, records[0][8]))
            else:
                continue
                print orf_id, ref, cds
                cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_start = %s, cdna_start = %s where id = %s", (records[0][4] - records[0][6] + 1, 1, records[0][8]))
        elif strand == "+":
            ref = str(ref_seq[chrom][records[0][4] - records[0][6]: records[0][4] - 1])
            cds = str(orf_seq[orf_id][:records[0][6] - 1])
            if ref != cds:
                #head[orf_id] = [ref, cds]
                if hamming_dist(ref, cds) / float(len(cds)) <= 0.25:
                    continue
                    print ref, cds, hamming_dist(ref, cds), len(cds)
                    cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_start = %s, cdna_start = %s where id = %s", (records[0][4] - records[0][6] + 1, 1, records[0][8]))
                #print ref, cds
                #continue
            else:
                continue
                #print orf_id, entrez, ref, cds
                print records[0][4], records[0][6], records[0][4] - records[0][6] + 1, 1, records[0][8]
                cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_start = %s, cdna_start = %s where id = %s", (records[0][4] - records[0][6] + 1, 1, records[0][8]))

cursor.execute("select orf_id, gi,  chrom, strand, ref_start, ref_end, cdna_start, cdna_end, id from horfeome_annotation.horfeome_cdna_mapping_hg38")
mapping_info = {}

for record in cursor:
    if (record[0],  record[2], record[3]) not in mapping_info:
        mapping_info[(record[0], record[2], record[3])] = []
    mapping_info[(record[0],  record[2] , record[3])].append(record)

for orf_id, chrom, strand in mapping_info:
    records = mapping_info[(orf_id, chrom, strand)]

    if len(orf_seq[orf_id]) > records[-1][7]:
        if strand == "+":
            ref = str(ref_seq[chrom][records[-1][5]: records[-1][5] + len(orf_seq[orf_id]) - records[-1][7]])
            cds = str(orf_seq[orf_id][records[-1][7]:])
            #print ref, cds
            if ref != cds:
                if len(ref) < 3:
                    ref = str(ref_seq[chrom][records[-1][5] + len(orf_seq[orf_id]) - records[-1][7] - 3: records[-1][5] + len(orf_seq[orf_id]) - records[-1][7]])
                    cds = str(orf_seq[orf_id][-3:])
                    if ref in ["TGA", "TAA", "TAG"] and cds in ["TGA", "TAA", "TAG"]:
                        print records[-1][5], records[-1][7], records[-1][5] + len(orf_seq[orf_id]) - records[-1][7], len(orf_seq[orf_id]), records[-1][-1]
                        cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_end = %s, cdna_end = %s where id = %s", (records[-1][5] + len(orf_seq[orf_id]) - records[-1][7], len(orf_seq[orf_id]), records[-1][-1]))
                else:
                    ref_last = ref[:-3]
                    cds_last = cds[:-3]
                    if ref_last == cds_last and ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"]: # there is no such case
                        print ref, cds
                        continue
                    elif ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"] and len(ref_last) <= 2:
                        print ref, cds
                        cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_end = %s, cdna_end = %s where id = %s", (records[-1][5] + len(orf_seq[orf_id]) - records[-1][7], len(orf_seq[orf_id]), records[-1][-1]))
                    elif len(ref_last)> 0 and hamming_dist(ref_last, cds_last) / float(len(cds_last)) <= 0.25 and ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"]:
                        #continue
                        print ref, cds
                        print records[-1][5], records[-1][7], records[-1][5] + len(orf_seq[orf_id]) - records[-1][7], len(orf_seq[orf_id]), records[-1][-1]
                    elif hamming_dist(ref, cds) / float(len(cds)) <= 0.25:
                        print ref, cds                                                    
            else:   
                print records[-1][5], records[-1][7], records[-1][5] + len(orf_seq[orf_id]) - records[-1][7], len(orf_seq[orf_id]), records[-1][-1]
                cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_end = %s, cdna_end = %s where id = %s", (records[-1][5] + len(orf_seq[orf_id]) - records[-1][7], len(orf_seq[orf_id]), records[-1][-1]))
        elif strand == "-":
            ref = rev_comp(str(ref_seq[chrom][records[-1][5]: records[-1][5] + len(orf_seq[orf_id]) - records[-1][7]]))
            cds = str(orf_seq[orf_id][:len(orf_seq[orf_id]) - records[-1][7]])
            if ref != cds:
#		print orf_id
#		print cds
#		print ref
                if hamming_dist(ref, cds) / float(len(cds)) <= 0.25:
                    print ref, cds
                    cursor.execute("update horfeome_annotation.horfeome_cdna_mapping_hg38 set ref_end = %s, cdna_end = %s where id = %s", (records[-1][5] + len(orf_seq[orf_id]) - records[-1][7], len(orf_seq[orf_id]), records[-1][-1]))
            else:
                print ref, cds, orf_id 
                print records[0][5], records[0][7], records[0][5] - records[0][7] + len(orf_seq[orf_id]), len(orf_seq[orf_id]), records[0][8]
