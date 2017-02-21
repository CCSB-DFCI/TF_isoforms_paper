# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 10:15:17 2016

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

con = MySQLdb.connect(host = "paros", user = "user", passwd = "getpass.getpass()")
cursor = con.cursor()
cursor.execute("select orf_id, gi, chrom, strand, ref_start, ref_end, cdna_start, cdna_end from horfeome_annotation.horfeome_cdna_mapping_hg38")
mapping_info = {}
for record in cursor:
    if (record[0], record[1]) not in mapping_info:
        mapping_info[(record[0], record[1])] = []
    mapping_info[(record[0], record[1])].append(record[2:])


mut_info = {}
orf_map_status = {}
insertion = {}
deletion = {}
substitution = {}
for orf_id, entrez_gene_id in mapping_info:
    flag = {"mapped": 1, "mapped_to_multi": 0, "tail_unmapped": 0, "splice_site_unmapped": 0, "insert_over_25": 0,
            "not_canon_splice": 0, "delete_over_25": 0, "with_snp": 0, "with_indel": 0}
    records = mapping_info[(orf_id, entrez_gene_id)]
    records.sort(key = lambda x: (x[4], -x[5]))
    ref_pos = [(i[4], i[5]) for i in records]
    if any(i[0]!=records[0][0] or i[1]!=records[0][1] for i in records):
        flag["mapped_to_multi"] = 1
    if records[0][4] > 1 or len(orf_seq[orf_id]) > records[-1][5]:
        flag["tail_unmapped"] = 1
        
    for i in xrange(1, len(records)):
        chrom, strand, ref_start, ref_end, cdna_start, cdna_end = records[i]
        if strand == "+":
            cdna_start = records[i - 1][5]
            cdna_end = records[i][4]
        else:
            cdna_end = len(orf_seq[orf_id]) - records[i - 1][5] + 1
            cdna_start = len(orf_seq[orf_id]) - records[i][4] + 1
        if records[i][4] - records[i - 1][5] > 1 and records[i][2] - records[i - 1][3] > 1: # unmapped near splice sites or indel substitution
            if records[i][2] - records[i - 1][3] >= 25:
                flag["splice_site_unmapped"] = 1 # splice site unmapped
            elif records[i][4] - records[i - 1][5] >= 25:
                flag["insert_over_25"] = 1
            else:
                flag["with_indel"] = 1
                if (orf_id, entrez_gene_id) not in substitution:
                    substitution[(orf_id, entrez_gene_id)] = []
                if strand == "+":
                    tmp = str(orf_seq[orf_id][records[i - 1][5]: records[i][4] - 1]).upper()
                else:
                    tmp = str(rev_comp(orf_seq[orf_id])[records[i - 1][5]: records[i][4] - 1]).upper()
                substitution[(orf_id, entrez_gene_id)].append([chrom, strand, records[i - 1][3], records[i][2], str(ref_seq[chrom][records[i - 1][3]: records[i][2] - 1]).upper(), cdna_start, cdna_end, tmp])
        elif records[i][4] - records[i - 1][5] > 1: # insertion
            if records[i][4] - records[i - 1][5] >= 25:
                flag["insert_over_25"] = 1
            else:
                flag["with_indel"] = 1
                if (orf_id, entrez_gene_id) not in insertion:
                    insertion[(orf_id, entrez_gene_id)] = []
                if strand == "+":
                    tmp = str(orf_seq[orf_id][records[i - 1][5]: records[i][4] - 1]).upper()
                else:
                    tmp = str(rev_comp(orf_seq[orf_id])[records[i - 1][5]: records[i][4] - 1]).upper()
                insertion[(orf_id, entrez_gene_id)].append([chrom, strand, records[i - 1][3], records[i][2], str(ref_seq[chrom][records[i - 1][3]: records[i][2] - 1]).upper(), cdna_start, cdna_end, tmp])
        elif records[i][2] - records[i - 1][3] > 1: # correct splice or deletion
            if (strand == "+" and str(ref_seq[chrom][records[i - 1][3]: records[i - 1][3]+2]) != "GT" or str(ref_seq[chrom][records[i][2]-3: records[i][2]-1]) != "AG") or (strand == "-" and str(ref_seq[chrom][records[i - 1][3]: records[i - 1][3]+2]) != "CT" or str(ref_seq[chrom][records[i][2]-3: records[i][2]-1]) != "AC"):
                flag["not_canon_splice"] = 1
                if records[i][2] - records[i - 1][3] > 50:
                    continue
                elif records[i][2] - records[i - 1][3] >= 25:
                    flag["delete_over_25"] = 1
                else:
                    flag["with_indel"] = 1
                    if (orf_id, entrez_gene_id) not in deletion:
                        deletion[(orf_id, entrez_gene_id)] = []
                    deletion[(orf_id, entrez_gene_id)].append([chrom, strand, records[i - 1][3], records[i][2], str(ref_seq[chrom][records[i - 1][3]: records[i][2] - 1]).upper(), cdna_start, cdna_end, str(orf_seq[orf_id][records[i - 1][5]: records[i][4] - 1]).upper()])
    
    if (orf_id, entrez_gene_id) not in mut_info:
        mut_info[(orf_id, entrez_gene_id)] = {}
    for i in records:
        chrom, strand, ref_start, ref_end, cdna_start, cdna_end = i
        orf_id, entrez, ref_start, ref_end, cdna_start, cdna_end = int(orf_id), int(entrez_gene_id), int(ref_start), int(ref_end), int(cdna_start), int(cdna_end)
        ref = ref_seq[chrom][ref_start - 1: ref_end].upper()
        if strand == "+":
            cds = orf_seq[orf_id][cdna_start - 1: cdna_end].upper()
        elif strand == "-":
            cds = rev_comp(orf_seq[orf_id])[cdna_start - 1: cdna_end].upper()
        for base in xrange(len(ref)):
	    if base<len(cds):	 	
            	if ref[base] != cds[base]:
                	#print ref[base], cds[base], len(orf_seq[orf_id]) - base - cdna_start + 1, base + ref_start              
                	if strand == "+":
                    		mut_info[(orf_id, entrez_gene_id)][base + cdna_start] = [ref[base], cds[base], base + ref_start, chrom, strand]
                    		flag["with_snp"] = 1
                	elif strand == "-":                 
                    		mut_info[(orf_id, entrez_gene_id)][len(orf_seq[orf_id]) - base - cdna_start + 1] = [ref[base], cds[base], base + ref_start, chrom, strand]
                    		flag["with_snp"] = 1
    orf_map_status[(orf_id, entrez_gene_id)] = flag
       
 
for (orf_id, entrez_id), flag in orf_map_status.items():
    cursor.execute("insert into horfeome_annotation.horfeome_cdna_mapping_status(orf_id, gi, mapped, mapped_to_multi, tail_unmapped, splice_site_unmapped, insert_over_25, not_canon_splice, delete_over_25, with_snp, with_indel ) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", (orf_id, entrez_id) + tuple([flag[i] for i in ["mapped", "mapped_to_multi", "tail_unmapped", "splice_site_unmapped", "insert_over_25", "not_canon_splice", "delete_over_25", "with_snp", "with_indel"]]))

mut_info_filter = {i: j for i, j in mut_info.iteritems() if len(j) > 0}

for (orf_id, entrez), v in mut_info_filter.items():
    for cdna_pos, value in sorted(v.items(), key = lambda x: x[1][2]):
        ref_base, cdna_base, ref_pos, chrom, strand = value
        cursor.execute("insert into horfeome_annotation.horfeome_cdna_snp_hg38 (orf_id, gi, chrom, strand, ref_pos, ref_base, cdna_pos, cdna_base) values (%s, %s, %s, %s, %s, %s, %s, %s)", (orf_id, entrez, chrom, strand, ref_pos, ref_base, cdna_pos, cdna_base))

for (orf_id, entrez), v in substitution.items():
    for i in v:
        cursor.execute("insert into horfeome_annotation.horfeome_cdna_indel_hg38 values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", (orf_id, entrez) + tuple(i))

for (orf_id, entrez), v in deletion.items():
    for i in v:
        cursor.execute("insert into horfeome_annotation.horfeome_cdna_indel_hg38 values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", (orf_id, entrez) + tuple(i))

for (orf_id, entrez), v in insertion.items():
    for i in v:
        cursor.execute("insert into horfeome_annotation.horfeome_cdna_indel_hg38 values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", (orf_id, entrez) + tuple(i))

"""
#deletion = {}
#insertion = {}
#tail = {}
#head = {}
#to_solve = {}

    chrom, strand = records[0][1], records[0][2]
    ## start from 1?
    if records[0][5] > 1:
        if strand == "-":
            ref = rev_comp(str(ref_seq[chrom][records[0][3] - records[0][5]: records[0][3] - 1]))
            cds = str(orf_seq[orf_id][-records[0][5] + 1:])
            if ref != cds:
                if len(ref) < 3:
                    ref = rev_comp(str(ref_seq[chrom][records[0][3] - records[0][5]: records[0][3] - records[0][5] + 3]))
                    cds = str(orf_seq[orf_id][-3:])
                    if ref in ["TGA", "TAA", "TAG"] and cds in ["TGA", "TAA", "TAG"]:
                        continue
                else:
                    ref_last = ref[:-3]
                    cds_last = cds[:-3]
                    if ref_last == cds_last and ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"]:
                        continue
                tail[orf_id] = [ref, cds]
        elif strand == "+":
            ref = str(ref_seq[chrom][records[0][3] - records[0][5]: records[0][3] - 1])
            cds = str(orf_seq[orf_id][:records[0][5] - 1])
            if ref != cds:
                head[orf_id] = [ref, cds]
    
    # any indel inside?
    for i in xrange(1, len(records)):
        if records[i][5] - records[i - 1][6] > 1 and records[i][3] - records[i - 1][4] > 1:
            ## TODO the mismatch near splice sites
            if orf_id not in to_solve:
                to_solve[orf_id] = []
            to_solve[orf_id].append(records[i])
        elif records[i][5] - records[i - 1][6] > 1: # insert
            if orf_id not in insertion:
                insertion[orf_id] = []
            insertion[orf_id].append([records[i - 1][6], records[i][5]])
        elif records[i][3] - records[i - 1][4] > 1 and records[i][3] - records[i - 1][4] < 100 and records[i][5] - records[i - 1][6] == 1: # TODO delete or splice
            if strand == "+":
                if str(ref_seq[chrom][records[i - 1][4]: records[i - 1][4]+2]) == "GT" and str(ref_seq[chrom][records[i][3]-3: records[i][3]-1]) == "AG":
                    continue
                else:
                    a += 1
                    if orf_id not in deletion:
                        deletion[orf_id] = {}
                    deletion[orf_id][(records[i - 1][6], records[i][5])] = str(ref_seq[chrom][records[i - 1][4]: records[i][3]])
            elif strand == "-":
                if str(ref_seq[chrom][records[i - 1][4]: records[i - 1][4]+2]) == "CT" and str(ref_seq[chrom][records[i][3]-3: records[i][3]-1]) == "AC":
                    continue
                else:
                    a += 1
                    if orf_id not in deletion:
                        deletion[orf_id] = {}
                    deletion[orf_id][(records[i - 1][6], records[i][5])] = str(ref_seq[chrom][records[i - 1][4]: records[i][3]])

    # stop at the last base?                    
    if len(orf_seq[orf_id]) > records[-1][6]:
        if strand == "+":
            ref = str(ref_seq[chrom][records[-1][4]: records[-1][4] + len(orf_seq[orf_id]) - records[-1][6]])
            cds = str(orf_seq[orf_id][records[-1][6]:])
            if ref != cds:
                if len(ref) < 3:
                    ref = str(ref_seq[chrom][records[-1][4] + len(orf_seq[orf_id]) - records[-1][6] - 3: records[-1][4] + len(orf_seq[orf_id]) - records[-1][6]])
                    cds = str(orf_seq[orf_id][-3:])
                    if ref in ["TGA", "TAA", "TAG"] and cds in ["TGA", "TAA", "TAG"]:
                        continue
                else:
                    ref_last = ref[:-3]
                    cds_last = cds[:-3]
                    if ref_last == cds_last and ref[-3:] in ["TGA", "TAA", "TAG"] and cds[-3:] in ["TGA", "TAA", "TAG"]:
                        continue
                tail[orf_id] = [ref, cds]
        elif strand == "-":
            ref = rev_comp(str(ref_seq[chrom][records[-1][4]: records[-1][4] + len(orf_seq[orf_id]) - records[-1][6]]))
            cds = str(orf_seq[orf_id][:len(orf_seq[orf_id]) - records[-1][6]])
            if ref != cds:
                head[orf_id] = [ref, cds]                

cursor.close()
"""
