# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 14:00:15 2016

@author: wentingbian
"""

import re
from Bio import SeqIO




#f = open("delta.axt")
#f = open("oc_seq.axt")
f = open("horf_full_seq.axt")
axt_records = {}
for i in f:
    index, ref, ref_start, ref_end, query, query_start, query_end, strand, score = i.strip().split(" ")
    ref = ref
    orf_id = int(re.findall("orf_id=(\d+)", query)[0])
    entrez_id = int(re.findall("gi=(\d+)", query)[0])
    ref_seq = f.next().strip()
    query_seq = f.next().strip()
    blank_line = f.next()
    #if orf_id in [3515,3905,4494,4731,4865,5290,6341,6474,8859,8920,8957,10266,11796,12774,13131] and entrez_id == 0: # I assign zero to orf_id without entrez gene id but these orf_id appeared multiple times with / without entrez gene id
    #    continue
    if (orf_id, entrez_id) not in axt_records:
        axt_records[(orf_id, entrez_id)] =  []
    axt_records[(orf_id, entrez_id)].append((index, ref, int(ref_start), int(ref_end), query, int(query_start), int(query_end), strand, int(score), ref_seq, query_seq))

orf_length = {}
orf_range = {}
for record in SeqIO.parse("horf_full_seq.fa", "fasta"):
    orf_id = int(re.findall("orf_id=(\d+)", record.id)[0])
    entrez = int(re.findall("gi=(\d+)", record.id)[0])
    try:
        cds_start = int(re.findall("cds_start=(\d+)", record.id)[0])
    except:
        cds_start = 0
    try:
        cds_end = int(re.findall("cds_end=(\d+)", record.id)[0])    
    except:
        cds_end = 0
    orf_length[(orf_id, entrez)] = len(record.seq)
    orf_range[(orf_id, entrez)] = (cds_start, cds_end)

            
def seperate_gap(ref_start, ref_end, query_start, query_end, ref_seq, query_seq):
    gaps = [(m.start(0)+1, m.end(0) - 1, "ref") for m in re.finditer(r"[atcg]-+[atcg]", ref_seq)]
    gaps.extend([(m.start(0) + 1, m.end(0) - 1, "query") for m in re.finditer(r"[atcg]-+[atcg]", query_seq)])
    results = []
    indel = []
    pos = 0
    ref_gap = 0
    query_gap = 0
    #print sorted(gaps)
    for gap in sorted(gaps):
        if gap[2] == "ref":
            indel.append([ref_start - ref_gap + gap[0] - 1, ref_start - ref_gap + gap[0], query_start - query_gap + gap[0] - 1, query_start - query_gap + gap[1], query_seq[gap[0]: gap[1]], "I"])
            results.append([ref_start + pos - ref_gap, ref_start - ref_gap + gap[0] - 1, query_start + pos - query_gap, query_start - query_gap + gap[0] - 1, ref_seq[pos: gap[0]], query_seq[pos : gap[0]]])
            pos = gap[1]
            ref_gap = ref_gap + gap[1] - gap[0]
        else:
            indel.append([ref_start - ref_gap + gap[0] - 1, ref_start - ref_gap + gap[1], query_start - query_gap + gap[0] - 1, query_start - query_gap + gap[0], ref_seq[gap[0]: gap[1]], "D"])
            results.append([ref_start + pos - ref_gap, ref_start - ref_gap + gap[0] - 1, query_start + pos - query_gap, query_start - query_gap + gap[0] - 1, ref_seq[pos : gap[0]], query_seq[pos : gap[0]]])
            pos = gap[1]
            query_gap = query_gap + gap[1] - gap[0]
    results.append([ref_start + pos - ref_gap, ref_end, query_start + pos - query_gap, query_end, ref_seq[pos: ], query_seq[pos: ]])
    #print results
    return results#, indel          

"""
## filter overlapping records
records_no_overlap = {}
for (orf_id, entrez), record in axt_records.items():
    record_dict = {}
    for i in record:
        if (i[1], i[7]) not in record_dict:
            record_dict[(i[1], i[7])] = []
        record_dict[(i[1], i[7])].append(i)
    rm_overlap = []
    for k, v in record_dict.items():
        v.sort(key = lambda x: (x[5], -x[6]))
        if len(v) == 1:
            rm_overlap.extend(v)
            continue
        tmp = [v[0]]
        minus_by = 1
        for i in xrange(1, len(v) - 1):
            if v[i][5] <= v[i - minus_by][6] and v[i][6] >= v[i + 1][5] and (v[i + 1][5] - v[i - minus_by][6] < 5):
                minus_by += 1
                continue
            elif v[i][5] >= v[i - minus_by][5] and v[i][6] <= v[i - minus_by][6]:
                minus_by += 1
                continue
            else:
                minus_by = 1
                tmp.append(v[i])
        if v[-1][6] > tmp[-1][6]:
            tmp.append(v[-1])
        rm_overlap.extend(tmp)
    records_no_overlap[(orf_id, entrez)] = rm_overlap            
## filter out additional chrom and strand
records_limit_pos = {}
for (orf_id, entrez), record in records_no_overlap.items():
    record_dict = {}
    for i in record:
        if (i[1], i[7]) not in record_dict:
            record_dict[(i[1], i[7])] = []
        record_dict[(i[1], i[7])].append(i)
    cover_info = {}
    for k, v in record_dict.items():
        base_covered = set()
        for i in v:
            base_covered = base_covered.union(range(i[5], i[6] + 1))
        cover_info[k] = orf_length[(orf_id, entrez)] - len(base_covered)
    min_missing = min(cover_info.values())
    best_pos = [k for k in record_dict if cover_info[k] == min_missing]
    tot_score = {}
    for k in best_pos:
        tot_score[k] = sum([i[8] for i in record_dict[k]])
    max_score = max(tot_score.values())
    best_score = [k for k in tot_score if tot_score[k] == max_score]
    tmp = []
    for k in best_score:
        tmp.extend(record_dict[k])
    records_limit_pos[(orf_id, entrez)] = tmp
"""
## deal with indel
records_no_indel = {}
for (orf_id, entrez), record in axt_records.items():
    records_no_indel[(orf_id, entrez)] = []
    for i in record:
        if "-" in i[9] or "-" in i[10]:
            if "-" in i[9]:
                print (orf_id, entrez)
            gaps = seperate_gap(i[2], i[3], i[5], i[6], i[9], i[10])
            for j in gaps:
                records_no_indel[(orf_id, entrez)].append([i[1], i[7]] + j[:4])
        else:
            records_no_indel[(orf_id, entrez)].append([i[1], i[7], i[2], i[3], i[5], i[6]])

## remove the UTR regions
records_final = {}
for (orf_id, entrez), record in records_no_indel.items():
    cds_start, cds_end = orf_range[(orf_id, entrez)]
    cds_length = orf_length[(orf_id, entrez)]
    records_final[(orf_id, entrez)] = []
    for i in record:
        chrom, strand, ref_start, ref_end, cdna_start, cdna_end = i
        if strand == "+":
            new_cdna_start = max(1, cdna_start - cds_start)
            new_ref_start =  ref_start - min(cdna_start - new_cdna_start - cds_start, 0)
            new_cdna_end = min(cds_length - cds_end - cds_start, cdna_end - cds_start)
            new_ref_end = ref_end - max(cdna_end - new_cdna_end - cds_start, 0)
        elif strand == "-":
            new_cdna_start = max(1, cdna_start - cds_end)
            new_ref_start =  ref_start - min(cdna_start - new_cdna_start - cds_end, 0)
            new_cdna_end = min(cds_length - cds_start - cds_end, cdna_end - cds_end)
            new_ref_end = ref_end - max(cdna_end - new_cdna_end - cds_end, 0)
        if new_cdna_start >= new_cdna_end:
            continue
        records_final[(orf_id, entrez)].append([chrom, strand, new_ref_start, new_ref_end, new_cdna_start, new_cdna_end])

#assert 0
import MySQLdb
import getpass
con = MySQLdb.connect(host = "paros", user = "user", passwd = getpass.getpass())
cursor = con.cursor()
for k, v in sorted(records_final.items(), key = lambda x: x[0][0]):
    for j in sorted(v, key = lambda x: (x[0], x[1], x[2])):
        cursor.execute("insert into horfeome_annotation.horfeome_cdna_mapping_hg38 (orf_id, entrez_gene_id, chrom, strand, ref_start, ref_end, cdna_start, cdna_end )values (%s, %s, %s, %s, %s, %s, %s, %s)", (k[0], k[1], j[0], j[1], j[2], j[3], j[4], j[5]))
