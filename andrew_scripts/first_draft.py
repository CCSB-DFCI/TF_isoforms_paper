"""
Strategy:
(1) Make a huge dictionary for all CDSs from "CDS.gtf". For each ENST, copy this CDS dictionary, and delete keys/values not in that ENST. Add remaining CDSs to a list.
To make this dictionary, I am going to have to open the chromosome fasta file, and concatenate the sequence.
If negative strand, do reverse complement of that concatenated string.



(2) Make a huge dictionary for all exons from "gencode.v25.annotation.gtf". Key is ENST, value is list of exon numbers. 
For each ENST, copy this exon dictionary, and delete keys/values not in that ENST. Add dictionary value (exon numbers) to a list.
(3) Make a huge dictionary for all ENSTs from "gencode.v25.annotation.gtf". Key is ENST, value is list: [ENSG, ENST, exon list from above, CDS dictionary from above].
(4) Make a huge dictionary for all ENSGs from "gencode.v25.annotation.gtf". For each ENSG, copy the ENST dictionary, and delete keys/values not in that ENSG. 
"""
# Task: Make a dictionary of all CDS's. d{CDS: [ENST, chr #, start, end, CDS sequence]}

c_hashComplements = {"A":"T", "C":"G", "G":"C", "T":"A"}

def reverse_complement( strDNA ):
    strRet = ""
    for s in reversed( strDNA ):
        strRet += c_hashComplements.get( s, s )
    return strRet

genome = open("GRCh38.p7.genome.fa").read()
chrom_header_seq = genome.split(">")[1:26]

CDS_dict = {}

CDS_file = open("CDS_toy.gtf", "r")
for line in CDS_file:
	chromosome = line.split("\t")[0]
	ENST = "ENST" + line.split("ENST")[1][:13]
	exon_num = line.split("exon_number ")[1].split(";")[0]
	start = int(line.split("\t")[3])
	end = int(line.split("\t")[4])+3
	# I added 3 to the end coordinate because the given coordinates do not include the stop codon.
	for i in chrom_header_seq:
		if i[:4] == chromosome:
			i = i[7:]
			if line.split("\t")[6] == "+":
				CDS = i[start:end+1]
				print CDS
			elif line.split("\t")[6] == "-":
				CDS = reverse_complement(i[start:end+1])
				print CDS
			#For some reason, the start of this CDS was 1150 bases upstream of reported start coordinate
			# And the end of the CDS was 1165 bases upstream of the reported end coordinate
	# temp_list = []
	# temp_list.append(ENST)
	# temp_list.append(chromosome)
	# temp_list.append(start)
	# temp_list.append(end)
	# temp_list.append(CDS)

	# CDS_dict[ENST+" CDS "+str(exon_num)] = temp_list
#print CDS_dict

# final_fa = open("checkCDS.txt", "w")
# for i in CDS_dict:
# 	final_fa.write(">" + i + "\n" + CDS_dict[i] + "\n")
# final_fa.close()





	