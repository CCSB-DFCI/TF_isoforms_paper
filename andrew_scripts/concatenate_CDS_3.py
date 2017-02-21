"""
Andrew Goldfarb
02/16/17
Task: From gtf "positive_strand_CDS.gtf", which has all the + strand CDS information for chromosome 22, collect
the exon coordinates to concatenate the exon sequences from the chromosome 22 sequence from file, 
"GRCh38.p2.genome_chr22_only.fa". 
"""
# (1) Read the CDS. Collect ENST Accessions, exon numbers, and [start, end] positions into lists:
ENST_accession = []
exon_number = []
start_end_list = []
positive_CDS_file = open("positive_strand_CDS.gtf", "r")
for line in positive_CDS_file:
	ENST = "ENST" + line.split("ENST")[1][:13]
	ENST_accession.append(ENST)
	exon_number.append(ENST + " exon " + line.split("exon_number ")[1].split(";")[0])
	start_end_list.append([int(line.split("\t")[3]), int(line.split("\t")[4])+3])
	# I added 3 to the end coordinate because the given coordinates do not include the stop codon.
ENST_accession = set(ENST_accession)

# (2) Make a dictionary that holds the start and end coordinates(value) for each ENST exon (key).
exon_coordinates_dict = {}
for i, j in zip(exon_number, start_end_list):
	exon_coordinates_dict[i] = j

# (3) Make a nested dictionary, where each ENST (key) holds a dictionary of the each exon and corresponding
# start and end coordinates. So: d[ENST][exon#][[start, end]]
import copy
ENST_exon_coordinates_dict = {}
for i in ENST_accession:
	temp_dict = copy.deepcopy(exon_coordinates_dict)
	for j in exon_coordinates_dict:
		if j[:17] != i:
			del temp_dict[j]
			ENST_exon_coordinates_dict[i] = temp_dict
	temp_dict = {}

chr22_string = open("GRCh38.p2.genome_chr22_only.fa").read()
chr22_string = chr22_string[8:]

# (4) Concatenate the exons together in the correct order
ENST_CDS_dict = {}
for i in ENST_exon_coordinates_dict:
	exon_numbers_list = []
	for index in str(ENST_exon_coordinates_dict[i]).split(str(i) + " exon ")[1:]:
		exon_numbers_list.append(int(index.split("'")[0]))
	exon_numbers_list = set(exon_numbers_list)
	# Example: set([2, 3, 4, 5, 6, 7, 8, 9])
	# Now concatenate the exons!
	CDS = ""
	for item in range(min(exon_numbers_list), max(exon_numbers_list)+1):
		start_pos = int(ENST_exon_coordinates_dict[str(i)][str(i) + " exon " + str(item)][0])
		end_pos = int(ENST_exon_coordinates_dict[str(i)][str(i) + " exon " + str(item)][1])
		CDS = CDS + chr22_string[start_pos:end_pos+1]
	ENST_CDS_dict[i] = CDS
# (5) Make a new file 
final_fa = open("concatenated_CDS.txt", "w")
for i in ENST_CDS_dict:
	final_fa.write(">" + i + "\n" + ENST_CDS_dict[i] + "\n")
final_fa.close()