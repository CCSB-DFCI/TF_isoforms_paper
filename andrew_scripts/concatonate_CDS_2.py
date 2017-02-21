"""
Andrew Goldfarb
02/15/17
Task: From gtf "positive_strand_CDS.gtf", which has all the + strand CDS information for chromosome 22, collect
the exon coordinates to concatenate the exon sequences from the chromosome 22 sequence in file, 
"GRCh38.p2.genome_chr22_only.fa". 
"""
# (1) Read the CDS, and collect the ENST Accessions, exon numbers, and start and end positions into lists:
ENST_accession = []
exon_number = []
start_end_list = []

positive_CDS_file = open("positive_strand_CDS.gtf", "r")
for line in positive_CDS_file:

	ENST = "ENST" + line.split("ENST")[1][:13]
	ENST_accession.append(ENST)

	exon = ENST + " exon " + line.split("exon_number ")[1].split(";")[0]
	exon_number.append(exon)

	start = int(line.split("\t")[3])
	end = int(line.split("\t")[4])
	start_end_list.append([start, end+3])
	# I added 8 and 11 because the coordinates for each CDS were off by those respective amounts.

# (2) Make a dictionary that holds the start and end coordinates(value) for each ENST exon (key).
exon_coordinates_dict = {}
for i, j in zip(exon_number, start_end_list):
	exon_coordinates_dict[i] = j

# (3) Make a nested dictionary, where each ENST (key) holds a dictionary of the each exon and corresponding
# start and end coordinates. So: d[ENST][exon#][[start, end]]
ENST_accession = set(ENST_accession)
ENST_exon_coordinates_dict = {}
for i in ENST_accession:
	temp_list_keys = []
	temp_list_values = []
	for j in exon_coordinates_dict:
		if j[:17] == i:
			temp_list_keys.append(j)
			temp_list_values.append(exon_coordinates_dict[j])
	temp_dict = {}
	for x,y in zip(temp_list_keys, temp_list_values):
		temp_dict[x] = y
	ENST_exon_coordinates_dict[i] = temp_dict

chr22_string = open("GRCh38.p2.genome_chr22_only.fa").read()
chr22_string = chr22_string[8:]
	# The first 8 indices are the header "chr22 22"
	
# (4) Concatenate the exons together in the correct order

ENST_CDS_dict = {}
for i in ENST_exon_coordinates_dict:
	ENST = str(i)
	exon_numbers_list = []
	# Want to make a list (exon_numbers_list) of all exon numbers in a transcript in order to get a 
	# to know the minimum and maximum exon number.
	a = str(ENST_exon_coordinates_dict[i]).split(ENST + " exon ")[1:]
	#print a
		#["3': ['20996886', '20997134'], '", "2': ['20996696', '20996801'], '", "1': ['20995996', '20996112']}"]
	for index in a:
		number = index.split("'")[0]
		exon_numbers_list.append(int(number))
	exon_numbers_list = set(exon_numbers_list)
	# Example: set([2, 3, 4, 5, 6, 7, 8, 9])
	# Now concatenate the exons!
	CDS = ""
	for item in range(min(exon_numbers_list), max(exon_numbers_list)+1):
		position = str(item)
		start_pos = int(ENST_exon_coordinates_dict[ENST][ENST + " exon " + position][0])
		#print start_pos
		end_pos = int(ENST_exon_coordinates_dict[ENST][ENST + " exon " + position][1])
		#print end_pos
		sequence = chr22_string[start_pos:end_pos+1]
		CDS = CDS + sequence
	ENST_CDS_dict[ENST] = CDS

final_fa = open("concatenated_CDS.txt", "w")
for i in ENST_CDS_dict:
	final_fa.write(">" + i + "\n" + ENST_CDS_dict[i] + "\n")
final_fa.close()