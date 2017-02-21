
# Create a new gif file containing only the + strands.
CDS_file = open("gencode.v25.CDS.gtf", "r")
positive_strand_lines = []
for line in CDS_file:
	split_line = line.split("\t")
	if "+" in split_line:
		positive_strand_lines.append(line)

final_gtf = open("positive_strand_CDS.gtf", "w")
for i in positive_strand_lines:
	i = i.strip()
	final_gtf.write(i + "\n")
final_gtf.close()


positive_CDS_string = open("pos_CDS.toy").read()
positive_CDS_file = open("pos_CDS.toy", "r")
ENST_header_list = []
for i in positive_CDS_string.split("ENST")[1:]:
	ENST_header = "ENST" + i[:13]
	ENST_header_list.append(ENST_header)
ENST_header_list = list(set(ENST_header_list))
"""
dict1 = {}
dict2 = {}
for ENST in ENST_header_list:
	print ENST
	for i in positive_CDS_file:
		print i
		if ENST in i:
			print i
			exon_number = i.split("exon_number ")[1][0]
			print exon_number
			start_position = i.split("\t")[3]
			print start_position
			end_position = i.split("\t")[4]
			print end_position
			dict2[ENST + " exon " + exon_number] = [start_position, end_position]
			print dict2
"""
dict1 = {}
dict2 = {}
a = positive_CDS_string.split("\n")
for ENST in ENST_header_list:
	#print ENST
	for i in a:
		if ENST in i:
			#print i
			exon_number = i.split("exon_number ")[1].split(";")[0]
			#print exon_number
			start_position = i.split("\t")[3]
			#print start_position
			end_position = i.split("\t")[4]
			#print end_position
			dict2[ENST + " exon " + exon_number] = [start_position, end_position]
		#dict1[ENST] = dict2
#print dict1


"""
open new CDS gtf file
list = split by "ENST"
dict1 = {}
for item in list:
	ENST_header = "ENST" + i[13]
	for line in gtf:
		if ENST_header in line:
			list = split by "exon_number "
			exon_number = int(list[0])
			list2 = split by "\t"
			starting position = int(list2[3])
			ending position = int(list2[4])
			dict2[exon_number] = [starting position, ending position]
	dict1[ENST_header] = dict2
variable = open ("chr22 file")
For ENST_header keys in dict1:
	lista = dict1[ENST_header].sort()  (want to sort by exon number)
	for i in lista:
		CDS = ""
		




Strategy:
	1. Open "gencode.v25.CDS.gtf" which already contains only CDSs
		Make a new gtf file containing only the + strand:
	
	for line in file:
		if "+" in line:
			append line to a list
	write a new file, printing each line from the list

	2. From the new gtf file with only + strands, (1) organize lines with same ENST into a common list, (2) get coordinates
	(3) concatonate into a string, (4) put the string into a dictionary:

	Open the new CDS gtf file
	list = split by "ENST"
	for item in list:
		ENST_header = "ENST" + i[:13]
		for line in gtf:
			if ENST_header in line:
				append that line to a common list
		CDS = ""
		variable = open ("chr22 file")
		for item in common list:
			list = split by "\t"
			starting position = int(list[3])
			ending position = int(list[4])
			in chr22 file, get sequence from string[starting position] to string[end position]
			CDS = CDS + string
		Dictionary[ENST_header] = CDS
	Write dictionary into a file

PROBLEM: Concatonation needs to take into account the correct order
Make nested dictionaries/lists. Can sort.
"""