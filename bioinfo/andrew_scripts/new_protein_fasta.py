"""
Andrew Goldfarb
02/07/2017
Given: "gencode.v25.annotation_fewchr22lines.gtf"
	A gtf of chromosome 22. Includes positions and identifiers of genes, exons, 
	transcripts, protein)
Given: "gencode.v25.pc_translations.fa"
	Fasta file of all proteins in the genome
Task: Create a fasta file of all the proteins in chromosome 22.
"""

#Make a list, header_and_sequence_list, where each index contains the header and sequence for 
#a protein in the genome. Contains 94359 proteins
translations_file = open("gencode.v25.pc_translations.fa").read()
header_and_sequence_list = translations_file.split(">")[1:]

#Change the header to be just the ENSP identifier
ENSP_alone_header = []
for pair in header_and_sequence_list:
	header = pair.split("\n")[0]
	sequence = pair.split("\n")[1]
	ENSP = header.split("|")[0]
	ENSP_holder = "".join(ENSP + "\n" + sequence)
	ENSP_alone_header.append(ENSP_holder)

#Make a list, ENSP_indentifiers_chr22, where each index is the ENSP identifier for all unique
#chr 22proteins; duplicates are removed.
annotation_file = open("gencode.v25.annotation_fewchr22lines.gtf").read()
temp_list = []
temp_list = annotation_file.split("ENSP")
ENSP_identifiers_chr22 = []
for i in temp_list:
	c = "ENSP" + i[:13]
	#because the first 12 indices contain the ENSP identifier
	ENSP_identifiers_chr22.append(c)
	#Remove duplicates by putting list in a set, then back into a list
ENSP_identifiers_chr22 = ENSP_identifiers_chr22[1:]
ENSP_identifiers_chr22 = list(set(ENSP_identifiers_chr22))

#Make a list, chr22_proteins, that contains the header and sequence of all chr22 proteins
chr22_proteins = []
for i in ENSP_alone_header:
	i = i.strip()
	for j in ENSP_identifiers_chr22:
		if j.find(i[0:17]) != -1:
			#because the first 16 indices are the ENSP identifier
			chr22_proteins.append(i)
#create a  file, chr22_proteins.txt, containing the chr22 proteins in fasta format
final_fasta = open("chr22_proteins.txt", "w")
for i in chr22_proteins:
	final_fasta.write(">" + i + "\n")
final_fasta.close()