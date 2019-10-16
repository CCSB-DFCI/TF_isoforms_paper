"""
Andrew Goldfarb
02/07/2017
Given: "gencode.v25.annotation_fewchr22lines.gtf"
	A gtf of chromosome 22. Includes positions and identifiers of genes, exons, 
	transcripts, protein)
Given: "gencode.v25.pc_transcripts.fa"
	Fasta file of all proteins in the genome
Task: Create a fasta file of all the transcripts in chromosome 22.
"""
#Make a list, header_and_sequence_list, where each index contains the header and sequence for 
#a transcript in the genome. Contains 94359 transcripts.
transcripts_file = open("gencode.v25.pc_transcripts.fa").read()
header_and_sequence_list = transcripts_file.split(">")[1:]

#Change the header to be just the ENST identifier
ENST_with_sequence = []
for pair in header_and_sequence_list:
	header = pair.split("\n")[0]
	ENST = header.split("|")[0]
	sequence = pair.split("\n")[1]
	ENST_holder = "".join(ENST + "\n" + sequence)
	ENST_with_sequence.append(ENST_holder)

#Make a list, ENST_identifiers_chr22, where each index is the ENST identifier for all unique
#chr22 transcripts; duplicates are removed.
annotation_file = open("gencode.v25.annotation_fewchr22lines.gtf").read()
temp_list = annotation_file.split("ENST")
ENST_identifiers_chr22 = []
for i in temp_list:
	c = "ENST" + i[:13]
	#because the first 12 indices contain the ENSP identifier
	ENST_identifiers_chr22.append(c)

	#Remove duplicates by putting list in a set, then back into a list
ENST_identifiers_chr22 = ENST_identifiers_chr22[1:]
ENST_identifiers_chr22 = list(set(ENST_identifiers_chr22))

#Make a list, chr22_transcripts, that contains the header and sequence of all chr22 transcripts
chr22_transcripts = []
for i in ENST_with_sequence:
	i = i.strip()
	for j in ENST_identifiers_chr22:
		if j.find(i[0:17]) != -1:
			#because the first 16 indices are the ENST identifier
			chr22_transcripts.append(i)
			
#create a  file, chr22_proteins.txt, containing the chr22 proteins in fasta format
final_fasta = open("chr22_transcripts.txt", "w")
for i in chr22_transcripts:
	final_fasta.write(">" + i + "\n")
final_fasta.close()