"""
Andrew Goldfarb
02/08/2017
Given: "gencode.v25.annotation_fewchr22lines.gtf"
	A gtf of chromosome 22. Includes positions and identifiers of genes, exons,
	transcripts, protein)
Task: Filter "gencode.v25.annotation_fewchr22lines.gtf" to a new gtf file that includes only
CDS features.
"""

annotation_file = open("gencode.v25.annotation_fewchr22lines.gtf", "r")

#Make a list, CDS_lines, that filters "gencode.v25.annotation_fewchr22lines.gtf" for CDS lines.
CDS_lines = []
for line in annotation_file:
	temp_list = line.split("\t")
	if "CDS" in temp_list:
		CDS_lines.append(line)

#Create a  file, gencode.v25.CDS.gtf, containing the CDS lines.
final_fasta = open("gencode.v25.CDS.gtf", "w")
for i in CDS_lines:
	i = i.strip()
	final_fasta.write(i + "\n")
final_fasta.close()
