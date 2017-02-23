"""
Andrew Goldfarb
02/16/17
Task: Filter "gencode.v25.annotation.gtf" so that only left with CDS lines. Final file is CDS.gtf
"""
gtf_file = open("gencode.v25.annotation.gtf", "r")
#Make a list, CDS_lines, that filters "gencode.v25.annotation.gtf" for CDS lines.
CDS_lines = []
for line in gtf_file:
	temp_list = line.split("\t")
	if "CDS" in temp_list:
		CDS_lines.append(line)
#Create a  file, gencode.v25.CDS.gtf, containing the CDS lines.
final_fasta = open("CDS.gtf", "w")
for i in CDS_lines:
	i = i.strip()
	final_fasta.write(i + "\n")
final_fasta.close()