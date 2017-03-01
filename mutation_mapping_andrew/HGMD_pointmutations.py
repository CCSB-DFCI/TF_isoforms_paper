"""
Andrew Goldfarb
02/21/17
Filters "HGMD_allmut.tsv" to include only missense/nonsense mutations.
"""

# **********GS_TODO: convert to function 'populate_variants_struc', format, docstrings
HGMD_file = open("HGMD_allmut.tsv", "r")

HGMD_lines = []
for line in HGMD_file:
	temp_list = line.split("\t")
	if temp_list[7] == "NULL" and temp_list[8] == "NULL" and temp_list[6] != "NULL":
		HGMD_lines.append(line)
# Keeping lines that have "NULL" in insertion and deletion information, but not amino acid change position.

#Create a  file, HGMD_pointmutations, containing the CDS lines.
final_fasta = open("HGMD_pointmutations.txt", "w")
for i in HGMD_lines:
	i = i.strip()
	final_fasta.write(i + "\n")
final_fasta.close()