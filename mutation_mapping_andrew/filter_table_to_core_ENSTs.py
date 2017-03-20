"""
Andrew Goldfarb
03/14/2017
Task: Previously created "table.txt", a file outputted from "data_structure.py", which summarizes the mutational mapping
data for all point/nonsense mutations in "HGMD_allmut.tsv" in "gencode.v25.annotation.gtf". Now, using a core set of
ENSTs that are reported to be protein coding "gencode25_list_core_pc_ensts.txt". Filtering "table.txt" for lines that
include these ENSTs.
"""
#Makes a list, "core_ENSTs" of the core ENSTs
core_ENSTs = []
for line in open("gencode25_list_core_pc_ensts.txt"):
	core_ENSTs.append(line.strip())

#Makes a list, "lines_w_core_ENSTs", of lines from "table.txt" that contain an ENST within the "core_ENSTs" list.
lines_w_core_ENSTs = []
for line in open("table.txt"):
	fields = line.split("\t")
	if fields[0] != "disease":
		ENST = fields[9].split(".")[0]
		if ENST in core_ENSTs:
			lines_w_core_ENSTs.append(line)

#Writes out the lines collected in the "lines_w_core_ENSTs" list into a new file "core_table.txt"
table_file = open("core_table.txt", "w")
table_file.write(("disease" + "\t" + "hg_gene_name" + "\t" + "gc25_gene_name" + "\t" + "ENSG" + "\t" + 
		"chromosome" + "\t" + "coordinate" + "\t" + "CDS_start_coord" + "\t" + "CDS_end_coord" + "\t" + "strand" + 
		"\t" + "ENST" + "\t" + "gc25_ref_nt" + "\t" + "hg_ref_nt" + "\t" + "hg_alt_nt" + "\t" + "gc25_ref_aa" + 
		"\t" + "hg_ref_aa" + "\t" + "gc25_alt_aa" + "\t" + "hg_alt_aa" + "\t" + "nt_mut_relative_position" + "\t" 
		+ "aa_mut_relative_position" + "\t" + "CDS_seq" + "\t" +
		"cds_seq_alt" + "\t" + "prot_seq_alt" + "\n"))
for i in lines_w_core_ENSTs:
	table_file.write(i)
