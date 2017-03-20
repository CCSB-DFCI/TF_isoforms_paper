"""
Andrew Goldfarb
03/14/2017
Task: From the "core_table.txt" file that contains the mutational mapping information of the gencode core ENSTs,
determine how many of the mutations occur on all isoforms of a gene, and how many mutations occur on only a subset
of isoforms for a gene.
"""

#Create a dictionary, ENSG_ENST_dict, where the key is the ENSG, and the value is a list of all ENSTs in that ENSG.
ENSG_ENST_dict = {}
for line in open("./data/gencode.v25.annotation.gtf"):
	if line[0] != "#":
		fields = line.strip().split("\t")
		chrom, annot_source, feature, start, end, dot1, strand, dot2, flag = fields

		if feature == "gene" and "protein_coding" in flag:
			words = flag.split(" ")
			ENSG = words[1][1:-2]
			if ENSG not in ENSG_ENST_dict:
				ENSG_ENST_dict[ENSG] = []

		elif feature == "transcript" and "protein_coding" in flag:
			words = flag.split(" ")
			ENST = words[3][1:-2]
			ENSG = words[1][1:-2]
			if ENSG in ENSG_ENST_dict:
				ensts = ENSG_ENST_dict[ENSG]
				if ENST not in ensts:
					ensts.append(ENST)
				else:
					print "error"
#print ENSG_ENST_dict["ENSG00000068078.17"]

# Create a dictionary, mutation_ENST_dict = {}, where the key is the mutation information, and value is list of all ENSTs.
mutation_ENST_dict = {}
for line in open("core_table.txt"):
	fields = line.strip().split("\t")
	if fields[0] != "disease":
		(disease, hg_gene_name, gc25_gene_name, ENSG, chromosome, coordinate, CDS_start_coord, CDS_end_coord, 
			strand, ENST, gc25_ref_nt, hg_ref_nt, hg_alt_nt, gc25_ref_aa, hg_ref_aa, gc25_alt_aa, hg_alt_aa, 
			nt_mut_relative_position, aa_mut_relative_position, CDS_seq, cds_seq_alt, prot_seq_alt) = fields
		mutation_identifier = disease + " " + coordinate + " " + hg_ref_nt + "to" + hg_alt_nt
		if mutation_identifier not in mutation_ENST_dict:
			mutation_ENST_dict[mutation_identifier] = [ENSG, []]
		if mutation_identifier in mutation_ENST_dict:
			ensts = mutation_ENST_dict[mutation_identifier][1]
			ensts.append(ENST)
#print mutation_ENST_dict["Hypercholesterolaemia 11105535 TtoA"]


#Organizes the mutations into one of two lists: (1) mutations that target all isoforms of a gene
#												(2) mutations that target a subset of isoforms in a gene
mutations_that_target_all_isoforms = []
mutations_in_subset_of_isoforms = []
for mutation in mutation_ENST_dict:
	ENSG = mutation_ENST_dict[mutation][0]
	ENSTs_on_mutation = mutation_ENST_dict[mutation][1]
	all_ENSTs = ENSG_ENST_dict[ENSG]
	if set(ENSTs_on_mutation) == set(all_ENSTs):
		mutations_that_target_all_isoforms.append(mutation)
	else:
		mutations_in_subset_of_isoforms.append(mutation)
print len(mutations_that_target_all_isoforms)
print len(mutations_in_subset_of_isoforms)



# To count how many lines in 'core_table.txt' have different gencode/HGMD (1) reported gene names, 
# (2) reported reference nucleotides at the mutation position, (3) reported reference amino acids at the mutation
#position, and (4) reported alternative amino acids at the mutation position
table = open("core_table.txt").read()
lines = table.split("\n")[1:]
total_lines = 0
lines_w_diff_gene_names = []
lines_w_diff_ref_nt = []
lines_w_diff_ref_AA = []
lines_w_diff_alt_AA = []

for line in lines:
	total_lines = total_lines + 1
	temp_list = line.split("\t")
	if temp_list[1] != temp_list[2]:
		lines_w_diff_gene_names.append(line)
	if temp_list[10] != temp_list[11]:
		lines_w_diff_ref_nt.append(line)
	if temp_list[13] != temp_list[14]:
		lines_w_diff_ref_AA.append(line)
	if temp_list[15] != temp_list[16]:
		lines_w_diff_alt_AA.append(line)

print total_lines
print len(lines_w_diff_gene_names)
print len(lines_w_diff_ref_nt)
print len(lines_w_diff_ref_AA)
print len(lines_w_diff_alt_AA)

