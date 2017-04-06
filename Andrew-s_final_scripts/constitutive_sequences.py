"""
Andrew Goldfarb
03/20/2017
Task: From all the transcripts of a gene, I want to determine the coordinates of nucleotides that are constitutive,
that is, which are present in all transcripts.
"""
import pickle
import os

# sublime text editor runs from active directory in finder
# make working directory one where script is located
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Makes a list, "core_ENSTs" of the core ENSTs
core_ENSTs = []
for line in open("gencode25_list_core_pc_ensts.txt"):
	core_ENSTs.append(line.strip())

#Dictionary, d, with data structure d[ENSG][ENST][cds#: abs_start, abs_end]
d_pickle = 'CDS_coords_dict.p'
if not os.path.exists(d_pickle):
	d = {}
	count = 0
	for line in open("./data/gencode.v25.annotation.gtf"):
		if line[0] != "#":
			count = count + 1
			if count%1000 == 0:
				print count
			fields = line.strip().split("\t")
			chrom, annot_source, feature, start, end, dot1, strand, dot2, flag = fields

			if feature == "gene" and "protein_coding" in flag:
				words = flag.split(" ")
				ENSG = words[1][1:-2]
				if ENSG not in d:
					d[ENSG] = {}

			elif feature == "transcript" and "protein_coding" in flag:
				words = flag.split(" ")
				ENST = words[3][1:-2]
				ENSG = words[1][1:-2]
				if ENSG in d:
					ensts = d[ENSG]
					if ENST not in ensts and ENST.split(".")[0] in core_ENSTs:
						d[ENSG][ENST] = {}

			elif feature == "CDS" and "protein_coding" in flag:
				words = flag.split(" ")
				CDS_index = int(words[17][0:-1])
				ENST = words[3][1:-2]
				ENSG = words[1][1:-2]
				if ENSG in d and ENST.split(".")[0] in core_ENSTs:
					cdss = d[ENSG][ENST]
					if CDS_index not in cdss:
						d[ENSG][ENST][CDS_index] = {"abs_start": start, "abs_end": end}
	pickle.dump(d, open('CDS_coords_dict.p', 'wb'))
else:
	d = pickle.load(open('CDS_coords_dict.p'))

#Dictionary, constitutive_ranges_in_gene, where key is the ENSG, and value is list of ranges of coordinates that are 
#constitutive sequences in all of that gene's isoforms
constitutive_pickle = 'constitutive_dict.p'
if not os.path.exists(constitutive_pickle):
	constitutive_ranges_in_gene = {}
	count_genes = 0
	for g, val_g in d.items():
		transcripts = d[g]
		number_of_transcripts = len(transcripts)
		list_of_all_CDS_ranges = []
		constitutive_sequences_in_gene = []
		for t, val_t in transcripts.items():
			cdss = d[g][t]
			for CDS,val_c in cdss.items():
				start_coord = int(d[g][t][CDS]["abs_start"])
				end_coord = int(d[g][t][CDS]["abs_end"])
				Range = set(range(start_coord, end_coord+1))
				list_of_all_CDS_ranges.append(Range)
		for i in list_of_all_CDS_ranges:
			count = 0
			common_range = set([])
			for j in list_of_all_CDS_ranges:
				constitutive = set([])
				if len(set(i) & set(j)) > 0:
					count = count + 1
					constitutive = set(i) & set(j)
					if len(common_range) > 0 and count > 1:
						common_range = common_range & constitutive
					elif len(common_range) == 0 and count == 1:
						common_range = constitutive
			if count == number_of_transcripts and len(common_range) > 0:
				common_range = list(common_range)
				constitutive_sequences_in_gene.append(str(min(common_range)) + "-" + str(max(common_range)))
		count_genes = count_genes + 1
		print count_genes
		constitutive_sequences_in_gene = list(set(constitutive_sequences_in_gene))
		constitutive_ranges_in_gene[g] = constitutive_sequences_in_gene
	pickle.dump(constitutive_ranges_in_gene, open('constitutive_dict.p', 'wb'))
else:
	constitutive_ranges_in_gene = pickle.load(open('constitutive_dict.p'))

# This below code is specific for HGMD mutation organization!
# Create a dictionary, mutation_dict, where the key is the mutation information, and value is ENSG and mutation coordinate.
# mutation_dict = {}
# for line in open("core_table.txt"):
# 	fields = line.strip().split("\t")
# 	if fields[0] != "disease":
# 		(disease, hg_gene_name, gc25_gene_name, ENSG, chromosome, coordinate, CDS_start_coord, CDS_end_coord, 
# 			strand, ENST, gc25_ref_nt, hg_ref_nt, hg_alt_nt, gc25_ref_aa, hg_ref_aa, gc25_alt_aa, hg_alt_aa, 
# 			nt_mut_relative_position, aa_mut_relative_position, CDS_seq, cds_seq_alt, prot_seq_alt) = fields
# 		mutation_identifier = disease + " " + coordinate + " " + hg_ref_aa + "to" + hg_alt_aa
# 		if mutation_identifier not in mutation_dict:
# 			mutation_dict[mutation_identifier] = [ENSG, coordinate]

# This below code is specific for dbSNP variant organization!
# Create a dictionary, mutation_dict, where the key is the mutation information, and value is ENSG and mutation coordinate.
mutation_dict = {}
for line in open("core_dbSNP_table.txt"):
	fields = line.strip().split("\t")
	if fields[0] != "SNP_ID":
		(SNP_ID, gc25_gene_name, ENSG, chromosome, coordinate, CDS_start_coord, CDS_end_coord, strand, ENST, gc25_ref_nt, 
			dbSNP_ref_nt, dbSNP_alt_nt, gc25_ref_aa, gc25_alt_aa, nt_mut_relative_position, aa_mut_relative_position,
			CDS_seq, cds_seq_alt, prot_seq_alt) = fields
		if SNP_ID not in mutation_dict:
			mutation_dict[SNP_ID] = [ENSG, coordinate]

#Determine how many mutations are in constitutive regions and how many are not in constitutive regions.
mutation_not_in_constitutive_region = []
mutation_in_constitutive_region = []
for mutation in mutation_dict:
	boolean = False
	ENSG = mutation_dict[mutation][0]
	coordinate = int(mutation_dict[mutation][1])
	if len(constitutive_ranges_in_gene[ENSG]) > 0:
		for i in constitutive_ranges_in_gene[ENSG]:
			start = int(i.split("-")[0])
			end = int(i.split("-")[1])
			if coordinate >= start and coordinate <= end:
				boolean = True
	if boolean == False:
		mutation_not_in_constitutive_region.append(mutation)
	else:
		mutation_in_constitutive_region.append(mutation)

#print len(mutation_not_in_constitutive_region)
#print len(mutation_in_constitutive_region)
