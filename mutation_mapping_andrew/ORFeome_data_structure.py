"""
Andrew Goldfarb
03/14/2017
Task: Create a data structure that organizes the "processed_orfeome_hg38.gtf".
PROBLEM - I'm concatenating exons, not CDSs! Translating this will be wrong!
PROBLEM??? - are the exons being concatenated in the correct order?
Do I even need to concatenate the exons and translate them? I don't know what Gloria is interested in.
	My alternative idea is to not translate or include any amino acid information. Just include the 
	sequence of the exon, with the mutation replacing the reference nt and in lowercase.
"""


import pickle
import os
import ORFeome_functions

# sublime text editor runs from active directory in finder
# make working directory one where script is located
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Create the data structure and fill it with information from "processed_orfeome_hg38.gtf"
d_pickle = 'ORFeome_gtf_dict.p'
if not os.path.exists(d_pickle):
	d = ORFeome_functions.populate_gene_annot_struc("./data/processed_orfeome_hg38.gtf")
	pickle.dump(d, open('ORFeome_gtf_dict.p', 'wb'))
else:
	d = pickle.load(open('ORFeome_gtf_dict.p'))
	print "test1"
print "test2"

# Create a dictionary, "genome_dict," where key is "chr#" and value is the sequence for that chromosome.
genome_pickle = 'ORFeome_genome_string.p'
if not os.path.exists(genome_pickle):
	genome_dict = ORFeome_functions.load_reference_genome("GRCh38.p7.genome.fa")
	pickle.dump(genome_dict, open('ORFeome_genome_string.p', 'wb'))
else:
	genome_dict = pickle.load(open('ORFeome_genome_string.p'))
	print "test3"
print "test4"


#Fill data structure with extracted "raw" CDS sequences
CDS_pickle = 'ORFeome_CDS_string.p'
if not os.path.exists(CDS_pickle):
	for t,val_t in d.items():
		chromosome = d[t][0]
		exons = d[t][3]
		for e,val_e in exons.items():
			start_coord = int(d[t][3][e][1]["start"])
			end_coord = int(d[t][3][e][1]["end"])
			d[t][3][e][0] = ORFeome_functions.extract_cds_sequence_from_genome(chromosome, start_coord, end_coord, genome_dict)
	pickle.dump(d, open('ORFeome_CDS_string.p', 'wb'))
	print "test5"
else:
	d = pickle.load(open('ORFeome_CDS_string.p'))
	print "test6"

#Create a dictionary, mutations_dict, that contains the mutation information from "HGMD_allmut.tsv".
mutation_pickle = 'ORFeome_mutation_dict.p'
if not os.path.exists(mutation_pickle):
	mutation_dict = ORFeome_functions.create_mutation_dictionary("./data/HGMD_allmut.tsv")
	pickle.dump(mutation_dict, open('ORFeome_mutation_dict.p', 'wb'))
	print "test7"
else:
	mutation_dict = pickle.load(open('ORFeome_mutation_dict.p'))
	print "test8"


#Write a file that contains information for mapping mutations
table_file = ORFeome_functions.create_table("ORFeome_table.txt")
count = 0
for m,val_m in mutation_dict.items():
	for t,val_t in d.items():
		if (mutation_dict[m]["chromosome"] == d[t][0][3:] and int(mutation_dict[m]["coordinate"]) >= 
			int(d[t][2]["start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[t][2]["end"])):
			exons = d[t][3]
			for e,val_e in exons.items():
				if (int(mutation_dict[m]["coordinate"]) >= int(d[t][3][e][1]["start"]) and 
					int(mutation_dict[m]["coordinate"]) <= int(d[t][3][e][1]["end"])):
					strand = d[t][1]
					#ref_full_CDS = ""
					#alt_full_CDS = ""
					#rel_start = 0
					#rel_end = 0
					for e, val_e in exons.items():
						if strand == "+":
							if (int(mutation_dict[m]["coordinate"]) >= int(d[t][3][e][1]["start"]) and 
								int(mutation_dict[m]["coordinate"]) <= int(d[t][3][e][1]["end"])):
								exon = e
								exon_seq = d[t][3][e][0]
								exon_start = d[t][3][e][1]["start"]
								exon_end = d[t][3][e][1]["end"]

								raw_cds_seq = list(d[t][3][e][0])
								difference = (int(int(mutation_dict[m]["coordinate"]) - 
									int(d[t][3][e][1]["start"])))
										
								ORF_ref_nt = raw_cds_seq[difference]
								#ref_full_CDS = ORFeome_functions.concatenate(ref_full_CDS, d[t][3][e][0])

								mutated_cds_seq = (ORFeome_functions.replace_ref_nt_w_alt_nt(raw_cds_seq, difference, 
									mutation_dict[m]["mut_nt"]))
								#alt_full_CDS = ORFeome_functions.concatenate(alt_full_CDS, mutated_cds_seq)
										
								#rel_start = ORFeome_functions.rel_start_func(rel_start, rel_end)
								#rel_end = ORFeome_functions.rel_end_func(rel_start, rel_end, d[t][3][e][0])
								#mutation_rel_position = ORFeome_functions.mutation_rel_position(rel_start, rel_end, difference)
										
							# else:
							# 	ref_full_CDS = ORFeome_functions.concatenate(ref_full_CDS, d[t][3][e][0])
							# 	alt_full_CDS = ORFeome_functions.concatenate(alt_full_CDS, d[t][3][e][0])
							# 	rel_start = ORFeome_functions.rel_start_func(rel_start, rel_end)
							# 	rel_end = ORFeome_functions.rel_end_func(rel_start, rel_end, d[t][3][e][0])

						if strand == "-":
							if (int(mutation_dict[m]["coordinate"]) >= int(d[t][3][e][1]["start"]) and 
								int(mutation_dict[m]["coordinate"]) <= int(d[t][3][e][1]["end"])):
								exon = e
								exon_seq = ORFeome_functions.reverse_complement(d[t][3][e][0])
								exon_start = d[t][3][e][1]["start"]
								exon_end = d[t][3][e][1]["end"]

								raw_cds_seq = list(d[t][3][e][0])
								difference = (int(int(mutation_dict[m]["coordinate"]) - 
									int(d[t][3][e][1]["start"])))

								ORF_ref_nt = ORFeome_functions.complement(raw_cds_seq[difference])
								#ref_full_CDS = ORFeome_functions.concatenate(ref_full_CDS, ORFeome_functions.reverse_complement(d[t][3][e][0]))

								mutated_cds_seq = (ORFeome_functions.reverse_complement(ORFeome_functions.replace_ref_nt_w_alt_nt(raw_cds_seq, difference, 
									ORFeome_functions.complement(mutation_dict[m]["mut_nt"]))))
								#alt_full_CDS = ORFeome_functions.concatenate(alt_full_CDS, mutated_cds_seq)

								#rel_start = ORFeome_functions.rel_start_func(rel_start, rel_end)
								#rel_end = ORFeome_functions.rel_end_func(rel_start, rel_end, d[t][3][e][0])
								#mutation_rel_position = ORFeome_functions.mutation_rel_position(rel_start, rel_end, (-(difference+1)))

							# else:
							# 	ref_full_CDS = ORFeome_functions.concatenate(ref_full_CDS, ORFeome_functions.reverse_complement(d[t][3][e][0]))
							# 	alt_full_CDS = ORFeome_functions.concatenate(alt_full_CDS, ORFeome_functions.reverse_complement(d[t][3][e][0]))
							# 	rel_start = ORFeome_functions.rel_start_func(rel_start, rel_end)
							# 	rel_end = ORFeome_functions.rel_end_func(rel_start, rel_end, d[t][3][e][0])

					# ref_translated = list(ORFeome_functions.translate_cds(ref_full_CDS.upper()))
					# alt_translated = list(ORFeome_functions.translate_cds(alt_full_CDS.upper()))
					# AA_rel_position = int((mutation_rel_position-1)/3) + 1
					# if (AA_rel_position -1) >= 0 and (AA_rel_position -1) < len(alt_translated):
					# 	ORF_ref_AA = ref_translated[AA_rel_position-1]
					# 	ORF_alt_AA = alt_translated[AA_rel_position-1]
					# 	alt_translated[AA_rel_position-1] = alt_translated[AA_rel_position-1].lower()
					# 	alt_translated = "".join(alt_translated) 

					table_file.write((mutation_dict[m]["disease"] + "\t" + mutation_dict[m]["gene"] + "\t" + 
						t + "\t" + str(mutation_dict[m]["chromosome"])  + "\t" + str(mutation_dict[m]["coordinate"]) + 
						"\t" + exon_start + "\t" + exon_end + "\t" + exon + "\t" + strand  + "\t" + ORF_ref_nt + "\t" + 
						mutation_dict[m]["ref_nt"] + "\t" + mutation_dict[m]["mut_nt"] + "\t" + mutation_dict[m]["ref_AA"] + 
						"\t" + mutation_dict[m]["mut_AA"] + "\t" + mutated_cds_seq + "\n"))
					count = count + 1
					if count%100 == 0:
						print count



						# table_file.write((mutation_dict[m]["disease"] + "\t" + mutation_dict[m]["gene"] + "\t" + 
						# 	t + "\t" + str(mutation_dict[m]["chromosome"])  + "\t" + str(mutation_dict[m]["coordinate"]) + 
						# 	"\t" + exon_start + "\t" + exon_end + "\t" + e + "\t" + strand  + "\t" + ORF_ref_nt + "\t" + 
						# 	mutation_dict[m]["ref_nt"] + "\t" + mutation_dict[m]["mut_nt"] + "\t" + ORF_ref_AA + 
						# 	"\t" + mutation_dict[m]["ref_AA"] + "\t" + ORF_alt_AA + "\t" + mutation_dict[m]["mut_AA"] + 
						# 	"\t" + str(mutation_rel_position) + "\t" + str(AA_rel_position) + "\t" + exon_seq + "\t" + 
						# 	alt_full_CDS + "\t" + alt_translated + "\n"))
						# count = count + 1
						# if count%100 == 0:
						# 	print count

