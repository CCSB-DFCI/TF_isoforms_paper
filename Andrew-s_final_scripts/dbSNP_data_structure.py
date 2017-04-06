"""
Andrew Goldfarb
03/13/2017
Task: Create a data structure that organizes the "gencode.v25.annotation.gtf" file by genes, transcripts, and CDS.
"""
import pickle
import os
import dbSNP_functions

# sublime text editor runs from active directory in finder
# make working directory one where script is located
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Create the data structure and fill it with information from "gencode.v25.annotation.gtf"
d_pickle = 'toy_gtf_dict.p'
if not os.path.exists(d_pickle):
	d = dbSNP_functions.populate_gene_annot_struc("./data/gencode.v25.annotation.gtf")
	pickle.dump(d, open('toy_gtf_dict.p', 'wb'))
else:
	d = pickle.load(open('toy_gtf_dict.p'))
	print "test1"
print "test2"

# Create a dictionary, "genome_dict," where key is "chr#" and value is the sequence for that chromosome.
genome_pickle = 'genome_string.p'
if not os.path.exists(genome_pickle):
	genome_dict = dbSNP_functions.load_reference_genome("GRCh38.p7.genome.fa")
	pickle.dump(genome_dict, open('genome_string.p', 'wb'))
else:
	genome_dict = pickle.load(open('genome_string.p'))
	print "test3"
print "test4"



#Fill data structure with extracted "raw" CDS sequences
CDS_pickle = 'CDS_string.p'
if not os.path.exists(CDS_pickle):
	for g,val_g in d.items():
		chromosome = d[g][1]
		transcripts = d[g][3]
		for t,val_t in transcripts.items():
			strand = d[g][3][t][0]
			CDSs = d[g][3][t][2]
			for CDS,val_c in CDSs.items():
				start_coord = int(d[g][3][t][2][CDS][1]["abs_start"])
				end_coord = int(d[g][3][t][2][CDS][1]["abs_end"])
				d[g][3][t][2][CDS][0] = (dbSNP_functions.extract_cds_sequence_from_genome(chromosome, start_coord, 
					end_coord, genome_dict))
	pickle.dump(d, open('CDS_string.p', 'wb'))
	print "test5"
else:
	d = pickle.load(open('CDS_string.p'))
	print "test6"


#Create a dictionary, mutations_dict, that contains the mutation information from "HGMD_allmut.tsv".
mutation_pickle = 'dbSNP_mutation_dict.p'
if not os.path.exists(mutation_pickle):
	mutation_dict = dbSNP_functions.populate_variant_struc("dbSNP_coding_variants.txt")
	pickle.dump(mutation_dict, open('dbSNP_mutation_dict.p', 'wb'))
	print "test7"
else:
	mutation_dict = pickle.load(open('dbSNP_mutation_dict.p'))
	print "test8"

#Write a file that contains information for mapping mutations
table_file = dbSNP_functions.create_table("dbSNP_table.txt")
count = 0
for m,val_m in mutation_dict.items():
	for g,val_g in d.items():
		if (mutation_dict[m]["chromosome"] == d[g][1][3:] and int(mutation_dict[m]["coordinate"]) >= 
			int(d[g][2]["start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][2]["end"])):
			transcripts = d[g][3]
			for t,val_t in transcripts.items():
				if (int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][1]["start"]) and 
					int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][1]["end"])):
					CDSs = d[g][3][t][2]
					for CDS, val_c in CDSs.items():
						if (int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][2][CDS][1]["abs_start"]) and 
							int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][2][CDS][1]["abs_end"])):
							strand = d[g][3][t][0]
							ref_full_CDS = ""
							alt_full_CDS = ""
							rel_start = 0
							rel_end = 0
							for CDS, val_c in CDSs.items():
								if strand == "+":
									if (int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][2][CDS][1]["abs_start"]) and 
										int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][2][CDS][1]["abs_end"])):
										
										CDS_seq = d[g][3][t][2][CDS][0]
										CDS_start = d[g][3][t][2][CDS][1]["abs_start"]
										CDS_end = d[g][3][t][2][CDS][1]["abs_end"]

										raw_cds_seq = list(d[g][3][t][2][CDS][0])
										difference = (int(int(mutation_dict[m]["coordinate"]) - 
											int(d[g][3][t][2][CDS][1]["abs_start"])))
										
										gencode_ref_nt = raw_cds_seq[difference]
										ref_full_CDS = dbSNP_functions.concatenate(ref_full_CDS, d[g][3][t][2][CDS][0])

										mutated_cds_seq = (dbSNP_functions.replace_ref_nt_w_alt_nt(raw_cds_seq, difference, 
											mutation_dict[m]["alt_nt"]))
										alt_full_CDS = dbSNP_functions.concatenate(alt_full_CDS, mutated_cds_seq)
										
										rel_start = dbSNP_functions.rel_start_func(rel_start, rel_end)
										rel_end = dbSNP_functions.rel_end_func(rel_start, rel_end, d[g][3][t][2][CDS][0])
										mutation_rel_position = dbSNP_functions.mutation_rel_position(rel_start, rel_end, difference)
										
									else:
										ref_full_CDS = dbSNP_functions.concatenate(ref_full_CDS, d[g][3][t][2][CDS][0])
										alt_full_CDS = dbSNP_functions.concatenate(alt_full_CDS, d[g][3][t][2][CDS][0])
										rel_start = dbSNP_functions.rel_start_func(rel_start, rel_end)
										rel_end = dbSNP_functions.rel_end_func(rel_start, rel_end, d[g][3][t][2][CDS][0])

								if strand == "-":
									if (int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][2][CDS][1]["abs_start"]) and 
										int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][2][CDS][1]["abs_end"])):
										
										CDS_seq = dbSNP_functions.reverse_complement(d[g][3][t][2][CDS][0])
										CDS_start = d[g][3][t][2][CDS][1]["abs_start"]
										CDS_end = d[g][3][t][2][CDS][1]["abs_end"]

										raw_cds_seq = list(d[g][3][t][2][CDS][0])
										difference = (int(int(mutation_dict[m]["coordinate"]) - 
											int(d[g][3][t][2][CDS][1]["abs_start"])))

										gencode_ref_nt = dbSNP_functions.complement(raw_cds_seq[difference])
										ref_full_CDS = dbSNP_functions.concatenate(ref_full_CDS, dbSNP_functions.reverse_complement(d[g][3][t][2][CDS][0]))

										mutated_cds_seq = (dbSNP_functions.reverse_complement(dbSNP_functions.replace_ref_nt_w_alt_nt(raw_cds_seq, difference, 
											mutation_dict[m]["alt_nt"])))
										alt_full_CDS = dbSNP_functions.concatenate(alt_full_CDS, mutated_cds_seq)

										rel_start = dbSNP_functions.rel_start_func(rel_start, rel_end)
										rel_end = dbSNP_functions.rel_end_func(rel_start, rel_end, d[g][3][t][2][CDS][0])
										mutation_rel_position = dbSNP_functions.mutation_rel_position(rel_start, rel_end, (-(difference+1)))

									else:
										ref_full_CDS = dbSNP_functions.concatenate(ref_full_CDS, dbSNP_functions.reverse_complement(d[g][3][t][2][CDS][0]))
										alt_full_CDS = dbSNP_functions.concatenate(alt_full_CDS, dbSNP_functions.reverse_complement(d[g][3][t][2][CDS][0]))
										rel_start = dbSNP_functions.rel_start_func(rel_start, rel_end)
										rel_end = dbSNP_functions.rel_end_func(rel_start, rel_end, d[g][3][t][2][CDS][0])

							ref_translated = list(dbSNP_functions.translate_cds(ref_full_CDS.upper()))
							alt_translated = list(dbSNP_functions.translate_cds(alt_full_CDS.upper()))
							AA_rel_position = int((mutation_rel_position-1)/3) + 1
							if (AA_rel_position -1) >= 0 and (AA_rel_position -1) < len(alt_translated):
								gencode_ref_AA = ref_translated[AA_rel_position-1]
								gencode_alt_AA = alt_translated[AA_rel_position-1]
								alt_translated[AA_rel_position-1] = alt_translated[AA_rel_position-1].lower()
								alt_translated = "".join(alt_translated) 

								if strand == "+":
									table_file.write((m + "\t" + d[g][0] + "\t" + g + "\t" + str(mutation_dict[m]["chromosome"]) + 
										"\t" + str(mutation_dict[m]["coordinate"]) + "\t" + CDS_start + "\t" + CDS_end + "\t" +
										strand + "\t" + t + "\t" + gencode_ref_nt + "\t" + mutation_dict[m]["ref_nt"] + "\t" + 
										mutation_dict[m]["alt_nt"] + "\t" + gencode_ref_AA + "\t" + gencode_alt_AA + "\t" + 
										str(mutation_rel_position) + "\t" + str(AA_rel_position) + "\t" + CDS_seq + "\t" + 
										alt_full_CDS + "\t" + alt_translated + "\n"))
								elif strand == "-":
									table_file.write((m + "\t" + d[g][0] + "\t" + g + "\t" + str(mutation_dict[m]["chromosome"]) + 
										"\t" + str(mutation_dict[m]["coordinate"]) + "\t" + CDS_start + "\t" + CDS_end + "\t" +
										strand + "\t" + t + "\t" + gencode_ref_nt + "\t" + dbSNP_functions.complement(mutation_dict[m]["ref_nt"]) + "\t" + 
										dbSNP_functions.complement(mutation_dict[m]["alt_nt"]) + "\t" + gencode_ref_AA + "\t" + gencode_alt_AA + "\t" + 
										str(mutation_rel_position) + "\t" + str(AA_rel_position) + "\t" + CDS_seq + "\t" + 
										alt_full_CDS + "\t" + alt_translated + "\n"))

								count = count + 1
								if count%100 == 0:
									print count



