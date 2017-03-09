"""
Andrew Goldfarb
02/23/2017
Task: Create a data structure that organizes the "gencode.v25.annotation.gtf" file by genes, transcripts, and CDS.
"""
import pickle
import os
import functions

# sublime text editor runs from active directory in finder
# make working directory one where script is located
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Create the data structure and fill it with information from "gencode.v25.annotation.gtf"
d_pickle = 'toy_gtf_dict.p'
if not os.path.exists(d_pickle):
	d = functions.populate_gene_annot_struc("./data/gencode.v25.annotation.gtf")
	pickle.dump(d, open('toy_gtf_dict.p', 'wb'))
else:
	d = pickle.load(open('toy_gtf_dict.p'))
	print "test1"
print "test2"

# Create a dictionary, "genome_dict," where key is "chr#" and value is the sequence for that chromosome.
genome_pickle = 'genome_string.p'
if not os.path.exists(genome_pickle):
	genome_dict = functions.load_reference_genome("GRCh38.p7.genome.fa")
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
				d[g][3][t][2][CDS][0] = functions.extract_cds_sequence_from_genome(chromosome, start_coord, end_coord, genome_dict)
	pickle.dump(d, open('CDS_string.p', 'wb'))
	print "test5"
else:
	d = pickle.load(open('CDS_string.p'))
	print "test6"


#Create a dictionary, mutations_dict, that contains the mutation information from "HGMD_allmut.tsv".
mutation_pickle = 'mutation_dict.p'
if not os.path.exists(mutation_pickle):
	mutation_dict = functions.create_mutation_dictionary("./data/HGMD_allmut.tsv")
	pickle.dump(mutation_dict, open('mutation_dict.p', 'wb'))
	print "test7"
else:
	mutation_dict = pickle.load(open('mutation_dict.p'))
	print "test8"


#Write a file that contains information for mapping mutations
check_file = open("Table.txt", "w")
check_file.write("disease" + "\t" + "hg_gene_name" + "\t" + "gc25_gene_name" + "\t" + "ENSG" + "\t" + "chromosome" + "\t" + "coordinate" + "\t" + "strand" + "\t" + "ENST" + "\t" + "gc25_ref_nt" + "\t" + "hg_ref_nt" + "\t" + "hg_alt_nt" + "\t" + "gc25_ref_aa" + "\t" + "hg_ref_aa" + "\t" + "gc25_alt_aa" + "\t" + "hg_alt_aa" + "\t" + "nt_mut_relative_position" + "\t" + "aa_mut_relative_position" + "\t" + "cds_seq_alt" + "\t" + "prot_seq_alt" + "\n")
for m,val_m in mutation_dict.items():
	for g,val_g in d.items():
		if mutation_dict[m]["chromosome"] == d[g][1][3:] and int(mutation_dict[m]["coordinate"]) >= int(d[g][2]["start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][2]["end"]):
			transcripts = d[g][3]
			for t,val_t in transcripts.items():
				if int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][1]["start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][1]["end"]):
					CDSs = d[g][3][t][2]
					for CDS, val_c in CDSs.items():
						if int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][2][CDS][1]["abs_start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][2][CDS][1]["abs_end"]):
							strand = d[g][3][t][0]
							ref_full_CDS = ""
							alt_full_CDS = ""
							rel_start = 0
							rel_end = 0
							for CDS, val_c in CDSs.items():
								if strand == "+":
									if int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][2][CDS][1]["abs_start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][2][CDS][1]["abs_end"]):
										raw_cds_seq = list(d[g][3][t][2][CDS][0])
										difference = int(int(mutation_dict[m]["coordinate"]) - int(d[g][3][t][2][CDS][1]["abs_start"]))

										gencode_ref_nt = raw_cds_seq[difference]
										ref_cds_seq = "".join(raw_cds_seq)
										ref_full_CDS = ref_full_CDS + ref_cds_seq

										raw_cds_seq[difference] = mutation_dict[m]["mut_nt"].lower()
										mutated_cds_seq = "".join(raw_cds_seq)
										alt_full_CDS = alt_full_CDS + mutated_cds_seq
										
										rel_start = rel_end + 1
										rel_end = rel_start + len(mutated_cds_seq) - 1
										mutation_rel_position = (range(rel_start, rel_end+1))[difference]
									else:
										ref_full_CDS = ref_full_CDS + d[g][3][t][2][CDS][0]
										alt_full_CDS = alt_full_CDS + d[g][3][t][2][CDS][0]
										rel_start = rel_end + 1
										rel_end = rel_start + len(d[g][3][t][2][CDS][0]) - 1
								if strand == "-":
									if int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][2][CDS][1]["abs_start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][2][CDS][1]["abs_end"]):
										raw_cds_seq = list(d[g][3][t][2][CDS][0])
										difference = int(int(mutation_dict[m]["coordinate"]) - int(d[g][3][t][2][CDS][1]["abs_start"]))

										gencode_ref_nt = functions.complement(raw_cds_seq[difference])	
										ref_cds_seq = "".join(functions.reverse_complement(raw_cds_seq))
										ref_full_CDS = ref_full_CDS + ref_cds_seq

										raw_cds_seq[difference] = functions.complement(mutation_dict[m]["mut_nt"]).lower()
										mutated_cds_seq = "".join(functions.reverse_complement(raw_cds_seq))
										alt_full_CDS = alt_full_CDS + mutated_cds_seq

										rel_start = rel_end + 1
										rel_end = rel_start + len(mutated_cds_seq) - 1
										mutation_rel_position = (range(rel_start, rel_end+1))[-(difference+1)]
									else:
										ref_full_CDS = ref_full_CDS + functions.reverse_complement(d[g][3][t][2][CDS][0])
										alt_full_CDS = alt_full_CDS + functions.reverse_complement(d[g][3][t][2][CDS][0])
										rel_start = rel_end + 1
										rel_end = rel_start + len(d[g][3][t][2][CDS][0]) - 1

							ref_translated = list(functions.translate_cds(ref_full_CDS.upper()))
							alt_translated = list(functions.translate_cds(alt_full_CDS.upper()))
							AA_rel_position = int((mutation_rel_position-1)/3) + 1
							if (AA_rel_position -1) >= 0 and (AA_rel_position -1) < len(alt_translated):
								gencode_ref_AA = ref_translated[AA_rel_position-1]
								gencode_alt_AA = alt_translated[AA_rel_position-1]
								alt_translated[AA_rel_position-1] = alt_translated[AA_rel_position-1].lower()
								alt_translated = "".join(alt_translated) 

								check_file.write(mutation_dict[m]["disease"] + "\t" + mutation_dict[m]["gene"] + "\t" + d[g][0] + "\t" + g + "\t" + str(mutation_dict[m]["chromosome"])  + "\t" + str(mutation_dict[m]["coordinate"]) + "\t" + strand  + "\t" + t  + "\t" + gencode_ref_nt + "\t" + mutation_dict[m]["ref_nt"] + "\t" + mutation_dict[m]["mut_nt"] + "\t" + gencode_ref_AA + "\t" + mutation_dict[m]["ref_AA"] + "\t" + gencode_alt_AA + "\t" + mutation_dict[m]["mut_AA"] + "\t" + str(mutation_rel_position) + "\t" + str(AA_rel_position) + "\t" + alt_full_CDS + "\t" + alt_translated + "\n")
							
							#When it works how I want, rewrite into functions
							



