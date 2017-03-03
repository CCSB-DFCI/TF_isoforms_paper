"""
Andrew Goldfarb
02/23/2017
Task: Create a data structure that organizes the "gencode.v25.annotation.gtf" file by genes, transcripts, and CDS.
Optimized for chromosome 22 toy file. 0.5 seconds
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


# #**********GS_TODO: convert to two functions 
# #'extract_cds_sequence_from_genome' (chromosome, start_coord, end_coord, genome_dict). Output: CDS sequence
# #'concatenate_cds_set'
# #'compute_relative_coords_for_cds_set'
CDS_pickle = 'CDS_string.p'
if not os.path.exists(CDS_pickle):
	for g,val_g in d.items():
		chromosome = d[g][1]
		transcripts = d[g][3]
		for t,val_t in transcripts.items():
			#rel_start = 0
			#rel_end = 0
			#full_CDS = ""
			strand = d[g][3][t][0]
			CDSs = d[g][3][t][2]
			for CDS,val_c in CDSs.items():
				start_coord = int(d[g][3][t][2][CDS][1]["abs_start"])
				end_coord = int(d[g][3][t][2][CDS][1]["abs_end"])
				#if chromosome in genome_dict and strand == "+":
				d[g][3][t][2][CDS][0] = functions.extract_cds_sequence_from_genome(chromosome, start_coord, end_coord, genome_dict)

					#chrom_sequence = genome_dict[chromosome]
					#d[g][3][t][3][CDS][0] = chrom_sequence[start_coord-1:end_coord]
					#full_CDS = full_CDS + d[g][3][t][3][CDS][0]
					#rel_start = rel_end + 1
					#rel_end = rel_start + len(chrom_sequence[start_coord-1:end_coord]) - 1
					#d[g][3][t][3][CDS][2]["rel_start"] = rel_start
					#d[g][3][t][3][CDS][2]["rel_end"] = rel_end
				# if chromosome in genome_dict and strand == "-":
				# 	chrom_sequence = genome_dict[chromosome]
				# 	d[g][3][t][3][CDS][0] = functions.reverse_complement(chrom_sequence[start_coord-1:end_coord])
				# 	full_CDS = full_CDS + d[g][3][t][3][CDS][0]
				# 	rel_start = rel_end + 1
				# 	rel_end = rel_start + len(chrom_sequence[start_coord-1:end_coord]) - 1
				# 	d[g][3][t][3][CDS][2]["rel_start"] = rel_start
				# 	d[g][3][t][3][CDS][2]["rel_end"] = rel_end
				#d[g][3][t][1] = full_CDS
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



for m,val_m in mutation_dict.items():
	for g,val_g in d.items():
		if mutation_dict[m]["chromosome"] == d[g][1][3:] and mutation_dict[m]["coordinate"] >= d[g][2]["start"] and mutation_dict[m]["coordinate"] <= d[g][2]["end"]:
			transcripts = d[g][3]
			for t,val_t in transcripts.items():
				if mutation_dict[m]["coordinate"] >= d[g][3][t][1]["start"] and mutation_dict[m]["coordinate"] <= d[g][3][t][1]["end"]:
					strand = d[g][3][t][0]
					CDSs = d[g][3][t][2]
					for CDS,val_c in CDSs.items():
						if mutation_dict[m]["coordinate"] >= d[g][3][t][2][CDS][1]["abs_start"] and mutation_dict[m]["coordinate"] <= d[g][3][t][2][CDS][1]["abs_end"]:
							raw_cds_seq = list(d[g][3][t][2][CDS][0])
							difference = int(int(mutation_dict[m]["coordinate"]) - int(d[g][3][t][2][CDS][1]["abs_start"]))
							#Is there a better strategy to get the position of the mutation?
							#IndexError: list index out of range
							print difference
							if raw_cds_seq[difference] == mutation_dict[m]["ref_nt"] and strand == "+":
								raw_cds_seq[difference] = mutation_dict[m]["mut_nt"]
								mutated_cds_seq = "".join(raw_cds_seq)
							if raw_cds_seq[difference] == mutation_dict[m]["ref_nt"] and strand == "-":
								raw_cds_seq[difference] = functions.complement(mutation_dict[m]["mut_nt"])
								mutated_cds_seq = "".join(raw_cds_seq)
								mutated_cds_seq = functions.reverse_complement(mutated_cds_seq)





