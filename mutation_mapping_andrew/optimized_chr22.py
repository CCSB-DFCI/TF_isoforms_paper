"""
Andrew Goldfarb
02/23/2017
Task: Optimized for chromosome 22 toy. Design is complete here, with all design to fill extract CDS sequence, 
concatenate the CDS sequences, and get the relative coordinates.
Next: Reorganize so that it can be done for the entire genome!
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
	d = functions.populate_gene_annot_struc("chr22_toy.gtf")
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
#print d
#{'ENSG00000100162.14': ['CENPM', 'chr22', {'start': '41938721', 'end': '41947164'}, {'ENST00000215980.9': ['-', {'start': '41938721', 'end': '41947164'}, {1: ['CAAGATGGTGGCCGTGTTCAGGCCGGGCAGCTTGTCCAGGGGCCTCAACACCGACAT', {'abs_end': '41947076', 'abs_start': '41947020'}], 2: ['ACCTTCAGCTCGGAGGCGCAGTCCTCTTTGAGCATCGAGTCCGCCAGCTGCTGCAGAAGAGCATCCTCCGTGCCCACCAG', {'abs_end': '41946496', 'abs_start': '41946417'}], 3: ['CTGTATTTGCTGTGAAGATTAACCACAAACACGATCAGGTCAATTCGGGGCCGATTCACACTGGAGGGCAAAGGGAGGGACTTTGCCAAGTGG', {'abs_end': '41946005', 'abs_start': '41945913'}], 4: ['CACCTGTGGCGAGGAAACACACCTTCCCCAAGAAGAAGCTGGCATCCACATGGCGCAGGGACTCCTCTGTGTTCTGGAGA', {'abs_end': '41945304', 'abs_start': '41945225'}], 5: ['CTCCAGGTCACAGTAGAGCAGGGGGCTTTGATAGGTGTGGGCCAGCTTCACCACGGTGTGCCGGTGAATGCTGCAGTGGCTCTCCCGCCCAG', {'abs_end': '41943701', 'abs_start': '41943610'}], 6: ['CAGGTCCTCCAGGGAGGGGCCCTCAGAGCTTCTCAGCAGGGACAGCAGGTTCAGAGCTGAGACACCGGGCACGTGGCCAGCACAGATCTGCAGCACGCGCACCAGGCGCTGCGCCATGGTGGCCCTAAAGCCTTCCAC', {'abs_end': '41939196', 'abs_start': '41939059'}]}]}], 'ENSG00000183066.14': ['WBP2NL', 'chr22', {'start': '41998725', 'end': '42058456'}, {'ENST00000412113.5': ['+', {'start': '41998725', 'end': '42032769'}, {1: ['ATGGCGGTGAATCAGAGCCACACCGAGAACCGCCGCGGAGCCCTCATCCCTAACGGTGAAAG', {'abs_end': '41998880', 'abs_start': '41998819'}], 2: ['TCTCTTGAAGCGGTCTCCGAATGTGGAGCTCTCCTTCCCACAGCGATCAGAAGGCTCAAATGTCTTTAGTGGTAGAAAGACAGGAACATTGTTTCTCACTTCATACCGG', {'abs_end': '42019419', 'abs_start': '42019311'}], 3: ['GTGATTTTCATAACTTCATGCTCCATCAGTGATCCCATGTTGTCTTTTATGATGCCATTTGATCTGATGACGAACCTCACTGTTGAACAACCAGTATTTGCTGCAAACTTCATTAAGGGAACTATTCAGGCAGCTCCATATG', {'abs_end': '42019803', 'abs_start': '42019662'}], 4: ['GTGGCTGGGAAGGACAAGCTACTTTTAAATTAGTCTTCAGAAATGGAGATGCCATTGAATTTGCCCAGTTGATGGTGAAAGCTGCCTCTGCTG', {'abs_end': '42020096', 'abs_start': '42020004'}], 5: ['TTATTGTCTATGGAGCCCCACCTGCAGGATATGGAGCCCCACCTCCCGGATACGGAGCCCCACCTGCAGGATATGGAGCCCAACCCGTAGGAAATGAAGGCCCGCCTGTGGGATACAGAGCCTCACCTGTGCGATATGGAGCCCCACCTCTTGGATACGGAGCCCCACCTGCAGGATATGGAGCCCCACCTCTAGGATATGGAGCCCCACCTCTTGGATATGGAACCCCACCTCTCGGATATGGAGCCCCACCTCTCGGATATGGAGCCCCACCTGCAGGAAATGAAGGCCCGCCTGCGGGATACAGAGCCTCACCTGCTGGATCAGGAGCCAGGCCTCAGGAATCTACAGCAGCCCAGGCTCCTGAAAACGAGGCTTCTCTTCCCTCTGCCTCCTCTTCTCAGGTCCATTCT', {'abs_end': '42027178', 'abs_start': '42026766'}]}]}]}

#Create a dictionary, mutations_dict, that contains the mutation information from "HGMD_allmut.tsv".
mutation_pickle = 'mutation_dict.p'
if not os.path.exists(mutation_pickle):
	mutation_dict = functions.create_mutation_dictionary("./data/HGMD_allmut.tsv")
	pickle.dump(mutation_dict, open('mutation_dict.p', 'wb'))
	print "test7"
else:
	mutation_dict = pickle.load(open('mutation_dict.p'))
	print "test8"
#print mutation_dict["A2M 9072739"]
#{'ref_nt': 'C', 'disease': 'Autism', 'ref_AA': 'Arg', 'coordinate': '9072739', 'mut_AA': 'Cys', 'gene': 'A2M', 'mut_nt': 'T', 'chromosome': '12'}


for m,val_m in mutation_dict.items():
	for g,val_g in d.items():
		if mutation_dict[m]["chromosome"] == d[g][1][3:] and int(mutation_dict[m]["coordinate"]) in range(int(d[g][2]["start"]), int(d[g][2]["end"])):
		#>= int(d[g][2]["start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][2]["end"]):
			transcripts = d[g][3]
			for t,val_t in transcripts.items():
				strand = d[g][3][t][0]
				CDSs = d[g][3][t][2]
				full_CDS = ""
				rel_start = 0
				rel_end = 0
				if int(mutation_dict[m]["coordinate"]) in range(int(d[g][3][t][1]["start"]), int(d[g][3][t][1]["end"])):
				#>= int(d[g][3][t][1]["start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][1]["end"]):
					for CDS,val_c in CDSs.items():
						if strand == "+":
							if int(mutation_dict[m]["coordinate"]) not in range(int(d[g][3][t][2][CDS][1]["abs_start"]), int(d[g][3][t][2][CDS][1]["abs_end"])):
								full_CDS = full_CDS + d[g][3][t][2][CDS][0]
								rel_start = rel_end + 1
								rel_end = rel_start + len(d[g][3][t][2][CDS][0]) - 1
							else:
							#int(mutation_dict[m]["coordinate"]) >= int(d[g][3][t][2][CDS][1]["abs_start"]) and int(mutation_dict[m]["coordinate"]) <= int(d[g][3][t][2][CDS][1]["abs_end"]):
								raw_cds_seq = list(d[g][3][t][2][CDS][0])
								difference = int(int(mutation_dict[m]["coordinate"]) - int(d[g][3][t][2][CDS][1]["abs_start"]))
								if raw_cds_seq[difference] == mutation_dict[m]["ref_nt"]:
									raw_cds_seq[difference] = mutation_dict[m]["mut_nt"]
									mutated_cds_seq = "".join(raw_cds_seq)
									full_CDS = full_CDS + mutated_cds_seq
									rel_start = rel_end + 1
									rel_end = rel_start + len(d[g][3][t][2][CDS][0]) - 1
									mutation_rel_position = (range(rel_start, rel_end))[difference]
						if strand == "-":
							if int(mutation_dict[m]["coordinate"]) not in range(int(d[g][3][t][2][CDS][1]["abs_start"]), int(d[g][3][t][2][CDS][1]["abs_end"])):
								full_CDS = full_CDS + functions.reverse_complement(d[g][3][t][2][CDS][0])
								rel_start = rel_end + 1
								rel_end = rel_start + len(d[g][3][t][2][CDS][0]) - 1
							else:
								raw_cds_seq = list(d[g][3][t][2][CDS][0])
								difference = int(int(mutation_dict[m]["coordinate"]) - int(d[g][3][t][2][CDS][1]["abs_start"]))
								if raw_cds_seq[difference] == functions.reverse_complement(mutation_dict[m]["ref_nt"]):
									raw_cds_seq[difference] = functions.complement(mutation_dict[m]["mut_nt"])
									mutated_cds_seq = "".join(raw_cds_seq)
									full_CDS = full_CDS + functions.reverse_complement(mutated_cds_seq)
									rel_start = rel_end + 1
									rel_end = rel_start + len(d[g][3][t][2][CDS][0]) - 1
									mutation_rel_position = (range(rel_start, rel_end))[difference]
				protein = functions.translate_cds(full_CDS)
				check_file = open("table.txt", "w")
				check_file.write(mutation_dict[m]["disease"] + "\t" + mutation_dict[m]["gene"] + "\t" + mutation_dict[m]["chromosome"] + "\t" + mutation_dict[m]["coordinate"] + "\t" + mutation_dict[m]["ref_nt"] + "\t" + mutation_dict[m]["mut_nt"] + "\t" + g + "\t" + t + "\t" + strand + "\t" + mutation_rel_position + "\t" + full_CDS + "\t" + protein + "\n")


							