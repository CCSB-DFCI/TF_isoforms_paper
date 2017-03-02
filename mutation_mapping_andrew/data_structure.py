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

#**********GS_TODO: convert to function 'populate_gene_annot_struc', add docstrings*
#*docstrings follow PEP247 format -> https://www.python.org/dev/peps/pep-0257/
d_pickle = 'toy_gtf_dict.p'
if not os.path.exists(d_pickle):
	# Create the data structure and fill it with information from "gencode.v25.annotation.gtf"
	d = {}
	for line in open("./data/gencode.v25.annotation.gtf"):
		if line[0] != "#":
			fields = line.strip().split("\t")
			chrom, annot_source, feature, start, end, dot1, strand, dot2, flag = fields

			if feature == "gene" and "protein_coding" in flag:
				words = flag.split(" ")
				ENSG = words[1][1:-2]
				gene_name = words[7][1:-2]
				if ENSG not in d:
					d[ENSG] = [gene_name, chrom, {"start": start, "end": end}, {}]

			elif feature == "transcript" and "protein_coding" in flag:
				words = flag.split(" ")
				ENST = words[3][1:-2]
				ENSG = words[1][1:-2]
				if ENSG in d:
					ensts = d[ENSG][3]
					if ENST not in ensts:
						d[ENSG][3][ENST] = [strand, "", {"start": start, "end": end}, {}]
					else:
						print "error"

			elif feature == "CDS" and "protein_coding" in flag:
				words = flag.split(" ")
				CDS_index = int(words[17][0:-1])
				ENST = words[3][1:-2]
				ENSG = words[1][1:-2]
				if ENSG in d:
					cdss = d[ENSG][3][ENST][3]
					if CDS_index not in cdss:
						d[ENSG][3][ENST][3][CDS_index] = ["", {"abs_start": start, "abs_end": end}, {}]

	pickle.dump(d, open('toy_gtf_dict.p', 'wb'))
else:
	d = pickle.load(open('toy_gtf_dict.p'))
	print "test1"

print "test2"

# **********GS_TODO: convert to function 'load_reference_genome', docstrings
genome_pickle = 'genome_string.p'
if not os.path.exists(genome_pickle):
	genome_dict = {}
	genome = open('GRCh38.p7.genome.fa').read()
	chrom_header_seq = genome.split(">")[1:26]
	for i in chrom_header_seq:
		chrom_key = i.split(" ")[0]
		i = i.split("\n")[1:]
		i = "\n".join(i)
		a = i.split("\n")
		seq = "".join(a)
		chrom_seq = seq
		genome_dict[chrom_key] = chrom_seq
		#instead of a list, organize into a dictionary. Key is "chr#", value is chromosome sequence
	pickle.dump(genome_dict, open('genome_string.p', 'wb'))
else:
	genome_dict = pickle.load(open('genome_string.p'))
	print "test3"
print "test4"
#**********GS_TODO: convert to two functions 
#'extract_cds_sequence_from_genome' (chromosome, start_coord, end_coord, genome_dict). Output: CDS sequence
#'concatenate_cds_set'
#'compute_relative_coords_for_cds_set'
CDS_pickle = 'CDS_string.p'
if not os.path.exists(CDS_pickle):
	for g,val_g in d.items():
		chromosome = val_g[1]
		transcripts = val_g[3]
		for t,val_t in transcripts.items():
			rel_start = 0
			rel_end = 0
			full_CDS = ""
			strand = val_t[0]
			CDSs = val_t[3]
			for CDS,val_c in CDSs.items():
				start_coord = int(val_c[1]["abs_start"])
				end_coord = int(val_c[1]["abs_end"])
				if chromosome in genome_dict and strand == "+":
					chrom_sequence = genome_dict[chromosome]
					val_c[0] = chrom_sequence[start_coord-1:end_coord]
					full_CDS = full_CDS + val_c[0]
					rel_start = rel_end + 1
					rel_end = rel_start + len(chrom_sequence[start_coord-1:end_coord]) - 1
					val_c[2]["rel_start"] = rel_start
					val_c[2]["rel_end"] = rel_end
				if chromosome in genome_dict and strand == "-":
					chrom_sequence = genome_dict[chromosome]
					val_c[0] = functions.reverse_complement(chrom_sequence[start_coord-1:end_coord])
					full_CDS = full_CDS + val_c[0]
					rel_start = rel_end + 1
					rel_end = rel_start + len(chrom_sequence[start_coord-1:end_coord]) - 1
					val_c[2]["rel_start"] = rel_start
					val_c[2]["rel_end"] = rel_end
				d[g][3][t][1] = full_CDS
	pickle.dump(d, open('CDS_string.p', 'wb'))
	print "test5"
else:
	d = pickle.load(open('CDS_string.p'))
	print "test6"
#print d["ENSG00000279457.3"]