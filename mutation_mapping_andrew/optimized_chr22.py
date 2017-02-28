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

d_pickle = 'toy_gtf_dict.p'
if not os.path.exists(d_pickle):
	# Create the data structure and fill it with information from "gencode.v25.annotation.gtf"
	d = {}
	# for line in open("./data/gencode.v25.annotation.gtf", "r"):
	for line in open("chr22_toy.gtf"):
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
#MAKE INTO A FUNCTION! INPUT IS THE DICTIONARY
#See which step is taking too long

genome_pickle = 'genome_string.p'
if not os.path.exists(genome_pickle):
	genome = open('GRCh38.p7.genome.fa').read()
	chrom_header_seq = genome.split(">")[1:26]
	for i in chrom_header_seq:
		if "chr22 " in i:
			chr22 = i
	pickle.dump(chr22, open('genome_string.p', 'wb'))
else:
	chr22 = pickle.load(open('genome_string.p'))
	print "test3"
 #lkajsdf
CDS_pickle = 'CDS_string_chr22.p'
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
				Chr22 = chr22
				start_coord = int(val_c[1]["abs_start"])
				end_coord = int(val_c[1]["abs_end"])
				#for i in chrom_header_seq:
					#if i.split(" ")[0] == chromosome and strand == "+":
						#i = i[7:]
						## TODO! (change! not always 7. split by \n then get first index out of list)
						#a = i.split("\n")
						#i = "".join(a)
						#val[0] = i[start_coord-1:end_coord]
			#if chr22.split(" ")[0] == chromosome and strand == "+":
				if strand == "+":
					Chr22 = Chr22.split("\n")[1:]
					Chr22 = "\n".join(Chr22)
					a = Chr22.split("\n")
					Chr22 = "".join(a)
					val_c[0] = Chr22[start_coord-1:end_coord]
					full_CDS = full_CDS + val_c[0]
					rel_start = rel_end + 1
					rel_end = rel_start + len(Chr22[start_coord-1:end_coord]) - 1
					val_c[2]["rel_start"] = rel_start
					val_c[2]["rel_end"] = rel_end
					Chr22 = chr22
			#elif chr22.split(" ")[0] == chromosome and strand == "-":
				if strand == "-":
					Chr22 = Chr22.split("\n")[1:]
					Chr22 = "\n".join(Chr22)
					a = Chr22.split("\n")
					Chr22 = "".join(a)
					val_c[0] = functions.reverse_complement(Chr22[start_coord-1:end_coord])
					full_CDS = full_CDS + val_c[0]
					rel_start = rel_end + 1
					rel_end = rel_start + len(Chr22[start_coord-1:end_coord]) - 1
					val_c[2]["rel_start"] = rel_start
					val_c[2]["rel_end"] = rel_end
					Chr22 = chr22
			d[g][3][t][1] = full_CDS
	pickle.dump(d, open('CDS_string_chr22.p', 'wb'))
	print "test4"
else:
	d = pickle.load(open('CDS_string_chr22.p'))
	print "test5"
print d