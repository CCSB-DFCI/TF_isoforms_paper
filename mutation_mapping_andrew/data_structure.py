"""
Andrew Goldfarb
02/23/2017
Task: Create a data structure that organizes the "gencode.v25.annotation.gtf" file by genes, transcripts, and CDS.
"""
#Andrew's comment
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
	for line in open("./toy.gtf"):
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
						d[ENSG][3][ENST] = [strand, {"start": start, "end": end}, {}]
					else:
						print "error"

			elif feature == "CDS" and "protein_coding" in flag:
				words = flag.split(" ")
				CDS_index = int(words[17][0:-1])
				ENST = words[3][1:-2]
				ENSG = words[1][1:-2]
				if ENSG in d:
					cdss = d[ENSG][3][ENST][2]
					if CDS_index not in cdss:
						d[ENSG][3][ENST][2][CDS_index] = ["", {"abs_start": start, "abs_end": end}, {}]

	pickle.dump(d, open('toy_gtf_dict.p', 'wb'))
else:
	d = pickle.load(open('toy_gtf_dict.p'))
	print "test1"


#CDS = "ATGAGCACAGGCCTGCGGTACAAGAGCAAGCTGGCGACCCCAGGTGAGGACAAGCAGGTAGACATTGACAAGCAGTACGTGGGCTTCGCCACACTGCCCAACCAGGTGCACCGCAAGTCGGTGAAGAAAGGCTTTGACTTCACACTCATGGTGGCTGGTGGTGAGTCAGGCCTGGGGAAGTCCACACTGGTCCACAGCCTCTTCCTGACAGACTTGTACAAGGACCGGAAGCTGCTCAGTGCTGAGGGTGAGCGCATCAGCCAGACGGTAGAGATTCTAAAACACACGGTGGACATTGAGGAGAAGGGAGTCAAGCTGAAGCTCACCATCGTGGACACGCCGGGATTCGGGGACGCTGTCAACAACACCGAGTGGTGCTGGAAGCCCATCACCGACTATGTGGACCAGCAGTTTGAGCAGTACTTCCGTGATGAGAGCGGCCTCAACCGAAAGAACATCCAAGACAACCGAGTGCACTGCTGCCTATACTTCATCTCCCCCTTCGGGCATGGGTGGCTGCGGCCAGTGGATGTGGGTTTCATGAAGGCATTGCATGAGAAGGTCAACATCGTGCCTCTCATCGCCAAAGCTGACTGTCTTGTCCCCAGTGAGATCCGGAAGCTGAAGGAGCGGGTGATCCGGGAGGAGATTGACAAGTTTGGGATCCATGTATACCAGTTCCCTGAGTGTGACTCGGACGAGGATGAGGACTTCAAGCAGCAGGACCGGGAACTGAAGGTGGAGAGCGCGCCCTTCGCCGTTATAGGCAGCAACACGGTGGTGGAGGCCAAGGGGCAGCGGGTCCGGGGCCGACTGTACCCCTGGGGGATCGTGGAGGGTGTGGAGAACCAGGCGCATTGCGACTTCGTGAAGCTGCGCAACATGCTCATCCGCACGCATATGCACGACCTCAAGGACGTGACGTGCGACGTGCACTACGAGAACTACCGCGCGCACTGCATCCAGCAGATGACCAGGTGCAAACTGACCCAGGACAGCCGCATGGAGAGCCCCATCCCGATCCTGCCGCTGCCCACCCCGGACGCCGAGACTGAGAAGCTTATCAGGATGAAGGATGAGGAAGTACTGAGGCGCATGCAGGAGATGCTGCAGAGGATGAAGCAGCAGATGCAGGACCAGTGA"
#print functions.reverse_complement(CDS)

print "test2"
#Extract CDS sequences from hg38 and add to data structure
#MAKE INTO A FUNCTION! INPUT IS THE DICTIONARY
#See which step is taking too long

genome_pickle = 'genome_string.p'
if not os.path.exists(genome_pickle):
	genome = open('GRCh38.p7.genome.fa').read()
	chrom_header_seq = genome.split(">")[1:26]
	pickle.dump(chrom_header_seq, open('genome_string.p', 'wb'))
else:
	chrom_header_seq = pickle.load(open('genome_string.p'))
	print "test3"
#TODO: make a toy of genome chromosome 22, and get 1 forward and 1 reverse from chr 22

CDS_pickle = 'CDS_string.p'
if not os.path.exists(CDS_pickle):
	for g,val in d.items():
		print g, val
		chromosome = val[1]
		transcripts = val[3]
		for t,val in transcripts.items():
			strand = val[0]
			CDSs = val[2]
			for CDS,val in CDSs.items():
				start_coord = int(val[1]["abs_start"])
				end_coord = int(val[1]["abs_end"])
				for i in chrom_header_seq:
					if i.split(" ")[0] == chromosome and strand == "+":
						i = i[7:]
						# TODO! (change! not always 7. split by \n then get first index out of list)
						a = i.split("\n")
						i = "".join(a)
						val[0] = i[start_coord-1:end_coord]
					elif i.split(" ")[0] == chromosome and strand == "-":
						i = i[7:]
						a = i.split("\n")
						i = "".join(a)
						val[0] = functions.reverse_complement(i[start_coord-1:end_coord])
else:
	d = pickle.load(open('CDS_string.p'))
	print "test4"
print d["ENSG00000164120.13"][3]["ENST00000542498.5"]