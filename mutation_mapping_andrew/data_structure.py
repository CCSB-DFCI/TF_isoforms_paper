"""
Andrew Goldfarb
02/23/2017
Task: Create a data structure that organizes the "gencode.v25.annotation.gtf" file by genes, transcripts, and CDS.
"""

d = {}
for line in open("gencode.v25.annotation.gtf", "r"):
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

#Function for getting reverse complement of a sequence
def reverse_complement(sequence):
    
    Complements = {"A":"T", "C":"G", "G":"C", "T":"A"}

    reverse_seq = ""
    for s in reversed(sequence):
        reverse_seq += Complements.get(s, s)
    return reverse_seq
#CDS = "ATGAGCACAGGCCTGCGGTACAAGAGCAAGCTGGCGACCCCAGGTGAGGACAAGCAGGTAGACATTGACAAGCAGTACGTGGGCTTCGCCACACTGCCCAACCAGGTGCACCGCAAGTCGGTGAAGAAAGGCTTTGACTTCACACTCATGGTGGCTGGTGGTGAGTCAGGCCTGGGGAAGTCCACACTGGTCCACAGCCTCTTCCTGACAGACTTGTACAAGGACCGGAAGCTGCTCAGTGCTGAGGGTGAGCGCATCAGCCAGACGGTAGAGATTCTAAAACACACGGTGGACATTGAGGAGAAGGGAGTCAAGCTGAAGCTCACCATCGTGGACACGCCGGGATTCGGGGACGCTGTCAACAACACCGAGTGGTGCTGGAAGCCCATCACCGACTATGTGGACCAGCAGTTTGAGCAGTACTTCCGTGATGAGAGCGGCCTCAACCGAAAGAACATCCAAGACAACCGAGTGCACTGCTGCCTATACTTCATCTCCCCCTTCGGGCATGGGTGGCTGCGGCCAGTGGATGTGGGTTTCATGAAGGCATTGCATGAGAAGGTCAACATCGTGCCTCTCATCGCCAAAGCTGACTGTCTTGTCCCCAGTGAGATCCGGAAGCTGAAGGAGCGGGTGATCCGGGAGGAGATTGACAAGTTTGGGATCCATGTATACCAGTTCCCTGAGTGTGACTCGGACGAGGATGAGGACTTCAAGCAGCAGGACCGGGAACTGAAGGTGGAGAGCGCGCCCTTCGCCGTTATAGGCAGCAACACGGTGGTGGAGGCCAAGGGGCAGCGGGTCCGGGGCCGACTGTACCCCTGGGGGATCGTGGAGGGTGTGGAGAACCAGGCGCATTGCGACTTCGTGAAGCTGCGCAACATGCTCATCCGCACGCATATGCACGACCTCAAGGACGTGACGTGCGACGTGCACTACGAGAACTACCGCGCGCACTGCATCCAGCAGATGACCAGGTGCAAACTGACCCAGGACAGCCGCATGGAGAGCCCCATCCCGATCCTGCCGCTGCCCACCCCGGACGCCGAGACTGAGAAGCTTATCAGGATGAAGGATGAGGAAGTACTGAGGCGCATGCAGGAGATGCTGCAGAGGATGAAGCAGCAGATGCAGGACCAGTGA"
#print reverse_complement(CDS)

#Extract CDS sequences from hg38 and add to data structure
genome = open("GRCh38.p7.genome.fa").read()
chrom_header_seq = genome.split(">")[1:26]

for g,val in d.items():
	chromosome = val[1]
	transcripts = val[3]
	for t,val in transcripts.items():
		strand = val[0]
		CDSs = val[2]
		for CDS,val in CDSs.items():
			start = int(val[1]["abs_start"])
			end = int(val[1]["abs_end"])
			for i in chrom_header_seq:
				if i.split(" ")[0] == chromosome and strand == "+":
					i = i[7:]
					val[0] = i[start:end+1]
				elif i.split(" ")[0] == chromosome and strand == "-":
					i = i[7:]
					val[0] = reverse_complement(i[start:end+1])
print d["ENSG00000188157.13"]