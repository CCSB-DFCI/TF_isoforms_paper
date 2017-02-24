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
			