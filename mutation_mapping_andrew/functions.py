"""
Andrew Goldfarb
03/01/2017
Listing all of the functions that I will use for mutational mapping
"""
#Function for translating a sequence
def translate_cds(sequence):

	"""Function: tranlsate a given nucleotide sequence into an amino acid sequence.

	Argument:
	sequence -- a nucleotide sequence

	Output:
	protein -- the translated amino acid sequence of the inputted sequence.

	For mutational mapping:
	The inputted sequence will be the concatenated CDSs for a given transcript.
	The output will be the amino acid sequence of that transcript's coding sequence."""


	genetic_code_dict = {
	    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
	    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
	    "TAT":"Y", "TAC":"Y", "TAA":None, "TAG":None,
	    "TGT":"C", "TGC":"C", "TGA":None, "TGG":"W",
	    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
	    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
	    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
	    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	   	"AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
	    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
	    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}	

	protein_seq = []
	for i in range(0, len(sequence), 3):
		codon = sequence[i:(i+3)]
		if genetic_code_dict.get(codon) != None and codon in genetic_code_dict:
			protein_seq.append(genetic_code_dict.get(codon))
			protein = "".join(protein_seq)
	return protein

#CDS = "ATGAGCACAGGCCTGCGGTACAAGAGCAAGCTGGCGACCCCAGGTGAGGACAAGCAGGTAGACATTGACAAGCAGTACGTGGGCTTCGCCACACTGCCCAACCAGGTGCACCGCAAGTCGGTGAAGAAAGGCTTTGACTTCACACTCATGGTGGCTGGTGGTGAGTCAGGCCTGGGGAAGTCCACACTGGTCCACAGCCTCTTCCTGACAGACTTGTACAAGGACCGGAAGCTGCTCAGTGCTGAGGGTGAGCGCATCAGCCAGACGGTAGAGATTCTAAAACACACGGTGGACATTGAGGAGAAGGGAGTCAAGCTGAAGCTCACCATCGTGGACACGCCGGGATTCGGGGACGCTGTCAACAACACCGAGTGGTGCTGGAAGCCCATCACCGACTATGTGGACCAGCAGTTTGAGCAGTACTTCCGTGATGAGAGCGGCCTCAACCGAAAGAACATCCAAGACAACCGAGTGCACTGCTGCCTATACTTCATCTCCCCCTTCGGGCATGGGTGGCTGCGGCCAGTGGATGTGGGTTTCATGAAGGCATTGCATGAGAAGGTCAACATCGTGCCTCTCATCGCCAAAGCTGACTGTCTTGTCCCCAGTGAGATCCGGAAGCTGAAGGAGCGGGTGATCCGGGAGGAGATTGACAAGTTTGGGATCCATGTATACCAGTTCCCTGAGTGTGACTCGGACGAGGATGAGGACTTCAAGCAGCAGGACCGGGAACTGAAGGTGGAGAGCGCGCCCTTCGCCGTTATAGGCAGCAACACGGTGGTGGAGGCCAAGGGGCAGCGGGTCCGGGGCCGACTGTACCCCTGGGGGATCGTGGAGGGTGTGGAGAACCAGGCGCATTGCGACTTCGTGAAGCTGCGCAACATGCTCATCCGCACGCATATGCACGACCTCAAGGACGTGACGTGCGACGTGCACTACGAGAACTACCGCGCGCACTGCATCCAGCAGATGACCAGGTGCAAACTGACCCAGGACAGCCGCATGGAGAGCCCCATCCCGATCCTGCCGCTGCCCACCCCGGACGCCGAGACTGAGAAGCTTATCAGGATGAAGGATGAGGAAGTACTGAGGCGCATGCAGGAGATGCTGCAGAGGATGAAGCAGCAGATGCAGGACCAGTGA"
#print translate_cds(CDS)




#Function for getting the reverse complement of a sequence
def reverse_complement(sequence):

	"""Function: generate the reverse complement for a given nucleotide sequence

	Argument:
	sequence -- a nucleotide sequence

	Output:
	reverse_seq -- the reverse complement of the inputted sequence

	For mutational mapping:
	The annotated genome "GRCh38.p7.genome.fa" provides only the positive-strand sequence.
	But CDS sequences that we want to extract are on the negative strand.
	For negative-stranded CDS's, we use the coordinates to get the positive-stranded complement,
	then use the reverse_complement function to get the sequence we are interested in."""
    

	complements = {"A":"T", "C":"G", "G":"C", "T":"A"}

	reverse_seq = ""
	for s in reversed(sequence):
		reverse_seq += complements.get(s, s)
	return reverse_seq
#CDS = "ATGAGCACAGGCCTGCGGTACAAGAGCAAGCTGGCGACCCCAGGTGAGGACAAGCAGGTAGACATTGACAAGCAGTACGTGGGCTTCGCCACACTGCCCAACCAGGTGCACCGCAAGTCGGTGAAGAAAGGCTTTGACTTCACACTCATGGTGGCTGGTGGTGAGTCAGGCCTGGGGAAGTCCACACTGGTCCACAGCCTCTTCCTGACAGACTTGTACAAGGACCGGAAGCTGCTCAGTGCTGAGGGTGAGCGCATCAGCCAGACGGTAGAGATTCTAAAACACACGGTGGACATTGAGGAGAAGGGAGTCAAGCTGAAGCTCACCATCGTGGACACGCCGGGATTCGGGGACGCTGTCAACAACACCGAGTGGTGCTGGAAGCCCATCACCGACTATGTGGACCAGCAGTTTGAGCAGTACTTCCGTGATGAGAGCGGCCTCAACCGAAAGAACATCCAAGACAACCGAGTGCACTGCTGCCTATACTTCATCTCCCCCTTCGGGCATGGGTGGCTGCGGCCAGTGGATGTGGGTTTCATGAAGGCATTGCATGAGAAGGTCAACATCGTGCCTCTCATCGCCAAAGCTGACTGTCTTGTCCCCAGTGAGATCCGGAAGCTGAAGGAGCGGGTGATCCGGGAGGAGATTGACAAGTTTGGGATCCATGTATACCAGTTCCCTGAGTGTGACTCGGACGAGGATGAGGACTTCAAGCAGCAGGACCGGGAACTGAAGGTGGAGAGCGCGCCCTTCGCCGTTATAGGCAGCAACACGGTGGTGGAGGCCAAGGGGCAGCGGGTCCGGGGCCGACTGTACCCCTGGGGGATCGTGGAGGGTGTGGAGAACCAGGCGCATTGCGACTTCGTGAAGCTGCGCAACATGCTCATCCGCACGCATATGCACGACCTCAAGGACGTGACGTGCGACGTGCACTACGAGAACTACCGCGCGCACTGCATCCAGCAGATGACCAGGTGCAAACTGACCCAGGACAGCCGCATGGAGAGCCCCATCCCGATCCTGCCGCTGCCCACCCCGGACGCCGAGACTGAGAAGCTTATCAGGATGAAGGATGAGGAAGTACTGAGGCGCATGCAGGAGATGCTGCAGAGGATGAAGCAGCAGATGCAGGACCAGTGA"
#print reverse_complement(CDS)




#Function to create and populate the data structure containing gene annotation information.
def populate_gene_annot_struc(gencode_file):

	"""Function: Create the data structure and fill it with information from "gencode.v25.annotation.gtf"

	Argument:
	gencode_file -- a file input that contains the information needed to populate the data structure.
					Function is specific for gencode format, such as "gencode.v25.annotation.gtf".

	Output:
	d -- the data structure populated with information from the input file.
		Some positions in the data structure are left unpopulated, to be populated by later functions.

	For mutational mapping:
	The file "gencode.v25.annotation.gtf" is inputted and read line-by-line to populate the data structure
	with relevant information. Data structure is organized the following way:
		d{gene:[gene_name, chromosome, {start: # , end: #}, {transcript: [strand, full_CDS_sequence, {start: # , end: #},
		{CDS: [CDS_sequence, {abs_start: #, abs_end: #}, {rel_start: #, rel_end: #}]}]}]}"""


	d = {}
	for line in open(gencode_file):
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
	return d
#print populate_gene_annot_struc("chr22_toy.gtf")




#Function to filter "HGMD_allmut.tsv" file for lines that contain missense/nonsense mutations.
def filter_for_coding_variants(mutation_file):

	"""Function: Filter the HGMD file "HGMD_allmut.tsv" for lines that include data for point/nonsense mutations only.

	Argument:
	mutation_file -- a file input that contains all mutational information associated with disease
					Function is specific to HGMD format, such as "HGMD_allmut.tsv"

	Output:
	lines_w_point_mutations -- a list of the lines from the input file that contain point/nonsense mutation data.
							Lines that contained mutations that aren't point/nonsense mutations are not included.

	For mutational mapping:
	The file "HGMD_allmut.tsv" is inputted and read line-by-line. Lines that contain point/nonsense mutations
	are added to a list. This list is the output."""


	lines_w_point_mutations = []
	for line in open(mutation_file):
		fields = line.strip().split("\t")
		if fields[7] == "NULL" and fields[8] == "NULL" and fields[6] != "NULL":
			lines_w_point_mutations.append(line)
	return lines_w_point_mutations
#filter_for_coding_variants("./data/HGMD_allmut.tsv")




#Function to create and populate the data structure containing HGMD mutation information.
#Keep lines under 80 characters!!! Use parentheses
def populate_variant_struc(lines_w_point_mutations):

	"""Function: Create and populate a data structure containing HGMD point mutation information.

	Argument:
	lines_w_point_mutations -- a list, where each index is a line from "HGMD_allmut.tsv" that has a
							missense/nonsense mutation. *** This input is the output of the 
							"filter_for_coding_variants" function!!!

	Output:
	mutations_dict -- a dictionary, containing the data structure that summarizes HGMD mutation information.

	For mutational mapping:
	First, the lines from "HGMD_allmut.tsv" containing missense/nonsense mutations are filtered into a list by
	the "filter_for_coding_variants" function, which outputs a list of these lines. In this function,
	"populate_variant_struc", the input is this list. The function goes through the list, extracting relevant
	information and organizing it into a dictionary, the output."""


	mutations_dict = {}
	for line in lines_w_point_mutations:
		fields = line.strip().split("\t")
		disease, gene, chrom, genename, gdbid, omimid, amino, deletion, insertion, codon, codonAff, descr, hgvs, hgvsAll, dbsnp, chromosome, startCoord, endCoord, tag, author, fullname, allname, vol, page, year, pmid, reftag, comments, acc_num, new_date, base = fields
		ref_nt = hgvs[-3]
		mut_nt = hgvs[-1]
		ref_AA = amino.split("-")[0]
		mut_AA = amino.split("-")[1]
		mutations_dict[gene + " " + startCoord] = [{"disease": disease, "gene": gene, "chromosome": chromosome, "coordinate":startCoord, "ref_nt": ref_nt, "mut_nt": mut_nt, "ref_AA": ref_AA, "mut_AA": mut_AA}]
	return mutations_dict
#print populate_variant_struc(filter_for_coding_variants("./data/HGMD_allmut.tsv"))["A2ML1 8851954"]



#Function: Load genome sequence data into a dictionary: key = "chr#", value = sequence
def load_reference_genome(genome_file):

	"""Function: Load genome sequence data into a dictionary, where key is "chr#" and
		value is the sequence.

		Argument:
		genome_file -- a file input in fasta format that contains the entire sequence of
					each chromosome. The format for this function is based on "GRCh38.p7.genome.fa".

		Output:
		genome_dict -- a dictionary that organizes the fasta file, where the key is "chr#" and the
					value is the sequence for that chromosome.

		For mutational mapping:
		The input file is the "GRCh38.p7.genome.fa" fasta file. Output is a dictionary containing 
		the sequence of chromosomes 1-22, X, Y, and M (mitochondrial), so 25 items total. """


	genome_dict = {}
	genome = open(genome_file).read()
	chrom_header_seq = genome.split(">")[1:26]
	for i in chrom_header_seq:
		chrom_key = i.split(" ")[0]
		i = i.split("\n")[1:]
		i = "\n".join(i)
		a = i.split("\n")
		seq = "".join(a)
		chrom_seq = seq
		genome_dict[chrom_key] = chrom_seq
	return genome_dict
#print len(load_reference_genome("GRCh38.p7.genome.fa"))


# **********GS_TODO: if time, add tester function, rather than now-commented out lines
# need to sit and discuss before staring