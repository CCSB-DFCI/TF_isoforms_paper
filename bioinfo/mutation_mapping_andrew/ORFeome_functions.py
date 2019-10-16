"""
Andrew Goldfarb
03/14/2017
Listing all of the functions that I will use for mutational mapping onto the ORFeome.
"""

def translate_cds(sequence):

	"""Function: tranlsate a given nucleotide sequence into an amino acid sequence.

	Argument:
	sequence -- a nucleotide sequence

	Output:
	protein -- the translated amino acid sequence of the inputted sequence.

	For mutational mapping:
	The inputted sequence will be the concatenated CDSs for a given transcript.
	The output will be the amino acid sequence of that transcript's coding sequence."""

	global protein

	genetic_code_dict = {
	    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
	    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
	    "TAT":"Y", "TAC":"Y", "TAA": "*", "TAG": "*",
	    "TGT":"C", "TGC":"C", "TGA": "*", "TGG":"W",
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
		if genetic_code_dict.get(codon) != "*" and codon in genetic_code_dict:
			protein_seq.append(genetic_code_dict.get(codon))
			protein = "".join(protein_seq)
		elif genetic_code_dict.get(codon) == "*":
			protein_seq.append("*")
			protein = "".join(protein_seq)
			break
	return protein



def reverse_complement(sequence):

	"""Function: generate the reverse complement for a given nucleotide sequence

	Argument:
	sequence -- a nucleotide sequence

	Output:
	reverse_seq -- the reverse complement of the inputted sequence

	For mutational mapping:
	The annotated genome "GRCh38.p7.genome.fa" provides only the positive-strand sequence.
	But some CDS sequences that we want to extract are on the negative strand.
	For negative-stranded CDS's, we use the coordinates to get the positive-stranded complement,
	then use the reverse_complement function to get the sequence we are interested in."""
    

	complements = {"A":"T", "C":"G", "G":"C", "T":"A", "a":"t", "c":"g", "g":"c", "t":"a"}

	reverse_seq = ""
	for s in reversed(sequence):
		reverse_seq += complements.get(s, s)
	return reverse_seq



def complement(sequence):

	"""Function: generate the complement for a given nucleotide sequence

	Argument:
	sequence -- a nucleotide sequence

	Output:
	complement_seq -- the complement of the inputted sequence

	For mutational mapping:
	The annotated genome "GRCh38.p7.genome.fa" provides only the positive-strand sequence.
	But some CDS sequences that we want to extract are on the negative strand.
	For reporting the reference nucleotide at the position of the mutation, the negative
	strand's reference nucleotide is really the complement of the nucleotide at the
	reported absolute position."""


	complements = {"A":"T", "C":"G", "G":"C", "T":"A", "a":"t", "c":"g", "g":"c", "t":"a"}
	complement_seq = ""
	for s in sequence:
		 complement_seq += complements.get(s, s)
	return complement_seq



def populate_gene_annot_struc(ORFeome_file):

	"""Function: Create the data structure and fill it with information from "processed_orfeome_hg38.gtf"

	Argument:
	ORFeome_file -- a file input that contains the information needed to populate the data structure.
					Function is specific for gencode format, such as "processed_orfeome_hg38.gtf".

	Output:
	d -- the data structure populated with information from the input file.
		Some positions in the data structure are left unpopulated, to be populated by later functions."""


	d = {}
	for line in open(ORFeome_file):
		fields = line.strip().split("\t")
		chrom, hORFEOME, feature, start, end, dot1, strand, dot2, flag = fields

		if feature == "transcript":
			exon_number = 0
			words = flag.split(" ")
			transcript_id = words[3][1:-2]
			if transcript_id not in d:
				d[transcript_id] = [chrom, strand, {"start": start, "end": end}, {}]

		elif feature == "exon":
			exon_number = exon_number + 1
			words = flag.split(" ")
			transcript_id = words[3][1:-2]
			if transcript_id in d:
				exons = d[transcript_id][3]
				exon_identifier = "exon_" + str(exon_number)
				if exon_identifier not in exons:
					d[transcript_id][3][exon_identifier] = ["", {"start": start, "end": end}]
	return d



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


def convert_3AA_to_1AA(amino_acid):

	"""Function: To convert an amino acid in the three-character notation to that in the one-character notation.

	Argument:
	amino_acid -- an amino acid in the three-character notation

	Output:
	amino_acid -- an amino acid in the one-character notation

	For mutational mapping:
	The file "HGMD_allmut.tsv" defines proteins in terms of their three-character amino acid notation, but I 
	want it in their one-character amino acid notation."""
	
	amino_acids = ({"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Cys":"C", "Glu":"E", "Gln":"Q", 
		"Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K", "Met":"M", "Phe":"F", "Pro":"P", 
		"Ser":"S", "Thr":"T", "Trp":"W", "Tyr":"Y", "Val":"V", "erm":"*"})
	#"erm" refers to "Term" for stop codon. "erm" because taking the last three indicies from descr in "HGMD_allmut.tsv"
	if amino_acid in amino_acids:
		amino_acid = amino_acids[amino_acid]
	return amino_acid



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
	the "filter_for_coding_variants" function, which outputs a list of these lines. In the current function,
	"populate_variant_struc", the input is this list. The function goes through the list, extracting relevant
	information and organizing it into a dictionary, the output."""


	mutations_dict = {}
	for line in lines_w_point_mutations:
		fields = line.strip().split("\t")
		(disease, gene, chrom, genename, gdbid, omimid, amino, deletion, insertion, codon, codonAff, descr, 
			hgvs, hgvsAll, dbsnp, chromosome, startCoord, endCoord, tag, author, fullname, allname, vol, page, 
			year, pmid, reftag, comments, acc_num, new_date, base) = fields
		ref_nt = hgvs[-3]
		mut_nt = hgvs[-1]
		#ref_AA = hgvsAll.split(" ")[-1][0]
		#mut_AA = hgvsAll.split(" ")[-1][-1]
		ref_AA = convert_3AA_to_1AA(descr[:3])
		mut_AA = convert_3AA_to_1AA(descr[-3:])
		#MAY HAVE TO CHANGE TO THESE BECAUSE NOT ALL "hgvsAll" ARE IN THE ASSUMED FORMAT!
		mutations_dict[gene + " " + startCoord] = ({"disease": disease, "gene": gene, "chromosome": chromosome, 
			"coordinate":startCoord, "ref_nt": ref_nt, "mut_nt": mut_nt, "ref_AA": ref_AA, "mut_AA": mut_AA})
	return mutations_dict



def create_mutation_dictionary(mutation_file):

	"""Function: Create and populate the HGMD data structure, given a set of filtered lines.
		Basically condenses the functions "filter_for_coding_variants" and "populate_variant_struc"
		into one.

	Argument:
		mutation_file -- a file input that contains all mutational information associated with disease
					Function is specific to HGMD format, such as "HGMD_allmut.tsv"

	Output:
	a dictionary, containing the data structure that summarizes HGMD mutation information.

	For mutational mapping:
	First, the lines from "HGMD_allmut.tsv" containing missense/nonsense mutations are filtered into a list by
	the "filter_for_coding_variants" function, which outputs a list of these lines. In the current function,
	"populate_variant_struc", the input is this list. The function goes through the list, extracting relevant
	information and organizing it into a dictionary, the output."""

	return populate_variant_struc(filter_for_coding_variants(mutation_file))



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



def extract_cds_sequence_from_genome(chromosome, start_coord, end_coord, genome_dict):

	"""Function: Extracts a sequence from GRCh38.p7.genome.fa dictionary (genome_dict) created 
			in the "load_reference_genome" function, given start and end coordinates, and the chromosome.

		Argument:
		chromosome -- the chromosome "key" for the genome_dict dictionary
		start_coord -- the start coordinate of a sequence to be extracted
		end_coord -- the end coordinate of a sequence to be extracted
		genome_dict -- a dictionary that organizes the fasta file, GRCh38.p7.genome.fa, where the key
		is "chr#" and the value is the sequence for that chromosome

		Output:
		extracted_sequence -- the extracted sequence with the defined coordinates of that chromosome

		For mutational mapping:
		This function is performed in order to fill the data structure created in the "populate_gene_annot_struc"
		function. The "raw" positive-stranded CDS sequences are extracted and put into the data structure. """

	chrom_sequence = genome_dict[chromosome]
	extracted_sequence = chrom_sequence[start_coord-1:end_coord]
	return extracted_sequence



def create_table(table_name):

	"""Function: Creates the tab-delimited table that will contain final information. Also, writes in the 
	"header" information that summarizes what is in the table. 

		Argument:
		table_name -- the name you wish to call the table. Put into quotes!

		Output:
		table_file -- the variable that the table is saved as

		For mutational mapping:
		Table header includes: (1) the disease the mutation is associated with, (2) the gene name, as reported
		by HGMD (3) the transcript identifier from "processed_orfeome_hg38.gtf", (4) the chromosome,
		from "processed_orfeome_hg38.gtf", (5) the mutation coordinate, from HGMD, (6) the start coordinate of 
		the exon where the mutation maps to, (7), the end coordinate of the exon where the mutation maps to
		(8) the strand ("+" or "-"), from "processed_orfeome_hg38.gtf", (9), ORFeome's reference nucleotide 
		at the position that HGMD reports a mutation, (10) the reference nucleotide as reported by HGMD, 
		(11) the alternative nucleotide, as reported by HGMD, (12) the reference amino acid, as translated by 
		my script, (13) the reference amino acid, reported by HGMD, (14) the alternative amino acid, as translated 
		by my script, (15) the alternative amino acid, as reported by HGMD, (16) the relative nucleotide position 
		of the mutation within the CDS, (17) the relative amino acid position of the mutation in the translated 
		CDS, (18) the sequence of the exon that the mutation maps to (19) the full-length coding sequence in 
		nucleotides, with the mutation designated in lower case, (20) the full-length translated coding sequence 
		in amino acids, with the altered amino acid designated in lower case."""

	table_file = open(table_name, "w")
	table_file.write(("disease" + "\t" + "hg_gene_name" + "\t" + "transcript_id" + "\t" + "chromosome" + "\t" + "coordinate" + 
		"\t" + "exon_start_coord" + "\t" + "exon_end_coord" + "\t" + "exon_number" + "\t" + "strand" + "\t" + "ORF_ref_nt" + 
		"\t" + "hg_ref_nt" + "\t" + "hg_alt_nt" + "\t" + "hg_ref_aa" + "\t" + "hg_alt_aa" + "\t" + "exon_seq_alt" + "\n"))
	return table_file

	

def concatenate(full_seq, seq_to_be_added):
	
	"""Function: Concatenates sequences. 

		Argument:
		full_seq -- the final sequence to be constructed
		seq_to_be_added -- the sequence to be added to construct a final sequence

		Output:
		full_seq -- the concatenated sequence

		For mutational mapping:
		This function is meant to concatenate individual CDS sequences together to construct a full-length
		coding sequence for a given transcript. """

	full_seq = full_seq + seq_to_be_added
	return full_seq



def rel_start_func(rel_start, rel_end):
		
	"""Function: To updated and return the relative start position of a string

		Argument:
		rel_start -- the relative start position
		rel_end -- the relative end position

		Output:
		rel_start -- the updated relative start position

		For mutational mapping:
		This function is performed to help get the relative position of the mutation in a full-length CDS. 
		As each CDS is concatenated the start position is 1 position after length of the previously
		concatenated sequence"""

	rel_start = rel_end + 1
	return rel_start



def rel_end_func(rel_start, rel_end, string):
		
	"""Function: To updated and return the relative end position of a string

		Argument:
		rel_start -- the relative start position
		rel_end -- the relative end position
		string -- the sequence being concatenated

		Output:
		rel_end -- the updated relative end position

		For mutational mapping:
		This function is performed to help get the relative position of the mutation in a full-length CDS. 
		As each CDS is concatenated, the end position is the previously determined relative start position,
		plus the length of the sequence that is added."""

	rel_end = rel_start + len(string) - 1
	return rel_end



def mutation_rel_position(rel_start, rel_end, difference):

	"""Function: To get the relative position of the mutation within the coding sequence.

		Argument:
		rel_start -- the relative start position
		rel_end -- the relative end position
		difference -- the difference between the mutation coordinate and relative start position. This value
			is equivalent to the index of the position of the mutation in the string sequence

		Output:
		mutation_rel_position -- a number, indicating the relative position within the CDS of where the mutation
		is positioned. This is assuming that the first index is 1, not 0.

		For mutational mapping:
		The relative position of the mutation within the CDS is printed out into the final table. We want to know this
		to help us be more confident where the mutation is when designing primers."""


	mutation_rel_position = (range(rel_start, rel_end+1))[difference]
	return mutation_rel_position
	


def replace_ref_nt_w_alt_nt(raw_cds_seq, difference, alt_nt):
	
	"""Function: To replace the reference nucleotide with the reported mutation in lower case.

		Argument:
		raw_cds_seq -- the "raw" extracted CDS from GRCh38.p7.genome.fa, with each nucleotide being an index
			in a list.
		difference -- the difference between the mutation coordinate and relative start position. This value
			is equivalent to the index of the position of the mutation in the string sequence
		alt_nt -- the nucleotide after the mutation, as reported by HGMD.

		Output:
		mutated_cds_seq -- the "raw" extracted CDS, but the index where the mutation is in is changed to
			the mutation nucleotide as reported by HGMD, but in lowercase.

		For mutational mapping:
		This functions changes the reference nucleotide to the alternative nucleotide, and makes it lowercase,
		so that it can be identified."""

	raw_cds_seq[difference] = alt_nt.lower()
	mutated_cds_seq = "".join(raw_cds_seq)
	return mutated_cds_seq

