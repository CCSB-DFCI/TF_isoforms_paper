"""
Andrew Goldfarb
03/16/2017
Task: (1) Determine how many unique HGMD mutations map to gencode's core ENST data set. 
	(2) Determine how many unique HGMD mutations map to the Vidal lab's ORFeome data set.
	(3) Determine how many of these HGMD mutations map (a) to both gencode's and the ORFeome set.
													   (b) to the gencode set but not ORFeome set.
													   (c) to the ORFeome set but not the gencode set.
"""




mutations_in_core_table = []
mutations_in_ORFeome = []
#Create a list, "mutations_in_core_table", that contains all unique mutations mapped to gencode's core ENSTs.
for line in open("core_table.txt"):
	fields = line.strip().split("\t")
	if fields[0] != "disease":
		(disease, hg_gene_name, gc25_gene_name, ENSG, chromosome, coordinate, CDS_start_coord, CDS_end_coord, 
			strand, ENST, gc25_ref_nt, hg_ref_nt, hg_alt_nt, gc25_ref_aa, hg_ref_aa, gc25_alt_aa, hg_alt_aa, 
			nt_mut_relative_position, aa_mut_relative_position, CDS_seq, cds_seq_alt, prot_seq_alt) = fields
		mutation_identifier = disease + " " + coordinate + " " + hg_ref_nt + "to" + hg_alt_nt
		mutations_in_core_table.append(mutation_identifier)
mutations_in_core_table = set(mutations_in_core_table)
print len(mutations_in_core_table)

#Create a list, "mutations_in_ORFeome", that contains all unique mutations mapped to the ORFeome.
for line in open("ORFeome_table.txt"):
	fields = line.strip().split("\t")
	if fields[0] != "disease":
		(disease, hg_gene_name, transcript_id, chromosome, coordinate, exon_start_coord, exon_end_coord, 
			exon_number, strand, ORF_ref_nt, hg_ref_nt, hg_alt_nt, hg_ref_aa, hg_alt_aa, exon_seq_alt) = fields
		mutation_identifier = disease + " " + coordinate + " " + hg_ref_nt + "to" + hg_alt_nt
		mutations_in_ORFeome.append(mutation_identifier)
mutations_in_ORFeome = set(mutations_in_ORFeome)
print len(mutations_in_ORFeome)

#Calculate how many mutations are mapped to gencode, but not the ORFeome.
mutations_in_core_AND_ORFeome = mutations_in_core_table & mutations_in_ORFeome
print len(mutations_in_core_AND_ORFeome)
mutations_in_core_NOT_ORFeome = mutations_in_core_table - mutations_in_ORFeome
print len(mutations_in_core_NOT_ORFeome)
mutations_in_ORFeome_NOT_core = mutations_in_ORFeome - mutations_in_core_table
print len(mutations_in_ORFeome_NOT_core)

