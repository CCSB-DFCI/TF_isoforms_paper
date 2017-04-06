#For a given ENSG, prints all reported disease mutations and dbSNPs, including which ENSTs each mutation maps to
import matplotlib.pyplot as plt
import pickle
import os
core_ENSTs = []
for line in open("gencode25_list_core_pc_ensts.txt"):
	core_ENSTs.append(line.strip())


d_pickle = 'variant_dict2.p'
if not os.path.exists(d_pickle):
	d = {}
	count = 0
	for line in open("./data/gencode.v25.annotation.gtf"):
		if line[0] != "#":
			count = count + 1
			if count%1000 == 0:
				print count
			fields = line.strip().split("\t")
			chrom, annot_source, feature, start, end, dot1, strand, dot2, flag = fields

			if feature == "gene" and "protein_coding" in flag:
				words = flag.split(" ")
				ENSG = words[1][1:-2]
				if ENSG not in d:
					d[ENSG] = [[], {}, {}]

			elif feature == "transcript" and "protein_coding" in flag:
				words = flag.split(" ")
				ENST = words[3][1:-2]
				ENSG = words[1][1:-2]
				if ENSG in d:
					if ENST not in d[ENSG][0] and ENST.split(".")[0] in core_ENSTs:
						d[ENSG][0].append(ENST)
	pickle.dump(d, open('variant_dict2.p', 'wb'))
else:
	d = pickle.load(open('variant_dict2.p'))

diseases = {}
for line in open("core_table.txt"):
	fields = line.strip().split("\t")
	if fields[0] != "disease":
		(disease, hg_gene_name, gc25_gene_name, ENSG, chromosome, coordinate, CDS_start_coord, CDS_end_coord, 
			strand, ENST, gc25_ref_nt, hg_ref_nt, hg_alt_nt, gc25_ref_aa, hg_ref_aa, gc25_alt_aa, hg_alt_aa, 
			nt_mut_relative_position, aa_mut_relative_position, CDS_seq, cds_seq_alt, prot_seq_alt) = fields
		mutation_identifier = disease + " " + coordinate + " " + hg_ref_nt + "to" + hg_alt_nt
		if disease not in diseases:
			diseases[disease] = []
			diseases[disease].append(ENSG)
		else:
			if ENSG not in diseases[disease]:
				diseases[disease].append(ENSG)
		if mutation_identifier not in d[ENSG][1]:
			d[ENSG][1][mutation_identifier] = []
			d[ENSG][1][mutation_identifier].append(ENST)
		else:
			if ENST not in d[ENSG][1][mutation_identifier]:
				d[ENSG][1][mutation_identifier].append(ENST)
print len(diseases) #7583 diseases
num_genes_per_disease = []
genes = []
for i in diseases:
	num_genes_per_disease.append(len(diseases[i]))
	for j in diseases[i]:
		if j not in genes:
			genes.append(j)
print len(genes)
num_genes_per_disease.sort()
print num_genes_per_disease

plt.hist(num_genes_per_disease, range=[0, 8], bins= 8, alpha=0.5)
plt.xlabel('genes per disease')
plt.ylabel('number of diseases')
plt.title('Number of genes per disease')
plt.show()
plt.savefig('genes_per_disease.pdf')

for line in open("core_dbSNP_table.txt"):
	fields = line.strip().split("\t")
	if fields[0] != "SNP_ID":
		(SNP_ID, gc25_gene_name, ENSG, chromosome, coordinate, CDS_start_coord, CDS_end_coord, strand, ENST, gc25_ref_nt, 
			dbSNP_ref_nt, dbSNP_alt_nt, gc25_ref_aa, gc25_alt_aa, nt_mut_relative_position, aa_mut_relative_position, CDS_seq,
			cds_seq_alt, prot_seq_alt) = fields
		if SNP_ID not in d[ENSG][2]:
			d[ENSG][2][SNP_ID] = []
			d[ENSG][2][SNP_ID].append(ENST)
		else:
			if ENST not in d[ENSG][2][SNP_ID]:
				d[ENSG][2][SNP_ID].append(ENST)

#print d['ENSG00000035403.16']
disease_mutations_per_gene = []
common_variants_per_gene = []
number_ENSTs_per_gene = []
for g, val_g in d.items():
	if len(d[g][1]) != 0:
		disease_mutations_per_gene.append(len(d[g][1]))
	if len(d[g][2])	!= 0:
		common_variants_per_gene.append(len(d[g][2]))
	if len(d[g][0]) != 0:
		number_ENSTs_per_gene.append(len(d[g][0]))
disease_mutations_per_gene.sort()
number_ENSTs_per_gene.sort()
#print number_ENSTs_per_gene #length is 19950, so every gene. But 169 0's, so 169 genes that have no core ENSTs.

plt.hist(number_ENSTs_per_gene, range=[0, 15], bins= 15, alpha=0.5)
plt.xlabel('isoforms per gene')
plt.ylabel('number of genes')
plt.title('Number of isoforms per gene')
plt.show()
plt.savefig('isoforms_per_gene.pdf')



# plt.bar(list_of_numbers_disease[:30], disease_proportions[:30], alpha=0.5, color='red')
# plt.bar(list_of_numbers_common[:30], common_proportions[:30], alpha=0.5, color='blue')
# plt.xlabel('mutations per gene')
# plt.ylabel('fraction of genes')
# plt.title('Number of mutations per gene')
# plt.show()
# plt.savefig('mutations_per_gene.pdf')







