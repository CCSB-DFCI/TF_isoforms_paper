# **********GS_TODO: convert to function 'filter_for_coding_variants', docstrings
d = {}
for line in open("./data/HGMD_allmut.tsv", "r"):
	fields = line.strip().split("\t")
	if fields[7] == "NULL" and fields[8] == "NULL" and fields[6] != "NULL":
		# **********GS_TODO: convert all trailing lines to <80 char with 
		# line continuation syntax
		disease, gene, chrom, genename, gdbid, omimid, amino, deletion, insertion, codon, codonAff, descr, hgvs, hgvsAll, dbsnp, chromosome, startCoord, endCoord, tag, author, fullname, allname, vol, page, year, pmid, reftag, comments, acc_num, new_date, base = fields
		ref_nt = hgvs[-3]
		mut_nt = hgvs[-1]
		ref_AA = amino.split("-")[0]
		mut_AA = amino.split("-")[1]
		d[gene + " " + startCoord] = [{"disease": disease, "gene": gene, "chromosome": chromosome, "coordinate":startCoord, "ref_nt": ref_nt, "mut_nt": mut_nt, "ref_AA": ref_AA, "mut_AA": mut_AA}]