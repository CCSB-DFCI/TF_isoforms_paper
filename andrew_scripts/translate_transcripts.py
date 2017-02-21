"""
Andrew Goldfarb
02/08/2017
Given: My previously made file, "transcripts_coding_sequence.txt", a Fasta file of each transcript's
coding sequence.
Task: Translate each coding sequence. Make a new fasta file containing each translated coding sequence.
"""

GeneticCode_dict = {
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

coding_nt_file = open("transcripts_coding_sequence.txt", "r")
header_sequence_dict = {}
proteins_with_wrong_codon = []
for line in coding_nt_file:
    if line[0] == ">":
        key = line
    if line[0] != ">":
        protein_seq = []
        sequence = line.strip()
        for i in range(0, len(sequence), 3):
            codon = sequence[i:(i+3)]
            if GeneticCode_dict.get(codon) != None and codon in GeneticCode_dict:
                protein_seq.append(GeneticCode_dict.get(codon))
        protein = "".join(protein_seq)
        header_sequence_dict[key] = protein

final_fasta = open("translated_coding_seq.txt", "w")
for i in header_sequence_dict:
    final_fasta.write(i + header_sequence_dict[i] + "\n")
final_fasta.close()