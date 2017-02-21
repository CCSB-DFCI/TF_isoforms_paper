"""
Andrew Goldfarb
02/08/2017
Given: gencode.v25.pc_transcripts.fa
	Fasta file of all transcripts in the genome. Sequences include UTRs and coding sequences. Coding
	sequence position specified by CDS: #-#
Task: For each transcript sequence, trim it down to just the coding sequence specified by the range
of numbers. Create a new fasta file of these trimmed transcript sequences.
"""
transcript_file = open("gencode.v25.pc_transcripts.fa", "r")
header_sequence_dict = {}
for line in transcript_file:
	if line[0] == ">":
		#example -> >ENST00000335137.3|ENSG00000186092.4|OTTHUMG00000001094.2|OTTHUMT00000003223.2|OR4F5-001|OR4F5|918|CDS:1-918|
		key = line
		temp_list = line.split("|")
		#example -> ['>ENST00000335137.3', 'ENSG00000186092.4', 'OTTHUMG00000001094.2', 'OTTHUMT00000003223.2', 'OR4F5-001', 'OR4F5', '918', 'CDS:1-918', '\n']
	for item in temp_list:
		if item[0:4] == "CDS:":
			a_list = item.split("-")
			#example -> ['CDS:1', '918']
			first_nucleotide = int(a_list[0][4:])
			last_nucleotide = int(a_list[1])
	sequence = transcript_file.next()[first_nucleotide -1 : last_nucleotide]
	header_sequence_dict[key] = sequence

final_fasta = open("transcripts_coding_sequence.txt", "w")
for i in header_sequence_dict:
	final_fasta.write(i + header_sequence_dict[i] + "\n")
final_fasta.close()