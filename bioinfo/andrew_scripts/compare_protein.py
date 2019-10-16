#How many of the 94359 amino acid sequences are identical? How many are different?

coding_seq_file = open("translated_coding_seq.txt").read()
temp_list = coding_seq_file.split(">")[1:]
coding_seq_dict = {}
for item in temp_list:
  a = item.split("\n")[0:2]
  ENST_identifier = a[0].split("|")[0]
  sequence = a[1]
  coding_seq_dict[ENST_identifier] = sequence
# From the translated transcripts, I now have a dictionary "coding_seq_dict" of their ENSTs and sequence

translations_file = open("gencode.v25.pc_translations.fa").read()
temp_list2 = translations_file.split(">")[1:]
translated_dict = {}
X_dict = {}
for i in temp_list2:
  b = i.split("\n")[0:2]
  ENST_identifier2 = b[0].split("|")[1]
  sequence2 = b[1]
  translated_dict[ENST_identifier2] = sequence2
  # From the reported translations.fa, I now have a dictionary "translated_dict" of their ENSTs and sequence

  #Make another dictionary, sub_dict, which holds sequences that begin with an "X". 7330 sequences begin with X
  if sequence2[0] == "X":
      X_dict[ENST_identifier2] = sequence
print len(X_dict)

# Quantify how many sequences are the same, and how many are different
matches_list = []
mismatches_list = []
for ENST in coding_seq_dict:
  my_seq = coding_seq_dict[ENST]
  GC_seq = translated_dict[ENST]
  if my_seq == GC_seq:
    matches_list.append(ENST)
  else:
    mismatches_list.append(ENST)
print len(matches_list)
print len(mismatches_list)
# 86943 in matches_list. 7416 in mismatches list.
# I noticed that many of the reported protein sequences in the mismatches list begin with an "X" and have very
# different sequence. That means that 7416 - 7330 = 86 mismatches are of a different type.

#Want to isolate which mismatched sequences don't begin with "X"
true_mismatches = []
for i in mismatches_list:
  if i not in X_dict:
    true_mismatches.append(i)
print true_mismatches

# The list, true_mismatches, holds the 86 ENSTs that are mismatches and whose sequences don't begin with "X"
# 13 of these are mitochondrial genes that have a different genetic code. Example: ENST00000361381.2
# The rest have one or more "U"s within the sequence. What are the U's? Example: ENST00000529546.5
print coding_seq_dict.get("ENST00000587648.5")
print translated_dict.get("ENST00000587648.5")
