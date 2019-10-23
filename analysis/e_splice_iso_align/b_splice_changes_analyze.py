# from the output of a_*, evaluate splice changes

from collections import defaultdict

data = defaultdict(list)  # match_chain -> alignment_info
for block in open('align_patterns.txt').read().split('\n\n'):
    cat = block.split('\n')[0]
    data[cat].append(block)

with open('b_ct_of_splice_change_categories.tsv', 'w') as ofile:
	for ct, cat in sorted([(len(blocks), cat) for cat, blocks in data.items()], reverse=True):
            ofile.write(cat + '\t' + str(ct) + '\n')
