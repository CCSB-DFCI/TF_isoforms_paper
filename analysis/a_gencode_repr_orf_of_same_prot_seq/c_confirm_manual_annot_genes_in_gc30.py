# book keeping on the genes in three lists - sachi, gloria, lambert
# check cases where genename is not in gencode30

import os
from isomodules import isofunc

data_dir = '/Users/gloriasheynkman/Documents/research_drive/files_ccsb/project_tf_isoforms/iso_master_take_3/data_input/'
path_gc_gtf = os.path.join(data_dir, 'gencode.v30.annotation.gtf')
path_tf_list_gloria = os.path.join(data_dir, 'human_tfs_annotation_full_and_stringent.tsv')
path_tf_list_sachi = os.path.join(data_dir, '1332_highconfidence_TFs_with_ensembl_updated_2017-08-08.txt')
path_tf_list_lambert = os.path.join(data_dir, 'Lambert_TFs_v_1.01.txt')

# %%

# get list of genenames in gencode30
gc_genes = set()
for line in open(path_gc_gtf):
    if line.startswith('#'):
        continue
    gene = line.split('gene_name "')[1].split('"')[0]
    gc_genes.add(gene)
print len(gc_genes)

# %%

tf_sachi, tf_gloria, tf_lambert = isofunc.load_tf_list(path_tf_list_sachi,
                                     path_tf_list_gloria, path_tf_list_lambert)
# 1332 in Sachi list, 1579 in Gloria list, 1639 in lambert list

# %%

# For genenames in the three tf lists (sachi, gloria, lambert), but not found
# in gencode30, I manually retrieved the updated name. Here check if the
# updated names are in gencode30.
new_gnames = []
for line in open('b_orig_genename_in_gencode30_man_annot.txt').readlines()[1:]:
    wds = line.rstrip().split('\t')
    if len(wds) == 6:
        new_gname = wds[5]
        new_gnames.append(new_gname)
for gname in new_gnames:
    if gname not in gc_genes:
        print gname
# confirmed - all manually updated genes *are* in gencode30
# *except* for the DUX* related genes that are only annotated in gloria list
