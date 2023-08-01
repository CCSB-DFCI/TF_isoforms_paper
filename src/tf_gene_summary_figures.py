
# coding: utf-8

# TODO:
# 
# - add AD and RD to cloned isoforms if not in reference isoform...

# In[1]:


import shutil
from pathlib import Path

from matplotlib import pyplot as plt
import pandas as pd

from data_loading import (load_isoform_and_paralog_y2h_data,
                          load_y1h_pdi_data,
                          load_m1h_activation_data,
                          load_valid_isoform_clones,
                          load_annotated_6k_collection)
from plotting import (y2h_ppi_per_tf_gene_plot,
                      y1h_pdi_per_tf_gene_plot,
                      m1h_activation_per_tf_gene_plot)


# In[2]:


shared_dir = Path('/Users/lukelambourne/Dropbox (Partners HealthCare)/TF_isoforms')
tf_webpage_dir = shared_dir / 'website'
shared_fig_dir = tf_webpage_dir / 'media'

y2h = load_isoform_and_paralog_y2h_data(add_missing_data=True)
y1h = load_y1h_pdi_data(add_missing_data=True)
m1h = load_m1h_activation_data(add_missing_data=True)
isoforms = load_valid_isoform_clones()
y2h = y2h.loc[y2h['ad_clone_acc'].isin(isoforms['clone_acc']).values, :]
y1h = y1h.loc[y1h['unique_acc'].isin(isoforms['clone_acc']).values, :]
m1h = m1h.loc[m1h['clone_acc'].isin(isoforms['clone_acc'].values), :]

tfs = load_annotated_6k_collection()


# In[12]:


with open('gene_summary_template.html', 'r') as f:
    template = f.read()
for tf in tfs.values():
    with open(tf_webpage_dir / 'pages/{}.html'.format(tf.name), 'w') as f:
        f.write(template.format(gene_name=tf.name,
                                ensembl_gene_id=tf.ensembl_gene_id,
                                uniprot_ac=tf.uniprot_ac,
                                tf_family=tf.tf_family))
shutil.copyfile('gene_summary.css', tf_webpage_dir / 'gene_summary.css')
tf_datalist = '\n'.join('        <option value="{}"/>'.format(name) for name in sorted(tfs.keys()))
with open(tf_webpage_dir / 'index.html', 'w') as f_index:
    with open('index_template.html', 'r') as f_index_template:
        f_index.write(f_index_template.read().format(tf_gene_name_list=tf_datalist))  


# In[42]:


len([orf for tf in tfs.values() for orf in tf.orfs])


# In[9]:


for tf in tfs.values():
    fig, ax = plt.subplots(1, 1)
    y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
    n_ppi_partners = ax.get_xlim()[1] + 0.5
    fig.set_size_inches(1 + 0.35 * n_ppi_partners, 1 + 0.35 * len(tf.cloned_isoforms))
    for fmt in ['.svg']:
        plt.savefig(shared_fig_dir / '{}_y2h-profile{}'.format(tf.name, fmt),
                    bbox_inches='tight')
    plt.close(plt.gcf())


# In[6]:


for tf in tfs.values():
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(6, 0.5 * len(tf.cloned_isoforms))
    m1h_activation_per_tf_gene_plot(tf.name, data=m1h)
    for fmt in ['.svg']:
        plt.savefig(shared_fig_dir / '{}_m1h-profile{}'.format(tf.name, fmt),
                    bbox_inches='tight')
    plt.close(plt.gcf())


# In[53]:


# looking at amino acid sequence differences between
# clone and gencode reference
# TODO: move somewhere else
from data_loading import load_annotated_gencode_tfs
tfs_gc = load_annotated_gencode_tfs()
from collections import Counter
count = 0
n_diff = []
for tf in tfs.values():
    for iso in tf.cloned_isoforms:
        if iso.is_novel_isoform():
            continue
        gc_aa_seq = tfs_gc[tf.name][iso.ensembl_transcript_names[0]].aa_seq
        if iso.aa_seq != gc_aa_seq:
            count += 1
            n_diff.append(sum(x != y for x, y in zip(iso.aa_seq, gc_aa_seq)))
            #print(iso.name, n_diff[-1])
print(count)
print(Counter(n_diff))


# In[3]:


# TODO: stop printing the No Y1H data available
for gene_name in isoforms['gene'].unique():
    if gene_name not in tfs:
        print('missing', gene_name)
        continue
    tf = tfs[gene_name]
    fig, ax = plt.subplots(1, 1)
    y1h_pdi_per_tf_gene_plot(gene_name, data=y1h)
    n_pdi_partners = ax.get_xlim()[1] + 0.5
    fig.set_size_inches(1 + 0.35 * n_pdi_partners, 1 + 0.35 * len(tf.cloned_isoforms))
    for fmt in ['.svg']:
        plt.savefig(shared_fig_dir / '{}_y1h-profile{}'.format(gene_name, fmt),
                    bbox_inches='tight')
    plt.close(plt.gcf())


# In[4]:


for tf in tfs.values():
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 0.65 * len(tf.cloned_isoforms))
    tf.protein_diagram(ax=ax)
    for fmt in ['.svg']:
        plt.savefig(shared_fig_dir / '{}_cloned-isoforms_protein-diagram{}'.format(tf.name, fmt),
                    bbox_inches='tight')
    plt.close(plt.gcf())


# In[19]:


for tf in tfs.values():
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 0.5 * len(tf.orfs))
    tf.exon_diagram(ax=ax, show_matched_transcripts=True)
    for fmt in ['.svg']:
        plt.savefig(shared_fig_dir / '{}_cloned-plus-ensembl-isoforms_exon-diagram{}'.format(tf.name, fmt),
                    bbox_inches='tight')
    plt.close(plt.gcf())

