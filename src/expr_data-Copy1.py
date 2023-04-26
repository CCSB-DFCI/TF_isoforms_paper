
# coding: utf-8

# ## Explore TF isoform expression data
# 
# TODO
# - try taking just 2 or 3 samples in the same data point
# - compare adult tissues in GTEx and development
# - can we do some statistics???

# In[215]:


import os
from itertools import combinations

import numpy as np
from scipy import stats
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd


# import ccsblib
# from ccsblib import ccsbplotlib as cplt

from data_loading import (load_isoform_and_paralog_y2h_data,
                          load_annotated_gencode_tfs,
                          load_y1h_pdi_data,
                          load_m1h_activation_data,
                          load_valid_isoform_clones,
                          load_seq_comparison_data,
                          load_gtex_gencode,
                          load_developmental_tissue_expression_gencode)


# In[216]:


tfs = load_annotated_gencode_tfs()

df_gtex, metadata_gtex, genes_gtex = load_gtex_gencode()

#TODO: move to data_loading
exclusion_list_gtex = {'Cells - Leukemia cell line (CML)',
                       'Cells - EBV-transformed lymphocytes',
                       'Cells - Cultured fibroblasts'}
df_gtex = df_gtex.loc[:, ~df_gtex.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex = metadata_gtex.loc[~metadata_gtex['body_site'].isin(exclusion_list_gtex), :]

means_gtex = df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean()

df_dev, metadata_dev, genes_dev = load_developmental_tissue_expression_gencode()

# TODO: move to data_loading
rename_dev_stage = {'8 week post conception,embryo': '08',
'11 week post conception,late embryo': '11',
'embryo,7 week post conception': '07',
'infant': 'infant',
'10 week post conception,late embryo': '10',
'young adult': 'young adult',
'13 week post conception,late embryo': '13',
'16 week post conception,late embryo': '16',
'4 week post conception,embryo': '04',
'neonate': 'neonate',
'19 week post conception,late embryo': '19',
'9 week post conception,late embryo': '09',
'adolescent': 'adolescent',
'5 week post conception,embryo': '05',
'embryo,6 week post conception': '06',
'12 week post conception,late embryo': '12',
'18 week post conception,late embryo': '18',
'toddler': 'toddler',
'elderly': 'elderly',
'middle adult': 'adult',
'school age child': 'child'}
metadata_dev['dev_stage'] = metadata_dev['Developmental_Stage'].map(rename_dev_stage)
means_dev = (df_dev.groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']), axis=1)
           .mean())
all_isos = {'|'.join(sorted(orf.ensembl_transcript_names))
            for tf in tfs.values() for orf in tf.orfs}
alt_isos = {'|'.join(sorted(orf.ensembl_transcript_names))
            for tf in tfs.values()
            for orf in tf.orfs
            if tf.has_MANE_select_isoform and not orf.is_MANE_select_transcript}
ref_isos = {'|'.join(sorted(orf.ensembl_transcript_names))
            for tf in tfs.values()
            for orf in tf.orfs
            if tf.has_MANE_select_isoform and orf.is_MANE_select_transcript}


# In[217]:


metadata_dev.shape


# In[218]:


metadata_gtex.shape


# In[219]:


# compare GTEx and dev datasets
df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']).value_counts().rename('# samples').value_counts().sort_index().to_frame()


# In[220]:


metadata_dev.groupby(['organism_part'])['dev_stage'].nunique().to_frame()


# In[221]:


df_gtex.columns.map(metadata_gtex['body_site']).value_counts().rename('# samples').to_frame()


# In[222]:


(means_gtex > 1).any(axis=1).value_counts()


# In[223]:


(means_gtex.loc[means_gtex.index.isin(alt_isos), :].sum(axis=1) >= 1).sum()


# In[224]:


# plot number of isoforms above 1 TPM

fig, axs = plt.subplots(2, 1, sharex=True)

n_iso = (means_gtex > 1).any(axis=1).groupby(genes_gtex).size()
x_max = n_iso.max()
xs = range(0, x_max + 1)
axs[0].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])

n_iso = (means_gtex > 1).any(axis=1).groupby(genes_gtex).sum()
axs[1].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])

axs[1].set_xticks(xs)
axs[1].set_xlabel('Unique protein isoforms per gene')
axs[0].text(x=7, y=400, s='All isoforms')
axs[1].text(x=7, y=400, s='≥ 1 TPM in at least one GTEx tissue')

def num2pct(y):
    return (y / n_iso.shape[0]) * 100

def pct2num(y):
    return (y / 100) * n_iso.shape[0]

for ax in axs:
    ax.set_ylim(0, 800)
    ax.set_yticks(range(0, 800, 100), minor=True)
    ax.set_ylabel('TF genes_gtex')
    for pos in ['top', 'right', 'bottom']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    pctax = ax.secondary_yaxis('right', functions=(num2pct, pct2num))
    pctax.set_ylabel('% of TF genes_gtex')
    pctax.set_yticks(range(0, 46, 5), minor=True)
fig.savefig('../figures/n-isoforms-per-gene_by-1TPM-cutoff_hist.pdf',
            bbox_inches='tight')


# In[225]:


# plot 2D
n_iso = (means_gtex > 1).any(axis=1).groupby(genes_gtex).sum().rename('tpm').to_frame()
n_iso['n'] = (means_gtex > 1).any(axis=1).groupby(genes_gtex).size()
xy_max = n_iso['n'].max()
pos = [(x, y) for x in range(xy_max + 1) for y in range(xy_max + 1)]
vals = [((n_iso['n'] == x) & (n_iso['tpm'] == y)).sum() for x, y in pos]
fig, ax = plt.subplots(1, 1)
ax.scatter(x=[x for x, _y in pos],
           y=[y for _x, y in pos],
           s=[v * 0.2 for v in vals])
ax.set_xticks(range(1, 26))
ax.set_yticks(range(0, 26))
ax.set_xlabel('Number of protein isoforms')
ax.set_ylabel('Number of protein isoforms\n≥ 1 TPM in at least one GTEx tissue')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
for pos in ['top', 'left']:
    ax.spines[pos].set_visible(False)
fig.savefig('../figures/n_isoforms-vs-n-gte1TPM_circle-plot.pdf',
            bbox_inches='tight')


# In[226]:


# fraction of alternative isoforms
sum(hasattr(tf.orfs[0], 'is_MANE_select_transcript') for tf in tfs.values())


# In[227]:


n_iso = (means_gtex > 1).any(axis=1).groupby(genes_gtex).sum()


# In[228]:


means_gtex.loc[~means_gtex.index.isin(all_isos), :]


# In[229]:


p = (means_gtex.loc[means_gtex.index.isin(alt_isos), :] >= 1).any(axis=1).sum()
n = means_gtex.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex.loc[means_gtex.index.isin(alt_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue')

p = (means_gtex.loc[means_gtex.index.isin(ref_isos), :] >= 1).any(axis=1).sum()
n = means_gtex.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex.loc[means_gtex.index.isin(ref_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = means_gtex.max(axis=1)

for isos, ax in zip([ref_isos, alt_isos], axs):
    ax.pie([(max_vals[max_vals.index.isin(isos)] < 1).sum(),
                ((max_vals[max_vals.index.isin(isos)] >= 1) &
                (max_vals[max_vals.index.isin(isos)] <= 5)).sum(),
                (max_vals[max_vals.index.isin(isos)] > 5).sum()
                ],
            labels=['< 1 TPM',
                    '1 - 5 TPM',
                    '> 5 TPM'],
            counterclock=False,
            startangle=90,
            autopct='%.0f%%')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
plt.subplots_adjust(wspace=0.4)
fig.savefig('../figures/GTEx-max-expression_by-reference-vs-alternative_pie.pdf',
            bbox_inches='tight')


# In[230]:


# Exlclude testis
cols = [c for c in means_gtex.columns if c != 'Testis']
p = (means_gtex.loc[means_gtex.index.isin(alt_isos), cols] >= 1).any(axis=1).sum()
n = means_gtex.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex.loc[means_gtex.index.isin(alt_isos), cols] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue')

p = (means_gtex.loc[means_gtex.index.isin(ref_isos), cols] >= 1).any(axis=1).sum()
n = means_gtex.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex.loc[means_gtex.index.isin(ref_isos), cols] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = means_gtex.loc[:, cols].max(axis=1)

for isos, ax in zip([ref_isos, alt_isos], axs):
    ax.pie([(max_vals[max_vals.index.isin(isos)] < 1).sum(),
                ((max_vals[max_vals.index.isin(isos)] >= 1) &
                (max_vals[max_vals.index.isin(isos)] <= 5)).sum(),
                (max_vals[max_vals.index.isin(isos)] > 5).sum()
                ],
            labels=['< 1 TPM',
                    '1 - 5 TPM',
                    '> 5 TPM'],
            counterclock=False,
            startangle=90,
            autopct='%.0f%%')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
plt.subplots_adjust(wspace=0.4)
fig.savefig('../figures/GTEx-max-expression_by-reference-vs-alternative_no-testis_pie.pdf',
            bbox_inches='tight')


# In[231]:


p = (means_dev.loc[means_dev.index.isin(alt_isos), :] >= 1).any(axis=1).sum()
n = means_dev.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue/dev stage')

p = (means_dev.loc[means_dev.index.isin(alt_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_dev.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue/dev stage')

p = (means_dev.loc[means_dev.index.isin(ref_isos), :] >= 1).any(axis=1).sum()
n = means_dev.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue/dev stage')

p = (means_dev.loc[means_dev.index.isin(ref_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_dev.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue/dev stage')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = means_dev.max(axis=1)

for isos, ax in zip([ref_isos, alt_isos], axs):
    ax.pie([(max_vals[max_vals.index.isin(isos)] < 1).sum(),
                ((max_vals[max_vals.index.isin(isos)] >= 1) &
                (max_vals[max_vals.index.isin(isos)] <= 5)).sum(),
                (max_vals[max_vals.index.isin(isos)] > 5).sum()
                ],
            labels=['< 1 TPM',
                    '1 - 5 TPM',
                    '> 5 TPM'],
            counterclock=False,
            startangle=90,
            autopct='%.0f%%')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
plt.subplots_adjust(wspace=0.4)
fig.savefig('../figures/developmental-max-expression_by-reference-vs-alternative_pie.pdf',
            bbox_inches='tight')


# In[232]:


# TODO: move this cell
#means_gtex_downsample = df_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean()


# In[233]:


# p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(alt_isos), :] >= 1).any(axis=1).sum()
# n = means_gtex_downsample.index.isin(alt_isos).sum()
# print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue')

# p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(alt_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
# n = means_gtex_downsample.index.isin(alt_isos).sum()
# print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue')

# p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(ref_isos), :] >= 1).any(axis=1).sum()
# n = means_gtex_downsample.index.isin(ref_isos).sum()
# print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue')

# p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(ref_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
# n = means_gtex_downsample.index.isin(ref_isos).sum()
# print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue')

# # pie chart
# fig, axs = plt.subplots(1, 2)
# max_vals = means_gtex_downsample.max(axis=1)

# for isos, ax in zip([ref_isos, alt_isos], axs):
#     ax.pie([(max_vals[max_vals.index.isin(isos)] < 1).sum(),
#                 ((max_vals[max_vals.index.isin(isos)] >= 1) &
#                 (max_vals[max_vals.index.isin(isos)] <= 5)).sum(),
#                 (max_vals[max_vals.index.isin(isos)] > 5).sum()
#                 ],
#             labels=['< 1 TPM',
#                     '1 - 5 TPM',
#                     '> 5 TPM'],
#             counterclock=False,
#             startangle=90,
#             autopct='%.0f%%')
# axs[0].set_title('MANE select isoforms')
# axs[1].set_title('Alternative isoforms')
# plt.subplots_adjust(wspace=0.4)
# fig.savefig('../figures/downsampled-GTEx-control-max-expression_by-reference-vs-alternative_pie.pdf',
#             bbox_inches='tight')


# In[234]:


# Apologies for the confusing code here
# getting the 3rd highest sample per isoform per tissue/dev timepoint

def third_highest(data):
    if data.shape[0] < 3:
        return np.nan
    return list(sorted(data, reverse=True))[2]

third_gtex = (df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
                     .apply(lambda x: x.apply(third_highest, axis=1)))
third_dev = (df_dev.groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']), axis=1)
                     .apply(lambda x: x.apply(third_highest, axis=1)))


# In[235]:


third_dev.notnull().sum(axis=1)


# In[236]:


# try requiring at least two samples in the same datapoint to reach the expression threshold
p = (third_gtex.loc[third_gtex.index.isin(alt_isos), :] >= 1).any(axis=1).sum()
n = third_gtex.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue')

p = (third_gtex.loc[third_gtex.index.isin(alt_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = third_gtex.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue')

p = (third_gtex.loc[third_gtex.index.isin(ref_isos), :] >= 1).any(axis=1).sum()
n = third_gtex.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue')

p = (third_gtex.loc[third_gtex.index.isin(ref_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = third_gtex.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = third_gtex.max(axis=1)

for isos, ax in zip([ref_isos, alt_isos], axs):
    ax.pie([(max_vals[max_vals.index.isin(isos)] < 1).sum(),
                ((max_vals[max_vals.index.isin(isos)] >= 1) &
                (max_vals[max_vals.index.isin(isos)] <= 5)).sum(),
                (max_vals[max_vals.index.isin(isos)] > 5).sum()
                ],
            labels=['< 1 TPM',
                    '1 - 5 TPM',
                    '> 5 TPM'],
            counterclock=False,
            startangle=90,
            autopct='%.0f%%')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
plt.subplots_adjust(wspace=0.4)
fig.savefig('../figures/GTEx-max-third-highest-sample-expression_by-reference-vs-alternative_pie.pdf',
            bbox_inches='tight')


# In[237]:


# try requiring at least two samples in the same datapoint to reach the expression threshold
p = (third_dev.loc[third_dev.index.isin(alt_isos), :] >= 1).any(axis=1).sum()
n = third_dev.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue')

p = (third_dev.loc[third_dev.index.isin(alt_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = third_dev.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue')

p = (third_dev.loc[third_dev.index.isin(ref_isos), :] >= 1).any(axis=1).sum()
n = third_dev.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue')

p = (third_dev.loc[third_dev.index.isin(ref_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = third_dev.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = third_dev.max(axis=1)

for isos, ax in zip([ref_isos, alt_isos], axs):
    ax.pie([(max_vals[max_vals.index.isin(isos)] < 1).sum(),
                ((max_vals[max_vals.index.isin(isos)] >= 1) &
                (max_vals[max_vals.index.isin(isos)] <= 5)).sum(),
                (max_vals[max_vals.index.isin(isos)] > 5).sum()
                ],
            labels=['< 1 TPM',
                    '1 - 5 TPM',
                    '> 5 TPM'],
            counterclock=False,
            startangle=90,
            autopct='%.0f%%')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
plt.subplots_adjust(wspace=0.4)
fig.savefig('../figures/developmental-max-third-highest-sample-expression_by-reference-vs-alternative_pie.pdf',
            bbox_inches='tight')


# In[238]:


# overlap between tissue and development
# need to deal with adult tissues in development data
from matplotlib_venn import venn2, venn3
tpm_thresholds = [1, 2, 5, 10]
fig, axs = plt.subplots(len(tpm_thresholds), 1)
fig.set_size_inches(w=2, h=2 * len(tpm_thresholds))
for tpm_threshold, ax in zip(tpm_thresholds, axs):
    a = set(means_gtex[(means_gtex.max(axis=1) >= np.log2(tpm_threshold + 1)) &
                    means_gtex.index.isin(alt_isos)].index.values)
    b = set(means_dev[(means_dev.max(axis=1) >= np.log2(tpm_threshold + 1)) &
                    means_dev.index.isin(alt_isos)].index.values)
    venn2([a, b], set_labels=['GTEx', 'Cardoso Moreira'], ax=ax)
    ax.set_title(f'TPM ≥ {tpm_threshold}')
plt.subplots_adjust(hspace=0.7)
fig.savefig('../figures/expressed-alt-isoform-overlap-by-TPM-threshold_Venn.pdf',
            bbox_inches='tight')


# In[239]:


metadata_dev['dev_stage'].value_counts()


# In[240]:


# remove adult tissues in development data
to_exclude = {'young adult', 'adult', 'elderly'}
df_dev.loc[:, df_dev.columns.map(~metadata_dev['dev_stage'].isin(to_exclude))]
means_dev_restricted = (df_dev.loc[:, df_dev.columns.map(~metadata_dev['dev_stage'].isin(to_exclude))].groupby(df_dev.loc[:, df_dev.columns.map(~metadata_dev['dev_stage'].isin(to_exclude))].columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']), axis=1)
           .mean())

tpm_thresholds = [1, 2, 5, 10]
fig, axs = plt.subplots(len(tpm_thresholds), 1)
fig.set_size_inches(w=2, h=2 * len(tpm_thresholds))
for tpm_threshold, ax in zip(tpm_thresholds, axs):
    a = set(means_gtex[(means_gtex.max(axis=1) >= np.log2(tpm_threshold + 1)) &
                    means_gtex.index.isin(alt_isos)].index.values)
    b = set(means_dev_restricted[(means_dev_restricted.max(axis=1) >= np.log2(tpm_threshold + 1)) &
                    means_dev_restricted.index.isin(alt_isos)].index.values)
    venn2([a, b], set_labels=['GTEx', 'Cardoso Moreira\nexcluding adult samples'], ax=ax)
    ax.set_title(f'TPM ≥ {tpm_threshold}')
plt.subplots_adjust(hspace=0.7)
fig.savefig('../figures/expressed-alt-isoform-overlap-restricted-by-TPM-threshold_Venn.pdf',
            bbox_inches='tight')


# In[241]:


(df_gtex.shape, df_dev.shape)


# In[242]:


n_samples_gtex = df_gtex.shape[1]
n_samples_gt1_gtex = (df_gtex >= 1).sum(axis=1)
n_samples_dev = df_dev.shape[1]
n_samples_gt1_dev = (df_dev >= 1).sum(axis=1)

fig, ax = plt.subplots(1, 1)
ax.plot([i / n_samples_gtex * 100 for i in range(1, n_samples_gtex + 1)],
        [(n_samples_gt1_gtex >= i).sum() for i in range(1, n_samples_gtex + 1)],
        label='GTEx')
ax.plot([i / n_samples_dev * 100 for i in range(1, n_samples_dev + 1)],
        [(n_samples_gt1_dev >= i).sum() for i in range(1, n_samples_dev + 1)],
        label='Cardoso Moreira')
ax.legend()
ax.set_ylabel('Isoforms with ≥ 1 TPM')
ax.set_xlabel('% of samples')
fig.savefig('../figures/n-isoforms-vs-pct-samples_GTEx-vs-development_line-plot.pdf',
            bbox_inches='tight')


# In[243]:


# is this difference just individual level variation?
(df_gtex.max(axis=1) > 1).value_counts()


# In[244]:


# is it possible to partition tissue and developmental stage specific isoforms?
# can I do this within the cardoso morosia data?


# In[245]:


(df_dev.max(axis=1) > 1).value_counts()


# In[246]:


# plot distribution of isoforms by TPM
fig, axs = plt.subplots(2, 1, sharex=True)
n_bins = 110
x_max = 11
axs[1].hist(means_gtex.max(axis=1)[means_gtex.index.isin(alt_isos)],
            bins=n_bins,
            range=(0, x_max))
axs[0].hist(means_gtex.max(axis=1)[means_gtex.index.isin(ref_isos)],
            bins=n_bins,
            range=(0, x_max))
axs[0].text(x=8, y=30, s='MANE select isoforms')
axs[1].text(x=8, y=120, s='Alternative isoforms')
for ax in axs:
    ax.set_ylabel('Number of isoforms')
    for loc in ['top', 'right']:
        ax.spines[loc].set_visible(False)
axs[1].set_xlabel('Expression, max across GTEx tissues – TPM')
x_ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
axs[1].set_xticks([np.log2(x + 1) for x in x_ticks])
axs[1].set_xticklabels(x_ticks)
fig.savefig('../figures/expression_GTEX_GENCODE-isoforms_by-reference-vs-alternative.pdf',
            bbox_inches='tight')


# In[247]:


# # plot distribution of isoforms by TPM
# fig, axs = plt.subplots(2, 1, sharex=True)
# n_bins = 110
# x_max = 11
# axs[1].hist(means_gtex_downsample.max(axis=1)[means_gtex_downsample.index.isin(alt_isos)],
#             bins=n_bins,
#             range=(0, x_max))
# axs[0].hist(means_gtex_downsample.max(axis=1)[means_gtex_downsample.index.isin(ref_isos)],
#             bins=n_bins,
#             range=(0, x_max))
# axs[0].text(x=8, y=axs[0].get_ylim()[1] * 0.95, s='MANE select isoforms')
# axs[1].text(x=8, y=axs[1].get_ylim()[1] * 0.95, s='Alternative isoforms')
# for ax in axs:
#     ax.set_ylabel('Number of isoforms')
#     for loc in ['top', 'right']:
#         ax.spines[loc].set_visible(False)
# axs[1].set_xlabel('Expression, max across GTEx tissues – TPM')
# x_ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
# axs[1].set_xticks([np.log2(x + 1) for x in x_ticks])
# axs[1].set_xticklabels(x_ticks)
# fig.savefig('../figures/expression_downsampled-GTEX-control_GENCODE-isoforms_by-reference-vs-alternative.pdf',
#             bbox_inches='tight')


# In[248]:


# include all genes

from data_loading import DATA_DIR


def load_gtex_gencode_all_genes():
    """
    NOTE: not summing up the transcripts with identical CDS
    
    """
    df = pd.read_csv(
        DATA_DIR / "processed/expression_2022-09-01/transcript.GTEx-GC30_Isoforms.txt",
        sep="\t",
    )
    metadata = pd.read_csv(
        DATA_DIR / "processed/gtex_2022/GTEx_SRARunTable.txt", sep="\t"
    )
    if metadata["Run"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    metadata = metadata.set_index("Run")
    df = df.set_index("UID")
    df = (df + 1.0).apply(np.log2)

    def extract_ensembl_gene_id(s):
        return s.split("|")[1].split(".")[0]

    genes = pd.Series(
        index=df.index,
        data=df.index.map(extract_ensembl_gene_id).values,
    )
    df.index = df.index.map(lambda x: x.split('|')[0].split('.')[0])
    genes.index = genes.index.map(lambda x: x.split('|')[0].split('.')[0])
    return df, metadata, genes

df_gtex_all, metadata_gtex_all, genes_gtex_all = load_gtex_gencode_all_genes()
df_gtex_all = df_gtex_all.loc[:, ~df_gtex_all.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex_all = metadata_gtex_all.loc[~metadata_gtex_all['body_site'].isin(exclusion_list_gtex), :]
means_gtex_all = df_gtex_all.groupby(df_gtex_all.columns.map(metadata_gtex_all['body_site']), axis=1).mean()

path_MANE_select=DATA_DIR / "external/MANE.GRCh38.v0.95.summary.txt"
mane = pd.read_csv(path_MANE_select, sep="\t")
mane_select = set(
    mane.loc[mane["MANE_status"] == "MANE Select", "Ensembl_nuc"]
    .str.slice(0, 15)
    .values
)

fig, axs = plt.subplots(2, 1, sharex=True)
n_bins = 110
x_max = 11
axs[0].hist(means_gtex_all.max(axis=1)[means_gtex_all.index.isin(mane_select)],
            bins=n_bins,
            range=(0, x_max))
axs[1].hist(means_gtex_all.max(axis=1)[~means_gtex_all.index.isin(mane_select)],
            bins=n_bins,
            range=(0, x_max))
axs[0].text(x=8, y=axs[0].get_ylim()[1], s='MANE select transcripts')
axs[1].text(x=8, y=axs[1].get_ylim()[1], s='Alternative transcripts')
for ax in axs:
    ax.set_ylabel('Number of transcripts')
    for loc in ['top', 'right']:
        ax.spines[loc].set_visible(False)
axs[1].set_xlabel('Expression, max across GTEx tissues – TPM')
x_ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
axs[1].set_xticks([np.log2(x + 1) for x in x_ticks])
axs[1].set_xticklabels(x_ticks)
fig.savefig('../figures/expression_GTEX_GENCODE-transcripts_by-reference-vs-alternative-for-all-genes-not-just-TFs.pdf',
            bbox_inches='tight')


# In[249]:


means_gtex_all.max(axis=1)[~means_gtex_all.index.isin(mane_select)]


# In[250]:


# plot distribution of isoforms by TPM
fig, axs = plt.subplots(2, 1, sharex=True)
n_bins = 110
x_max = 11
axs[1].hist(means_dev.max(axis=1)[means_dev.index.isin(alt_isos)],
            bins=n_bins,
            range=(0, x_max))
axs[0].hist(means_dev.max(axis=1)[means_dev.index.isin(ref_isos)],
            bins=n_bins,
            range=(0, x_max))
axs[0].text(x=8, y=axs[0].get_ylim()[1] * 0.95, s='MANE select isoforms')
axs[1].text(x=8, y=axs[1].get_ylim()[1] * 0.95, s='Alternative isoforms')
for ax in axs:
    ax.set_ylabel('Number of isoforms')
    for loc in ['top', 'right']:
        ax.spines[loc].set_visible(False)
axs[1].set_xlabel('Expression, max across Cardoso Moreira tissues/timepoint – TPM')
x_ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
axs[1].set_xticks([np.log2(x + 1) for x in x_ticks])
axs[1].set_xticklabels(x_ticks)
fig.savefig('../figures/expression_development_GENCODE-isoforms_by-reference-vs-alternative.pdf',
            bbox_inches='tight')


# In[251]:


def load_developmental_tissue_expression_gencode_all_genes():
    """
    Cardoso-Moreira et al. Nature


    """
    metadata = pd.read_csv(DATA_DIR / "processed/Cardoso-Moreira_et_al/metadata.txt")
    df = pd.read_csv(
        DATA_DIR
        / "processed/expression_2022-09-01/transcript.Cardoso-Moreira-et-al-Nature-2019-GC30_Isoforms.txt",
        sep="\t",
    )
    if metadata["Run"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    metadata = metadata.set_index("Run")
    df = df.set_index("UID")
    df = (df + 1.0).apply(np.log2)

    def extract_ensembl_gene_id(s):
        return s.split("|")[1].split(".")[0]

    genes = pd.Series(
        index=df.index,
        data=df.index.map(extract_ensembl_gene_id).values,
    )
   
    genes = genes[~genes.index.duplicated(keep="first")]

    # the file has ERR2598062.fastq.gz instead of ERR2598062
    df.columns = df.columns.str.slice(0, -len(".fastq.gz"))

    if not df.columns.isin(metadata.index).all():
        raise UserWarning("Missing metadata")
    metadata = metadata.loc[metadata.index.isin(df.columns), :]
    df.index = df.index.map(lambda x: x.split('|')[0].split('.')[0])
    genes.index = genes.index.map(lambda x: x.split('|')[0].split('.')[0])
    return df, metadata, genes


df_dev_all, metadata_dev_all, genes_dev_all = load_developmental_tissue_expression_gencode_all_genes()
metadata_dev_all['dev_stage'] = metadata_dev_all['Developmental_Stage'].map(rename_dev_stage)
means_dev_all = (df_dev_all.groupby(df_dev_all.columns.map(metadata_dev_all['organism_part'] + ' ' + metadata_dev_all['dev_stage']), axis=1)
           .mean())

fig, axs = plt.subplots(2, 1, sharex=True)
n_bins = 110
x_max = 11
axs[0].hist(means_dev_all.max(axis=1)[means_dev_all.index.isin(mane_select)],
            bins=n_bins,
            range=(0, x_max))
axs[1].hist(means_dev_all.max(axis=1)[~means_dev_all.index.isin(mane_select)],
            bins=n_bins,
            range=(0, x_max))
axs[0].text(x=8, y=axs[0].get_ylim()[1] * 0.95, s='MANE select isoforms')
axs[1].text(x=8, y=axs[1].get_ylim()[1] * 0.95, s='Alternative isoforms')
for ax in axs:
    ax.set_ylabel('Number of isoforms')
    for loc in ['top', 'right']:
        ax.spines[loc].set_visible(False)
axs[1].set_xlabel('Expression, max across Cardoso Moreira tissues/timepoint – TPM')
x_ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
axs[1].set_xticks([np.log2(x + 1) for x in x_ticks])
axs[1].set_xticklabels(x_ticks)
fig.savefig('../figures/expression_development_GENCODE-transcripts_by-reference-vs-alternative-for-all-genes-not-just-TFs.pdf',
            bbox_inches='tight')


# In[252]:


(
    means_gtex.loc[means_gtex.index.isin(ref_isos)].max(axis=1).mean(),
    means_dev.loc[means_dev.index.isin(ref_isos)].max(axis=1).mean(),
    means_gtex_all.loc[means_gtex_all.index.isin(mane_select)].max(axis=1).mean(),
    means_dev_all.loc[means_dev_all.index.isin(mane_select)].max(axis=1).mean(),
)


# In[253]:


# GAPDH and actin
(
    means_gtex_all.groupby(genes_gtex_all).sum().max(axis=1)['ENSG00000111640'],
    means_dev_all.groupby(genes_dev_all).sum().max(axis=1)['ENSG00000111640'],
    means_gtex_all.groupby(genes_gtex_all).sum().max(axis=1)['ENSG00000075624'],
    means_dev_all.groupby(genes_dev_all).sum().max(axis=1)['ENSG00000075624']

    )


# In[254]:


# number of isoforms vs gene expression, publications, and exons 
tpm_per_gene = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .sum()
                .groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
                .mean())
gn = tpm_per_gene.max(axis=1).rename('TPM - gene-level, max across GTEx tissues').to_frame()
gn['n_isoforms'] = gn.index.map(genes_gtex.value_counts())


# In[255]:


#from ccsblib.huri import load_number_publications_per_gene, load_id_map

# n_pub = load_number_publications_per_gene().to_frame().reset_index()
# ensembl_to_hgnc = (load_id_map('ensembl_gene_id', 'hgnc_symbol')
#                     .drop_duplicates('ensembl_gene_id')
#                     .set_index('ensembl_gene_id')
#                     ['hgnc_symbol'])
# n_pub['hgnc_symbol'] = n_pub['ensembl_gene_id'].map(lambda x: ensembl_to_hgnc.get(x, np.nan))
# n_pub = n_pub.groupby('hgnc_symbol').sum()['n_pubmed_ids']
# gn['Number of publications'] = gn.index.map(n_pub)


# In[256]:


# TODO: change to get reference isoform
gn['Number of exons in reference isoform'] = gn.index.map({name: len(tf.orfs[0].exons) for name, tf in tfs.items()})


# In[257]:


gn.head()


# In[258]:


# log-scale?
fig, ax = plt.subplots(1, 1)
x_col = 'TPM - gene-level, max across GTEx tissues'
y_col = 'n_isoforms'
x = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), x_col].values
y = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), y_col].values
ax.scatter(x, y, alpha=0.05)
ax.set_ylabel('Number of unique protein isoforms')
ax.set_xlabel(x_col)
r = stats.pearsonr(x, y)[0]
rho = stats.spearmanr(x, y)[0]
ax.text(x=x.max() * 0.9, y=20, s=f'r = {r:.2f}', ha='right')
ax.text(x=x.max() * 0.9, y=18, s=f'rho = {rho:.2f}', ha='right')


# In[259]:


import seaborn as sns

# fig, axs = plt.subplots(3, 1, sharex=True)
# fig.set_size_inches(w=8, h=12)

# y_col = 'n_isoforms'
# x = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), x_col].values
# y = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), y_col].values

# ax = axs[0]
# x_col = 'TPM - gene-level, max across GTEx tissues'
# sns.boxplot(data=gn, 
#                x='n_isoforms', 
#                y=x_col,
#             ax=ax,
#             color='C0')
# ax.set_yscale('log')
# x = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), x_col].values
# y = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), y_col].values
# r = stats.pearsonr(np.log2(x + 1), y)[0]
# rho = stats.spearmanr(x, y)[0]
# ax.text(y=x.max() * 0.9, x=18, s=f'Pearson\'s r = {r:.2f}', ha='right')
# ax.text(y=x.max() * 0.5, x=18, s=f'Spearman\'s rho = {rho:.2f}', ha='right')

# ax = axs[1]
# x_col = 'Number of publications'
# sns.boxplot(data=gn, 
#                x='n_isoforms', 
#                y=x_col,
#                 ax=ax,
#                 color='C0')
# ax.set_yscale('log')
# x = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), x_col].values
# y = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), y_col].values
# r = stats.pearsonr(np.log2(x + 1), y)[0]
# rho = stats.spearmanr(x, y)[0]
# ax.text(y=x.max() * 0.9, x=18, s=f'Pearson\'s r = {r:.2f}', ha='right')
# ax.text(y=x.max() * 0.5, x=18, s=f'Spearman\'s rho = {rho:.2f}', ha='right')

# ax = axs[2]
# x_col = 'Number of exons in reference isoform'
# sns.boxplot(data=gn,
#                x='n_isoforms', 
#                y=x_col,
#              ax=ax,
#              color='C0')
# x = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), x_col].values
# y = gn.loc[gn[x_col].notnull() & gn[y_col].notnull(), y_col].values
# r = stats.pearsonr(x, y)[0]
# rho = stats.spearmanr(x, y)[0]
# ax.text(y=x.max() * 0.9, x=18, s=f'Pearson\'s r = {r:.2f}', ha='right')
# ax.text(y=x.max() * 0.8, x=18, s=f'Spearman\'s rho = {rho:.2f}', ha='right')
# fig.savefig('../figures/n-isoforms-vs-expression-publication-n-exons_boxplot.pdf',
#             bbox_inches='tight')


# In[260]:


from data_loading import load_tf_families

fam = load_tf_families()
gn['family'] = gn.index.map(fam)
gn['is_nuclear_receptor'] = (gn['family'] == 'Nuclear receptor')


# In[261]:


gn['family'].value_counts().head(10)


# In[262]:


# fig, ax = plt.subplots(1, 1)
# sns.boxplot(data=gn,
#             x='is_nuclear_receptor',
#             y='Number of publications',
#             ax=ax)
# ax.set_yscale('log')


# In[263]:


fig, ax = plt.subplots(1, 1)
sns.boxplot(data=gn,
            x='is_nuclear_receptor',
            y='Number of exons in reference isoform',
            ax=ax)


# In[264]:


fig, ax = plt.subplots(1, 1)
sns.boxplot(data=gn,
            x='is_nuclear_receptor',
            y='TPM - gene-level, max across GTEx tissues',
            ax=ax)
ax.set_yscale('log')


# In[265]:


# percentage of alternative isoform
# plot distribution of fraction of gene expression for ref and alt

# has to be fraction where isoform TPM is at least 1, right (fill na with 0)

per_gene_gtex = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .transform('sum'))
f_gtex = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
        .mean())
f_gtex = f_gtex * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex = f_gtex * 100
fig, axs = plt.subplots(2, 1, sharex=True)
n_bins=100
axs[0].hist(f_gtex.max(axis=1).loc[f_gtex.index.isin(ref_isos)].values,
            range=(0, 100),
            bins=n_bins)
axs[1].hist(f_gtex.max(axis=1).loc[f_gtex.index.isin(alt_isos)].values,
            range=(0, 100),
            bins=n_bins)
axs[1].set_xlabel('% of gene expression, GTEx, maximum across tissues')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
for ax in axs:
    ax.set_ylabel('Number of isoforms')
    for pos in ['top', 'right']:
        ax.spines[pos].set_visible(False)
plt.subplots_adjust(hspace=0.6)
fig.savefig('../figures/expression-fraction-GTEx-max-across-tissues_ref-vs-alt_hist.pdf',
            bbox_inches='tight')


# In[266]:


# percentage of alternative isoform
# plot distribution of fraction of gene expression for ref and alt

# has to be fraction where isoform TPM is at least 1, right (fill na with 0)

per_gene_dev = ((2 ** df_dev - 1)
                .groupby(genes_dev)
                .transform('sum'))
f_dev = (((2 ** df_dev - 1) / per_gene_dev)
        .groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']),
         axis=1)
        .mean())
f_dev = f_dev * ((per_gene_dev.groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']),
                                             axis=1)
                                             .mean() >= 1)
                                         .applymap(lambda x: {False: np.nan, True: 1}[x]))  # only count fractions if gene TPM is >= 1

f_dev = f_dev * 100
fig, axs = plt.subplots(2, 1, sharex=True)
n_bins=100
axs[0].hist(f_dev.max(axis=1).loc[f_dev.index.isin(ref_isos)].values,
            range=(0, 100),
            bins=n_bins)
axs[1].hist(f_dev.max(axis=1).loc[f_dev.index.isin(alt_isos)].values,
            range=(0, 100),
            bins=n_bins)
axs[1].set_xlabel('% of gene expression, Cardoso Moreira, maximum across tissues')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
for ax in axs:
    ax.set_ylabel('Number of isoforms')
    for pos in ['top', 'right']:
        ax.spines[pos].set_visible(False)
plt.subplots_adjust(hspace=0.6)
fig.savefig('../figures/expression-fraction-development-max-across-tissues_ref-vs-alt_hist.pdf',
            bbox_inches='tight')


# In[267]:


TPM_THRESHOLD = 2
per_gene_dev = ((2 ** df_dev - 1)
                .groupby(genes_dev)
                .transform('sum'))
f_dev = (((2 ** df_dev - 1) / per_gene_dev)
        .groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']),
         axis=1)
        .mean())
f_dev = f_dev * ((per_gene_dev.groupby(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']),
                                             axis=1)
                                             .mean() >= TPM_THRESHOLD)
                                         .applymap(lambda x: {False: np.nan, True: 1}[x]))  # only count fractions if gene TPM is >= 1
f_dev = f_dev * 100

per_gene_gtex = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .transform('sum'))
f_gtex = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
        .mean())
f_gtex = f_gtex * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean() >= TPM_THRESHOLD).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1
f_gtex = f_gtex * 100

# alternative isoform > 50% in dev with corresponding ref iso > 50% in either
# plus some magnitude of change?
putative_switching_alt_isos = set(f_dev.loc[f_dev.index.isin(alt_isos) &
                                           (f_dev.max(axis=1) > 50), :].index.values)
print(len(putative_switching_alt_isos))
valid_ref_isos = (set(f_dev.loc[f_dev.index.isin(ref_isos) &
                               (f_dev.max(axis=1) > 50), :].index.values)
                        .union(
                  set(f_gtex.loc[f_gtex.index.isin(ref_isos) &
                               (f_gtex.max(axis=1) > 50), :].index.values)
                        ))
if not (genes_gtex == genes_dev).all():
        raise UserWarning()
genes = genes_gtex
genes_with_valid_ref_isos = {genes[iso] for iso in valid_ref_isos}
switching_alt_isos = {iso for iso in putative_switching_alt_isos if genes[iso] in genes_with_valid_ref_isos}
print(len(switching_alt_isos))


# ### note to luke - this block of code seems to finish correctly only sometimes, something about the .sample() step

# In[268]:


# percentage of alternative isoform
# plot distribution of fraction of gene expression for ref and alt

# has to be fraction where isoform TPM is at least 1, right (fill na with 0)

per_gene_gtex = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .transform('sum'))
f_gtex = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1)
        .mean())
f_gtex = f_gtex * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex = f_gtex * 100
fig, axs = plt.subplots(2, 1, sharex=True)
n_bins=100
axs[0].hist(f_gtex.max(axis=1).loc[f_gtex.index.isin(ref_isos)].values,
            range=(0, 100),
            bins=n_bins)
axs[1].hist(f_gtex.max(axis=1).loc[f_gtex.index.isin(alt_isos)].values,
            range=(0, 100),
            bins=n_bins)
axs[1].set_xlabel('% of gene expression, downsampled GTEx, maximum across dummy tissues')
axs[0].set_title('MANE select isoforms')
axs[1].set_title('Alternative isoforms')
for ax in axs:
    ax.set_ylabel('Number of isoforms')
    for pos in ['top', 'right']:
        ax.spines[pos].set_visible(False)
plt.subplots_adjust(hspace=0.6)
fig.savefig('../figures/expression-fraction-downsampled-GTEx-control-max-across-tissues_ref-vs-alt_hist.pdf',
            bbox_inches='tight')


# In[269]:


# TODO: diagonal line
# Comparing adult samples between GTEx and develpment
import statsmodels.api as sm
paired_tissues = [('liver adult', 'Liver'),
                  ('heart young adult', 'Heart - Atrial Appendage'),
                  ('testis adult', 'Testis')]
for cm_tissue, gtex_tissue in paired_tissues:
    fig, axs = plt.subplots(1, 2)
    fig.set_size_inches(w=8, h=4)
    paired = pd.merge(means_dev.loc[:, [cm_tissue]],
                      means_gtex.loc[:, [gtex_tissue]],
            how='inner',
            left_index=True, right_index=True)
    (paired.loc[paired.index.isin(ref_isos), :]
    .plot.scatter(x=gtex_tissue, 
                        y=cm_tissue, 
                        alpha=0.2,
                        ax=axs[0]))
    upper = 8
    #upper = paired.max().max()
    r = paired.loc[paired.index.isin(ref_isos), :].corr().loc[cm_tissue, gtex_tissue]

    x = paired.loc[paired.index.isin(ref_isos), gtex_tissue].values
    y = paired.loc[paired.index.isin(ref_isos), cm_tissue].values
    intercept, slope = sm.OLS(y, sm.add_constant(x)).fit().params
    axs[0].plot(range(upper + 1), [intercept + slope * x for x in range(upper + 1)], '--', color='grey', zorder=-1)

    axs[0].set_title('{} – MANE select isoforms – R^2 = {:.2f}'.format(gtex_tissue.split()[0], r**2),
                    fontsize=10)
    (paired.loc[paired.index.isin(alt_isos), :]
    .plot.scatter(x=gtex_tissue, 
                        y=cm_tissue, 
                        alpha=0.2,
                        ax=axs[1]))
    r = paired.loc[paired.index.isin(alt_isos), :].corr().loc[cm_tissue, gtex_tissue]
    x = paired.loc[paired.index.isin(alt_isos), gtex_tissue].values
    y = paired.loc[paired.index.isin(alt_isos), cm_tissue].values
    intercept, slope = sm.OLS(y, sm.add_constant(x)).fit().params
    axs[1].plot(range(upper + 1), [intercept + slope * x for x in range(upper + 1)], '--', color='grey', zorder=-1)
    axs[1].set_title('{} alternative isoforms – R^2 = {:.2f}'.format(gtex_tissue.split()[0], r**2),
                     fontsize=10)

    for ax in axs:
        ax.set_ylim(0, upper + 0.5)
        ax.set_xlim(0, upper + 0.5)
        ax.set_yticks(range(0, int(upper+1)))
        ax.set_xticks(range(0, int(upper+1)))
        ax.plot([0, upper], [0, upper], '-', color='black')
        ax.set_xlabel('GTEx – Mean log2(TPM + 1)')
        ax.set_ylabel('Developmental dataset – Mean log2(TPM + 1)')
    plt.savefig('../figures/GTEx-vs-dev_{}_GENCODE.pdf'.format(gtex_tissue), bbox_inches='tight')


# In[270]:


# calculate tissue-specificity score for each isoform
def calculate_tiss_spec(x):
    # assumes input of a row of only the tissue expression data
    med = x.median() 
    iqr = x.quantile(0.75) - x.quantile(0.25)
    for tiss in x.index:
        score = (x[tiss] - med) / float(iqr)
        x[tiss] = score
    x['max_tip'] = x.max()
    return x


# In[271]:


# dot / box plot
df = df_gtex.T.copy()
isoform = 'AEBP2-207'
df['body_site'] = df.index.map(metadata_gtex['body_site'])
sns.boxplot(data=df, x='body_site', y=isoform)
sns.stripplot(data=df, x='body_site', y=isoform)


# ### kaia's new plots

# In[272]:


ref_alt_map = pd.DataFrame([ref_isos]).T
ref_alt_map.columns = ["ref"]
ref_alt_map["gene"] = ref_alt_map["ref"].str.split("|", expand=True)[0].str[:-4]

alt_isos_df = pd.DataFrame([alt_isos]).T
alt_isos_df.columns = ["alt"]
alt_isos_df["gene"] = alt_isos_df["alt"].str.split("|", expand=True)[0].str[:-4]

ref_alt_map = ref_alt_map.merge(alt_isos_df, on="gene", how="left")
print(len(ref_alt_map))
ref_alt_map_nonan = ref_alt_map[~pd.isnull(ref_alt_map["alt"])]
print(len(ref_alt_map_nonan))
ref_alt_map_nonan.head()


# In[273]:


ref_alt_map_nonan[ref_alt_map_nonan["gene"] == "NKX2-5"]


# In[274]:


means_dev["max_dev"] = means_dev.max(axis=1)

# max out anything above 11 (2000 tpm) to make plots more readable, as luke did above
means_dev[means_dev["max_dev"] > 11] = 11

print(means_dev["max_dev"].max())
print(means_dev["max_dev"].min())
means_dev_ri = means_dev.reset_index()
means_dev_ri["UID_rep"] = means_dev_ri["UID"].str.replace("_", "|")
means_dev_ri.head()


# In[275]:


means_gtex["max_gtex"] = means_gtex.max(axis=1)

# max out anything above 11 (2000 tpm) to make plots more readable, as luke did above
means_gtex[means_gtex["max_gtex"] > 11] = 11

print(means_gtex["max_gtex"].max())
print(means_gtex["max_gtex"].min())
means_gtex_ri = means_gtex.reset_index()
means_gtex_ri["UID_rep"] = means_gtex_ri["UID"].str.replace("_", "|")
means_gtex_ri.head()


# In[276]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(means_gtex_ri[["UID_rep", "max_gtex"]], left_on="ref", 
                                            right_on="UID_rep", how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(means_gtex_ri[["UID_rep", "max_gtex"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(means_dev_ri[["UID_rep", "max_dev"]], left_on="ref", 
                                            right_on="UID_rep", how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(means_dev_ri[["UID_rep", "max_dev"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="inner")
print(len(ref_alt_map_nonan))

ref_alt_map_nonan.head()


# In[277]:


ref_alt_map_nonan.sort_values(by="max_gtex_ref", ascending=False).head()


# In[278]:


fig = plt.figure(figsize=(3.5, 3))

ax = sns.histplot(data=ref_alt_map_nonan, x="max_gtex_ref", y="max_gtex_alt",
                  bins=30, cbar=True, cbar_kws={"label": "# isoform pairs"})

ax.set_xlim((-0.3, 11.5))
ax.set_ylim((-0.3, 11.5))
ax.set_xlabel("maximum expression of reference")
ax.set_ylabel("maximum expression of alternative")
ax.set_title("GTEx dataset\n(n=%s ref/alt pairs)" % len(ref_alt_map_nonan))

ticks = [0, 1, 2, 5, 10, 50, 100, 500, 2000]
ticklabels = [0, 1, 2, 5, 10, 50, 100, 500, "2000+"]
ax.set_xticks([np.log2(x + 1) for x in ticks])
ax.set_xticklabels(ticklabels)
ax.set_yticks([np.log2(y + 1) for y in ticks])
ax.set_yticklabels(ticklabels)

# find num where ref > alt
ra = len(ref_alt_map_nonan[ref_alt_map_nonan["max_gtex_ref"] > ref_alt_map_nonan["max_gtex_alt"]])
ax.text(8.5, 7, "%s\n(%s%%)" % (ra, round(ra/len(ref_alt_map_nonan), 2)*100), ha="left", va="bottom")

# find num where alt > ref
ar = len(ref_alt_map_nonan[ref_alt_map_nonan["max_gtex_ref"] < ref_alt_map_nonan["max_gtex_alt"]])
ax.text(9, 10.5, "%s\n(%s%%)" % (ar, round(ar/len(ref_alt_map_nonan), 2)*100), ha="right", va="top")

ax.plot([-0.3,11.5], [-0.3, 11.5], color="black", linestyle="dashed")

fig.savefig('../figures/expression-scatter-ref_v_alt-gtex.pdf',
            bbox_inches='tight')


# In[279]:


fig = plt.figure(figsize=(3.5, 3))

ax = sns.histplot(data=ref_alt_map_nonan, x="max_dev_ref", y="max_dev_alt",
                  bins=30, cbar=True, cbar_kws={"label": "# isoform pairs"})

ax.set_xlim((-0.3, 11.5))
ax.set_ylim((-0.3, 11.5))
ax.set_xlabel("maximum expression of reference")
ax.set_ylabel("maximum expression of alternative")
ax.set_title("developmental dataset\n(n=%s ref/alt pairs)" % len(ref_alt_map_nonan))

ticks = [0, 1, 2, 5, 10, 50, 100, 500, 2000]
ticklabels = [0, 1, 2, 5, 10, 50, 100, 500, "2000+"]
ax.set_xticks([np.log2(x + 1) for x in ticks])
ax.set_xticklabels(ticklabels)
ax.set_yticks([np.log2(y + 1) for y in ticks])
ax.set_yticklabels(ticklabels)

# find num where ref > alt
ra = len(ref_alt_map_nonan[ref_alt_map_nonan["max_dev_ref"] > ref_alt_map_nonan["max_dev_alt"]])
ax.text(8.5, 7, "%s\n(%s%%)" % (ra, round(ra/len(ref_alt_map_nonan), 2)*100), ha="left", va="bottom")

# find num where alt > ref
ar = len(ref_alt_map_nonan[ref_alt_map_nonan["max_dev_ref"] < ref_alt_map_nonan["max_dev_alt"]])
ax.text(9, 10.5, "%s\n(%s%%)" % (ar, round(ar/len(ref_alt_map_nonan), 2)*100), ha="right", va="top")

ax.plot([-0.3,11.5], [-0.3, 11.5], color="black", linestyle="dashed")

fig.savefig('../figures/expression-scatter-ref_v_alt-dev.pdf',
            bbox_inches='tight')


# In[280]:


# putting NA as 0!
f_dev["max_ratio_dev"] = f_dev.max(axis=1)
f_dev["max_ratio_dev"].fillna(0, inplace=True)
f_dev["min_ratio_dev"] = f_dev.min(axis=1)
f_dev["min_ratio_dev"].fillna(0, inplace=True)
print(f_dev["max_ratio_dev"].max())
print(f_dev["max_ratio_dev"].min())
f_dev_ri = f_dev.reset_index()
f_dev_ri["UID_rep"] = f_dev_ri["UID"].str.replace("_", "|")
f_dev_ri.head()


# In[281]:


# putting NA as 0!
f_gtex["max_ratio_gtex"] = f_gtex.max(axis=1)
f_gtex["max_ratio_gtex"].fillna(0, inplace=True)
f_gtex["min_ratio_gtex"] = f_gtex.min(axis=1)
f_gtex["min_ratio_gtex"].fillna(0, inplace=True)
print(f_gtex["max_ratio_gtex"].max())
print(f_gtex["max_ratio_gtex"].min())
f_gtex_ri = f_gtex.reset_index()
f_gtex_ri["UID_rep"] = f_gtex_ri["UID"].str.replace("_", "|")
f_gtex_ri.head()


# In[282]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(f_gtex_ri[["UID_rep", "max_ratio_gtex", "min_ratio_gtex"]], left_on="ref", 
                                            right_on="UID_rep", how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(f_gtex_ri[["UID_rep", "max_ratio_gtex", "min_ratio_gtex"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(f_dev_ri[["UID_rep", "max_ratio_dev", "min_ratio_dev"]], left_on="ref", 
                                            right_on="UID_rep", how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(f_dev_ri[["UID_rep", "max_ratio_dev", "min_ratio_dev"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="inner")
print(len(ref_alt_map_nonan))
ref_alt_map_nonan.head()


# In[283]:


fig = plt.figure(figsize=(3.5, 3))

ax = sns.histplot(data=ref_alt_map_nonan, x="max_ratio_gtex_ref", y="max_ratio_gtex_alt",
                  bins=30, cbar=True)
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xlabel("maximum expression of reference\nfraction of gene tpm")
ax.set_ylabel("maximum expression of alternative\nfraction of gene tpm")
ax.set_title("GTEx dataset\n(n=%s ref/alt pairs)" % len(ref_alt_map_nonan))

# find num where ref > alt
ra = len(ref_alt_map_nonan[ref_alt_map_nonan["max_ratio_gtex_ref"] > ref_alt_map_nonan["max_ratio_gtex_alt"]])
ax.text(90, 0, "%s\n(%s%%)" % (ra, round(ra/len(ref_alt_map_nonan), 2)*100), ha="right", va="bottom")

# find num where alt > ref
ar = len(ref_alt_map_nonan[ref_alt_map_nonan["max_ratio_gtex_ref"] < ref_alt_map_nonan["max_ratio_gtex_alt"]])
ax.text(0, 99, "%s\n(%s%%)" % (ar, round(ar/len(ref_alt_map_nonan), 2)*100), ha="left", va="top")

ax.plot([-2,102], [-2, 102], color="black", linestyle="dashed")

fig.savefig('../figures/expression-ratio-scatter-ref_v_alt-gtex.pdf',
            bbox_inches='tight')


# In[284]:


fig = plt.figure(figsize=(3.5, 3))

ax = sns.histplot(data=ref_alt_map_nonan, x="max_ratio_dev_ref", y="max_ratio_dev_alt",
                  bins=30, cbar=True)
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xlabel("maximum expression of reference\nfraction of gene tpm")
ax.set_ylabel("maximum expression of alternative\nfraction of gene tpm")
ax.set_title("developmental dataset\n(n=%s ref/alt pairs)" % len(ref_alt_map_nonan))

# find num where ref > alt
ra = len(ref_alt_map_nonan[ref_alt_map_nonan["max_ratio_dev_ref"] > ref_alt_map_nonan["max_ratio_dev_alt"]])
ax.text(90, 0, "%s\n(%s%%)" % (ra, round(ra/len(ref_alt_map_nonan), 2)*100), ha="right", va="bottom")

# find num where alt > ref
ar = len(ref_alt_map_nonan[ref_alt_map_nonan["max_ratio_dev_ref"] < ref_alt_map_nonan["max_ratio_dev_alt"]])
ax.text(0, 99, "%s\n(%s%%)" % (ar, round(ar/len(ref_alt_map_nonan), 2)*100), ha="left", va="top")

ax.plot([-2,102], [-2, 102], color="black", linestyle="dashed")

fig.savefig('../figures/expression-ratio-scatter-ref_v_alt-dev.pdf',
            bbox_inches='tight')


# In[285]:


ref_alt_map_nonan["mm_dev_alt_delta"] = ref_alt_map_nonan["max_ratio_dev_alt"] - ref_alt_map_nonan["min_ratio_dev_alt"]
ref_alt_map_nonan["mm_dev_ref_delta"] = ref_alt_map_nonan["max_ratio_dev_ref"] - ref_alt_map_nonan["min_ratio_dev_ref"]
ref_alt_map_nonan["mm_gtex_alt_delta"] = ref_alt_map_nonan["max_ratio_gtex_alt"] - ref_alt_map_nonan["min_ratio_gtex_alt"]
ref_alt_map_nonan["mm_gtex_ref_delta"] = ref_alt_map_nonan["max_ratio_gtex_ref"] - ref_alt_map_nonan["min_ratio_gtex_ref"]
ref_alt_map_nonan.sample(5)


# In[345]:


fig = plt.figure(figsize=(3.5, 2))

dd_ref = ref_alt_map_nonan[["ref", "mm_dev_ref_delta"]].drop_duplicates()
dd_alt = ref_alt_map_nonan[["alt", "mm_dev_alt_delta"]].drop_duplicates()

ax = sns.distplot(dd_alt["mm_dev_alt_delta"], color=sns.color_palette("deep")[3], kde=False,
                  label="alternative (n=%s)" % len(dd_alt))
sns.distplot(dd_ref["mm_dev_ref_delta"], color=sns.color_palette("deep")[0],
             kde=False, label="reference (n=%s)" % len(dd_ref), ax=ax)

ax.set_xlabel("maximum ∆ ratio across samples")
ax.set_ylabel("isoform count")
ax.set_title("developmental data")
ax.set_xlim((-1, 101))

ax.axvline(x=80, linestyle="dashed", color="black")

ref_o80 = len(dd_ref[dd_ref["mm_dev_ref_delta"] > 80])
alt_o80 = len(dd_alt[dd_alt["mm_dev_alt_delta"] > 80])

ax.text(101, 323, "%s\n(%s%%)" % (alt_o80, round(alt_o80/len(dd_alt)*100, 2)), ha="right", va="top", 
        color=sns.color_palette("deep")[3])
ax.text(101, 250, "%s\n(%s%%)" % (ref_o80, round(ref_o80/len(dd_ref)*100, 1)), ha="right", va="top", 
        color=sns.color_palette("deep")[0])

plt.legend(loc=2, bbox_to_anchor=(0.12, 1))
fig.savefig('../figures/expression-histogram-max_delta-dev.pdf',
            bbox_inches='tight')


# In[346]:


fig = plt.figure(figsize=(3.5, 2))

dd_ref = ref_alt_map_nonan[["ref", "mm_gtex_ref_delta"]].drop_duplicates()
dd_alt = ref_alt_map_nonan[["alt", "mm_gtex_alt_delta"]].drop_duplicates()

ax = sns.distplot(dd_alt["mm_gtex_alt_delta"], color=sns.color_palette("deep")[3], kde=False,
                  label="alternative (n=%s)" % len(dd_alt))
sns.distplot(dd_ref["mm_gtex_ref_delta"], color=sns.color_palette("deep")[0],
             kde=False, label="reference (n=%s)" % len(dd_ref), ax=ax)

ax.set_xlabel("maximum ∆ ratio across samples")
ax.set_ylabel("isoform count")
ax.set_title("GTEx data (down-sampled)")
ax.set_xlim((-1, 101))

ax.axvline(x=80, linestyle="dashed", color="black")

ref_o80 = len(dd_ref[dd_ref["mm_gtex_ref_delta"] > 80])
alt_o80 = len(dd_alt[dd_alt["mm_gtex_alt_delta"] > 80])

ax.text(101, 375, "%s\n(%s%%)" % (alt_o80, round(alt_o80/len(dd_alt)*100, 1)), ha="right", va="top", 
        color=sns.color_palette("deep")[3])
ax.text(101, 290, "%s\n(%s%%)" % (ref_o80, round(ref_o80/len(dd_ref)*100, 1)), ha="right", va="top", 
        color=sns.color_palette("deep")[0])

plt.legend(loc=2, bbox_to_anchor=(0.07, 1))
fig.savefig('../figures/expression-histogram-max_delta-gtex.pdf',
            bbox_inches='tight')


# ### # make domain figure - move this into domain notebook at some point

# In[323]:


# loop through ref/alt pairs above and calculate total num AAs inserted/deleted/frameshifted
tot_ins = []
tot_perc_ins = []
tot_dd = []
tot_perc_dd = []
tot_f = []
tot_perc_f = []

tot_ins_dom = []
tot_perc_ins_dom = []
tot_dd_dom = []
tot_perc_dd_dom = []
tot_f_dom = []
tot_perc_f_dom = []

for i, row in ref_alt_map_nonan.iterrows():
    ref = row.ref.split("|")[0]
    alt = row.alt.split("|")[0]
    gene = ref[:-4]
    
    # manual fixes
    if gene == "AC092072.1":
        gene = "ZNF223"
    if gene == "AC008554.1":
        gene = "ZNF737"
    if gene == "AC073611.1":
        gene = "SP7"
    if gene == "AC118549.1":
        gene = "ZZZ3"
    if gene == "ZUP1":
        gene = "ZUFSP"
    if gene == "AC139768.1":
        gene = "POU6F1"
    if gene == "PHF19":
        gene = "PHF19 "
    #print("gene: %s | ref: %s | alt: %s" % (gene, ref, alt))
    
    pp_str = tfs[gene].pairwise_changes_relative_to_reference(ref, alt)
    aa_ftr = tfs[gene].aa_feature_disruption(ref)
    if len(aa_ftr) == 0:
        ins_dom = 0
        perc_ins_dom = 0
        dd_dom = 0
        perc_dd_dom = 0
        f_dom = 0
        perc_f_dom = 0
    else:
        aa_ftr_alt = aa_ftr[aa_ftr["alt_iso"] == alt]
        aa_ftr_alt_dom_grp = aa_ftr_alt.groupby("alt_iso")[["deletion",
                                                            "insertion",
                                                            "frameshift"]].agg("sum").reset_index()
        
        ins_dom = aa_ftr_alt_dom_grp.insertion.iloc[0]
        perc_ins_dom = ins_dom/len(pp_str)*100
        dd_dom = aa_ftr_alt_dom_grp.deletion.iloc[0]
        perc_dd_dom = dd_dom/len(pp_str)*100
        f_dom = aa_ftr_alt_dom_grp.frameshift.iloc[0]
        perc_f_dom = f_dom/len(pp_str)*100
    
    
    ins = pp_str.count("I")
    perc_ins = ins/len(pp_str)*100
    dd = pp_str.count("D")
    perc_dd = dd/len(pp_str)*100
    f = pp_str.count("F")
    f += pp_str.count("f")
    perc_f = f/len(pp_str)*100
    
#     print("# insertions: %s (%s%%) | # deletions: %s (%s%%) | # frameshifts: %s (%s%%)" % (ins, perc_ins,
#                                                                                            dd, perc_dd,
#                                                                                            f, perc_f))
    
    tot_ins.append(ins)
    tot_perc_ins.append(perc_ins)
    tot_dd.append(dd)
    tot_perc_dd.append(perc_dd)
    tot_f.append(f)
    tot_perc_f.append(perc_f)
    
    tot_ins_dom.append(ins_dom)
    tot_perc_ins_dom.append(perc_ins_dom)
    tot_dd_dom.append(dd_dom)
    tot_perc_dd_dom.append(perc_dd_dom)
    tot_f_dom.append(f_dom)
    tot_perc_f_dom.append(perc_f_dom)

ref_alt_map_nonan["n_ins"] = tot_ins
ref_alt_map_nonan["perc_ins"] = tot_perc_ins
ref_alt_map_nonan["n_dd"] = tot_dd
ref_alt_map_nonan["perc_dd"] = tot_perc_dd
ref_alt_map_nonan["n_f"] = tot_f
ref_alt_map_nonan["perc_f"] = tot_perc_f

ref_alt_map_nonan["n_ins_dom"] = tot_ins_dom
ref_alt_map_nonan["perc_ins_dom"] = tot_perc_ins_dom
ref_alt_map_nonan["n_dd_dom"] = tot_dd_dom
ref_alt_map_nonan["perc_dd_dom"] = tot_perc_dd_dom
ref_alt_map_nonan["n_f_dom"] = tot_f_dom
ref_alt_map_nonan["perc_f_dom"] = tot_perc_f_dom
ref_alt_map_nonan.sample(5)


# In[324]:


def mimic_r_boxplot(ax):
    for i, patch in enumerate(ax.artists):
        r, g, b, a = patch.get_facecolor()
        col = (r, g, b, 1)
        patch.set_facecolor((r, g, b, .5))
        patch.set_edgecolor((r, g, b, 1))

        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        line_order = ["lower", "upper", "whisker_1", "whisker_2", "med", "fliers"]
        for j in range(i*6,i*6+6):
            elem = line_order[j%6]
            line = ax.lines[j]
            if "whisker" in elem:
                line.set_visible(False)
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)
            if "fliers" in elem:
                line.set_alpha(0.5)


# In[327]:


to_plot = pd.melt(ref_alt_map_nonan, id_vars=["ref", "gene", "alt"], value_vars=["n_ins", "perc_ins",
                                                                                 "n_dd", "perc_dd",
                                                                                 "n_f", "perc_f",
                                                                                 "n_ins_dom", "perc_ins_dom",
                                                                                 "n_dd_dom", "perc_dd_dom",
                                                                                 "n_f_dom", "perc_f_dom"])
to_plot["dom_cat"] = to_plot["variable"].str.contains("dom")
to_plot["n_or_perc"] = to_plot["variable"].str.split("_", expand=True)[0]
to_plot["type"] = to_plot["variable"].str.split("_", expand=True)[1]
to_plot.sample(5)


# In[353]:


fig = plt.figure(figsize=(3.6, 3))
ax = sns.boxplot(data=to_plot[to_plot["n_or_perc"] == "perc"], 
                 x="type", y="value", hue="dom_cat", order=["dd", "ins", "f"],
                 palette=sns.color_palette("deep"), fliersize=5, notch=True,
                 flierprops={"marker": "o"})
mimic_r_boxplot(ax)
ax.set_xlabel("")
ax.set_xticklabels(["deletions", "insertions", "frameshift"], rotation=30, ha="right", va="top")
ax.set_ylabel("% of amino acids affected")
ax.set_title("sequence changes in alternative TF isoforms\ncompared to reference TF isoforms")
handles, labels = ax.get_legend_handles_labels()
labels = ["in entire protein", "in annotated domains"]
ax.legend(handles, labels)
fig.savefig('../figures/domain-overall-boxplot.pdf',
            bbox_inches='tight')


# In[355]:


to_plot[to_plot["n_or_perc"] == "perc"].groupby(["type", "dom_cat"]).agg("median")


# In[365]:


len(ref_alt_map_nonan[ref_alt_map_nonan["perc_f_dom"] > 0])


# In[357]:


len(ref_alt_map_nonan[ref_alt_map_nonan["perc_f"] >= 10])


# In[358]:


len(ref_alt_map_nonan)


# In[366]:


214/2305

