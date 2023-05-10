
# coding: utf-8

# ## Explore TF isoform expression data
# 
# kaia cleaning up luke's original code

# In[1]:


import os
import itertools
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
                          load_developmental_tissue_expression_gencode,
                          load_tf_families)


# In[2]:


np.random.seed(2023)


# In[3]:


PAPER_PRESET = {"style": "ticks", "font": "Helvetica", "context": "paper", 
                "rc": {"font.size":7,"axes.titlesize":7,
                       "axes.labelsize":7, 'axes.linewidth':0.5,
                       "legend.fontsize":6, "xtick.labelsize":6,
                       "ytick.labelsize":6, "xtick.major.size": 3.0,
                       "ytick.major.size": 3.0, "axes.edgecolor": "black",
                       "xtick.major.pad": 3.0, "ytick.major.pad": 3.0}}
PAPER_FONTSIZE = 7


# In[4]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# ## 1. load Gencode TFs + GTEx + Dev RNA-seq

# In[5]:


tfs = load_annotated_gencode_tfs()

df_gtex, metadata_gtex, genes_gtex = load_gtex_gencode()

exclusion_list_gtex = {'Cells - Leukemia cell line (CML)',
                       'Cells - EBV-transformed lymphocytes',
                       'Cells - Cultured fibroblasts'}

df_gtex = df_gtex.loc[:, ~df_gtex.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex = metadata_gtex.loc[~metadata_gtex['body_site'].isin(exclusion_list_gtex), :]

means_gtex = df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean()

df_dev, metadata_dev, genes_dev = load_developmental_tissue_expression_gencode()

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


# In[6]:


metadata_dev.shape


# In[7]:


metadata_gtex.shape


# In[8]:


len(all_isos)


# In[9]:


len(ref_isos)


# In[10]:


len(alt_isos)


# In[11]:


(means_gtex > 1).any(axis=1).value_counts()


# In[12]:


(means_gtex.loc[means_gtex.index.isin(alt_isos), :].sum(axis=1) >= 1).sum()


# ## 2. cumulative plot: % of samples with isoforms ≥ 1 TPM in both datasets

# In[13]:


n_samples_gtex = df_gtex.shape[1]
n_samples_gt1_gtex = (df_gtex >= 1).sum(axis=1)
n_samples_dev = df_dev.shape[1]
n_samples_gt1_dev = (df_dev >= 1).sum(axis=1)

fig, ax = plt.subplots(1, 1, figsize=(2, 2))
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


# ## 3. isoforms per family

# In[14]:


# number of isoforms vs gene expression, publications, and exons 
tpm_per_gene = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .sum()
                .groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
                .mean())
gn = tpm_per_gene.max(axis=1).rename('TPM - gene-level, max across GTEx tissues').to_frame()
gn['n_isoforms'] = gn.index.map(genes_gtex.value_counts())


# In[15]:


fam = load_tf_families()
gn['family'] = gn.index.map(fam)
gn['is_nuclear_receptor'] = (gn['family'] == 'Nuclear receptor')
gn.head()


# In[16]:


len(gn)


# In[17]:


len(gn[gn["n_isoforms"] > 1])


# In[18]:


gn.n_isoforms.mean()


# In[19]:


gn.sort_values(by="n_isoforms", ascending=False).head()


# In[20]:


fam_members = pd.DataFrame(gn['family'].value_counts()).reset_index()
fam_members_ov12 = fam_members[fam_members["family"] >= 12]
fam_members_ov12


# In[21]:


# plot distribution of isoforms by TPM
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(3, 3))

df = gn[gn["family"].isin(fam_members_ov12["index"])]
df_grp = df.groupby("family")["n_isoforms"].agg("median").reset_index()
x_order = list(df_grp.sort_values(by="n_isoforms", ascending=False)["family"])

sns.boxplot(data=gn[gn["family"].isin(fam_members_ov12["index"])], 
            x="family", y="n_isoforms", order=x_order, ax=axs[0])
sns.boxplot(data=gn[gn["family"].isin(fam_members_ov12["index"])], 
            x="family", y="TPM - gene-level, max across GTEx tissues", order=x_order, ax=axs[1])

axs[0].set_xlabel("")
axs[1].set_xlabel("")
axs[1].set_yscale("log")
axs[1].set_ylabel("max gene expression")
axs[1].set_xticklabels(x_order, rotation=30, ha="right", va="top")

plt.tight_layout()


# ## 4. downsample GTEx
# 
# GTEx has more samples per condition than Dev, but Dev has more conditions

# In[22]:


# conditions (body sites): gtex
len(metadata_gtex['body_site'].value_counts())


# In[23]:


# samples per body site: gtex
metadata_gtex['body_site'].value_counts()


# In[24]:


# conditions (body sites): dev
metadata_dev['body_site'] = metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']
len(metadata_dev['body_site'].value_counts())


# In[25]:


# samples per body site: dev
metadata_dev['body_site'].value_counts()


# ### loop through GTEx tissues and pick the # of samples by randomly matching to a dev dataset
# this is inherently unstable when sampling w/o replacement as will end up with times where there are more samps in the dev that you're randomly matching to than the gtex (rare but happens)

# In[26]:


# loop through gtex tissues
# pick number of samples according to dev dataset
# loop again
# make fake metadata file
n_samples_dev = df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']).value_counts().values
np.random.shuffle(n_samples_dev)
gtex_tissues = metadata_gtex['body_site'].value_counts().index.values

metadata_gtex_dummy = {}
for i, (n_samples, tissue) in enumerate(zip(n_samples_dev, itertools.cycle(gtex_tissues))):
    metadata_gtex_dummy[tissue + '_' + str(i)] = (metadata_gtex.loc[(metadata_gtex['body_site'] == tissue)
                                                                    & ~metadata_gtex.index.isin({s for samples in metadata_gtex_dummy.values() for s in samples}),
                                                                    :]
                                                        .sample(n_samples).index.values)

# TODO: check it is sampling with replacement and ends up same size as dev   
# NOTE: this block of code is unstable depending on seed

metadata_gtex_dummy = (pd.Series({v: k for k, vs in metadata_gtex_dummy.items() for v in vs}, name='body_site')
                         .to_frame())

if metadata_dev.shape[0] != metadata_gtex_dummy.shape[0]:
    raise UserWarning('Problem with downsampling code')
if sorted(n_samples_dev) != sorted(metadata_gtex_dummy.groupby('body_site').size().values):
    raise UserWarning('Problem with downsampling code')
if metadata_gtex_dummy.index.duplicated().any():
    raise UserWarning('Unexpected duplicates')


# In[27]:


metadata_gtex_dummy.shape


# In[28]:


len(metadata_gtex_dummy.body_site.unique())


# In[29]:


len(metadata_gtex_dummy.body_site.str.split("_", expand=True)[0].unique())


# In[30]:


metadata_dev.shape


# In[31]:


len(df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']).unique())


# In[112]:


df_dev.columns.map(metadata_dev['organism_part'] + ' ' + metadata_dev['dev_stage']).unique()


# In[119]:


tmp = metadata_dev.groupby(["organism_part", "dev_stage"])["BioSample"].agg("count").reset_index()
tmp.sort_values(by="BioSample")


# #### this dataframe is now the same shape as the dev data in both # of samples and # of "sites"
# gets to the same # of "sites" by re-sampling among GTEx tissues

# In[140]:


# write this file so we can load it in the DN section later
metadata_gtex_dummy.to_csv("../data/processed/metadata_gtex_dummy.csv")


# In[32]:


means_gtex_downsample = df_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean()


# ## 5. histograms: isoforms per gene + thresholded on expression

# ### GTEx: all

# In[33]:


# plot number of isoforms above 1 TPM

fig, axs = plt.subplots(2, 1, sharex=False, figsize=(3.5, 2.2))

n_iso = (means_gtex > 1).any(axis=1).groupby(genes_gtex).size()
x_max = n_iso.max()
xs = range(0, x_max + 1)
axs[0].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs], color="slategrey")

# label n
for h, x in zip([n_iso.value_counts().to_dict().get(x, 0) for x in xs], xs):
    if h == 0:
        continue
    axs[0].text(x, h, " %s" % h, rotation=90, fontsize=fontsize-2, ha="center", va="bottom",
                color="slategrey")

n_iso = (means_gtex > 1).any(axis=1).groupby(genes_gtex).sum()
axs[1].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])

axs[0].set_xticks(xs)
axs[0].set_xlabel("Unique annotated protein isoforms per gene")
axs[0].tick_params(axis='x', labelsize=fontsize-2)
axs[1].set_xticks(xs)
axs[1].tick_params(axis='x', labelsize=fontsize-2)
axs[1].set_xlabel('Unique protein isoforms per gene')
#axs[0].text(x=7, y=400, s='All isoforms')
axs[1].text(x=7, y=400, s='≥ 1 TPM in ≥ 1 tissue')

def num2pct(y):
    return (y / n_iso.shape[0]) * 100

def pct2num(y):
    return (y / 100) * n_iso.shape[0]

for ax in axs:
    ax.set_ylim(0, 800)
    ax.set_yticks(range(0, 800, 100), minor=True)
    ax.set_ylabel('TF genes')
    for pos in ['top', 'right', 'bottom']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    pctax = ax.secondary_yaxis('right', functions=(num2pct, pct2num))
    pctax.set_ylabel('% of TF genes')
    pctax.set_yticks(range(0, 46, 5), minor=True)
fig.savefig('../figures/n-isoforms-per-gene_by-1TPM-cutoff_hist-GTEx.pdf',
            bbox_inches='tight')


# ### GTEx: downsample

# In[34]:


# plot number of isoforms above 1 TPM

fig, axs = plt.subplots(2, 1, sharex=False, figsize=(3.3, 2.2))

n_iso = (means_gtex_downsample > 1).any(axis=1).groupby(genes_gtex).size()
x_max = n_iso.max()
xs = range(0, x_max + 1)
axs[0].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])

n_iso = (means_gtex_downsample > 1).any(axis=1).groupby(genes_gtex).sum()
axs[1].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])

axs[1].set_xticks(xs)
axs[1].set_xlabel('Unique protein isoforms per gene')
axs[0].text(x=7, y=400, s='All isoforms')
axs[1].text(x=7, y=400, s='≥ 1 TPM in at least one GTEx down-sampled tissue')

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
fig.savefig('../figures/n-isoforms-per-gene_by-1TPM-cutoff_hist-GTEx_downsamp.pdf',
            bbox_inches='tight')


# ### Dev

# In[35]:


# plot number of isoforms above 1 TPM

fig, axs = plt.subplots(2, 1, sharex=False, figsize=(3.3, 2.2))

n_iso = (means_dev > 1).any(axis=1).groupby(genes_dev).size()
x_max = n_iso.max()
xs = range(0, x_max + 1)
axs[0].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])

n_iso = (means_dev > 1).any(axis=1).groupby(genes_dev).sum()
axs[1].bar(x=xs, height=[n_iso.value_counts().to_dict().get(x, 0) for x in xs])

axs[1].set_xticks(xs)
axs[1].set_xlabel('Unique protein isoforms per gene')
axs[0].text(x=7, y=400, s='All isoforms')
axs[1].text(x=7, y=400, s='≥ 1 TPM in at least one dev tissue')

for ax in axs:
    ax.set_ylim(0, 800)
    ax.set_yticks(range(0, 800, 100), minor=True)
    ax.set_ylabel('TF genes_dev')
    for pos in ['top', 'right', 'bottom']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_tick_params(length=0)
    pctax = ax.secondary_yaxis('right', functions=(num2pct, pct2num))
    pctax.set_ylabel('% of TF genes_dev')
    pctax.set_yticks(range(0, 46, 5), minor=True)
fig.savefig('../figures/n-isoforms-per-gene_by-1TPM-cutoff_hist-GTEx_dev.pdf',
            bbox_inches='tight')


# ## 6. pie charts: num isoforms in each expression bucket

# ### GTEx: all

# In[36]:


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


# #### GTEx all: exclude testis

# In[37]:


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


# ### GTEx: downsample

# In[38]:


p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(alt_isos), :] >= 1).any(axis=1).sum()
n = means_gtex_downsample.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(alt_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex_downsample.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue')

p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(ref_isos), :] >= 1).any(axis=1).sum()
n = means_gtex_downsample.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(ref_isos), :] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex_downsample.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = means_gtex_downsample.max(axis=1)

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
fig.savefig('../figures/downsampled-GTEx-control-max-expression_by-reference-vs-alternative_pie.pdf',
            bbox_inches='tight')


# In[39]:


# Exlclude testis
cols = [c for c in means_gtex_downsample.columns if 'Testis' not in c]
p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(alt_isos), cols] >= 1).any(axis=1).sum()
n = means_gtex_downsample.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(alt_isos), cols] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex_downsample.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue')

p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(ref_isos), cols] >= 1).any(axis=1).sum()
n = means_gtex_downsample.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue')

p = (means_gtex_downsample.loc[means_gtex_downsample.index.isin(ref_isos), cols] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_gtex_downsample.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = means_gtex_downsample.loc[:, cols].max(axis=1)

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
fig.savefig('../figures/downsampled-GTEx-max-expression_by-reference-vs-alternative_no-testis_pie.pdf',
            bbox_inches='tight')


# ### Dev

# In[40]:


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


# #### Dev: no testis

# In[41]:


cols = [c for c in means_dev.columns if 'testis' not in c]
p = (means_dev.loc[means_dev.index.isin(alt_isos), cols] >= 1).any(axis=1).sum()
n = means_dev.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 1 TPM in at least one tissue/dev stage')

p = (means_dev.loc[means_dev.index.isin(alt_isos), cols] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_dev.index.isin(alt_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) alternative isoforms ≥ 5 TPM in at least one tissue/dev stage')

p = (means_dev.loc[means_dev.index.isin(ref_isos), cols] >= 1).any(axis=1).sum()
n = means_dev.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 1 TPM in at least one tissue/dev stage')

p = (means_dev.loc[means_dev.index.isin(ref_isos), cols] >= np.log2(5 + 1)).any(axis=1).sum()
n = means_dev.index.isin(ref_isos).sum()
print(f'{p} out of {n} ({p/n:.0%}) MANE select isoforms ≥ 5 TPM in at least one tissue/dev stage')

# pie chart
fig, axs = plt.subplots(1, 2)
max_vals = means_dev.loc[:, cols].max(axis=1)

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
fig.savefig('../figures/developmental-max-expression_by-reference-vs-alternative_no-testis_pie.pdf',
            bbox_inches='tight')


# ## 7. distribution of TPMs across isoforms

# ### GTEx: all

# In[42]:


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


# ### GTEx: downsample

# In[43]:


# plot distribution of isoforms by TPM
fig, axs = plt.subplots(2, 1, sharex=True)
n_bins = 110
x_max = 11
axs[1].hist(means_gtex_downsample.max(axis=1)[means_gtex_downsample.index.isin(alt_isos)],
            bins=n_bins,
            range=(0, x_max))
axs[0].hist(means_gtex_downsample.max(axis=1)[means_gtex_downsample.index.isin(ref_isos)],
            bins=n_bins,
            range=(0, x_max))
axs[0].text(x=8, y=axs[0].get_ylim()[1] * 0.95, s='MANE select isoforms')
axs[1].text(x=8, y=axs[1].get_ylim()[1] * 0.95, s='Alternative isoforms')
for ax in axs:
    ax.set_ylabel('Number of isoforms')
    for loc in ['top', 'right']:
        ax.spines[loc].set_visible(False)
axs[1].set_xlabel('Expression, max across GTEx tissues – TPM')
x_ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
axs[1].set_xticks([np.log2(x + 1) for x in x_ticks])
axs[1].set_xticklabels(x_ticks)
fig.savefig('../figures/expression_downsampled-GTEX-control_GENCODE-isoforms_by-reference-vs-alternative.pdf',
            bbox_inches='tight')


# ### Dev

# In[44]:


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


# ## 8. distribution of ratios across isoforms

# ### GTEx: all

# In[45]:


# percentage of alternative isoform
# plot distribution of fraction of gene expression for ref and alt

# fraction where gene tpm > 1

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


# ### GTEx: downsample

# In[46]:


genes_gtex


# In[47]:


# percentage of alternative isoform
# plot distribution of fraction of gene expression for ref and alt

# has to be fraction where isoform TPM is at least 1, right (fill na with 0)

per_gene_gtex = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .transform('sum'))
f_gtex_downsample = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1)
        .mean())
f_gtex_downsample = f_gtex_downsample * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex_downsample = f_gtex_downsample * 100
fig, axs = plt.subplots(2, 1, sharex=True)
n_bins=100
axs[0].hist(f_gtex_downsample.max(axis=1).loc[f_gtex_downsample.index.isin(ref_isos)].values,
            range=(0, 100),
            bins=n_bins)
axs[1].hist(f_gtex_downsample.max(axis=1).loc[f_gtex_downsample.index.isin(alt_isos)].values,
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


# ### Dev

# In[48]:


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


# #### aside: calculating isoforms that switch

# In[49]:


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


# ## 9. ref v alt 2D heatmaps: max expression

# In[50]:


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


# In[51]:


ref_alt_map_nonan[ref_alt_map_nonan["gene"] == "NKX2-5"]


# ### GTEx: all

# In[52]:


means_gtex["max_gtex"] = means_gtex.max(axis=1)
means_gtex["min_gtex"] = means_gtex.min(axis=1)

# max out anything above 11 (2000 tpm) to make plots more readable, as luke did above
means_gtex[means_gtex["max_gtex"] > 11] = 11

print(means_gtex["max_gtex"].max())
print(means_gtex["max_gtex"].min())
means_gtex_ri = means_gtex.reset_index()
means_gtex_ri["UID_rep"] = means_gtex_ri["UID"].str.replace("_", "|")


# In[53]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(means_gtex_ri[["UID_rep", "max_gtex", "min_gtex"]], left_on="ref", 
                                            right_on="UID_rep", how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(means_gtex_ri[["UID_rep", "max_gtex", "min_gtex"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="inner")


# In[54]:


fig = plt.figure(figsize=(2, 1.5))

ax = sns.histplot(data=ref_alt_map_nonan, x="max_gtex_ref", y="max_gtex_alt",
                  bins=30, cbar=True, cbar_kws={"label": "# isoform pairs",
                                                "ticks": [0, 5, 10, 15, 20, 25, 30]}, cmap="rocket_r",
                  vmin=0, vmax=30)

ax.set_xlim((-0.3, 11.5))
ax.set_ylim((-0.3, 11.5))
ax.set_xlabel("max expression of ref")
ax.set_ylabel("max expression of alt")
ax.set_title("GTEx dataset\n(n=%s ref/alt pairs)" % len(ref_alt_map_nonan))

ticks = [0, 1, 5, 10, 100, 400, 2000]
ticklabels = [0, 1, 5, 10, 100, 400, "2000+"]
ax.set_xticks([np.log2(x + 1) for x in ticks])
ax.set_xticklabels(ticklabels)
ax.set_yticks([np.log2(y + 1) for y in ticks])
ax.set_yticklabels(ticklabels)
ax.tick_params(axis='x', labelsize=fontsize-2)
ax.tick_params(axis='y', labelsize=fontsize-2)

cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "5", "10", "15", "20", "25", "30+"])

# find num where ref > alt
ra = len(ref_alt_map_nonan[ref_alt_map_nonan["max_gtex_ref"] > ref_alt_map_nonan["max_gtex_alt"]])
print(ra)
#ax.text(8.5, 7, "%s\n(%s%%)" % (ra, round(ra/len(ref_alt_map_nonan), 2)*100), ha="left", va="bottom")

# find num where alt > ref
ar = len(ref_alt_map_nonan[ref_alt_map_nonan["max_gtex_ref"] < ref_alt_map_nonan["max_gtex_alt"]])
print(ar)
#ax.text(9, 10.5, "%s\n(%s%%)" % (ar, round(ar/len(ref_alt_map_nonan), 2)*100), ha="right", va="top")

ax.plot([-0.3,11.5], [-0.3, 11.5], color="black", linestyle="dashed")

fig.savefig('../figures/expression-scatter-ref_v_alt-gtex.pdf',
            bbox_inches='tight')


# ### GTEx: downsampled

# In[55]:


means_gtex_downsample["max_gtex_downsample"] = means_gtex_downsample.max(axis=1)
means_gtex_downsample["min_gtex_downsample"] = means_gtex_downsample.min(axis=1)

# max out anything above 11 (2000 tpm) to make plots more readable, as luke did above
means_gtex_downsample[means_gtex_downsample["max_gtex_downsample"] > 11] = 11

print(means_gtex_downsample["max_gtex_downsample"].max())
print(means_gtex_downsample["max_gtex_downsample"].min())
means_gtex_downsample_ri = means_gtex_downsample.reset_index()
means_gtex_downsample_ri["UID_rep"] = means_gtex_downsample_ri["UID"].str.replace("_", "|")


# In[56]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(means_gtex_downsample_ri[["UID_rep", "max_gtex_downsample",
                                                                      "min_gtex_downsample"]], 
                                            left_on="ref", right_on="UID_rep", how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(means_gtex_downsample_ri[["UID_rep", "max_gtex_downsample",
                                                                      "min_gtex_downsample"]], 
                                            left_on="alt", right_on="UID_rep", suffixes=("_ref", "_alt"), how="inner")


# In[57]:


fig = plt.figure(figsize=(2, 1.5))

ax = sns.histplot(data=ref_alt_map_nonan, x="max_gtex_downsample_ref", y="max_gtex_downsample_alt",
                  bins=30, cbar=True, cbar_kws={"label": "# isoform pairs",
                                                "ticks": [0, 5, 10, 15, 20, 25, 30]}, cmap="rocket_r",
                  vmin=0, vmax=30)

ax.set_xlim((-0.3, 11.5))
ax.set_ylim((-0.3, 11.5))
ax.set_xlabel("max tpm of ref")
ax.set_ylabel("max tpm of alt")
ax.set_title("(n=%s ref/alt pairs)" % len(ref_alt_map_nonan))

ticks = [0, 1, 5, 10, 100, 400, 2000]
ticklabels = [0, 1, 5, 10, 100, 400, "2000+"]
ax.set_xticks([np.log2(x + 1) for x in ticks])
ax.set_xticklabels(ticklabels)
ax.set_yticks([np.log2(y + 1) for y in ticks])
ax.set_yticklabels(ticklabels)
ax.tick_params(axis='x', labelsize=fontsize-2)
ax.tick_params(axis='y', labelsize=fontsize-2)

cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "5", "10", "15", "20", "25", "30+"])

# find num where ref > alt
ra = len(ref_alt_map_nonan[ref_alt_map_nonan["max_gtex_downsample_ref"] > ref_alt_map_nonan["max_gtex_downsample_alt"]])
print(ra)
#ax.text(8.5, 7, "%s\n(%s%%)" % (ra, round(ra/len(ref_alt_map_nonan), 2)*100), ha="left", va="bottom")

# find num where alt > ref
ar = len(ref_alt_map_nonan[ref_alt_map_nonan["max_gtex_downsample_ref"] < ref_alt_map_nonan["max_gtex_downsample_alt"]])
print(ar)
#ax.text(9, 10.5, "%s\n(%s%%)" % (ar, round(ar/len(ref_alt_map_nonan), 2)*100), ha="right", va="top")

ax.plot([-0.3,11.5], [-0.3, 11.5], color="black", linestyle="dashed")

fig.savefig('../figures/expression-scatter-ref_v_alt-gtex-downsample.pdf',
            bbox_inches='tight')


# ### Dev

# In[58]:


means_dev["max_dev"] = means_dev.max(axis=1)
means_dev["min_dev"] = means_dev.min(axis=1)

# max out anything above 11 (2000 tpm) to make plots more readable, as luke did above
means_dev[means_dev["max_dev"] > 11] = 11

print(means_dev["max_dev"].max())
print(means_dev["max_dev"].min())
means_dev_ri = means_dev.reset_index()
means_dev_ri["UID_rep"] = means_dev_ri["UID"].str.replace("_", "|")


# In[59]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(means_dev_ri[["UID_rep", "max_dev", "min_dev"]], left_on="ref", 
                                            right_on="UID_rep", how="inner")
ref_alt_map_nonan = ref_alt_map_nonan.merge(means_dev_ri[["UID_rep", "max_dev", "min_dev"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="inner")


# In[60]:


fig = plt.figure(figsize=(2, 1.5))

ax = sns.histplot(data=ref_alt_map_nonan, x="max_dev_ref", y="max_dev_alt",
                  bins=30, cbar=True, cbar_kws={"label": "# isoform pairs",
                                                "ticks": [0, 5, 10, 15, 20, 25, 30]}, cmap="rocket_r",
                  vmin=0, vmax=30)

ax.set_xlim((-0.3, 11.5))
ax.set_ylim((-0.3, 11.5))
ax.set_xlabel("max tpm of ref")
ax.set_ylabel("max tpm of alt")
ax.set_title("n=%s ref/alt pairs" % len(ref_alt_map_nonan))

ticks = [0, 1, 5, 10, 100, 400, 2000]
ticklabels = [0, 1, 5, 10, 100, 400, "2000+"]
ax.set_xticks([np.log2(x + 1) for x in ticks])
ax.set_xticklabels(ticklabels)
ax.set_yticks([np.log2(y + 1) for y in ticks])
ax.set_yticklabels(ticklabels)
ax.tick_params(axis='x', labelsize=fontsize-2)
ax.tick_params(axis='y', labelsize=fontsize-2)

cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "5", "10", "15", "20", "25", "30+"])

# find num where ref > alt
ra = len(ref_alt_map_nonan[ref_alt_map_nonan["max_dev_ref"] > ref_alt_map_nonan["max_dev_alt"]])
print(ra)
#ax.text(8.5, 7, "%s\n(%s%%)" % (ra, round(ra/len(ref_alt_map_nonan), 2)*100), ha="left", va="bottom")

# find num where alt > ref
ar = len(ref_alt_map_nonan[ref_alt_map_nonan["max_dev_ref"] < ref_alt_map_nonan["max_dev_alt"]])
print(ar)
#ax.text(9, 10.5, "%s\n(%s%%)" % (ar, round(ar/len(ref_alt_map_nonan), 2)*100), ha="right", va="top")

ax.plot([-0.3,11.5], [-0.3, 11.5], color="black", linestyle="dashed")

fig.savefig('../figures/expression-scatter-ref_v_alt-dev.pdf',
            bbox_inches='tight')


# ## 10. per isoform: max v min ratio
# 
# removing NaNs - not counting anything where *gene* expression < 1

# ### GTEx: all

# In[61]:


print(len(f_gtex))
f_gtex["max_ratio_gtex"] = f_gtex.max(axis=1)
f_gtex["min_ratio_gtex"] = f_gtex.min(axis=1)
f_gtex_nonan = f_gtex[(~pd.isnull(f_gtex["max_ratio_gtex"])) & (~pd.isnull(f_gtex["min_ratio_gtex"]))]
print(len(f_gtex_nonan))

f_gtex_ri = f_gtex_nonan.reset_index()
f_gtex_ri["UID_rep"] = f_gtex_ri["UID"].str.replace("_", "|")


# In[62]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(f_gtex_ri[["UID_rep", "max_ratio_gtex", "min_ratio_gtex"]], left_on="ref", 
                                            right_on="UID_rep", how="left")
ref_alt_map_nonan = ref_alt_map_nonan.merge(f_gtex_ri[["UID_rep", "max_ratio_gtex", "min_ratio_gtex"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="left")


# In[63]:


fig = plt.figure(figsize=(2, 1.5))

df = ref_alt_map_nonan[["ref", "min_ratio_gtex_ref", "max_ratio_gtex_ref"]].drop_duplicates()
df = df[(~pd.isnull(df["min_ratio_gtex_ref"])) & (~pd.isnull(df["max_ratio_gtex_ref"]))]

ax = sns.histplot(data=df, x="min_ratio_gtex_ref", y="max_ratio_gtex_ref",
                  bins=30, cbar=True, cmap="rocket_r", vmin=0, vmax=120, cbar_kws={"label": "# isoform pairs",
                                                                                   "ticks": [0, 20, 40, 60,
                                                                                             80, 100, 120]})
cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "20", "40", "60", "80", "100", "120+"])
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("max isoform fraction")
ax.set_ylabel("min isoform fraction")
ax.set_title("n=%s ref isoforms" % len(df))

fig.savefig('../figures/expression-ratio-scatter-ref-gtex.pdf',
            bbox_inches='tight')


# In[64]:


fig = plt.figure(figsize=(2, 1.5))

df = ref_alt_map_nonan[["alt", "min_ratio_gtex_alt", "max_ratio_gtex_alt"]].drop_duplicates()
df = df[(~pd.isnull(df["min_ratio_gtex_alt"])) & (~pd.isnull(df["max_ratio_gtex_alt"]))]

ax = sns.histplot(data=ref_alt_map_nonan, x="min_ratio_gtex_alt", y="max_ratio_gtex_alt",
                  bins=30, cbar=True, cmap="rocket_r", vmin=0, vmax=120, cbar_kws={"label": "# isoform pairs",
                                                                                   "ticks": [0, 20, 40, 60,
                                                                                             80, 100, 120]})
cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "20", "40", "60", "80", "100", "120+"])
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("min isoform fraction")
ax.set_ylabel("max isoform fraction")
ax.set_title("n=%s alt isoforms" % len(df))
fig.savefig('../figures/expression-ratio-scatter-alt-gtex.pdf',
            bbox_inches='tight')


# ### GTEx: downsample

# In[65]:


print(len(f_gtex_downsample))
f_gtex_downsample["max_ratio_gtex_downsample"] = f_gtex_downsample.max(axis=1)
f_gtex_downsample["min_ratio_gtex_downsample"] = f_gtex_downsample.min(axis=1)
f_gtex_downsample_nonan = f_gtex_downsample[(~pd.isnull(f_gtex_downsample["max_ratio_gtex_downsample"])) & 
                                            (~pd.isnull(f_gtex_downsample["min_ratio_gtex_downsample"]))]
print(len(f_gtex_downsample_nonan))

f_gtex_downsample_ri = f_gtex_downsample_nonan.reset_index()
f_gtex_downsample_ri["UID_rep"] = f_gtex_downsample_ri["UID"].str.replace("_", "|")


# In[66]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(f_gtex_downsample_ri[["UID_rep", "max_ratio_gtex_downsample", 
                                                       "min_ratio_gtex_downsample"]], left_on="ref", 
                                            right_on="UID_rep", how="left")
ref_alt_map_nonan = ref_alt_map_nonan.merge(f_gtex_downsample_ri[["UID_rep", "max_ratio_gtex_downsample", 
                                                                  "min_ratio_gtex_downsample"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="left")


# In[67]:


fig = plt.figure(figsize=(2, 1.5))

df = ref_alt_map_nonan[["ref", "min_ratio_gtex_downsample_ref", "max_ratio_gtex_downsample_ref"]].drop_duplicates()
df = df[(~pd.isnull(df["min_ratio_gtex_downsample_ref"])) & (~pd.isnull(df["max_ratio_gtex_downsample_ref"]))]

n_switches = df[(df["min_ratio_gtex_downsample_ref"] < 10) & (df["max_ratio_gtex_downsample_ref"] > 90)]
n_off = df[(df["min_ratio_gtex_downsample_ref"] < 10) & (df["max_ratio_gtex_downsample_ref"] < 10)]
print(len(n_switches))
p_switches_ref_gtex_ds = len(n_switches)/len(df)
p_off_ref_gtex_ds = len(n_off)/len(df)

ax = sns.histplot(data=df, x="min_ratio_gtex_downsample_ref", y="max_ratio_gtex_downsample_ref",
                  bins=30, cbar=True, cmap="mako_r", vmin=0, vmax=120, cbar_kws={"label": "# isoform pairs",
                                                                                   "ticks": [0, 20, 40, 60,
                                                                                             80, 100, 120]})
cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "20", "40", "60", "80", "100", "120+"])
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("min isoform fraction")
ax.set_ylabel("max isoform fraction")
ax.set_title("n=%s ref isoforms" % len(df))
fig.savefig('../figures/expression-ratio-scatter-ref-gtex-downsample.pdf',
            bbox_inches='tight')


# In[68]:


fig = plt.figure(figsize=(2, 1.5))

df = ref_alt_map_nonan[["alt", "min_ratio_gtex_downsample_alt", "max_ratio_gtex_downsample_alt"]].drop_duplicates()
df = df[(~pd.isnull(df["min_ratio_gtex_downsample_alt"])) & (~pd.isnull(df["max_ratio_gtex_downsample_alt"]))]

n_switches = df[(df["min_ratio_gtex_downsample_alt"] < 10) & (df["max_ratio_gtex_downsample_alt"] > 90)]
n_off = df[(df["min_ratio_gtex_downsample_alt"] < 10) & (df["max_ratio_gtex_downsample_alt"] < 10)]
print(len(n_switches))
p_switches_alt_gtex_ds = len(n_switches)/len(df)
p_off_alt_gtex_ds = len(n_off)/len(df)

ax = sns.histplot(data=ref_alt_map_nonan, x="min_ratio_gtex_downsample_alt", y="max_ratio_gtex_downsample_alt",
                  bins=30, cbar=True, cmap="mako_r", vmin=0, vmax=120, cbar_kws={"label": "# isoform pairs",
                                                                                   "ticks": [0, 20, 40, 60,
                                                                                             80, 100, 120]})
cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "20", "40", "60", "80", "100", "120+"])
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("min isoform fraction")
ax.set_ylabel("max isoform fraction")
ax.set_title("n=%s alt isoforms" % len(df))

# add lines to distinguish events
ax.plot([10, 0], [10, 10], linestyle="dotted", color="black")
ax.plot([10, 10], [0, 10], linestyle="dotted", color="black")
ax.plot([0, 10], [90, 90], linestyle="dotted", color="black")
ax.plot([10, 10], [90, 100], linestyle="dotted", color="black")
ax.text(10, 5, " low", ha="left", va="center", fontstyle="italic", color="slategrey")
ax.text(10, 95, " switch", ha="left", va="center", fontstyle="italic", color=sns.color_palette("mako")[1])

fig.savefig('../figures/expression-ratio-scatter-alt-gtex-downsample.pdf',
            bbox_inches='tight')


# ### Dev

# In[69]:


print(len(f_dev))
f_dev["max_ratio_dev"] = f_dev.max(axis=1)
f_dev["min_ratio_dev"] = f_dev.min(axis=1)
f_dev_nonan = f_dev[(~pd.isnull(f_dev["max_ratio_dev"])) & (~pd.isnull(f_dev["min_ratio_dev"]))]
print(len(f_dev_nonan))

f_dev_ri = f_dev_nonan.reset_index()
f_dev_ri["UID_rep"] = f_dev_ri["UID"].str.replace("_", "|")


# In[70]:


ref_alt_map_nonan = ref_alt_map_nonan.merge(f_dev_ri[["UID_rep", "max_ratio_dev", "min_ratio_dev"]], left_on="ref", 
                                            right_on="UID_rep", how="left")
ref_alt_map_nonan = ref_alt_map_nonan.merge(f_dev_ri[["UID_rep", "max_ratio_dev", "min_ratio_dev"]], left_on="alt", 
                                            right_on="UID_rep", suffixes=("_ref", "_alt"), how="left")


# In[71]:


fig = plt.figure(figsize=(2, 1.5))

df = ref_alt_map_nonan[["ref", "min_ratio_dev_ref", "max_ratio_dev_ref"]].drop_duplicates()
df = df[(~pd.isnull(df["min_ratio_dev_ref"])) & (~pd.isnull(df["max_ratio_dev_ref"]))]

n_switches = df[(df["min_ratio_dev_ref"] < 10) & (df["max_ratio_dev_ref"] > 90)]
n_off = df[(df["min_ratio_dev_ref"] < 10) & (df["max_ratio_dev_ref"] < 10)]
print(len(n_switches))
p_switches_ref_dev = len(n_switches)/len(df)
p_off_ref_dev = len(n_off)/len(df)

ax = sns.histplot(data=ref_alt_map_nonan, x="min_ratio_dev_ref", y="max_ratio_dev_ref",
                  bins=30, cbar=True, cmap="mako_r", vmin=0, vmax=120, cbar_kws={"label": "# isoform pairs",
                                                                                   "ticks": [0, 20, 40, 60,
                                                                                             80, 100, 120]})
cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "20", "40", "60", "80", "100", "120+"])
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("min isoform fraction")
ax.set_ylabel("max isoform fraction")
ax.set_title("n=%s ref isoforms" % len(df))

# # add lines to distinguish events
# ax.plot([20, 0], [20, 20], linestyle="dotted", color="black")
# ax.plot([20, 20], [0, 20], linestyle="dotted", color="black")
# ax.plot([0, 20], [70, 70], linestyle="dotted", color="black")
# ax.plot([20, 20], [70, 100], linestyle="dotted", color="black")
# ax.text(20, 10, " low", ha="left", va="center", fontstyle="italic", color="slategrey")
# ax.text(20, 85, " switch", ha="left", va="center", fontstyle="italic", color=sns.color_palette("mako")[0])

fig.savefig('../figures/expression-ratio-scatter-ref-dev.pdf',
            bbox_inches='tight')


# In[72]:


fig = plt.figure(figsize=(2, 1.5))

df = ref_alt_map_nonan[["alt", "min_ratio_dev_alt", "max_ratio_dev_alt"]].drop_duplicates()
df = df[(~pd.isnull(df["min_ratio_dev_alt"])) & (~pd.isnull(df["max_ratio_dev_alt"]))]

n_switches = df[(df["min_ratio_dev_alt"] < 10) & (df["max_ratio_dev_alt"] > 90)]
n_off = df[(df["min_ratio_dev_alt"] < 10) & (df["max_ratio_dev_alt"] < 10)]
print(len(n_switches))
p_switches_alt_dev = len(n_switches)/len(df)
p_off_alt_dev = len(n_off)/len(df)

ax = sns.histplot(data=ref_alt_map_nonan, x="min_ratio_dev_alt", y="max_ratio_dev_alt",
                  bins=30, cbar=True, cmap="mako_r", vmin=0, vmax=120, cbar_kws={"label": "# isoform pairs",
                                                                                   "ticks": [0, 20, 40, 60,
                                                                                             80, 100, 120]})
cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["0", "20", "40", "60", "80", "100", "120+"])
ax.set_xlim((-2, 102))
ax.set_ylim((-2, 102))
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("min isoform fraction")
ax.set_ylabel("max isoform fraction")
ax.set_title("n=%s alt isoforms" % len(df))

# add lines to distinguish events
ax.plot([10, 0], [10, 10], linestyle="dotted", color="black")
ax.plot([10, 10], [0, 10], linestyle="dotted", color="black")
ax.plot([0, 10], [90, 90], linestyle="dotted", color="black")
ax.plot([10, 10], [90, 100], linestyle="dotted", color="black")
ax.text(10, 5, " low", ha="left", va="center", fontstyle="italic", color="slategrey")
ax.text(10, 95, " switch", ha="left", va="center", fontstyle="italic", color=sns.color_palette("mako")[1])

fig.savefig('../figures/expression-ratio-scatter-alt-dev.pdf',
            bbox_inches='tight')


# In[73]:


bar = pd.DataFrame.from_dict({"gtex_ds_ref": {"switch": p_switches_ref_gtex_ds*100, "low": p_off_ref_gtex_ds*100},
                              "gtex_ds_alt": {"switch": p_switches_alt_gtex_ds*100, "low": p_off_alt_gtex_ds*100},
                              "dev_ref": {"switch": p_switches_ref_dev*100, "low": p_off_ref_dev*100},
                              "dev_alt": {"switch": p_switches_alt_dev*100, "low": p_off_alt_dev*100}}, 
                             orient="index").reset_index()
bar["shift"] = 100-(bar["switch"]+bar["low"])
bar = bar[["index", "low", "switch", "shift"]]
bar


# In[74]:


palette = {"low": "lightgrey",
           "switch": sns.color_palette("mako")[1],
           "shift": sns.color_palette("mako")[5]}
palette


# In[75]:


gtex_bar = bar[bar["index"].str.contains("gtex")]
ax = gtex_bar.plot.bar(x="index", stacked=True, color=palette.values(), figsize=(1.4, 1.4))
ax.set_ylabel("% of isoforms")
ax.set_xlabel("")
#ax.set_title("GTEx")

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))
ax.set_xticklabels(["ref", "alt"], ha="right", va="top", rotation=30)

plt.savefig('../figures/expression-switch-bar-gtex.pdf',
            bbox_inches='tight')


# In[76]:


dev_bar = bar[bar["index"].str.contains("dev")]
ax = dev_bar.plot.bar(x="index", stacked=True, color=palette.values(), figsize=(1.4, 1.4))
ax.set_ylabel("% of isoforms")
ax.set_xlabel("")
#ax.set_title("dev")

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))
ax.set_xticklabels(["ref", "alt"], ha="right", va="top", rotation=30)

plt.savefig('../figures/expression-switch-bar-dev.pdf',
            bbox_inches='tight')


# ### example plot: TF gene whose isoform ratios change across tissues

# In[77]:


tmp = ref_alt_map_nonan
tmp["mm_gtex_ds_ref"] = tmp["max_ratio_gtex_downsample_ref"]-tmp["min_ratio_gtex_downsample_ref"]
tmp["mm_gtex_ds_alt"] = tmp["max_ratio_gtex_downsample_alt"]-tmp["min_ratio_gtex_downsample_alt"]
tmp["mm_dev_ref"] = tmp["max_ratio_dev_ref"]-tmp["min_ratio_dev_ref"]
tmp["mm_dev_alt"] = tmp["max_ratio_dev_alt"]-tmp["min_ratio_dev_alt"]
tmp["dg_ref"] = tmp["mm_dev_ref"]-tmp["mm_gtex_ds_ref"]
tmp["dg_alt"] = tmp["mm_dev_alt"]-tmp["mm_gtex_ds_alt"]
#tmp.sort_values(by="dg_alt", ascending=False).head(30)


# In[78]:


def developmental_tissue_expression_plot(gene_name, figsize, ylim, means, cols, fig_suffix):
    n_isos = len(means.loc[genes == gene_name])
    palette = sns.color_palette("Spectral", as_cmap=False, n_colors=n_isos)
    fig, axes = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(figsize)
    ### bar chart ###
    (means.loc[genes == gene_name, cols]
          .T
          .plot.bar(ax=axes[0],
                    legend=False,
                    width=0.7,
                    color=list(palette)))
    ### percentages ###
    raw_means = 2 ** means.loc[genes == gene_name, cols] - 1.
    (raw_means.div(raw_means.sum(axis=0))
              .T.plot.bar(ax=axes[1], 
                          stacked=True,
                          legend=False,
                          color=list(palette)))
    axes[0].set_ylabel('log2(tpm + 1)\n')
    axes[0].set_ylim(ylim)
    axes[1].set_ylabel('percent')
    axes[1].set_yticklabels(['{:.0%}'.format(t) for t in axes[1].get_yticks()])
    axes[1].legend(loc='lower left', bbox_to_anchor=(1, 0))
    axes[0].axhline(y=1, color='black', linewidth=0.5, linestyle="dashed")
    plt.subplots_adjust(hspace=0.25)
    plt.savefig('../figures/expression_' + gene_name + '_' + fig_suffix + '.pdf',
                bbox_inches='tight')


# In[79]:


liver_cols = [x for x in means_dev.columns if "liver" in x]
kidney_cols = [x for x in means_dev.columns if "kidney" in x]
developmental_tissue_expression_plot("HIF1A", (6, 1.75), (0, 6), means_dev, liver_cols + kidney_cols, 
                                     "means_dev_liver_kidney")


# In[80]:


liver_cols = [x for x in means_gtex.columns if "Liver" in x]
kidney_cols = [x for x in means_gtex.columns if "Kidney" in x]
developmental_tissue_expression_plot("HIF1A", (0.35, 1.75), (0, 6), means_gtex, liver_cols + kidney_cols, 
                                     "means_gtex_liver_kidney")


# In[81]:


tmp_nn = tmp[(~pd.isnull(tmp["min_ratio_dev_alt"])) & (~pd.isnull(tmp["max_ratio_dev_alt"]))]
tmp_srt = tmp_nn[(tmp_nn["max_ratio_dev_alt"] > 90) & (tmp_nn["min_ratio_dev_alt"] < 10)]
tmp_srt[["gene", "ref", "alt", "max_ratio_dev_alt", 
         "min_ratio_dev_alt", "max_dev_alt", "min_dev_alt"]].sort_values(by="max_ratio_dev_alt", ascending=False).head(20)


# In[138]:


notestis_cols = [x for x in means_dev.columns if "testis" not in x]
notestis_cols = [x for x in notestis_cols if "_dev" not in x]
notestis_cols = [x for x in notestis_cols if "ovary" not in x]
developmental_tissue_expression_plot("PKNOX1", (5, 1.75), (0, 6), means_dev, notestis_cols, 
                                     "means_dev_notestis")


# In[83]:


heart_cols = [x for x in means_dev.columns if "heart" in x]
ovary_cols = [x for x in means_dev.columns if "ovary" in x]
developmental_tissue_expression_plot("HEY2", (4, 1.75), (0, 6), means_dev, heart_cols + ovary_cols, 
                                     "means_dev_heart_ovary")


# In[84]:


heart_cols = [x for x in means_gtex.columns if "Heart" in x]
ovary_cols = [x for x in means_gtex.columns if "Ovary" in x]
developmental_tissue_expression_plot("HEY2", (0.5, 1.75), (0, 6), means_gtex, heart_cols + ovary_cols, 
                                     "means_gtex_heart_ovary")


# In[85]:


brain_cols = [x for x in means_dev.columns if "brain" in x]
developmental_tissue_expression_plot("NFKB1", (4, 1.75), (0, 6), means_dev, brain_cols, 
                                     "means_dev_brain")


# In[86]:


brain_cols = [x for x in means_gtex.columns if "Brain" in x]
developmental_tissue_expression_plot("NFKB1", (1.75, 1.75), (0, 6), means_gtex, brain_cols, 
                                     "means_gtex_brain")


# In[87]:


fig, ax = plt.subplots(figsize=(5, 0.5))

tfs["HEY2"].exon_diagram(ax=ax)

fig.savefig("../figures/HEY1_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[88]:


fig, ax = plt.subplots(figsize=(7.5, 1))

tfs["NFKB1"].exon_diagram(ax=ax)

fig.savefig("../figures/NFKB1_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[89]:


fig, ax = plt.subplots(figsize=(6, 1.5))

tfs["HIF1A"].exon_diagram(ax=ax)


# In[136]:


fig, ax = plt.subplots(figsize=(5, 0.75))

tfs["CREB5"].exon_diagram(ax=ax)

fig.savefig("../figures/CREB5_exon_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[137]:


fig, ax = plt.subplots(figsize=(5, 0.75))

tfs["CREB5"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)

fig.savefig("../figures/CREB5_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[90]:


fig, ax = plt.subplots(figsize=(3, 1))

tfs["HEY2"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)

fig.savefig("../figures/HEY1_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[91]:


fig, ax = plt.subplots(figsize=(7, 1.5))

tfs["NFKB1"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)

fig.savefig("../figures/NFKB1_protein_diagram.pdf", bbox_inches="tight", dpi="figure")


# In[92]:


fig, ax = plt.subplots(figsize=(6, 1.5))

tfs["HIF1A"].protein_diagram(only_cloned_isoforms=False, draw_legend=False, ax=ax)


# In[124]:


ss_alt = len(ref_alt_map_nonan[(ref_alt_map_nonan["max_ratio_gtex_downsample_alt"] >= 10)].gene.unique())
ss_alt


# In[125]:


ss_alt = len(ref_alt_map_nonan[(ref_alt_map_nonan["max_ratio_dev_alt"] >= 10)].gene.unique())
ss_alt


# In[126]:


len(ref_alt_map_nonan.gene.unique())


# In[128]:


864/909


# ### # make domain figure - move this into domain notebook at some point

# In[93]:


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

tot_ins_eff = []
tot_perc_ins_eff = []
tot_dd_eff = []
tot_perc_dd_eff = []
tot_f_eff = []
tot_perc_f_eff = []

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
        
        ins_eff = 0
        perc_ins_eff = 0
        dd_eff = 0
        perc_dd_eff = 0
        f_eff = 0
        perc_f_eff = 0
    else:
        aa_ftr_alt = aa_ftr[aa_ftr["alt_iso"] == alt]
        
        # separate pfam and effector domains
        pfam = aa_ftr_alt[aa_ftr_alt["category"] == "Pfam_domain"]
        eff = aa_ftr_alt[aa_ftr_alt["category"] == "effector_domain"]
        
        if len(pfam) > 0:
            pfam_grp = pfam.groupby("alt_iso")[["deletion", "insertion", "frameshift"]].agg("sum").reset_index()
        
            ins_dom = pfam_grp.insertion.iloc[0]
            perc_ins_dom = ins_dom/len(pp_str)*100
            dd_dom = pfam_grp.deletion.iloc[0]
            perc_dd_dom = dd_dom/len(pp_str)*100
            f_dom = pfam_grp.frameshift.iloc[0]
            perc_f_dom = f_dom/len(pp_str)*100
        else:
            ins_dom = 0
            perc_ins_dom = 0
            dd_dom = 0
            perc_dd_dom = 0
            f_dom = 0
            perc_f_dom = 0
            
        if len(eff) > 0:
            eff_grp = eff.groupby("alt_iso")[["deletion", "insertion", "frameshift"]].agg("sum").reset_index()
        
            ins_eff = eff_grp.insertion.iloc[0]
            perc_ins_eff = ins_eff/len(pp_str)*100
            dd_eff = eff_grp.deletion.iloc[0]
            perc_dd_eff = dd_eff/len(pp_str)*100
            f_eff = eff_grp.frameshift.iloc[0]
            perc_f_eff = f_eff/len(pp_str)*100
        else:
            ins_eff = 0
            perc_ins_eff = 0
            dd_eff = 0
            perc_dd_eff = 0
            f_eff = 0
            perc_f_eff = 0
        
        
    
    
    ins = pp_str.count("I")
    perc_ins = ins/len(pp_str)*100
    dd = pp_str.count("D")
    perc_dd = dd/len(pp_str)*100
    f = pp_str.count("F")
    f += pp_str.count("f")
    perc_f = f/len(pp_str)*100
    
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
    
    tot_ins_eff.append(ins_eff)
    tot_perc_ins_eff.append(perc_ins_eff)
    tot_dd_eff.append(dd_eff)
    tot_perc_dd_eff.append(perc_dd_eff)
    tot_f_eff.append(f_eff)
    tot_perc_f_eff.append(perc_f_eff)

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

ref_alt_map_nonan["n_ins_eff"] = tot_ins_eff
ref_alt_map_nonan["perc_ins_eff"] = tot_perc_ins_eff
ref_alt_map_nonan["n_dd_eff"] = tot_dd_eff
ref_alt_map_nonan["perc_dd_eff"] = tot_perc_dd_eff
ref_alt_map_nonan["n_f_eff"] = tot_f_eff
ref_alt_map_nonan["perc_f_eff"] = tot_perc_f_eff

ref_alt_map_nonan.sample(5)


# In[94]:


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


# In[95]:


def comp_cat(row):
    if "dom" in row.variable:
        return "pfam"
    elif "eff" in row.variable:
        return "effector"
    else:
        return "total"


# In[96]:


to_plot = pd.melt(ref_alt_map_nonan, id_vars=["ref", "gene", "alt"], value_vars=["n_ins", "perc_ins",
                                                                                 "n_dd", "perc_dd",
                                                                                 "n_f", "perc_f",
                                                                                 "n_ins_dom", "perc_ins_dom",
                                                                                 "n_dd_dom", "perc_dd_dom",
                                                                                 "n_f_dom", "perc_f_dom",
                                                                                 "n_ins_eff", "perc_ins_eff",
                                                                                 "n_dd_eff", "perc_dd_eff",
                                                                                 "n_f_eff", "perc_f_eff"])
to_plot["n_or_perc"] = to_plot["variable"].str.split("_", expand=True)[0]
to_plot["type"] = to_plot["variable"].str.split("_", expand=True)[1]
to_plot["dom_cat"] = to_plot.apply(comp_cat, axis=1)
to_plot.sample(5)


# In[97]:


fig = plt.figure(figsize=(2.3, 1.25))
ax = sns.boxplot(data=to_plot[to_plot["n_or_perc"] == "perc"], 
                 x="type", y="value", hue="dom_cat", order=["dd", "ins", "f"],
                 palette=sns.color_palette("Set2"), fliersize=5, notch=True,
                 flierprops={"marker": "."})
mimic_r_boxplot(ax)


ax.set_xlabel("")
ax.set_xticklabels(["deletions", "insertions", "frameshift"], rotation=30, ha="right", va="top")
ax.set_ylabel("% AA affected")
ax.set_title("alt v. ref TF isoforms")
handles, labels = ax.get_legend_handles_labels()
labels = ["all", "Pfam", "effector"]
ax.legend(handles, labels, loc=2, bbox_to_anchor=(1.01, 1))
fig.savefig('../figures/domain-overall-boxplot.pdf',
            bbox_inches='tight')


# In[98]:


to_plot[to_plot["n_or_perc"] == "perc"].groupby(["type", "dom_cat"]).agg("median")


# In[99]:


ref = "ESR1-201|ESR1-202|ESR1-207|ESR1-208".split("|")[0]
alt = "ESR1-210".split("|")[0]

pp_str = tfs["ESR1"].pairwise_changes_relative_to_reference(ref, alt)
aa_ftr = tfs["ESR1"].aa_feature_disruption(ref)


# In[100]:


pp_str


# In[101]:


aa_ftr[aa_ftr["alt_iso"] == alt]


# In[102]:


len(ref_alt_map_nonan[ref_alt_map_nonan["perc_f_dom"] > 0])


# In[103]:


len(ref_alt_map_nonan[ref_alt_map_nonan["perc_ins"] >= 10])


# In[104]:


len(ref_alt_map_nonan[ref_alt_map_nonan["perc_f"] >= 10])


# In[105]:


len(ref_alt_map_nonan)


# In[106]:


214/2305


# In[110]:


len(ref_alt_map_nonan[(ref_alt_map_nonan["perc_dd_eff"] > 0) |
                      (ref_alt_map_nonan["perc_ins_eff"] > 0) |
                      (ref_alt_map_nonan["perc_f_eff"] > 0)])


# In[111]:


684/2305

