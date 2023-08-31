
# coding: utf-8

# ## TODO
# 
# - Add disordered regions
# - Are there N/C terminal regexes that are not matching correctly?
# - Do I require each isoform has at least one PPI?
# - vizualize motifs

# In[1]:


# TODO: # remove KRTAPS / LCE (only 2% of PPIs)
import re

import numpy as np
import pandas as pd
import seaborn as sns

from data_loading import (load_isoform_and_paralog_y2h_data,
                          load_valid_isoform_clones,
                          load_annotated_6k_collection,
                          load_human_tf_db,
                          load_tf_families,
                          load_pfam_domains_horfeome)


pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 255)


# In[2]:


# looking at Lambert et al. binding mode data, to see if it can help
tfdb = load_human_tf_db()
print(tfdb['Binding mode'].value_counts())
print()
print(tfdb.loc[tfdb['Binding mode'] == 'Obligate heteromer', 'DBD'].value_counts())


# In[3]:


def load_3did_DDIs():
    fpath = '../data/external/3did_flat_2022-05.txt'
    domain_pairs = []
    for line in open(fpath, 'r'):
        if line.startswith('#=ID'):
            domain_pairs.append((line.split()[3][1:8], line.split()[4][:7]))
    df = pd.DataFrame(data=domain_pairs, columns=['pfam_a', 'pfam_b'])
    if df.duplicated().any():
        raise UserWarning('unexpected duplicates')
    return df


def load_cofactors():
    cof = pd.read_csv('../data/external/AnimalTFDB3_Homo_sapiens_TF_cofactors.txt',
                    sep='\t')
    if cof['Symbol'].duplicated().any():
        raise UserWarning('unexpected duplicates')
    return cof


# In[4]:


tfs = load_annotated_6k_collection()
ppi = load_isoform_and_paralog_y2h_data()
ppi = ppi.loc[ppi['category'] == 'tf_isoform_ppis', :]
ppi['gene_level_pair'] = ppi['ad_gene_symbol'] + '_' + ppi['db_gene_symbol']
ppi['ad_iso_id'] = ppi['ad_clone_acc'].apply(lambda x: x.split('|')[0] + '-' + x.split('|')[1].split('/')[0])
# dropping gene with insertion relative to reference geneome
ppi = ppi.loc[~(ppi['ad_gene_symbol'] == 'PCGF6'), :]
iso = load_valid_isoform_clones()
iso = iso.set_index('clone_acc')

# TODO: add in pfam AC for ZF for ZF array
tf_pfam_domains = {tf.name: {dom.accession for iso in tf.cloned_isoforms
                             for dom in iso.aa_seq_features 
                             if dom.category == 'Pfam_domain'}
                   for tf in tfs.values()}
pfam = load_pfam_domains_horfeome()

# kaia added this line as types weren't auto-converting in pandas version whatever
pfam['orf_id'] = pfam['orf_id'].astype(int)

partner_domains = pfam.groupby('orf_id')['pfam_accession'].apply(set)
ppi['partner_domains'] = ppi['db_orf_id'].map(partner_domains)
ppi['tf_domains'] = ppi['ad_gene_symbol'].map(tf_pfam_domains)

# check no-overlap with domain motif
ddi = load_3did_DDIs()

def matching_DDIs(row):
    if (
        pd.isnull(row['partner_domains']) or
        pd.isnull(row['tf_domains'])
        ):
        return np.nan
    if (
        len(row['partner_domains']) == 0 or
        len(row['tf_domains']) == 0
    ):
        return np.nan
    matches = ddi.loc[(ddi['pfam_a'].isin(row['partner_domains']) &
                    ddi['pfam_b'].isin(row['tf_domains'])) |
                   (ddi['pfam_a'].isin(row['tf_domains']) &
                    ddi['pfam_b'].isin(row['partner_domains']))].values
    if len(matches) == 0:
        return np.nan
    return frozenset((a, b) for a, b in matches)

ppi['matching_DDI'] = ppi.apply(matching_DDIs, axis=1)

ggi = ppi.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'db_orf_id', 'matching_DDI']].drop_duplicates()

df = pd.DataFrame(data={tuple(sorted(p)) for ps in ggi.loc[ggi['matching_DDI'].notnull(), 'matching_DDI'].unique() for p in ps},
                  columns=['pfam_a', 'pfam_b'])
df['name_a'] = df['pfam_a'].map(pfam.loc[:, ['pfam_accession', 'domain_name']]
                                    .drop_duplicates()
                                    .set_index('pfam_accession')['domain_name'])
df['name_b'] = df['pfam_b'].map(pfam.loc[:, ['pfam_accession', 'domain_name']]
                                    .drop_duplicates()
                                    .set_index('pfam_accession')['domain_name'])
df['description_a'] = df['pfam_a'].map(pfam.loc[:, ['pfam_accession', 'domain_description']]
                                    .drop_duplicates()
                                    .set_index('pfam_accession')['domain_description'])
df['description_b'] = df['pfam_b'].map(pfam.loc[:, ['pfam_accession', 'domain_description']]
                                    .drop_duplicates()
                                    .set_index('pfam_accession')['domain_description'])
df.sort_values(['name_a', 'name_b']).to_csv('../data/internal/DDI_for_manual_annotation.tsv', index=False, sep='\t')

# filter and consolidate DDIs

ddi_annot = pd.read_csv('../data/internal/DDI_manual_annotation.tsv', sep='\t')

valid_ddi_pairs = {frozenset(p) for p in ddi_annot.loc[ddi_annot['to_use'], ['pfam_a', 'pfam_b']].values}
def filter_ddi(pairs):
    if pd.isnull(pairs):
        return np.nan
    passed = frozenset(p for p in pairs if frozenset(p) in valid_ddi_pairs)
    if len(passed) == 0:
        return np.nan
    return passed

ggi['matching_DDI_filtered'] = ggi['matching_DDI'].apply(filter_ddi)

ddi_to_merge = {
    ('PF00170', 'PF07716'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF00170', 'PF03131'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF07716', 'PF07716'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF03131', 'PF03131'): ('PF00170', 'PF00170'),  # bZIP-bZIP
    ('PF00046', 'PF05920'): ('PF00046',	'PF00046'),  # homeobox-homeobox
    ('PF05920', 'PF05920'): ('PF00046',	'PF00046'),  # homeobox-homeobox
    ('PF00989', 'PF14598'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF00989', 'PF08447'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF00989', 'PF13426'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF14598', 'PF14598'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF08447', 'PF14598'): ('PF00989', 'PF00989'),  # PAS-PAS
    ('PF08447', 'PF08447'): ('PF00989', 'PF00989'),  # PAS-PAS
}
ggi['matching_DDI_filtered'] = ggi['matching_DDI_filtered'].apply(lambda x: frozenset(ddi_to_merge.get(tuple(sorted(p)), tuple(sorted(p))) for p in x) if pd.notnull(x) else x)

pfam_name = (pfam.loc[:, ['pfam_accession', 'domain_name']]
                 .drop_duplicates()
                 .set_index('pfam_accession')
                 ['domain_name']
                 .to_dict())
pfam_name['PF17725'] = 'YBD'  # one missing from that table because it was added more recently
ggi['DDI'] = ggi['matching_DDI_filtered'].apply(lambda x: ' and '.join(pfam_name[a] + '|' + pfam_name[b] for a, b in x) if pd.notnull(x) else x)

tfdb = load_human_tf_db()
ggi['is_tf_tf'] = ggi['db_gene_symbol'].isin(tfdb['HGNC symbol'].values)
fam = load_tf_families()
ggi['is_same_family_tf'] = (ggi['ad_gene_symbol'].map(fam) == ggi['db_gene_symbol'].map(fam))

dlm = pd.read_csv('../data/external/elm_interaction_domains.tsv',
                  sep='\t')
elm = pd.read_csv('../data/external/elm_classes.tsv', header=0, comment='#', sep='\t')
orfs_with_slim_dom = pfam.loc[pfam['pfam_accession'].isin(dlm['Interaction Domain Id']),
                              'orf_id'].unique()
partner_slim_binding_doms = {orf_id: frozenset([d for d in pfam.loc[pfam['orf_id'] == orf_id, 'pfam_accession'].unique() 
                                                if d in dlm['Interaction Domain Id'].unique()]) 
                             for orf_id in orfs_with_slim_dom}
ggi['slim_binding_domains_in_partner'] = ggi['db_orf_id'].map(partner_slim_binding_doms)


# In[5]:


print(len(orfs_with_slim_dom), pfam['orf_id'].nunique())
ggi['putative_DMI'] = ggi['db_orf_id'].isin(orfs_with_slim_dom)
print(ggi['putative_DMI'].sum(), ggi['putative_DMI'].mean(), (ggi['DDI'].isnull() & ggi['putative_DMI']).sum())

candidates = (ggi.loc[ggi['putative_DMI'], :]
                                       .groupby(['db_gene_symbol', 'slim_binding_domains_in_partner'])
                                       .size()
                                       .to_frame())
candidates['names'] = (ggi.loc[ggi['putative_DMI'], :]
                                       .groupby(['db_gene_symbol', 'slim_binding_domains_in_partner'])
                                       ['ad_gene_symbol']
                                       .apply(lambda x: ' '.join(sorted(x))))
                                       
print('top slim binding partners:', candidates.sort_values(0, ascending=False).head(25))

print('')


# In[6]:


ggi.loc[ggi['putative_DMI'] & ggi['DDI'].notnull(), 'DDI'].value_counts()


# In[7]:


dlm.loc[dlm['Interaction Domain Id'].isin(['PF07531', 'PF01753', 'PF08788'])]


# In[8]:


# TFs with slim binding domains
slim_binding_doms = set(dlm['Interaction Domain Id'].unique())
tfs_with_slim_binding_doms = {tf_gene_symbol: frozenset(d for d in doms if d in slim_binding_doms) for tf_gene_symbol, doms in tf_pfam_domains.items()
                              if any(d in slim_binding_doms for d in doms)}
print(len(tfs_with_slim_binding_doms), 'TFs in TFiso1.0 with SLiM binding domain(s)')
ggi['slim_binding_domains_in_tf'] = ggi['ad_gene_symbol'].map(tfs_with_slim_binding_doms)
print(ggi['slim_binding_domains_in_tf'].notnull().sum(), 'gene-level PPI where TF has slim bindig domain')


# In[9]:


ggi['slim_binding_domains_in_tf'].value_counts()


# In[10]:


# STAT3 has SH2 domain -- but conserved in all isoforms
# RBPJ has BTD domain -- but conserved in all isoforms


# In[11]:


ggi.loc[ggi['slim_binding_domains_in_tf'] == frozenset(['PF00104']), 'ad_gene_symbol'].value_counts()


# In[12]:


ggi.loc[ggi['slim_binding_domains_in_tf'] == frozenset(['PF00104']), 'db_gene_symbol'].value_counts()


# In[13]:


ggi.loc[ggi['db_gene_symbol'] == 'NR0B1']


# In[14]:


dlm.loc[dlm['Interaction Domain Id'] == 'PF09270']


# In[15]:


n_ggi = ggi.shape[0]
n_ddi = ggi['matching_DDI'].notnull().sum()
print(n_ggi, 'gene-level PPIs')
print(f'{n_ddi} ({n_ddi/n_ggi:.0%}) mapped to DDIs')
ggi['matching_DDI'].value_counts()


# In[16]:


# check the DDIs that were removed
#ggi.loc[ggi['matching_DDI'].notnull() & ggi['matching_DDI_filtered'].isnull()]


# In[17]:


n_ggi = ggi.shape[0]
n_ddi = ggi['matching_DDI_filtered'].notnull().sum()
print(f'{n_ddi} ({n_ddi/n_ggi:.0%}) mapped to DDIs, after filtering')


# In[18]:


# names
from matplotlib import pyplot as plt
fig, ax = plt.subplots(1, 1)
ggi['DDI'].value_counts().plot.pie(ax=ax)
ax.set_ylabel('')
fig.savefig('../figures/DDI-types_pie.pdf',
            bbox_inches='tight')


# In[19]:


(ggi['is_tf_tf'].mean(),
 ggi['is_same_family_tf'].mean(),
ggi.loc[ggi['DDI'].notnull(), 'is_tf_tf'].mean(),
ggi.loc[ggi['DDI'].notnull(), 'is_same_family_tf'].mean())


# In[20]:


# domain removal
# for each alt iso, for each DDI, calc fraction of domain
# removed and fraction of PPIs retained
dom = pd.concat([g.aa_feature_disruption(g.cloned_reference_isoform.name) for g in tfs.values()])


# In[21]:


df = ppi.loc[ppi['matching_DDI'].notnull(), :].copy()
ref_isos = {tf.cloned_reference_isoform.clone_acc for tf in tfs.values()}
positive_in_ref = df.loc[df['ad_clone_acc'].isin(ref_isos) & 
                         (df['Y2H_result'] == True), 'gene_level_pair'].unique()
df = df.loc[df['gene_level_pair'].isin(positive_in_ref) &
            ~df['ad_clone_acc'].isin(ref_isos), :]
df['matching_DDI'] = df['matching_DDI'].apply(filter_ddi)
df['matching_DDI'] = df['matching_DDI'].apply(lambda x: frozenset(ddi_to_merge.get(tuple(sorted(p)), tuple(sorted(p))) for p in x) if pd.notnull(x) else x)
df = df.loc[df['matching_DDI'].notnull(), :]

merged_domains = {'PF07716': 'PF00170',  # bZIP
 'PF03131': 'PF00170',  # bZIP
 'PF05920': 'PF00046',  # homeobox
'PF14598': 'PF00989',  # PAS
'PF08447': 'PF00989',  # PAS
'PF13426': 'PF00989',  # PAS
}
df['tf_domains'] = df['tf_domains'].apply(lambda x: {merged_domains.get(d, d) for d in x})
dom['accession'] = dom['accession'].apply(lambda x: merged_domains.get(x, x))


def pick_the_one_domain(row):
    ddi_domains = {x for a, b in row['matching_DDI'] for x in [a, b]}
    ds = {d for d in row['tf_domains'] if d in ddi_domains}
    if len(ds) == 0:
        print(row)
        raise UserWarning('something wrong')
    return ds

df['tf_domains'] = df.apply(pick_the_one_domain, axis=1)

def filter_for_domain_in_cloned_reference_isoform(row):
    in_ref = {d.accession for d in tfs[row['ad_gene_symbol']].cloned_reference_isoform.aa_seq_features}
    return {d for d in row['tf_domains'] if d in in_ref}


df['tf_domains'] = df.apply(filter_for_domain_in_cloned_reference_isoform, axis=1)
df = df.loc[df['tf_domains'].map(lambda x: len(x) > 0), :]

def fraction_of_DDI_domains_removed(row):
    ds = dom.loc[(dom['alt_iso'] == row['ad_iso_id']) 
                  & dom['accession'].isin(row['tf_domains']), :]
    if ds.shape[0] == 0:
        print(row)
        raise UserWarning('something wrong')
    return ds[['deletion', 'frameshift']].sum().sum() / ds['length'].sum()

def insertion_in_DDI_domains(row):
    ds = dom.loc[(dom['alt_iso'] == row['ad_iso_id']) 
                  & dom['accession'].isin(row['tf_domains']), :]
    return ds['insertion'].sum()

df['fraction_of_DDI_domains_removed'] = df.apply(fraction_of_DDI_domains_removed, axis=1)
df['insertion_in_DDI_domains'] = df.apply(insertion_in_DDI_domains, axis=1)

df['tf_domains'] = df['tf_domains'].apply(frozenset)


# In[22]:


data = (df.loc[df['Y2H_result'].notnull(), :]
          .groupby(['ad_iso_id', 'tf_domains'])
          ['fraction_of_DDI_domains_removed']
          .mean()
          .to_frame())

# kaia added these lines bc type conversion wasn't working in my version of pandas
nonan = df.loc[df['Y2H_result'].notnull(), :]
nonan['Y2H_result'] = nonan['Y2H_result'].astype(int)

data['Y2H_result_mean'] = nonan.groupby(['ad_iso_id', 'tf_domains'])['Y2H_result'].mean()
data['insertion_in_DDI_domains'] = nonan.groupby(['ad_iso_id', 'tf_domains'])['insertion_in_DDI_domains'].mean()


# In[23]:


data.loc[data['insertion_in_DDI_domains'] > 0]


# In[24]:


data.reset_index()[data.reset_index()['ad_iso_id'].duplicated(keep=False)]


# In[25]:


# add distance from domain
# TODO move to isolib
def n_aa_change_from_feature(gene, ref_iso_name, alt_iso_name, domain_start, domain_end):
    algn = gene.pairwise_changes_relative_to_reference(ref_iso_name, alt_iso_name)

    def _coords_transform_aa_seq_to_alignment(i, alignment):
        if i > len(alignment.replace("I", "")):
            raise ValueError("position is not in isoform AA sequence")
        aa_seq_indices = [
            "" if c == "I" else len(alignment[:j].replace("I", ""))
            for j, c in enumerate(alignment)
        ]
        return aa_seq_indices.index(i)
    
    start = _coords_transform_aa_seq_to_alignment(domain_start, algn)
    end = _coords_transform_aa_seq_to_alignment(domain_end - 1, algn) + 1

    if not all(x == 'M' for x in algn[start:end]):
        return 0  # change is within the domain
    
    big_number = 9999999999999999999999999
    c_dist = big_number
    n_dist = big_number
    for i, l in enumerate(reversed(algn[:start])):
        if l != 'M':
            c_dist = i + 1
            break
    for i, l in enumerate(algn[end:]):
        if l != 'M':
            n_dist = i + 1
            break
    if c_dist == big_number and n_dist == big_number:
        raise UserWarning('problem calculating distance')
    return min([c_dist, n_dist])


def n_aa_to_all_features(self, ref_iso_name):
    results = []
    ref_iso = self._orf_dict[ref_iso_name]
    row = {"gene": self.name, "ref_iso": ref_iso_name}
    for aa_feature in ref_iso.aa_seq_features:
        for alt_iso_name, alt_iso in self._orf_dict.items():
            if alt_iso_name == ref_iso_name:
                continue
            row.update(
                {
                    "alt_iso": alt_iso_name,
                    "accession": aa_feature.accession,
                    "category": aa_feature.category,
                    "start_in_ref_iso": aa_feature.start,
                    "end_in_ref_iso": aa_feature.end,
                    "length": aa_feature.end - aa_feature.start,
                }
            )
            row.update({"n_aa_change_to_domain": n_aa_change_from_feature(self, ref_iso_name, alt_iso_name, aa_feature.start, aa_feature.end)})
            results.append(row.copy())
    results = pd.DataFrame(results)
    return results

dist = pd.concat([n_aa_to_all_features(g, g.cloned_reference_isoform.name) for g in tfs.values()])

def get_dist(row):
    alt_iso, dom_accessions = row.name
    return dist.loc[(dist['alt_iso'] == alt_iso) & (dist['accession'].isin(dom_accessions)),
                    'n_aa_change_to_domain'].min()

data['domain_n_aa_to_change'] = data.apply(get_dist, axis=1)


# In[26]:


from matplotlib import pyplot as plt
COLOR_PURPLE = (155 / 255, 97 / 255, 153 / 255)


gs_kw = dict(width_ratios=[1, 1, 2])

fig, axs = plt.subplots(1, 3, sharey=True, gridspec_kw=gs_kw)
fig.set_size_inches(w=6, h=2.5)

point_size = 7

axs[0].set_title('Full loss of domain',
fontsize=10)
sns.swarmplot(data=data.loc[data['fraction_of_DDI_domains_removed'] == 1, :],
              y='Y2H_result_mean', 
              x='fraction_of_DDI_domains_removed',
              size=point_size,
         #     order=[
         #            'Full loss\nof DBD',
         #            ],
            clip_on=False,
              ax=axs[0],
              color=COLOR_PURPLE,
              alpha=1)

axs[1].set_title('Partial loss of domain',
fontsize=10)
partial_loss = (data['fraction_of_DDI_domains_removed'] > 0) & (data['fraction_of_DDI_domains_removed'] < 1)
axs[1].scatter(data.loc[partial_loss, 'fraction_of_DDI_domains_removed'].values,
               data.loc[partial_loss, 'Y2H_result_mean'].values,
           alpha=1,
           s=point_size**2/1.5,  # I don't think there should be a divide by anything here....
            color=COLOR_PURPLE,
           clip_on=False)
axs[1].set_xlabel('Proportion missing')
axs[1].set_xlim(1, 0)
axs[1].set_xticks([0.99, 0.5, 0.01])
axs[1].set_xticklabels([f'{x:.0%}' for x in axs[1].get_xticks()])
#axs[1].set_xticks(range(10, 91, 10), minor=True)



axs[2].set_title('Full domain in\nalternative isoform', fontsize=10)
axs[2].scatter(data.loc[(data['fraction_of_DDI_domains_removed'] == 0), 'domain_n_aa_to_change'].values,
               data.loc[(data['fraction_of_DDI_domains_removed'] == 0), 'Y2H_result_mean'].values,
           alpha=1,
           s=point_size**2/1.5,  # I don't think there should be a divide by anything here....
            color=COLOR_PURPLE,
           clip_on=False)
axs[2].set_xlabel('Distance of alternative\nsequence from domain\n(number of AA)')

for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim(0, 1)
for ax in axs[1:]:
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_tick_params(which='both', length=0)
for i in [0]:
    axs[i].set_xlabel('')
    axs[i].set_ylabel('')
    axs[i].spines['bottom'].set_visible(False)
    axs[i].xaxis.set_tick_params(length=0)
    axs[i].set_xticks([])
axs[0].set_yticks([0, 0.25, 0.5, 0.75, 1])
axs[0].set_yticks(np.linspace(0, 1, 21), minor=True)
axs[0].set_yticklabels(['{:.0%}'.format(y) for y in axs[0].get_yticks()])
axs[0].set_ylabel('Fraction of domain-domain mediated\nPPIs with alternative isoform')
fig.savefig('../figures/PPI_vs_domain_removal.pdf', bbox_inches='tight')


# In[27]:


data.loc[(data['fraction_of_DDI_domains_removed'] > 0) & (data['Y2H_result_mean'] > 0)]


# ### Cases of partial domain removal without full loss of PPIs
# 
# - E2F3-3, coiled-coil marked box domain binds DP TFs, binds DP1, DP2, not DP3. All novel isoforms. 
# - FOS-2. Retains all the bzip partners despite missing some of the domain. Maybe the basic region and might lose DNA binding (no data). 
# - LHX9-2 loses one amino acid of LIM domain but doesn't lose any binding partners
# - SMAD4-1. Loses some of MH2 domain. Doens't lose partners. Novel isoform.

# In[28]:


# kaia added these lines because type conversion wasn't working
nonan = df.loc[df['Y2H_result'].notnull() & (df['insertion_in_DDI_domains'] > 0), :]
nonan['Y2H_result'] = nonan['Y2H_result'].astype(int)

x = nonan.groupby(['ad_iso_id', 'tf_domains'])['insertion_in_DDI_domains'].mean().values
y = nonan.groupby(['ad_iso_id', 'tf_domains'])['Y2H_result'].mean().values

fig, ax = plt.subplots(1, 1)
ax.scatter(x, y, alpha=0.5, color='purple')
ax.set_ylabel('Fraction of alternative\nisoforms retaining PPI')
ax.set_xlabel('# aa insertion within DDI domain')
fig.savefig('../figures/PPI_vs_domain_insertion.pdf', bbox_inches='tight')


# In[29]:


# break down DDI vs DMI etc etc.


# In[30]:


(ggi.loc[ggi['is_tf_tf'], 
        'is_same_family_tf']
.value_counts()
.plot.pie())


# In[31]:


cof = load_cofactors()
ggi['is_tf_cf'] = ggi['db_gene_symbol'].isin(cof['Symbol'].values)


# In[32]:


(~ggi['is_tf_tf'] & ggi['is_tf_cf']).mean()


# In[33]:


from data_loading import load_ppi_partner_categories

cats = load_ppi_partner_categories()


# In[34]:


# this is by partner, plot by PPI instead
cats['category'].value_counts().plot.pie()


# In[35]:


cats.loc[cats['category'] == 'nuclear transport']


# In[36]:


ggi.loc[ggi['db_gene_symbol'].isin(
        cats.loc[cats['category'] == 'nuclear transport',
                 'partner'].values)]


# In[37]:


cats.head()


# In[38]:


fig, ax = plt.subplots(1, 1)
(ggi['db_gene_symbol']
.map(cats.drop_duplicates('partner').set_index('partner')['category'])
.value_counts()
.plot.pie(ax=ax))
ax.set_ylabel('')
fig.savefig('../figures/PPIs-gene-level-manual-categories_pie.pdf',
            bbox_inches='tight')


# In[39]:


ggi['db_gene_symbol'].map(cats.drop_duplicates('partner').set_index('partner')['category']).value_counts()


# In[40]:


ggi['db_gene_symbol'].map(cats.drop_duplicates('partner').set_index('partner')['category']).value_counts().sum()


# In[44]:


print(len(ggi))
ggi = ggi.merge(cats, left_on='db_gene_symbol', right_on='partner').drop_duplicates()
print(len(ggi))
ggi.cofactor_type.value_counts()


# In[45]:


def categorize_PPI_partner(row):
    if row['is_tf_tf']:
        return 'TF'
    elif row['is_tf_cf']:
        return 'cofactor'
    elif row['category'] == "signaling":
        return 'signaling'
    else:
        return 'other'

ggi['partner_category'] = ggi.apply(categorize_PPI_partner, axis=1)
ggi['partner_category'].value_counts()[['TF', 'cofactor', 'signaling', 'other']].plot.pie()


# In[46]:


len(ggi)


# In[47]:


len(ggi.db_gene_symbol.unique())


# In[48]:


ggi['partner_category'].value_counts()


# In[50]:


ys = np.array([len(ggi[ggi["partner_category"] == "TF"]), len(ggi[ggi["partner_category"] == "cofactor"]),
               len(ggi[ggi["partner_category"] == "signaling"]),
              len(ggi[ggi["partner_category"] == "other"])])
labels = ["TF", "cofactor", "signaling", "other"]
colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], sns.color_palette("Set2")[5], "darkgray"]

fig, ax = plt.subplots(figsize=(2.0, 2.0), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.0f%%', startangle=-120, 
                    explode=(0.02, 0.2, 0.05, 0.05))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")

fig.savefig("../figures/PPIs-gene-level-manual-categories_simplified.pdf", dpi="figure", bbox_inches="tight")


# In[51]:


ggi[ggi["partner_category"] != "cofactor"].cofactor_type.value_counts()


# In[52]:


cofacs = ggi[ggi["partner_category"] == "cofactor"]

ys = np.array([len(cofacs[cofacs["cofactor_type"] == "corepressor"]), 
               len(cofacs[cofacs["cofactor_type"] == "coactivator"]),
               len(cofacs[cofacs["cofactor_type"] == "both"]),
               len(cofacs[cofacs["cofactor_type"] == "unknown"])])
labels = ["coactivator", "corepressor", "both", "unknown"]
colors = [sns.color_palette("Set2")[3], sns.color_palette("Set2")[0], sns.color_palette("Set2")[4], "darkgray"]

fig, ax = plt.subplots(figsize=(2.25, 2.25), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.0f%%', startangle=-90, 
                    explode=(0.05, 0.05, 0.05, 0.75))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")

fig.savefig("../figures/PPIs-gene-level-manual-categories_cofactors.pdf", dpi="figure", bbox_inches="tight")


# In[94]:


# RNA binding proteins
rbp = pd.read_csv('../data/external/RBPDB_v1.3.1_proteins_human_2012-11-21.tdt',
                 sep='\t',
                 header=None)
# BUG: there are 10 genes missing names
rbp = set(rbp.loc[rbp.loc[:, 4].notnull(), 4].unique())
ggi['is_tf_rbp'] = ggi['db_gene_symbol'].isin(rbp)
(~ggi['is_tf_tf'] & ~ggi['is_tf_cf'] & ggi['is_tf_rbp']).mean()


# In[95]:


# kinase / phosphatase
from Bio import SeqIO

kinase_genes = {r.description.split()[2][len('gene='):] for r in SeqIO.parse('../data/external/human_kinases.fa', 'fasta')}
ggi['is_tf_kinase'] = ggi['db_gene_symbol'].isin(kinase_genes)
(~ggi['is_tf_tf'] & ~ggi['is_tf_cf'] & ggi['is_tf_kinase']).mean()


# In[96]:


# relation of domain loss to DDI PPI loss


# In[97]:


pair = ('PF00023', 'PF16179')
pair = ('PF12796', 'PF16179')
pair = ('PF00010',	'PF00989')
pair = ('PF00505',	'PF00505')

pair = ('PF00989',	'PF00989')
pair = ('PF00989',	'PF14598')
#pair = ('PF00989',	'PF08447')
#pair = ('PF08447',	'PF08447')

# NPAS2	ARNTL
# NPAS2	ARNTL2
# ARNT2	SIM1
# HIF1A	ARNTL

pair = ('PF00170',	'PF00170')

ggi.loc[ggi['matching_DDI'].apply(lambda x: not pd.isnull(x) and pair in x)]


# In[98]:


pair = ('PF03166',	'PF03166')
ggi.loc[ggi['matching_DDI'].apply(lambda x: not pd.isnull(x) and pair in x)]


# In[99]:


pair = ('PF03166',	'PF08782')
ggi.loc[ggi['matching_DDI'].apply(lambda x: not pd.isnull(x) and pair in x)]


# In[100]:


pair = ('PF07714',	'PF00017')
ggi.loc[ggi['matching_DDI'].apply(lambda x: not pd.isnull(x) and pair in x)]


# In[101]:


ggi.loc[ggi['matching_DDI'].apply(lambda x: not pd.isnull(x) and ('PF00023', 'PF16179') in x)]


# In[102]:


# proportion of gene-level PPIs where the partner has a slim binding domain
# or the TF has a slim binding domain
ggi.head()


# In[103]:


print(ppi.loc[:, ['ad_gene_symbol', 'db_gene_symbol']].drop_duplicates().shape[0],
      'gene level PPIs')
print(ppi.loc[ppi['db_orf_id'].isin(orfs_with_slim_dom),
              ['ad_gene_symbol', 'db_gene_symbol']].drop_duplicates().shape[0],
      'gene level PPIs involving motif binding domain')


# In[104]:


dom_slim_ppi = set(ppi.loc[ppi['db_orf_id'].isin(orfs_with_slim_dom),
                       'gene_level_pair'].unique())
pos_ppi = set(ppi.loc[(ppi['Y2H_result'] == True), 'gene_level_pair'].unique())
neg_ppi = set(ppi.loc[(ppi['Y2H_result'] == False), 'gene_level_pair'].unique())
diff_slim_ppi = dom_slim_ppi.intersection(pos_ppi.intersection(neg_ppi))
print(len(diff_slim_ppi), 
      'gene-level PPIs where partner has SLiM binding domain and some isoforms bind and some dont')
print('involving {} TF genes'.format(len({p.split('_')[0] for p in diff_slim_ppi})))


# In[105]:


def match_elm_motifs(orfs, elm):
    """Take the ELM motifs file and match the regex to given aa seqs
    Args:
        orfs (Series): amino acid sequences indexed by ORF ID
    Returns:
        DataFrame: one row for each seperate match of a motif to an ORF
    """
    motifs = []
    for elmIdx, elmRow in elm.iterrows():
        # inserting ?: to make groups non-capture
        motifMatches = orfs.str.findall(elmRow['Regex'].replace('(', '(?:'))
        for orfID, matches in zip(motifMatches.index, motifMatches.values):
            for match in matches:
                start = orfs.loc[orfID].find(match)
                # switch from 0 to 1-based indexing
                motifs.append((orfID,
                               elmRow['Accession'],
                               elmRow['ELMIdentifier'],
                               start + 1,
                               start + len(match)))
    motifs = pd.DataFrame(data=motifs,
                          columns=('orf_id', 'region_id', 'ELM_ID', 'start', 'end'))
    motifs['source'] = 'ELM_prediction'
    motifs['type'] = 'SLiM'
    return motifs


slims = match_elm_motifs(iso['aa_seq'], elm)


# In[106]:


slims


# In[107]:


def isoform_specific_regions(gene, subset=None):
    """The name is a bit misleading because it's not specific to one isoform but just
       not common to all isoforms.

    Returns: dict(frozenset: list(str)): isoform IDs and list of contiguous AA sequences
                                         that map to them only 

    """
    algn = gene.genomic_alignment_of_aa_seqs(subset=subset)
    subset_prev = None
    isr = {}
    len_algn = len(list(algn.values())[0])
    for i in range(len_algn):
        subset = frozenset({k for k, v in algn.items() if v[i] != '-'})  # the isoforms that have an aa at that genomic position
        if subset_prev is None:
            if (len(subset) < len(algn)) and (len(subset) > 0):
                start = i
                subset_prev = subset
        else:
            if subset != subset_prev:
                if (len(subset_prev) < len(algn)) and (len(subset_prev) > 0):
                    subseq = (algn[list(subset_prev)[0]][start:i], start, i)
                    isr[subset_prev] = isr.get(subset_prev, []) + [subseq]
                start = i
                subset_prev = subset
            elif  i == (len_algn - 1):
                if (len(subset_prev) < len(algn)) and (len(subset_prev) > 0):
                    subseq = (algn[list(subset_prev)[0]][start:], start, i + 1)
                    isr[subset_prev] = isr.get(subset_prev, []) + [subseq]
                start = i
                subset_prev = subset
    merged = {}
    for iso_subset, subseqs in isr.items():
        merged[iso_subset] = []
        prev_end = np.inf
        prev_subseq = ''
        for subseq, start, end in subseqs:
            if start <= prev_end + 2:
                prev_subseq += subseq
                prev_end = end
            else:
                if prev_subseq != '':
                    merged[iso_subset].append(prev_subseq)
                prev_subseq = subseq
                prev_end = end
        merged[iso_subset].append(prev_subseq)
    return merged


def ppi_tf_gene(data, gene_name):
    tf = data.loc[(data['category'] == 'tf_isoform_ppis') &
                (data['ad_gene_symbol'] == gene_name),
                ['ad_clone_name', 'db_gene_symbol', 'Y2H_result']].copy()
    tf = tf.pivot(index='ad_clone_name',
                  columns='db_gene_symbol',
                  values='Y2H_result')
    return tf


def ppi_linked_isoform_specific_regions(ppi_data, gene):
    """
    For now, do not use cases where there are missing values
    """
    isr = isoform_specific_regions(gene, subset={iso.name for iso in gene.cloned_isoforms})
    ppi = ppi_tf_gene(ppi_data, gene.name)
    ppi_isr = {}
    ppi_iso = {partner: set(ppi.index[ppi[partner]])
               for partner in ppi.columns
               if ppi[partner].notnull().all()}
    for partner, ppi_iso_subset in ppi_iso.items():
        for isr_subset, aa_seqs in isr.items():
            if ppi_iso_subset == isr_subset:
                ppi_isr[partner] = (isr_subset, aa_seqs)
    return ppi_isr


# In[108]:


slim_binding_domains = pd.merge(pfam,
                                dlm,
                                how='inner',
                                left_on='pfam_accession',
                                right_on='Interaction Domain Id')
slim_binding_domains = pd.merge(slim_binding_domains,
                                ppi,
                                how='inner',
                                left_on='orf_id',
                                right_on='db_orf_id')
slim_ppis = slim_binding_domains.loc[:, ['ad_gene_symbol', 
                                        'db_gene_symbol',
                                        'pfam_accession',
                                        'domain_name',
                                        'domain_description',
                                        'ELM identifier']].drop_duplicates()
ppi_isr = {}
for gene_name in slim_ppis['ad_gene_symbol'].unique():
    if gene_name not in tfs:
        print(gene_name, 'missing')
        continue
    ppi_isr[gene_name] = ppi_linked_isoform_specific_regions(ppi, tfs[gene_name])
slim_ppis = slim_ppis.loc[slim_ppis.apply(lambda x: x['db_gene_symbol'] in ppi_isr[x['ad_gene_symbol']],
                              axis=1),
              :]
slim_ppis['aa_seq_isr'] = slim_ppis.apply(lambda x: ppi_isr[x['ad_gene_symbol']][x['db_gene_symbol']][1], axis=1)


# In[109]:


def isr_contains_slim(row):
    # inserting ?: to make groups non-capture
    if row['ELM identifier'] not in elm['ELMIdentifier'].values:
        #raise UserWarning('Missing ELM entry for: ', row['ELM identifier'])
        print('Missing ELM entry for: ', row['ELM identifier'])
        return False
    regex = elm.loc[elm['ELMIdentifier'] == row['ELM identifier'], 'Regex'].values[0].replace('(', '(?:')
    return any(bool(re.search(regex, aa_seq)) for aa_seq in row['aa_seq_isr'])


slim_ppis['slim_match'] = slim_ppis.apply(isr_contains_slim, axis=1)


# In[110]:


slim_ppis['slim_match'].sum()


# In[111]:


slim_ppis.loc[slim_ppis['slim_match'] & (slim_ppis['db_gene_symbol'] == 'TXK')].sort_values('ad_gene_symbol')


# In[112]:


# use disordered info
slim_ppis


# In[113]:


slim_ppis.loc[slim_ppis['slim_match'], 'ELM identifier'].value_counts()


# In[114]:


slim_ppis.loc[slim_ppis['slim_match'] & (slim_ppis['ELM identifier'] == 'LIG_SH3_3')]


# In[115]:


slim_ppis.loc[slim_ppis['ad_gene_symbol'] == 'PBX1']


# In[116]:


# count inside region, count outside region...
slim_ppis['aa_seq_ref_iso'] = slim_ppis['ad_gene_symbol'].apply(lambda x: tfs[x].orfs[0].aa_seq)
slim_ppis['regex'] = slim_ppis['ELM identifier'].map(elm[['ELMIdentifier', 'Regex']]
                                                        .drop_duplicates()
                                                        .set_index('ELMIdentifier')
                                                        ['Regex'].str.replace('(', '(?:'))
# TODO: fix missing regex.
#       Difficult because they don't appear on the website
#       or in the list of (e.g. MOD_PLK). I should email
#       the ELM people.
slim_ppis = slim_ppis.dropna(subset=['regex'])
slim_ppis['match_count_ref_iso'] = slim_ppis.apply(lambda x: len(re.findall(x['regex'], x['aa_seq_ref_iso'])),
                axis=1)
slim_ppis['match_count_isr'] = slim_ppis.apply(lambda x: sum(len(re.findall(x['regex'], s)) for s in x['aa_seq_isr']),
                axis=1)
slim_ppis['aa_seq_len_ref_iso'] = slim_ppis['aa_seq_ref_iso'].str.len()
slim_ppis['isr_len'] = slim_ppis['aa_seq_isr'].apply(lambda x: sum(len(s) for s in x))
slim_ppis['motif_probability'] = slim_ppis['ELM identifier'].map(elm.set_index('ELMIdentifier')['Probability'])


# In[117]:


slim_ppis.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_isr'].sum().head(30)


# In[118]:


# per ISR / domain, how many get a match
(slim_ppis.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_isr'].sum() > 0).value_counts()


# In[119]:


slim_ppis.head()


# In[120]:


# for each gene-gene pair, pick random sequence on longest isoform of length of ISR and look for match
gg = slim_ppis.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'aa_seq_ref_iso', 'isr_len']].drop_duplicates().copy().set_index(['ad_gene_symbol', 'db_gene_symbol'])

def rand_subseq(row):
    rand_start = np.random.randint((len(row['aa_seq_ref_iso']) - row['isr_len']) + 1)
    return row['aa_seq_ref_iso'][rand_start:rand_start + row['isr_len']]

null_dist = []
n_samples = 100
for __ in range(n_samples):
    rand_seq = gg.apply(rand_subseq, axis=1)
    rand_seq = pd.merge(rand_seq.rename('rand_aa_seq').reset_index(),
                        slim_ppis.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'regex']],
                        how='inner')
    rand_seq['match_count_rand_aa_seq'] = rand_seq.apply(lambda x: len(re.findall(x['regex'], x['rand_aa_seq'])),
                    axis=1)
    null_dist.append((rand_seq.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_rand_aa_seq'].sum() > 0).sum())
real = (slim_ppis.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_isr'].sum() > 0).sum()
n = slim_ppis[['ad_gene_symbol', 'db_gene_symbol']].drop_duplicates().shape[0]
print('real =', real, 'out of', n)
print('P =', len([n for n in null_dist if n >= real]) / n_samples)


# In[121]:


# trying same thing using most specific regex
slim_ppis_spec = slim_ppis.sort_values('motif_probability', ascending=True).drop_duplicates(['ad_gene_symbol', 'db_gene_symbol', 'pfam_accession']).copy()


# In[122]:


gg = slim_ppis_spec.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'aa_seq_ref_iso', 'isr_len']].drop_duplicates().copy().set_index(['ad_gene_symbol', 'db_gene_symbol'])

def rand_subseq(row):
    rand_start = np.random.randint((len(row['aa_seq_ref_iso']) - row['isr_len']) + 1)
    return row['aa_seq_ref_iso'][rand_start:rand_start + row['isr_len']]

null_dist = []
n_samples = 100
for __ in range(n_samples):
    rand_seq = gg.apply(rand_subseq, axis=1)
    rand_seq = pd.merge(rand_seq.rename('rand_aa_seq').reset_index(),
                        slim_ppis_spec.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'regex']],
                        how='inner')
    rand_seq['match_count_rand_aa_seq'] = rand_seq.apply(lambda x: len(re.findall(x['regex'], x['rand_aa_seq'])),
                    axis=1)
    null_dist.append((rand_seq.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_rand_aa_seq'].sum() > 0).sum())
real = (slim_ppis_spec.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_isr'].sum() > 0).sum()
n = slim_ppis_spec[['ad_gene_symbol', 'db_gene_symbol']].drop_duplicates().shape[0]
print('real = ', real, 'out of', n)
print('P =', len([n for n in null_dist if n >= real]) / n_samples)


# In[123]:


slim_ppis['isr_frac'] = (slim_ppis['isr_len'] / slim_ppis['aa_seq_len_ref_iso'])


# In[124]:


(slim_ppis[['ad_gene_symbol', 'db_gene_symbol', 'isr_frac']].drop_duplicates()['isr_frac'] < 0.2).sum()


# In[125]:


from matplotlib import pyplot as plt
(slim_ppis.loc[: ,['ad_gene_symbol', 'db_gene_symbol', 'isr_frac']].drop_duplicates()['isr_frac'] * 100).plot.hist()
plt.xlabel('% of reference isoform in PPI specific region')
plt.ylabel('Number of gene-level PPIs')
plt.savefig('../figures/fraction_isr.pdf', bbox_inches='tight')


# In[126]:


slim_ppis_spec = slim_ppis.loc[slim_ppis['isr_frac'] <= 0.1, :].copy()
gg = slim_ppis_spec.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'aa_seq_ref_iso', 'isr_len']].drop_duplicates().copy().set_index(['ad_gene_symbol', 'db_gene_symbol'])

def rand_subseq(row):
    rand_start = np.random.randint((len(row['aa_seq_ref_iso']) - row['isr_len']) + 1)
    return row['aa_seq_ref_iso'][rand_start:rand_start + row['isr_len']]

null_dist = []
n_samples = 100
for __ in range(n_samples):
    rand_seq = gg.apply(rand_subseq, axis=1)
    rand_seq = pd.merge(rand_seq.rename('rand_aa_seq').reset_index(),
                        slim_ppis_spec.loc[:, ['ad_gene_symbol', 'db_gene_symbol', 'regex']],
                        how='inner')
    rand_seq['match_count_rand_aa_seq'] = rand_seq.apply(lambda x: len(re.findall(x['regex'], x['rand_aa_seq'])),
                    axis=1)
    null_dist.append((rand_seq.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_rand_aa_seq'].sum() > 0).sum())
real = (slim_ppis_spec.groupby(['ad_gene_symbol', 'db_gene_symbol'])['match_count_isr'].sum() > 0).sum()
n = slim_ppis_spec[['ad_gene_symbol', 'db_gene_symbol']].drop_duplicates().shape[0]
print('real = ', real, 'out of', n)
print('P =', len([n for n in null_dist if n >= real]) / n_samples)


# In[127]:


print(slim_ppis['match_count_ref_iso'].sum())
print(slim_ppis['match_count_isr'].sum())
print(slim_ppis['match_count_ref_iso'].sum() / slim_ppis['aa_seq_len_ref_iso'].sum())
print(slim_ppis['match_count_isr'].sum() / slim_ppis['isr_len'].sum())


# In[128]:


res = slim_ppis.sort_values('motif_probability', ascending=True).drop_duplicates(['ad_gene_symbol', 'db_gene_symbol', 'pfam_accession'])
print(res['match_count_ref_iso'].sum() / res['aa_seq_len_ref_iso'].sum())
print(res['match_count_isr'].sum() / res['isr_len'].sum())
print((res['match_count_ref_iso'] - res['match_count_isr']).sum() / (res['aa_seq_len_ref_iso'] - res['isr_len']).sum())
print((slim_ppis.loc[slim_ppis['slim_match'],
               ['ad_gene_symbol', 'db_gene_symbol']]
          .drop_duplicates()
          .shape[0]))


# In[129]:


slim_ppis.columns


# In[130]:


(slim_ppis.loc[slim_ppis['slim_match'], ['ad_gene_symbol', 'db_gene_symbol', 'domain_description', 'ELM identifier', 'regex', 'motif_probability']].sort_values('motif_probability')
.drop_duplicates(['ad_gene_symbol', 'db_gene_symbol', 'domain_description'])).to_csv('../output/slim_matches.csv', index=False)


# In[131]:


match_pairs = (slim_ppis.loc[slim_ppis['slim_match'],
               ['ad_gene_symbol', 'db_gene_symbol']]
                .drop_duplicates())
match_pairs.head()


# In[132]:


print(match_pairs.shape[0], 'gene-level PPIs corresponding to motif in the isoform regions')
print('out of {} possible'.format(slim_ppis[['ad_gene_symbol', 'db_gene_symbol']].drop_duplicates().shape[0]))


# In[133]:


match_pairs['ad_gene_symbol'].value_counts().head()


# In[134]:


slim_ppis.loc[slim_ppis['ad_gene_symbol'] == 'FOXP2', :]


# In[135]:


match_pairs['db_gene_symbol'].value_counts().head()


# In[136]:


slim_ppis.loc[slim_ppis['slim_match'] & (slim_ppis['db_gene_symbol'] == 'PIN1'), :]


# Pin1 Links the Activities of c-Abl
# and p300 in Regulating p73 Function, Mol. Cell, 2004
# 
# It has been showing that TP73 isoforms lacking the PIN1 interacting motif have reduced transcriptional activity.
# 
# Also the mouse FOS has the PIN1 motif binding shown.

# In[137]:


ppi.loc[ppi['db_gene_symbol'] == 'PIN1',
        'ad_gene_symbol'].nunique()
ppi.loc[ppi['db_gene_symbol'] == 'PIN1',
        'ad_gene_symbol'].unique()


# In[138]:


# smallest region
slim_ppis['isr_length'] = slim_ppis['aa_seq_isr'].apply(lambda x: min(len(s) for s in x))
(slim_ppis.loc[slim_ppis['slim_match'], :]
          .sort_values('isr_length')
          .drop_duplicates(['ad_gene_symbol', 'db_gene_symbol', 'pfam_accession'])
          .head(30))


# In[139]:


# most specific regex
(slim_ppis.loc[slim_ppis['slim_match'], :]
          .sort_values('motif_probability')
          .drop_duplicates(['ad_gene_symbol', 'db_gene_symbol', 'pfam_accession'])
          .head(30))


# In[140]:


slim_ppis.loc[slim_ppis['slim_match'] & (slim_ppis['ad_gene_symbol'] == 'STAT3'), :]

