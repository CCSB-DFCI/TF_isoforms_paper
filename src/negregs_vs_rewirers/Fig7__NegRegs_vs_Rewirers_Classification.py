#!/usr/bin/env python
# coding: utf-8

# # Figure 7: Classifying TFs into Neg. Regs. vs. Rewirers
# 
# - classification based on assay results

# In[1]:


import warnings
warnings.filterwarnings('ignore')


# In[2]:


import matplotlib as mpl
import matplotlib.pyplot as plt
import met_brewer
import pandas as pd
import numpy as np
import seaborn as sns
import sys
import upsetplot

from Bio.Seq import Seq
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
import statsmodels.stats.multitest as smt
from upsetplot import plot

# import utils
sys.path.append("../")
sys.path.append("../data_loading")

import plotting
from plotting import PAPER_PRESET, PAPER_FONTSIZE, nice_boxplot, mimic_r_boxplot, annotate_pval

from data_loading import (load_annotated_TFiso1_collection,
                          load_developmental_tissue_expression_remapped,
                          load_gtex_remapped,
                          load_ref_vs_alt_isoforms_table,
                          load_condensate_data)

get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'svg'")
mpl.rcParams['figure.autolayout'] = False


# In[3]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# In[4]:


np.random.seed(2023)


# ## functions

# In[5]:


def calculate_tau(df):
    array = df.values
    
    ## will return NaN as tau for every row that has any NaNs
    array_max = np.max(array, axis=1)
    tmp = array.T / array_max
    tmp = 1 - tmp.T
    nonan_taus = np.sum(tmp, axis=1) / (array.shape[1])
    
    ## will ignore NaNs and compute on the rest of the values
    array_max = np.nanmax(array, axis=1)
    tmp = array.T / array_max
    tmp = 1 - tmp.T
    nan_taus = np.nansum(tmp, axis=1) / np.count_nonzero(~np.isnan(array), axis=1)
    
    
    return nonan_taus, nan_taus, array_max


# In[6]:


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


# In[7]:


dn_pal = {"ref": sns.color_palette("Set2")[0],
       "similar": sns.color_palette("Set2")[0],
       "rewire": sns.color_palette("Set2")[2],
       "DN": sns.color_palette("Set2")[1],
       "NA": "lightgray",
       "likely": "darkgray",
          "combination": sns.color_palette("Set2")[5],
          "assessed": sns.color_palette("Set2")[5]}


# ## variables

# In[8]:


joung_orf_f = "../../data/external/joung_files/Joung_ORF_lib.txt"
joung_data_f = "../../data/external/joung_files/Joung_ORF_scores.txt"
joung_cells_f = "../../data/external/joung_files/Figure3B_celltype_mapping.csv"

joung_down_map_batch_f = "../../data/external/joung_files/subsample_mapping_batch.txt"
joung_down_map_TF_f = "../../data/external/joung_files/subsample_mapping_TF.txt"
joung_down_map_louvain_f = "../../data/external/joung_files/subsample_mapping_louvain.txt"


# ## 1. import data

# In[9]:


pairs = load_ref_vs_alt_isoforms_table()

# RORC-1 alt iso is causing an error - filter out here - there's no data for it?
pairs = pairs[pairs["clone_acc_alt"] != "RORC|1/6|05F11"]

pairs['ref_iso'] = pairs['clone_acc_ref'].apply(lambda x: x.split('|')[0] + '-' + x.split('|')[1].split('/')[0])
pairs['alt_iso'] = pairs['clone_acc_alt'].apply(lambda x: x.split('|')[0] + '-' + x.split('|')[1].split('/')[0])


# In[10]:


len(pairs[pd.isnull(pairs["dbd_pct_lost"])])


# In[11]:


pairs[pairs['alt_iso'] == 'TBX5-2'].iloc[0]


# In[12]:


joung_orf = pd.read_table(joung_orf_f)
joung_orf["Name"] = joung_orf["Name"].str.strip()

joung_data = pd.read_table(joung_data_f)
joung_data["Name"] = joung_data["TF ORF"].str.split("-", expand=True)[0].str.strip()

joung_cells = pd.read_table(joung_cells_f, sep=",")

joung_down_map_batch = pd.read_table(joung_down_map_batch_f, index_col=0)
print(len(joung_down_map_batch))
joung_down_map_TF = pd.read_table(joung_down_map_TF_f, index_col=0)
print(len(joung_down_map_TF))
joung_down_map_louvain = pd.read_table(joung_down_map_louvain_f, index_col=0)
print(len(joung_down_map_louvain))

joung_down_map = joung_down_map_batch.join(joung_down_map_TF).join(joung_down_map_louvain)
print(len(joung_down_map))


# In[13]:


df_gtex, metadata_gtex, genes_gtex = load_gtex_remapped()

exclusion_list_gtex = {'Cells - Leukemia cell line (CML)',
                       'Cells - EBV-transformed lymphocytes',
                       'Cells - Cultured fibroblasts'}

df_gtex = df_gtex.loc[:, ~df_gtex.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex = metadata_gtex.loc[~metadata_gtex['body_site'].isin(exclusion_list_gtex), :]

means_gtex = df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean()


# In[14]:


metadata_gtex_dummy = pd.read_table("../../data/processed/metadata_gtex_dummy.csv", sep=",", index_col=0)


# In[15]:


df_dev, metadata_dev, genes_dev = load_developmental_tissue_expression_remapped()

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


# ## 2. categorize based on assay data

# ### M1H

# In[16]:


def m1h_cat(row):
    
    # ref is activator
    if row.activation_ref >= 1:
        if row.activation_fold_change_log2 <= -1 and row.activation_alt <= 1 and row.activation_alt >= -1:
            return "activation loss"
        elif not pd.isnull(row.activation_fold_change_log2):
            
            # only consider iso to be rewired if foldchange > 2x
            if np.abs(row.activation_fold_change_log2) >= 1:
                return "rewire"
            else:
                return "similar"

        else:
            return "NA"
    
    # ref is repressor
    elif row.activation_ref <= -1:
        if row.activation_fold_change_log2 >= 1 and row.activation_alt <= 1 and row.activation_alt >= -1:
            return "repression loss"
        elif not pd.isnull(row.activation_fold_change_log2):
            
            # only consider iso to be rewired if foldchange > 2x
            if np.abs(row.activation_fold_change_log2) >= 1:
                return "rewire"
            else:
                return "similar"
        else:
            return "NA"
        
    # no ref data so can't make conclusions
    elif pd.isnull(row.activation_ref):
        return "NA"
    
    # ref is middling so can be GoF
    else:
        if row.activation_fold_change_log2 >= 1:
            return "activation GoF"
        elif row.activation_fold_change_log2 <= -1:
            return "repression GoF"
        
        # if both isoforms are middling, consider similar
        elif not pd.isnull(row.activation_fold_change_log2):
            return "similar"
        
        else:
            return "NA"
        
pairs["m1h_cat"] = pairs.apply(m1h_cat, axis=1)
pairs.m1h_cat.value_counts()


# ### Y1H

# In[17]:


def y1h_cat(row):
    if row.n_positive_PDI_ref_filtered > 0:
        if row.n_positive_PDI_alt_filtered == 0:
            return "PDI loss"
        elif row.n_shared_PDI == row.n_PDI_successfully_tested_in_ref_and_alt:
            return "no PDI change"
        elif pd.isnull(row.n_positive_PDI_alt_filtered):
            return "NA"
        else:
            return "PDI rewire"
    elif row.n_positive_PDI_ref_filtered == 0:
        if row.n_positive_PDI_alt_filtered > 0:
            return "PDI gain"
        else:
            return "NA"
    else:
        return "NA"
    
pairs["y1h_cat"] = pairs.apply(y1h_cat, axis=1)
pairs.y1h_cat.value_counts()


# ### Y2H

# In[18]:


def y2h_cat(row):
    if row.dimer_ppi == "loses all" or row.tf_cofactor_ppi == "loses all" or row.tf_signalling_ppi == "loses all":
        n = []
        if row.dimer_ppi == "loses all":
            n.append("dimer")
        if row.tf_cofactor_ppi == "loses all":
            n.append("cofactor")
        if row.tf_signalling_ppi == "loses all":
            n.append("signalling")
        s = ",".join(n)
        s = "PPI loss: %s" % s
        return s
    
    elif row.n_positive_PPI_ref > 0 and row.n_positive_PPI_alt == 0:
        return "PPI loss: all"
    
    elif row.dimer_ppi == "retains all" and row.tf_cofactor_ppi == "retains all" and row.tf_signalling_ppi == "retains all":
        if row.other_than_dimer_ppi == "retains all":
            return "no PPI change (all PPIs)"
        else:
            return "no PPI change (important PPIs)"
    
    elif pd.isnull(row.dimer_ppi) and pd.isnull(row.tf_cofactor_ppi) and pd.isnull(row.tf_signalling_ppi) and pd.isnull(row.other_than_dimer_ppi) and pd.isnull(row.tf_tf_ppi):
        return "NA"
    
    else:
        
        # all PPIs retained but some above categories null so missed in if statements
        if row.PPI_jaccard == 1:
            return "no PPI change (all PPIs)"
        else:
            return "PPI rewire"
    
pairs["y2h_cat"] = pairs.apply(y2h_cat, axis=1)
pairs.y2h_cat.value_counts()


# ## localization

# In[19]:


con_pairs, con_df = load_condensate_data()


# In[20]:


pairs = pairs.merge(con_pairs[["gene_symbol", "clone_acc_ref", "clone_acc_alt",
                               'condensates_observed_HEK_ref',
       'condensates_observed_U2OS_ref', 'HEK_Condensate_ref',
       'U2OS_Condensate_ref', 'localization_HEK_ref', 'localization_U2OS_ref',
       'condensates_observed_HEK_alt', 'condensates_observed_U2OS_alt',
       'HEK_Condensate_alt', 'U2OS_Condensate_alt', 'localization_HEK_alt',
       'localization_U2OS_alt', 'condensate_cat_HEK', 'condensate_cat_U2OS',
       'condensate_cat_merged_HEK', 'condensate_cat_merged_U2OS',
       'condensate_cat_only_HEK', 'condensate_cat_only_U2OS',
       'localization_cat_HEK', 'localization_cat_U2OS',
       'condensate_cat_only_detailed_HEK', 'condensate_cat_only_detailed_U2OS',
       'condensate_or_loc_change_HEK', 'condensate_or_loc_change_U2OS',
       'combined_cat_HEK', 'combined_cat_U2OS',
       'condensate_or_loc_change_both', 'condensate_cat_detailed_HEK',
       'condensate_cat_detailed_U2OS']], how="left", on=["gene_symbol", "clone_acc_ref", "clone_acc_alt"])


# In[21]:


pairs.localization_HEK_ref.value_counts()


# In[22]:


def loc_cat_HEK(row):
    if row.localization_HEK_ref == "nucleus" or row.localization_HEK_ref == "both":
        if row.localization_HEK_alt == "cytoplasm":
            return "localization loss"
        elif row.localization_HEK_alt == "nucleus" or row.localization_HEK_alt == "both":
            return "no localization change"
        elif not pd.isnull(row.localization_HEK_alt):
            return "localization change"
        else:
            return "NA"
        
    # if ref is exclusively cytoplasm, call loc change but not loc loss
    elif row.localization_HEK_ref == "cytoplasm":
        if row.localization_HEK_alt == "cytoplasm":
            return "no localization change"
        elif row.localization_HEK_alt == "nucleus" or row.localization_HEK_alt == "both":
            return "localization change"
        elif not pd.isnull(row.localization_HEK_alt):
            return "localization change"
        else:
            return "NA"
    
    else:
        return "NA"
    
def loc_cat_U2OS(row):
    if row.localization_U2OS_ref == "nucleus" or row.localization_U2OS_ref == "both":
        if row.localization_U2OS_alt == "cytoplasm":
            return "localization loss"
        elif row.localization_U2OS_alt == "nucleus" or row.localization_U2OS_alt == "both":
            return "no localization change"
        elif not pd.isnull(row.localization_U2OS_alt):
            return "localization change"
        else:
            return "NA"
        
    # if ref is exclusively cytoplasm, call loc change but not loc loss
    elif row.localization_U2OS_ref == "cytoplasm":
        if row.localization_U2OS_alt == "cytoplasm":
            return "no localization change"
        elif row.localization_U2OS_alt == "nucleus" or row.localization_U2OS_alt == "both":
            return "localization change"
        elif not pd.isnull(row.localization_U2OS_alt):
            return "localization change"
        else:
            return "NA"
    
    else:
        return "NA"
        
def loc_cat(row):
    if row.loc_cat_HEK == "localization loss" or row.loc_cat_U2OS == "localization loss":
        return "localization loss"
    elif row.loc_cat_HEK == "localization change" or row.loc_cat_U2OS == "localization change":
        return "localization change"
    elif row.loc_cat_HEK == "no localization change" and row.loc_cat_U2OS == "no localization change":
        return "no localization change"
    else:
        return "NA"
    
pairs["loc_cat_HEK"] = pairs.apply(loc_cat_HEK, axis=1)
pairs["loc_cat_U2OS"] = pairs.apply(loc_cat_U2OS, axis=1)
pairs["loc_cat"] = pairs.apply(loc_cat, axis=1)
pairs.loc_cat.value_counts()


# In[23]:


len(pairs[pairs["loc_cat"] != "NA"])


# ## 3. categorize negative regulators
# include any whose DBD loss is >= 10%

# In[24]:


def dn_cat(row, dbd_thresh=10):
    
    # if activity loss
    if "loss" in row.m1h_cat:
        if "loss" in row.y1h_cat and "loss: all" in row.y2h_cat:
            return "likely nf"
        else:
            if row.y1h_cat == "NA" and row.y2h_cat == "NA":
                return "NA"
            else:
                n = ["activ"]
                if "loss" in row.y1h_cat:
                    n.append("PDIs")
                if "loss" in row.y2h_cat:
                    n.append("PPIs")
                if row.dbd_pct_lost >= dbd_thresh:
                    n.append("DBD loss")
                s = ",".join(n)
                s = "DN (%s)" % s
                return s
    
    # otherwise, if no evidence of m1h activity
    elif row.activation_alt <= 1 and row.activation_alt >= -1:
        if "loss" in row.y1h_cat and "loss: all" in row.y2h_cat:
            return "likely nf"
        else:
            if row.y1h_cat == "NA" and row.y2h_cat == "NA":
                return "NA"
            else:
                n = []
                if "loss" in row.y1h_cat:
                    n.append("PDIs")
                if "loss" in row.y2h_cat:
                    n.append("PPIs")
                if row.dbd_pct_lost >= dbd_thresh:
                    n.append("DBD loss")
                
                if len(n) > 0:
                    s = ",".join(n)
                    s = "DN (%s)" % s
                    
                else:
                    
                    # if m1h category is similar, and y1h/y2h are also no change, consider similar
                    if row.m1h_cat == "similar":
                        if row.y1h_cat == "no PDI change" and "no PPI change" in row.y2h_cat:
                            s = "similar"
                        elif row.y1h_cat == "no PDI change" and row.y2h_cat == "NA":
                            s = "similar"
                        elif row.y1h_cat == "NA" and "no PPI change" in row.y2h_cat:
                            s = "similar"
                        else:
                            s = "rewire"
                    else:
                        s = "rewire"
                
                return s
    
    # otherwise, if no m1h data
    elif pd.isnull(row.activation_alt):
        if "loss" in row.y1h_cat and "loss: all" in row.y2h_cat:
            return "NA"
        else:
            if row.y1h_cat == "NA" and row.y2h_cat == "NA":
                return "NA"
            elif row.y1h_cat != "NA" and row.y2h_cat != "NA":
                n = []
                if "loss" in row.y1h_cat:
                    n.append("PDIs")
                if "loss" in row.y2h_cat:
                    n.append("PPIs")
                if row.dbd_pct_lost >= dbd_thresh:
                    n.append("DBD loss")
                
                if len(n) > 0:
                    s = ",".join(n)
                    s = "DN (%s)" % s
                else:
                    
                    # if both y1h and y2h are similar, consider similar
                    if row.y1h_cat == "no PDI change" and "no PPI change" in row.y2h_cat:
                        s = "similar"
                    else:
                        s = "rewire"
                
                return s
            
            else:
                return "NA"
            
    # otherwise, if evidence of m1h functionality
    else:
        if row.y1h_cat == "NA" and row.y2h_cat == "NA":
            return "NA"
        else:
            n = []
            if "loss" in row.y1h_cat:
                n.append("PDIs")
            if "loss" in row.y2h_cat:
                n.append("PPIs")
            if row.dbd_pct_lost >= dbd_thresh:
                n.append("DBD loss")

            if len(n) > 0:
                s = ",".join(n)
                s = "DN (%s)" % s
            else:
                
                # if m1h category is similar, and y1h/y2h are also no change, consider similar
                if row.m1h_cat == "similar":
                    if row.y1h_cat == "no PDI change" and "no PPI change" in row.y2h_cat:
                        s = "similar"
                    elif row.y1h_cat == "no PDI change" and row.y2h_cat == "NA":
                        s = "similar"
                    elif row.y1h_cat == "NA" and "no PPI change" in row.y2h_cat:
                        s = "similar"
                    else:
                        s = "rewire"
                else:
                    s = "rewire"

            return s
            
pairs["dn_cat"] = pairs.apply(dn_cat, axis=1)
pairs.dn_cat.value_counts()


# In[25]:


pairs[pairs['dn_cat'] == 'DN (activ,PDIs,PPIs,DBD loss)']


# In[26]:


# now update anything based on localization calls
# but only if we were initially able to categorize it as DN/sim/rw
# (NAs stay NAs as loc is not evidence of function on its own)
def dn_loc(row):
    if row.dn_cat == "NA":
        return "NA"
    elif row.dn_cat == "likely nf":
        return "likely nf"
    else:
        
        # loc loss = DN
        if row.loc_cat == "localization loss":
            if row.dn_cat == "rewire" or row.dn_cat == "similar":
                return "DN (loc loss)"
            else:
                s = row.dn_cat[:-1]
                s += ", loc loss)"
                return s
        
        # loc change = rw
        elif row.loc_cat == "localization change":
            if row.dn_cat == "rewire" or row.dn_cat == "similar":
                return "rewire"
            else:
                return row.dn_cat
        
        else:
            return row.dn_cat
        
pairs["dn_cat_plus_loc"] = pairs.apply(dn_loc, axis=1)
pairs["dn_cat_plus_loc"].value_counts()


# In[27]:


pairs["dn_cat"] = pairs["dn_cat_plus_loc"]


# In[28]:


pairs[pairs["dn_cat"] == "DN (loc loss)"][["gene_symbol", "clone_acc_ref", "clone_acc_alt",
                                           "m1h_cat", "y1h_cat", "y2h_cat", "dbd_pct_lost", "loc_cat_HEK", 
                                           "loc_cat_U2OS", "loc_cat", "dn_cat"]]


# In[29]:


# double check that rewirers aren't all the same across axes


# In[30]:


to_plot = pairs[pairs["dn_cat"].isin(["rewire", "similar"])][["gene_symbol", "clone_acc_ref", "clone_acc_alt",
                                              "PPI_jaccard", "PDI_jaccard", "dbd_pct_lost",
                                              "activation_fold_change_log2", "loc_cat", "dn_cat"]]

# make cols that are easier to visualize on plot
to_plot["fc_abs_activ"] = 2**-np.abs(to_plot["activation_fold_change_log2"])
to_plot["1m_dbd_pct"] = (100-to_plot["dbd_pct_lost"])/100

# comment this out if we delete loc stuff
def loc_quant(row):
    if row.loc_cat == "NA":
        return np.nan
    elif row.loc_cat == "localization loss":
        return 0
    elif row.loc_cat == "no localization change":
        return 1
    else:
        return 0.5
to_plot["loc_quant"] = to_plot.apply(loc_quant, axis=1)


# In[31]:


isos = to_plot.copy()
isos = isos.sort_values(by="dn_cat")
isos = isos.reset_index().reset_index()

theoretical_iso = {"level_0": [isos.level_0.max()+1], "index": [0], "gene_symbol": ["theoretical"],
                   "clone_acc_ref": ["theoretical"], "clone_acc_alt": ["theoretical"],
                   "PPI_jaccard": [1.0], "PDI_jaccard": [1.0], "dbd_pct_lost": [0.0],
                   "activation_fold_change_log2": [np.log2(1)], "fc_abs_activ": [1.0],
                   "1m_dbd_pct": [1.0], "loc_cat": ["no localization change"], "loc_quant": [1.0], 
                   "dn_cat": ["theoretical"]}

isos = pd.concat([isos, pd.DataFrame.from_dict(theoretical_iso)])
isos


# In[32]:


columns_to_plot = ["level_0", "PPI_jaccard", "PDI_jaccard", "fc_abs_activ", "1m_dbd_pct", "loc_quant"]
column_titles = ["", "PPI\n(jaccard)", "PDI\n(jaccard)",
                 "activ.\n(% of ref.)", "DBD\n(% of ref.)", "∆ localization"]

# separate into similar v rewires
isos_sim = isos[isos["dn_cat"].isin(["similar", "theoretical"])]
isos_rw = isos[isos["dn_cat"].isin(["rewire"])]

colors = sns.color_palette("husl", n_colors=len(isos_rw))


# ### similar first

# In[33]:


# Creating the figure and axes (subplots) aligned in a single row without shared y-axis
gs_kw = dict(width_ratios=[0.75, 1, 1, 1, 1, 1])
fig, axs = plt.subplots(1, 6, sharey=False, gridspec_kw=gs_kw)
fig.set_size_inches(w=4, h=4)
df = isos_sim

# Adjusting the space between subplots for better alignment of lines
fig.subplots_adjust(wspace=0.05)

# zorder of axis labels
plt.rcParams["axes.axisbelow"] = False

# Re-plotting each column in a separate subplot
for i, col in enumerate(columns_to_plot):
    axs[i].clear()  # Clear previous axes to avoid duplication
    axs[i].set_title(column_titles[i])  # Setting the title of the subplot to the column name
    axs[i].set_xticks([])  # Hide x-axis labels

    # Setting different limits for y-axis based on column data
    axs[i].set_ylim(np.nanmin(df[col]) - 0.1, np.nanmax(df[col]) + 0.1)
    
    j = 0
    black_js = []
    for name, cat in zip(df['clone_acc_alt'], df['dn_cat']):
        if name == "theoretical":
            color = "black"
            black_js.append(j)
        elif cat == "similar":
            color = "black"
            black_js.append(j)
        else:
            color = colors[j]
        
        # Mark the value with a dot
        axs[i].plot(0, df[col].iloc[j], color=color, marker='o', 
                    markersize=6)
        j += 1
    
    axs[i].spines['top'].set_visible(False)
    axs[i].spines['bottom'].set_visible(False)
    
    if i == 0:
        axs[i].spines['right'].set_visible(False)
    elif i == 5:
        axs[i].yaxis.tick_right()
        axs[i].yaxis.set_label_position("right")
        axs[i].set_ylabel("Percent similarity")
    else:
        axs[i].spines['right'].set_visible(False)
        axs[i].yaxis.set_tick_params(labelleft=False)
        axs[i].set_yticks([])

# Remove previous lines
fig.lines = []

# set axis limits
tickpos = list(df["level_0"])
axs[0].set_ylim((tickpos[0]-0.5), (tickpos[-1]+0.5))
axs[1].set_ylim((-0.05, 1.05))
axs[2].set_ylim((-0.05, 1.05))
axs[3].set_ylim((-0.05, 1.05))
axs[4].set_ylim((-0.05, 1.05))
axs[5].set_ylim((-0.05, 1.05))


# Data preparation for plotting
data_lines = []
for _, row in df.iterrows():
    data_lines.append([row[columns_to_plot[0]], row[columns_to_plot[1]], 
                       row[columns_to_plot[2]], row[columns_to_plot[3]],
                       row[columns_to_plot[4]], row[columns_to_plot[5]]])

# Re-connecting dots across subplots with lines
for j, data_line in enumerate(data_lines):
    if j in black_js:
        color = "black"
    else:
        color = colors[j]
        
    skip_flag = False
    
    # Adding lines between subplots; adjusting coordinates for subplot borders
    for k in range(len(axs)-1):
        # We use figure coordinates to draw lines between subplots
        transFigure = fig.transFigure.inverted()
        
        
        if not skip_flag:
            pt1_k = k
            pt2_k = k+1
            pt1 = data_line[pt1_k]
            pt2 = data_line[pt2_k]
        
        # if we prev skipped due to na, keep the same pt1 and use current pt2
        else:
            pt1_k = k-1
            pt2_k = k+1
            pt1 = data_line[pt1_k]
            pt2 = data_line[pt2_k]
            
        if np.isnan(pt2):
            #print("skip flag true")
            skip_flag = True
            continue
        else:
            skip_flag = False
        
        #print("pt1: %s | pt2: %s | ax1: %s | ax2: %s" % (pt1, pt2, pt1_k, pt2_k))
        
        # Get coordinates of the points in figure coordinate system
        coord1 = transFigure.transform(axs[pt1_k].transData.transform([0, pt1]))
        coord2 = transFigure.transform(axs[pt2_k].transData.transform([0, pt2]))

        # Calculate the space adjustment based on the subplot spacing
        #space_adjustment = 0.015 + 0.23 * k

        # Adding space adjustment to x-coordinate for accurate alignment
        line = plt.Line2D((coord1[0], coord2[0]), 
                          (coord1[1], coord2[1]), 
                          transform=fig.transFigure, color=color, linestyle="dashed")
        fig.add_artist(line)
        


# relabel ticks
ticklabels = list(isos["clone_acc_alt"])
axs[0].set_yticks(tickpos)
#_ = axs[0].set_yticklabels(ticklabels)


# In[34]:


# Creating the figure and axes (subplots) aligned in a single row without shared y-axis
gs_kw = dict(width_ratios=[0.75, 1, 1, 1, 1, 1])
fig, axs = plt.subplots(1, 6, sharey=False, gridspec_kw=gs_kw)
fig.set_size_inches(w=4, h=10)
df = isos_rw

# Adjusting the space between subplots for better alignment of lines
fig.subplots_adjust(wspace=0.05)

# zorder of axis labels
plt.rcParams["axes.axisbelow"] = False

# Re-plotting each column in a separate subplot
for i, col in enumerate(columns_to_plot):
    axs[i].clear()  # Clear previous axes to avoid duplication
    axs[i].set_title(column_titles[i])  # Setting the title of the subplot to the column name
    axs[i].set_xticks([])  # Hide x-axis labels

    # Setting different limits for y-axis based on column data
    axs[i].set_ylim(np.nanmin(df[col]) - 0.1, np.nanmax(df[col]) + 0.1)
    
    j = 0
    black_js = []
    for name, cat in zip(df['clone_acc_alt'], df['dn_cat']):
        if name == "theoretical":
            color = "black"
            black_js.append(j)
        elif cat == "similar":
            color = "black"
            black_js.append(j)
        else:
            color = colors[j]
        
        # Mark the value with a dot
        axs[i].plot(0, df[col].iloc[j], color=color, marker='o', 
                    markersize=6)
        j += 1
    
    axs[i].spines['top'].set_visible(False)
    axs[i].spines['bottom'].set_visible(False)
    
    if i == 0:
        axs[i].spines['right'].set_visible(False)
    elif i == 5:
        axs[i].yaxis.tick_right()
        axs[i].yaxis.set_label_position("right")
        axs[i].set_ylabel("Percent similarity")
    else:
        axs[i].spines['right'].set_visible(False)
        axs[i].yaxis.set_tick_params(labelleft=False)
        axs[i].set_yticks([])

# Remove previous lines
fig.lines = []

# set axis limits
tickpos = list(df["level_0"])
axs[0].set_ylim((tickpos[0]-0.5), (tickpos[-1]+0.5))
axs[1].set_ylim((-0.05, 1.05))
axs[2].set_ylim((-0.05, 1.05))
axs[3].set_ylim((-0.05, 1.05))
axs[4].set_ylim((-0.05, 1.05))
axs[5].set_ylim((-0.05, 1.05))


# Data preparation for plotting
data_lines = []
for _, row in df.iterrows():
    data_lines.append([row[columns_to_plot[0]], row[columns_to_plot[1]], 
                       row[columns_to_plot[2]], row[columns_to_plot[3]],
                       row[columns_to_plot[4]], row[columns_to_plot[5]]])

# Re-connecting dots across subplots with lines
for j, data_line in enumerate(data_lines):
    if j in black_js:
        color = "black"
    else:
        color = colors[j]
        
    skip_flag = False
    
    # Adding lines between subplots; adjusting coordinates for subplot borders
    for k in range(len(axs)-1):
        # We use figure coordinates to draw lines between subplots
        transFigure = fig.transFigure.inverted()
        
        
        if not skip_flag:
            pt1_k = k
            pt2_k = k+1
            pt1 = data_line[pt1_k]
            pt2 = data_line[pt2_k]
        
        # if we prev skipped due to na, keep the same pt1 and use current pt2
        else:
            pt1_k = k-1
            pt2_k = k+1
            pt1 = data_line[pt1_k]
            pt2 = data_line[pt2_k]
            
        if np.isnan(pt2):
            #print("skip flag true")
            skip_flag = True
            continue
        else:
            skip_flag = False
        
        #print("pt1: %s | pt2: %s | ax1: %s | ax2: %s" % (pt1, pt2, pt1_k, pt2_k))
        
        # Get coordinates of the points in figure coordinate system
        coord1 = transFigure.transform(axs[pt1_k].transData.transform([0, pt1]))
        coord2 = transFigure.transform(axs[pt2_k].transData.transform([0, pt2]))

        # Calculate the space adjustment based on the subplot spacing
        #space_adjustment = 0.015 + 0.23 * k

        # Adding space adjustment to x-coordinate for accurate alignment
        line = plt.Line2D((coord1[0], coord2[0]), 
                          (coord1[1], coord2[1]), 
                          transform=fig.transFigure, color=color, linestyle="dashed")
        fig.add_artist(line)
        


# relabel ticks
ticklabels = list(isos["clone_acc_alt"])
axs[0].set_yticks(tickpos)
#_ = axs[0].set_yticklabels(ticklabels)


# ## 4. summary plots of DN categorization

# In[35]:


pairs["dn_short"] = pairs["dn_cat"].str.split(" ", expand=True)[0]
pairs.dn_short.value_counts()


# In[36]:


pairs[pairs["dn_short"] == "likely"]


# In[37]:


pairs[pairs["dn_cat"] == "DN (loc loss)"]


# In[38]:


pairs[pairs["gene_symbol"] == "FOXP2"][["gene_symbol", "clone_acc_alt", "clone_acc_ref",
                                        "m1h_cat", "y1h_cat", "y2h_cat", "loc_cat", "dn_cat_plus_loc"]]


# In[39]:


def mech_bool(row, mech_col):
    if "DN" in row.dn_cat:
        if mech_col in row.dn_cat:
            return True
        else:
            return False
    else:
        return np.nan
    
pairs["dn_ppi"] = pairs.apply(mech_bool, mech_col="PPIs", axis=1)
pairs["dn_pdi"] = pairs.apply(mech_bool, mech_col="PDIs", axis=1)
pairs["dn_activ"] = pairs.apply(mech_bool, mech_col="activ", axis=1)
pairs["dn_dbd"] = pairs.apply(mech_bool, mech_col="DBD loss", axis=1)
pairs["dn_loc"] = pairs.apply(mech_bool, mech_col="loc loss", axis=1)
pairs[pairs["dn_short"] == "DN"].sample(5)


# In[40]:


fig = plt.figure(figsize=(1.5, 1.75))

ax = sns.countplot(data=pairs, x="dn_short", palette=sns.color_palette("Set2"),
                   order=["DN", "rewire", "similar", "NA"])
ax.set_xticklabels(["putative DN", "putative re-wirer", "similar", "NA"], ha="right", va="top", rotation=30)
ax.set_xlabel("")
ax.set_ylabel("Number of alternative isoforms")

fig.savefig("../../figures/fig7/DN_countplot.pdf", dpi="figure", bbox_inches="tight")


# In[41]:


ppis = list(set(list(pairs[pairs["dn_ppi"] == True]["clone_acc_alt"])))
pdis = list(set(list(pairs[pairs["dn_pdi"] == True]["clone_acc_alt"])))
activ = list(set(list(pairs[pairs["dn_activ"] == True]["clone_acc_alt"])))
dbd = list(set(list(pairs[pairs["dn_dbd"] == True]["clone_acc_alt"])))
loc = list(set(list(pairs[pairs["dn_loc"] == True]["clone_acc_alt"])))

contents = {"loss of PPIs": ppis, "loss of PDIs": pdis, "loss of activity": activ, "loss of DBD": dbd,
            "loss of localization": loc}
contents = upsetplot.from_contents(contents)

all_dn = set(ppis).union(set(pdis)).union(set(activ)).union(set(dbd))
print(len(all_dn))

fig = plt.figure(figsize=(3, 2))
d = plot(contents, fig=fig, sort_by="cardinality", show_counts=True, element_size=12, 
     intersection_plot_elements=4, totals_plot_elements=3)
d["intersections"].set_ylabel("Number of\nalternative isoforms")
d["intersections"].grid(False)
d["totals"].grid(False)

fig.savefig("../../figures/fig7/DN_negreg_upset.pdf", dpi="figure", bbox_inches="tight")


# In[42]:


rw = pairs[pairs["dn_cat"] == "rewire"]
ppis = list(set(list(rw[rw["y2h_cat"] == "PPI rewire"]["clone_acc_alt"])))
pdis = list(set(list(rw[rw["y1h_cat"].str.contains("PDI")]["clone_acc_alt"])))
activ = list(set(list(rw[rw["m1h_cat"] != "NA"]["clone_acc_alt"])))
loc = list(set(list(rw[rw["loc_cat"] == "localization change"]["clone_acc_alt"])))

contents = {"change in PPIs": ppis, "change in PDIs": pdis, "change in activity": activ, 
            "change in localization": loc}
contents = upsetplot.from_contents(contents)

fig = plt.figure(figsize=(3, 2))
d = plot(contents, fig=fig, sort_by="cardinality", show_counts=True, element_size=12, 
         intersection_plot_elements=4, totals_plot_elements=3)
d["intersections"].set_ylabel("Number of\nalternative isoforms")
d["intersections"].grid(False)
d["totals"].grid(False)

fig.savefig("../../figures/fig7/DN_rewire_upset.pdf", dpi="figure", bbox_inches="tight")


# In[43]:


ppis = list(set(list(pairs[pairs["y2h_cat"] != "NA"]["clone_acc_alt"])))
pdis = list(set(list(pairs[pairs["y1h_cat"] != "NA"]["clone_acc_alt"])))
activ = list(set(list(pairs[pairs["m1h_cat"] != "NA"]["clone_acc_alt"])))
loc = list(set(list(pairs[pairs["loc_cat"] != "NA"]["clone_acc_alt"])))

contents = {"assessed PPIs": ppis, "assessed PDIs": pdis, "assessed activity": activ, "assessed localization": loc}
contents = upsetplot.from_contents(contents)

all_as = set(ppis).union(set(pdis)).union(set(activ))
print(len(all_as))

fig = plt.figure(figsize=(3, 2))
d = plot(contents, fig=fig, sort_by="cardinality", show_counts=True, element_size=12, 
         intersection_plot_elements=4, totals_plot_elements=3)
d["intersections"].set_ylabel("Number of\nalternative isoforms")
d["intersections"].grid(False)
d["totals"].grid(False)

fig.savefig("../../figures/fig7/DN_pairs_assessed_upset.pdf", dpi="figure", bbox_inches="tight")


# In[44]:


y = np.array([len(pairs[pairs["dn_short"] == "rewire"]), 
              len(pairs[pairs["dn_short"] == "DN"]),
              len(pairs[pairs["dn_short"] == "similar"]),
              len(pairs[pairs["dn_short"] == "NA"]), 
              len(pairs[pairs["dn_short"] == "likely"])])
labels = ["rewire", "negative regulator", "similar", "NA", "likely non-functional"]
colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], sns.color_palette("Set2")[0],
          "lightgray", "darkgray"]

fig = plt.figure(figsize=(1.75, 1.75))
ws, ls, ns = plt.pie(y, labels=labels, colors=colors, autopct='%1.1f%%', startangle=-45, explode=(0, 0.1, 0, 0, 0))
for w in ws:
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
ns[4].set_color("white")

fig.savefig("../../figures/fig7/dn_pie.incl_NA.pdf", dpi="figure", bbox_inches="tight")


# In[45]:


ys = np.array([len(pairs[pairs["dn_short"] == "similar"]),
               len(pairs[pairs["dn_short"] == "rewire"]), 
               len(pairs[pairs["dn_short"] == "DN"]),
               len(pairs[pairs["dn_short"] == "likely"])])
labels = ["similar", "rewirer", "negative regulator", "likely non-functional"]
colors = [sns.color_palette("Set2")[0], sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], "darkgray"]

fig, ax = plt.subplots(figsize=(2.0, 2.0), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.1f%%', startangle=90, 
                    explode=(0.05, 0.05, 0.05, 0.15))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")
ns[3].set_text("")

fig.savefig("../../figures/fig7/dn_pie.no_NA.pdf", dpi="figure", bbox_inches="tight")


# In[46]:


outer_ys = np.array([len(pairs[(pairs["dn_short"] == "similar")]),
                     len(pairs[(pairs["dn_short"] == "rewire")]), 
                     len(pairs[(pairs["dn_short"] == "DN")])])
outer_labels = ["similar\n%0.0f%%" % round((outer_ys[0]/np.sum(outer_ys)*100)),
                "rewirer\n%0.0f%%" % round((outer_ys[1]/np.sum(outer_ys)*100)), 
                "negative regulator\n%0.0f%%" % round((outer_ys[2]/np.sum(outer_ys)*100))]
outer_colors = [sns.color_palette("Set2")[0], sns.color_palette("Set2")[2], sns.color_palette("Set2")[1]]

inner_ys = np.array([len(pairs[(pairs["dn_short"] == "similar") & (pairs["is_alt_novel_isoform"])]), 
                     len(pairs[(pairs["dn_short"] == "similar") & (~pairs["is_alt_novel_isoform"])]),
                     len(pairs[(pairs["dn_short"] == "rewire") & (pairs["is_alt_novel_isoform"])]), 
                     len(pairs[(pairs["dn_short"] == "rewire") & (~pairs["is_alt_novel_isoform"])]), 
                     len(pairs[(pairs["dn_short"] == "DN") & (pairs["is_alt_novel_isoform"])]),
                     len(pairs[(pairs["dn_short"] == "DN") & (~pairs["is_alt_novel_isoform"])])])
inner_colors = [sns.light_palette(sns.color_palette("Set2")[0])[3], 
                sns.light_palette(sns.color_palette("Set2")[0])[1], 
                sns.light_palette(sns.color_palette("Set2")[2])[3], 
                sns.light_palette(sns.color_palette("Set2")[2])[1], 
                sns.light_palette(sns.color_palette("Set2")[1])[3], 
                sns.light_palette(sns.color_palette("Set2")[1])[1]]
hatches = ['++', '', '++', '', '++', '']


fig, ax = plt.subplots(figsize=(2.2, 2.2), subplot_kw=dict(aspect="equal"))

o_ws, o_ls = ax.pie(outer_ys, colors=outer_colors, labels=outer_labels,
                    startangle=90, radius=1, wedgeprops=dict(width=0.3, edgecolor='w'))
i_ws, i_ls, i_ns = ax.pie(inner_ys, colors=inner_colors, autopct='%0.0f%%', 
                          startangle=90, radius=0.7, 
                          wedgeprops=dict(width=0.4, edgecolor='w'),
                          textprops={'fontsize': 5}, pctdistance=0.8)

for i, w in enumerate(i_ws):
    w.set(hatch=hatches[i])
    
ax.set_title("alternative TF isoform categories\n(%s reference-alternative pairs)" % (np.sum(outer_ys)))

fig.savefig("../../figures/fig7/dn_pie.novel_nested.pdf", dpi="figure", bbox_inches="tight")


# In[47]:


outer_ys = np.array([len(pairs[(pairs["dn_short"] == "similar")]),
                     len(pairs[(pairs["dn_short"] == "rewire")]), 
                     len(pairs[(pairs["dn_short"] == "DN")])])
outer_labels = ["similar\n%0.0f%%" % round((outer_ys[0]/np.sum(outer_ys)*100)),
                "rewirer\n%0.0f%%" % round((outer_ys[1]/np.sum(outer_ys)*100)), 
                "negative\nregulator\n%0.0f%%" % round((outer_ys[2]/np.sum(outer_ys)*100))]
outer_colors = [sns.color_palette("Set2")[0], sns.color_palette("Set2")[2], sns.color_palette("Set2")[1]]

inner_ys = np.array([len(pairs[(pairs["dn_short"] == "similar") & (pairs["is_alt_novel_isoform"])]), 
                     len(pairs[(pairs["dn_short"] == "similar") & (~pairs["is_alt_novel_isoform"])]),
                     len(pairs[(pairs["dn_short"] == "rewire") & (pairs["is_alt_novel_isoform"])]), 
                     len(pairs[(pairs["dn_short"] == "rewire") & (~pairs["is_alt_novel_isoform"])]), 
                     len(pairs[(pairs["dn_short"] == "DN") & (pairs["is_alt_novel_isoform"])]),
                     len(pairs[(pairs["dn_short"] == "DN") & (~pairs["is_alt_novel_isoform"])])])
inner_colors = [sns.light_palette(sns.color_palette("Set2")[0])[3], 
                sns.light_palette(sns.color_palette("Set2")[0])[1], 
                sns.light_palette(sns.color_palette("Set2")[2])[3], 
                sns.light_palette(sns.color_palette("Set2")[2])[1], 
                sns.light_palette(sns.color_palette("Set2")[1])[3], 
                sns.light_palette(sns.color_palette("Set2")[1])[1]]
hatches = ['++', '', '++', '', '++', '']


fig, ax = plt.subplots(figsize=(1.2, 1.2), subplot_kw=dict(aspect="equal"))

o_ws, o_ls = ax.pie(outer_ys, colors=outer_colors, labels=outer_labels,
                    startangle=90, radius=1, wedgeprops=dict(width=0.3, edgecolor='w'))
i_ws, i_ns = ax.pie(inner_ys, colors=inner_colors, 
                          startangle=90, radius=0.7, 
                          wedgeprops=dict(width=0.3, edgecolor='w'))

for i, w in enumerate(i_ws):
    w.set(hatch=hatches[i])
    
ax.set_title("alternative TF isoform categories\n(%s reference-alternative pairs)" % (np.sum(outer_ys)))

fig.savefig("../../figures/fig7/dn_pie.novel_nested.no_labels.pdf", dpi="figure", bbox_inches="tight")


# In[48]:


# create df that contains info about why alt isos were sorted into a particular category
# for DNs, include loss of DBD
tmp = pairs[pairs["dn_short"] == "DN"]
dn_pdi_change = len(tmp[tmp["dn_pdi"] == True])
dn_ppi_change = len(tmp[tmp["dn_ppi"] == True])
dn_activ_change = len(tmp[tmp["dn_activ"] == True])
dn_dbd_change = len(tmp[tmp["dn_dbd"] == True])
dn_loc_change = len(tmp[tmp["dn_loc"] == True])
tot_dn = len(tmp)

# for rewirers, include change in DBD (all will be <10% else they'd be DNs)
tmp = pairs[pairs["dn_short"] == "rewire"]
rw_pdi_change = len(tmp[~tmp["y1h_cat"].isin(["NA", "no PDI change"])])
rw_ppi_change = len(tmp[~tmp["y2h_cat"].isin(["NA", "no PPI change (important PPIs)"])])
rw_activ_change = len(tmp[~tmp["m1h_cat"].isin(["NA", "similar"])])
rw_dbd_change = len(tmp[tmp["dbd_pct_lost"] > 0])
rw_loc_change = len(tmp[tmp["loc_cat"] == "localization change"])
tot_rw = len(tmp)

# then include all assessed among categorized
tmp = pairs[pairs["dn_short"] != "NA"]
pdi_assessed = len(tmp[tmp["y1h_cat"] != "NA"])
ppi_assessed = len(tmp[tmp["y2h_cat"] != "NA"])
activ_assessed = len(tmp[tmp["m1h_cat"] != "NA"])
dbd_assessed = len(tmp[~pd.isnull(tmp["dbd_pct_lost"])])
loc_assessed = len(tmp[tmp["loc_cat"] != "NA"])
tot_assessed = len(tmp)

df = pd.DataFrame.from_dict({"DN": {"pdi": dn_pdi_change/tot_dn*100, 
                                    "ppi": dn_ppi_change/tot_dn*100,
                                    "activ": dn_activ_change/tot_dn*100, 
                                    "dbd": dn_dbd_change/tot_dn*100,
                                    "loc": dn_loc_change/tot_dn*100},
                             "rewire": {"pdi": rw_pdi_change/tot_rw*100,
                                        "ppi": rw_ppi_change/tot_rw*100,
                                        "activ": rw_activ_change/tot_rw*100,
                                        "dbd": rw_dbd_change/tot_rw*100,
                                        "loc": rw_loc_change/tot_rw*100},
                             "assessed": {"pdi": pdi_assessed/tot_assessed*100, 
                                    "ppi": ppi_assessed/tot_assessed*100,
                                    "activ": activ_assessed/tot_assessed*100,
                                          "dbd": dbd_assessed/tot_assessed*100,
                                          "loc": loc_assessed/tot_assessed*100}})

df = df.reset_index()
df


# In[49]:


colors = met_brewer.met_brew(name="Hokusai3", n=5, brew_type="discrete")
sns.palplot(colors)


# In[50]:


to_plot = pd.melt(df, id_vars="index")
to_plot


# In[51]:


fig = plt.figure(figsize=(2, 1))

ax = sns.barplot(data=to_plot, x="variable", y="value", hue="index", palette=colors,
                 order=["DN", "rewire"], edgecolor="black", linewidth=0.5,
                 hue_order=["pdi", "ppi", "activ", "dbd", "loc"])

ax.set_xlabel("")
ax.set_xticklabels(["negative regulators\n\n(loss of function)", "rewirers\n\n(change in function)"])

# # Define some hatches
# hatches = ['', 'xxxx', '', '', 'xxxx', '', '', 'xxxx', '', '', 'xxxx', '', '', 'xxxx', '']
# fills = ['', '', 'white', '', '', 'white', '', '', 'white', '', '', 'white', '', '', 'white']
# mpl.rcParams['hatch.linewidth'] = 0.5

# # Loop over the bars
# for i, bar in enumerate(ax.patches):
#     # Set a different hatch for each bar
#     bar.set_hatch(hatches[i])
#     fill = fills[i]
#     if fill == "white":
#         bar.set_alpha(0.5)

handles, previous_labels = ax.get_legend_handles_labels()
plt.legend(loc=2, bbox_to_anchor=(1.01, 1), frameon=False,
           handles=handles, labels=["PDIs", "PPIs", "activity", "DBD", "localization"])

ax.set_ylim((0, 100))
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
ax.set_ylabel("Percentage of\nalternative isoforms")

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/dn_function_unstacked_bar.pdf", dpi="figure", bbox_inches="tight")


# In[52]:


genes_w_dn = pairs[pairs["dn_short"] == "DN"][["gene_symbol", "family"]].drop_duplicates()
genes_w_rw = pairs[pairs["dn_short"] == "rewire"][["gene_symbol", "family"]].drop_duplicates()
genes_w_sim = pairs[pairs["dn_short"] == "similar"][["gene_symbol", "family"]].drop_duplicates()
tot_genes = pairs[["gene_symbol", "family"]].drop_duplicates()

tot_genes_per_f = tot_genes.groupby("family")["gene_symbol"].agg("count").reset_index()
dn_genes_per_f = genes_w_dn.groupby("family")["gene_symbol"].agg("count").reset_index()
rw_genes_per_f = genes_w_rw.groupby("family")["gene_symbol"].agg("count").reset_index()
sim_genes_per_f = genes_w_sim.groupby("family")["gene_symbol"].agg("count").reset_index()

family_cats = tot_genes_per_f.merge(dn_genes_per_f, 
                                    on="family", how="left")
family_cats = family_cats.merge(rw_genes_per_f, 
                                on="family", how="left").merge(sim_genes_per_f,
                                                               on="family",
                                                               how="left", suffixes=("_a", "_b"))
family_cats.fillna(0, inplace=True)
family_cats.columns = ["family", "tot", "dn", "rw", "sim"]

family_cats["tot_p"] = family_cats["tot"]/family_cats["tot"].sum(axis=0)*100
family_cats["dn_p"] = family_cats["dn"]/family_cats["dn"].sum(axis=0)*100
family_cats["rw_p"] = family_cats["rw"]/family_cats["rw"].sum(axis=0)*100
family_cats["sim_p"] = family_cats["sim"]/family_cats["sim"].sum(axis=0)*100
family_cats.sort_values(by="tot", ascending=False).head(11)
family_cats


# In[53]:


colors = met_brewer.met_brew(name="Renoir", n=11, brew_type="discrete")
sns.palplot(colors)


# In[54]:


fig, ax = plt.subplots(figsize=(1.5, 1.75))

xs = ["total", "negative regulator", "rewirer", "similar"]

b = np.zeros(4)
c = 0
for i, row in family_cats.sort_values(by="tot", ascending=False).head(11).iterrows():
    y = list(row[["tot_p", "dn_p", "rw_p", "sim_p"]])
    ax.bar(xs, y, bottom=b, label=row.family, color=colors[c])
    b = np.add(b, y)
    c += 1

ax.bar(xs, np.subtract([100, 100, 100, 100], b), bottom=b, label="smaller families", color="lightgray")
ax.set_ylabel("Percentage of categorized isoforms")
ax.set_xticklabels(xs, rotation=30, ha="right", va="top")

# add legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/dn_families_stacked_bar.pdf", dpi="figure", bbox_inches="tight")


# ## 5. are there differences in condensate formation between categories?

# In[55]:


def condensate_cat_either(row):
    if row.condensate_cat_only_HEK == "Difference" or row.condensate_cat_only_U2OS == "Difference":
        return "Difference"
    elif row.condensate_cat_only_HEK == "No difference" and row.condensate_cat_only_U2OS == "No difference":
        return "No difference"
    elif row.condensate_cat_only_HEK == "No difference" or row.condensate_cat_only_U2OS == "No difference":
        return "No difference"
    else:
        return np.nan
    
pairs["condensate_cat_either"] = pairs.apply(condensate_cat_either, axis=1)
pairs.condensate_cat_either.value_counts()


# In[56]:


def condensate_cat_both(row):
    if row.condensate_cat_only_HEK == "Difference" and row.condensate_cat_only_U2OS == "Difference":
        return "Difference"
    elif row.condensate_cat_only_HEK == "No difference" and row.condensate_cat_only_U2OS == "No difference":
        return "No difference"
    else:
        return np.nan
    
pairs["condensate_cat_both"] = pairs.apply(condensate_cat_both, axis=1)
pairs.condensate_cat_both.value_counts()


# In[57]:


from scipy.stats import fisher_exact


# In[58]:


fe = np.zeros((2, 2))

tmp_nonan = pairs[~pd.isnull(pairs["condensate_cat_either"])]

fe[0, 0] = len(tmp_nonan[(tmp_nonan["dn_short"] == "DN") & 
                    (tmp_nonan["condensate_cat_either"] == "Difference")].clone_acc_alt.unique())
fe[1, 0] = len(tmp_nonan[(tmp_nonan["dn_short"] != "DN") & 
                    (tmp_nonan["condensate_cat_either"] == "Difference")].clone_acc_alt.unique())
fe[0, 1] = len(tmp_nonan[(tmp_nonan["dn_short"] == "DN") &
                    (tmp_nonan["condensate_cat_either"] == "No difference")].clone_acc_alt.unique())
fe[1, 1] = len(tmp_nonan[(tmp_nonan["dn_short"] != "DN") & 
                    (tmp_nonan["condensate_cat_either"] == "No difference")].clone_acc_alt.unique())
fe


# In[59]:


p = fisher_exact(fe)[1]
print(p)


# In[60]:


tots = pd.DataFrame(tmp_nonan.dn_short.value_counts())
sig = pd.DataFrame(tmp_nonan[tmp_nonan["condensate_cat_either"] == "Difference"].dn_short.value_counts())
con_st = tots.join(sig, lsuffix="_tot", rsuffix="_diff")
con_st = con_st.loc[["DN", "rewire", "similar", "NA"]]
con_st = con_st/con_st.sum(axis=0)
con_st.fillna(0, inplace=True)
con_st


# In[61]:


fig, ax = plt.subplots(figsize=(0.4, 1.5))

xs = ["total", "diff. in condensates"]
y1 = list(con_st[["count_tot", "count_diff"]].loc["NA"])
y2 = list(con_st[["count_tot", "count_diff"]].loc["similar"])
b2 = np.add(y1, y2)
y3 = list(con_st[["count_tot", "count_diff"]].loc["rewire"])
b3 = np.add(b2, y3)
y4 = list(con_st[["count_tot", "count_diff"]].loc["DN"])

ax.bar(xs, y1, color=dn_pal["NA"], label="NA", edgecolor="black", linewidth=0.5)
ax.bar(xs, y2, bottom=y1, color=dn_pal["similar"], label="similar", edgecolor="black", linewidth=0.5)
ax.bar(xs, y3, bottom=b2, color=dn_pal["rewire"], label="rewire", edgecolor="black", linewidth=0.5)
ax.bar(xs, y4, bottom=b3, color=dn_pal["DN"], label="neg. reg.", edgecolor="black", linewidth=0.5)

# annotate pval
annotate_pval(ax, 0, 1, 1.025, 0, 1.025, p, fontsize-1)

# add legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)

ax.set_ylabel("Percentage of\nalternative isoforms")
ax.set_xticklabels(xs, rotation=30, ha="right", va="top")
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_yticklabels(["%s%%" % x for x in [0, 20, 40, 60, 80, 100]])
ax.set_ylim(0, 1)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_tick_params(length=0)

fig.savefig("../../figures/fig7/dn_stacked_bar_condensates.either.pdf", dpi="figure", bbox_inches="tight")


# In[62]:


fe = np.zeros((2, 2))

tmp_nonan = pairs[~pd.isnull(pairs["condensate_cat_both"])]

fe[0, 0] = len(tmp_nonan[(tmp_nonan["dn_short"] == "DN") & 
                    (tmp_nonan["condensate_cat_both"] == "Difference")].clone_acc_alt.unique())
fe[1, 0] = len(tmp_nonan[(tmp_nonan["dn_short"] != "DN") & 
                    (tmp_nonan["condensate_cat_both"] == "Difference")].clone_acc_alt.unique())
fe[0, 1] = len(tmp_nonan[(tmp_nonan["dn_short"] == "DN") &
                    (tmp_nonan["condensate_cat_both"] == "No difference")].clone_acc_alt.unique())
fe[1, 1] = len(tmp_nonan[(tmp_nonan["dn_short"] != "DN") & 
                    (tmp_nonan["condensate_cat_both"] == "No difference")].clone_acc_alt.unique())
fe


# In[63]:


p = fisher_exact(fe)[1]
print(p)


# In[64]:


tots = pd.DataFrame(tmp_nonan.dn_short.value_counts())
sig = pd.DataFrame(tmp_nonan[tmp_nonan["condensate_cat_both"] == "Difference"].dn_short.value_counts())
con_st = tots.join(sig, lsuffix="_tot", rsuffix="_diff")
con_st = con_st.loc[["DN", "rewire", "similar", "NA"]]
con_st = con_st/con_st.sum(axis=0)
con_st.fillna(0, inplace=True)
con_st


# In[65]:


fig, ax = plt.subplots(figsize=(0.4, 1.5))

xs = ["total", "diff. in condensates"]
y1 = list(con_st[["count_tot", "count_diff"]].loc["NA"])
y2 = list(con_st[["count_tot", "count_diff"]].loc["similar"])
b2 = np.add(y1, y2)
y3 = list(con_st[["count_tot", "count_diff"]].loc["rewire"])
b3 = np.add(b2, y3)
y4 = list(con_st[["count_tot", "count_diff"]].loc["DN"])

ax.bar(xs, y1, color=dn_pal["NA"], label="NA", edgecolor="black", linewidth=0.5)
ax.bar(xs, y2, bottom=y1, color=dn_pal["similar"], label="similar", edgecolor="black", linewidth=0.5)
ax.bar(xs, y3, bottom=b2, color=dn_pal["rewire"], label="rewire", edgecolor="black", linewidth=0.5)
ax.bar(xs, y4, bottom=b3, color=dn_pal["DN"], label="neg. reg.", edgecolor="black", linewidth=0.5)

# annotate pval
annotate_pval(ax, 0, 1, 1.025, 0, 1.025, p, fontsize-1)

# add legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)

ax.set_ylabel("Percentage of\nalternative isoforms")
ax.set_xticklabels(xs, rotation=30, ha="right", va="top")
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_yticklabels(["%s%%" % x for x in [0, 20, 40, 60, 80, 100]])
ax.set_ylim(0, 1)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_tick_params(length=0)

fig.savefig("../../figures/fig7/dn_stacked_bar_condensates.both.pdf", dpi="figure", bbox_inches="tight")


# ## 6. plot expression profiles of isoform categories

# In[66]:


# use same downsample as prev figs
means_gtex_downsample = df_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean()


# In[67]:


# calculate expression ratios - dev
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


# In[68]:


# calculate expression ratios - gtex
per_gene_gtex = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .transform('sum'))
f_gtex = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
        .mean())
f_gtex = f_gtex * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex = f_gtex * 100


# In[69]:


# calculate expression ratios -gtex downsampled
per_gene_gtex_ds = ((2 ** df_gtex.loc[:,metadata_gtex_dummy.index] - 1)
                   .groupby(genes_gtex)
                   .transform('sum'))

f_gtex_downsample = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1)
        .mean())
f_gtex_downsample = f_gtex_downsample * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex_downsample = f_gtex_downsample * 100


# In[70]:


# calculate gene-level tissue specificities
gene_dev_nonan_taus, gene_dev_nan_taus, gene_dev_array_max = calculate_tau(per_gene_dev.drop_duplicates())
gene_gtex_nonan_taus, gene_gtex_nan_taus, gene_gtex_array_max = calculate_tau(per_gene_gtex.drop_duplicates())
gene_gtex_ds_nonan_taus, gene_gtex_ds_nan_taus, gene_gtex_ds_array_max = calculate_tau(per_gene_gtex_ds.drop_duplicates())


# In[71]:


gene_taus = pd.DataFrame()
gene_taus["UID"] = per_gene_dev.drop_duplicates().index
gene_taus["dev_tau"] = gene_dev_nan_taus
gene_taus["gtex_tau"] = gene_gtex_nan_taus
gene_taus["gtex_ds_tau"] = gene_gtex_ds_nan_taus
gene_taus["gene_name"] = gene_taus["UID"].str.split("|", expand=True)[0]
gene_taus.sample(5)


# In[72]:


# join w pairs table
dev_ratios = f_dev.reset_index()
dev_ratios["clone_acc"] = dev_ratios["UID"].str.split(" ", expand=True)[0].astype(str)
dev_ratios = dev_ratios[dev_ratios["clone_acc"] != "noclone"]
len(dev_ratios)


# In[73]:


gtex_ratios = f_gtex.reset_index()
gtex_ratios["clone_acc"] = gtex_ratios["UID"].str.split(" ", expand=True)[0].astype(str)
gtex_ratios = gtex_ratios[gtex_ratios["clone_acc"] != "noclone"]
len(gtex_ratios)


# In[74]:


gtex_ds_ratios = f_gtex_downsample.reset_index()
gtex_ds_ratios["clone_acc"] = gtex_ds_ratios["UID"].str.split(" ", expand=True)[0].astype(str)
gtex_ds_ratios = gtex_ds_ratios[gtex_ds_ratios["clone_acc"] != "noclone"]
len(gtex_ds_ratios)


# In[75]:


dn_ref = pairs[["gene_symbol", "family", "clone_acc_ref", "is_ref_novel_isoform", "is_MANE_select_isoform_cloned",
             "dn_short"]].drop_duplicates()
dn_ref.columns = ["gene_name", "family", "tf1p0_id", "is_novel", "is_MANE_select", "dn_cat"]
dn_ref["dn_cat"] = "ref"
dn_ref["iso_status"] = "ref"


# In[76]:


dn_alt = pairs[["gene_symbol", "family", "clone_acc_alt", "is_alt_novel_isoform", "is_MANE_select_isoform_cloned",
             "dn_short"]].drop_duplicates()
dn_alt.columns = ["gene_name", "family", "tf1p0_id", "is_novel", "is_MANE_select", "dn_cat"]
dn_alt["is_MANE_select"] = False # assuming none of the alts are the MANE select
dn_alt["iso_status"] = "alt"


# In[77]:


dn_cats = pd.concat([dn_ref, dn_alt]).drop_duplicates()


# In[78]:


dev_ratios = dev_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id")
gtex_ratios = gtex_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id")
gtex_ds_ratios = gtex_ds_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id")
print(len(dev_ratios))
print(len(gtex_ratios))
print(len(gtex_ds_ratios))


# In[79]:


dn_cats = dn_cats.merge(gene_taus, on="gene_name")
print(len(dn_cats))
dn_cats.head()


# In[80]:


ref_expr = dn_cats.groupby(["gene_name", "family", "dn_cat", "dev_tau",
                            "gtex_tau", "gtex_ds_tau"])["tf1p0_id"].agg("count").reset_index()
ref_expr = ref_expr.pivot(index="gene_name",
                          columns="dn_cat", values="tf1p0_id")
ref_expr.fillna(0, inplace=True)


# In[81]:


def categorize_gene(row):
    if row.DN > 0 and row.rewire == 0 and row.similar == 0:
        return "DN"
    elif row.rewire > 0 and row.DN == 0 and row.similar == 0:
        return "rewire"
    elif row.similar > 0 and row.DN == 0 and row.rewire == 0:
        return "similar"
    elif row.NA > 0:
        return "NA"
    else:
        return "combination"
    
ref_expr["gene_cat"] = ref_expr.apply(categorize_gene, axis=1)
ref_expr.reset_index(inplace=True)
ref_expr = ref_expr.merge(dn_cats[["gene_name", "family", "dev_tau", "gtex_tau", "gtex_ds_tau"]],
                          on="gene_name").drop_duplicates()
print(len(ref_expr))
ref_expr.sample(5)


# In[82]:


fig, ax = nice_boxplot(ref_expr, "dev_tau", "gene_cat", dn_pal, ["rewire", "DN", "similar", "combination", "NA"], 
            [1.01, 1.045, 1.08, 0.99], 0.35, "", ["rewirer", "negative regulator", "similar", "combination", "NA"], 
            "gene-level tissue specificity (tau)", False, (0.3, 1.23), 
            "developmental gene expression\nclassified TF genes")

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticklabels(["rewirer", "negative regulator", "similar", "combination", "NA"], rotation=30, 
                   ha="right", va="top")
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xlabel("")
ax.set_ylabel("Gene-level tissue specificity (tau)", position=(0, 0.38))

# manually set left axis so it stops at 1.0
ax.set_ylim((0.45, 1.2))
ax.spines['left'].set_visible(False)
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
axes_to_data = ax.transAxes + ax.transData.inverted()
left_spine_in_data_coords = axes_to_data.transform((0, 0))
ax.plot([left_spine_in_data_coords[0], left_spine_in_data_coords[0]], [0.45, 1],
         color=ax.spines['bottom'].get_edgecolor(), linewidth=ax.spines['bottom'].get_linewidth())


fig.savefig("../../figures/fig7/DN_DevTau_Gene_Boxplot.pdf", dpi="figure", bbox_inches="tight")


# In[83]:


fig, ax = nice_boxplot(ref_expr, "gtex_ds_tau", "gene_cat", dn_pal, ["rewire", "DN", "similar", "combination", "NA"], 
            [1.01, 1.045, 1.08, 0.999], 0.35, "", ["rewirer", "negative regulator", "similar", "combination", "NA"], 
            "gene-level tissue specificity (tau)", False, (0.3, 1.23), 
            "developmental gene expression\nclassified TF genes")

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticklabels(["rewirer", "negative regulator", "similar", "combination", "NA"], rotation=30, 
                   ha="right", va="top")
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xlabel("")
ax.set_ylabel("Gene-level tissue specificity (tau)", position=(0, 0.38))

# manually set left axis so it stops at 1.0
ax.set_ylim((0.5, 1.2))
ax.spines['left'].set_visible(False)
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
axes_to_data = ax.transAxes + ax.transData.inverted()
left_spine_in_data_coords = axes_to_data.transform((0, 0))
ax.plot([left_spine_in_data_coords[0], left_spine_in_data_coords[0]], [0.5, 1],
         color=ax.spines['bottom'].get_edgecolor(), linewidth=ax.spines['bottom'].get_linewidth())

fig.savefig("../../figures/fig7/DN_GTExDsTau_Gene_Boxplot.pdf", dpi="figure", bbox_inches="tight")


# ## 7. make supplemental tables

# In[84]:


# supp table: DN classifications
supp_negregs = pairs[["gene_symbol", "Ensembl_gene_ID", "family", "ref_iso",
                      "alt_iso", "dbd_pct_lost", "m1h_cat", "y1h_cat", "y2h_cat", "loc_cat", "dn_cat", "dn_short"]]
supp_negregs.dn_short.value_counts()


# In[85]:


supp_negregs.columns = ["gene_symbol", "Ensembl_gene_ID", "family", "reference_isoform",
                        "alternative_isoform", "DBD_pct_lost_in_alt",
                        "M1H_activation_category", "PDI_category", "PPI_category", "localization_category",
                        "detailed_alt_iso_classification",
                        "alt_iso_classification"]
supp_negregs = supp_negregs[["gene_symbol", "Ensembl_gene_ID", "family", "reference_isoform",
                             "alternative_isoform", "DBD_pct_lost_in_alt",
                             "PDI_category", "PPI_category", "M1H_activation_category", "localization_category",
                             "alt_iso_classification", "detailed_alt_iso_classification"]]

supp_negregs["alt_iso_classification"].replace("DN", "negative regulator", inplace=True)
supp_negregs["detailed_alt_iso_classification"] = supp_negregs["detailed_alt_iso_classification"].str.replace("DN", "negative regulator")

supp_negregs["alt_iso_classification"].replace("rewire", "rewirer", inplace=True)
supp_negregs["detailed_alt_iso_classification"].replace("rewire", "rewirer", inplace=True)

supp_negregs["alt_iso_classification"].replace("similar", "similar to ref.", inplace=True)
supp_negregs["detailed_alt_iso_classification"].replace("similar", "similar to ref.", inplace=True)

supp_negregs["alt_iso_classification"].replace("likely", "likely non-functional", inplace=True)
supp_negregs["detailed_alt_iso_classification"].replace("likely", "likely non-functional", inplace=True)

supp_negregs[supp_negregs["alt_iso_classification"] == "negative regulator"].sample(5)


# In[86]:


supp_negregs.alt_iso_classification.value_counts()


# In[87]:


supp_negregs.to_csv("../../supp/SuppTable_NegRegs.txt", sep="\t", index=False)

