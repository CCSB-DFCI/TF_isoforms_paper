#!/usr/bin/env python
# coding: utf-8

# # Figure 7: Classifying TFs into Neg. Regs. vs. Rewirers

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

# import utils
sys.path.append("../")
sys.path.append("../data_loading")

import plotting
from plotting import PAPER_PRESET, PAPER_FONTSIZE, nice_boxplot, mimic_r_boxplot, annotate_pval
from isoform_pairwise_metrics import load_ref_vs_alt_isoforms_table

from data_loading import (load_annotated_TFiso1_collection,
                          load_developmental_tissue_expression_remapped,
                          load_gtex_remapped)

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
          "combination": sns.color_palette("Set2")[5]}


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


# In[11]:


df_gtex, metadata_gtex, genes_gtex = load_gtex_remapped()

exclusion_list_gtex = {'Cells - Leukemia cell line (CML)',
                       'Cells - EBV-transformed lymphocytes',
                       'Cells - Cultured fibroblasts'}

df_gtex = df_gtex.loc[:, ~df_gtex.columns.map(metadata_gtex['body_site']).isin(exclusion_list_gtex)]
metadata_gtex = metadata_gtex.loc[~metadata_gtex['body_site'].isin(exclusion_list_gtex), :]

means_gtex = df_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean()


# In[12]:


metadata_gtex_dummy = pd.read_table("../../data/processed/metadata_gtex_dummy.csv", sep=",", index_col=0)


# In[13]:


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

# In[14]:


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

# In[15]:


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

# In[16]:


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


# ## 3. categorize negative regulators
# include any whose DBD loss is >= 10%?

# In[17]:


def dn_cat(row):
    
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
                if row.dbd_pct_lost >= 10:
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
                if row.dbd_pct_lost >= 50:
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
                if row.dbd_pct_lost >= 10:
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
            if row.dbd_pct_lost >= 10:
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


# In[18]:


# double check that rewirers aren't all the same across axes


# In[19]:


to_plot = pairs[pairs["dn_cat"].isin(["rewire", "similar"])][["gene_symbol", "clone_acc_ref", "clone_acc_alt",
                                              "PPI_jaccard", "PDI_jaccard", "dbd_pct_lost",
                                              "activation_fold_change_log2", "dn_cat"]]

# make cols that are easier to visualize on plot
to_plot["fc_abs_activ"] = 2**-np.abs(to_plot["activation_fold_change_log2"])
to_plot["1m_dbd_pct"] = (100-to_plot["dbd_pct_lost"])/100


# In[20]:


isos = to_plot.copy()
isos = isos.sort_values(by="dn_cat")
isos = isos.reset_index().reset_index()

theoretical_iso = {"level_0": [isos.level_0.max()+1], "index": [0], "gene_symbol": ["theoretical"],
                   "clone_acc_ref": ["theoretical"], "clone_acc_alt": ["theoretical"],
                   "PPI_jaccard": [1.0], "PDI_jaccard": [1.0], "dbd_pct_lost": [0.0],
                   "activation_fold_change_log2": [np.log2(1)], "fc_abs_activ": [1.0],
                   "1m_dbd_pct": [1.0], "dn_cat": ["theoretical"]}

isos = isos.append(pd.DataFrame.from_dict(theoretical_iso))
isos


# In[21]:


columns_to_plot = ["level_0", "PPI_jaccard", "PDI_jaccard", "fc_abs_activ", "1m_dbd_pct"]
column_titles = ["", "PPI\n(jaccard)", "PDI\n(jaccard)",
                 "activ.\n(% of ref.)", "DBD\n(% of ref.)"]

# separate into similar v rewires
isos_sim = isos[isos["dn_cat"].isin(["similar", "theoretical"])]
isos_rw = isos[isos["dn_cat"].isin(["rewire"])]

colors = sns.color_palette("husl", n_colors=len(isos_rw))


# ### similar first

# In[22]:


# Creating the figure and axes (subplots) aligned in a single row without shared y-axis
gs_kw = dict(width_ratios=[0.75, 1, 1, 1, 1])
fig, axs = plt.subplots(1, 5, sharey=False, gridspec_kw=gs_kw)
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
    elif i == 4:
        axs[i].yaxis.tick_right()
        axs[i].yaxis.set_label_position("right")
        axs[i].set_ylabel("percent similarity")
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


# Data preparation for plotting
data_lines = []
for _, row in df.iterrows():
    data_lines.append([row[columns_to_plot[0]], row[columns_to_plot[1]], 
                       row[columns_to_plot[2]], row[columns_to_plot[3]],
                       row[columns_to_plot[4]]])

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
_ = axs[0].set_yticklabels(ticklabels)


# In[23]:


# Creating the figure and axes (subplots) aligned in a single row without shared y-axis
gs_kw = dict(width_ratios=[0.75, 1, 1, 1, 1])
fig, axs = plt.subplots(1, 5, sharey=False, gridspec_kw=gs_kw)
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
    elif i == 4:
        axs[i].yaxis.tick_right()
        axs[i].yaxis.set_label_position("right")
        axs[i].set_ylabel("percent similarity")
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


# Data preparation for plotting
data_lines = []
for _, row in df.iterrows():
    data_lines.append([row[columns_to_plot[0]], row[columns_to_plot[1]], 
                       row[columns_to_plot[2]], row[columns_to_plot[3]],
                       row[columns_to_plot[4]]])

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
_ = axs[0].set_yticklabels(ticklabels)


# ## 4. summary plots of DN categorization

# In[24]:


pairs["dn_short"] = pairs["dn_cat"].str.split(" ", expand=True)[0]
pairs.dn_short.value_counts()


# In[25]:


pairs[pairs["dn_short"] == "likely"]


# In[26]:


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
pairs[pairs["dn_short"] == "DN"].sample(5)


# In[27]:


fig = plt.figure(figsize=(1.5, 1.75))

ax = sns.countplot(data=pairs, x="dn_short", palette=sns.color_palette("Set2"),
                   order=["DN", "rewire", "similar", "NA"])
ax.set_xticklabels(["putative DN", "putative re-wirer", "similar", "NA"], ha="right", va="top", rotation=30)
ax.set_xlabel("")
ax.set_ylabel("count of alternative TF isoforms")

fig.savefig("../../figures/fig7/DN_countplot.pdf", dpi="figure", bbox_inches="tight")


# In[28]:


from upsetplot import plot


# In[29]:


ppis = list(set(list(pairs[pairs["dn_ppi"] == True]["clone_acc_alt"])))
pdis = list(set(list(pairs[pairs["dn_pdi"] == True]["clone_acc_alt"])))
activ = list(set(list(pairs[pairs["dn_activ"] == True]["clone_acc_alt"])))
dbd = list(set(list(pairs[pairs["dn_dbd"] == True]["clone_acc_alt"])))

contents = {"loss of PPIs": ppis, "loss of PDIs": pdis, "loss of activity": activ, "loss of DBD": dbd}
contents = upsetplot.from_contents(contents)

all_dn = set(ppis).union(set(pdis)).union(set(activ)).union(set(dbd))
print(len(all_dn))

fig = plt.figure(figsize=(3, 2))
d = plot(contents, fig=fig, sort_by="cardinality", show_counts=True, element_size=12, 
     intersection_plot_elements=4, totals_plot_elements=3)
d["intersections"].set_ylabel("# isoforms")
d["intersections"].grid(False)
d["totals"].grid(False)

fig.savefig("../../figures/fig7/DN_negreg_upset.pdf", dpi="figure", bbox_inches="tight")


# In[30]:


rw = pairs[pairs["dn_cat"] == "rewire"]
ppis = list(set(list(rw[rw["y2h_cat"] == "PPI rewire"]["clone_acc_alt"])))
pdis = list(set(list(rw[rw["y1h_cat"].str.contains("PDI")]["clone_acc_alt"])))
activ = list(set(list(rw[rw["m1h_cat"] != "NA"]["clone_acc_alt"])))

contents = {"change in PPIs": ppis, "change in PDIs": pdis, "change in activity": activ}
contents = upsetplot.from_contents(contents)

fig = plt.figure(figsize=(3, 2))
d = plot(contents, fig=fig, sort_by="cardinality", show_counts=True, element_size=12, 
         intersection_plot_elements=4, totals_plot_elements=3)
d["intersections"].set_ylabel("# isoforms")
d["intersections"].grid(False)
d["totals"].grid(False)

fig.savefig("../../figures/fig7/DN_rewire_upset.pdf", dpi="figure", bbox_inches="tight")


# In[31]:


ppis = list(set(list(pairs[pairs["y2h_cat"] != "NA"]["clone_acc_alt"])))
pdis = list(set(list(pairs[pairs["y1h_cat"] != "NA"]["clone_acc_alt"])))
activ = list(set(list(pairs[pairs["m1h_cat"] != "NA"]["clone_acc_alt"])))

contents = {"assessed PPIs": ppis, "assessed PDIs": pdis, "assessed activity": activ}
contents = upsetplot.from_contents(contents)

all_as = set(ppis).union(set(pdis)).union(set(activ))
print(len(all_as))

fig = plt.figure(figsize=(3, 2))
d = plot(contents, fig=fig, sort_by="cardinality", show_counts=True, element_size=12, 
         intersection_plot_elements=4, totals_plot_elements=3)
d["intersections"].set_ylabel("# isoforms")
d["intersections"].grid(False)
d["totals"].grid(False)

fig.savefig("../../figures/fig7/DN_pairs_assessed_upset.pdf", dpi="figure", bbox_inches="tight")


# In[32]:


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


# In[33]:


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


# In[34]:


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


# In[35]:


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


# In[36]:


# create df for stacked bar chart
delta_pdis = pairs[~pairs["y1h_cat"].isin(["NA", "no PDI change"])]
pdis_vc = pd.DataFrame(delta_pdis.dn_short.value_counts()).reset_index()

delta_ppis = pairs[~pairs["y2h_cat"].isin(["NA", "no PPI change (important PPIs)", "no PPI change (all PPIs)"])]
ppis_vc = pd.DataFrame(delta_ppis.dn_short.value_counts()).reset_index()

delta_activ = pairs[~pairs["m1h_cat"].isin(["NA", "similar"])]
activ_vc = pd.DataFrame(delta_activ.dn_short.value_counts()).reset_index()

mrg = pdis_vc.merge(ppis_vc, on="index", how="outer").merge(activ_vc, on="index", how="outer")
mrg.fillna(0, inplace=True)
mrg.columns = ["index", "PDIs", "PPIs", "activity"]

to_plot = pd.melt(mrg, id_vars="index")
to_plot.sample(5)


# In[37]:


# make stacked barchart situation
tmp = pairs[pairs["dn_short"] == "DN"]
dn_pdi_change = len(tmp[tmp["dn_pdi"] == True])
dn_ppi_change = len(tmp[tmp["dn_ppi"] == True])
dn_activ_change = len(tmp[tmp["dn_activ"] == True])
dn_dbd_change = len(tmp[tmp["dn_dbd"] == True])
tot_dn = dn_pdi_change + dn_ppi_change + dn_activ_change + dn_dbd_change

tmp = pairs[pairs["dn_short"] == "rewire"]
rw_pdi_change = len(tmp[~tmp["y1h_cat"].isin(["NA", "no PDI change"])])
rw_ppi_change = len(tmp[~tmp["y2h_cat"].isin(["NA", "no PPI change (important PPIs)"])])
rw_activ_change = len(tmp[tmp["m1h_cat"] != "NA"])
rw_dbd_change = len(tmp[tmp["dbd_pct_lost"] > 0])
tot_rw = rw_pdi_change + rw_ppi_change + rw_activ_change + rw_dbd_change

df = pd.DataFrame.from_dict({"DN": {"pdi_change": dn_pdi_change/tot_dn*100, 
                                    "ppi_change": dn_ppi_change/tot_dn*100,
                                    "activ_change": dn_activ_change/tot_dn*100, 
                                    "dbd_change": dn_dbd_change/tot_dn*100},
                             "rewire": {"pdi_change": rw_pdi_change/tot_rw*100,
                                        "ppi_change": rw_ppi_change/tot_rw*100,
                                        "activ_change": rw_activ_change/tot_rw*100,
                                        "dbd_change": rw_dbd_change/tot_rw*100}})
df["DN_cumsum"] = np.cumsum(df["DN"])
df["rw_cumsum"] = np.cumsum(df["rewire"])
df


# In[38]:


colors = met_brewer.met_brew(name="Hokusai3", n=4, brew_type="discrete")
sns.palplot(colors)


# In[39]:


fig, ax = plt.subplots(figsize=(0.85, 1.5))

xs = ["negative regulator", "rewirer"]
y1 = list(df[["DN", "rewire"]].loc["pdi_change"])
y2 = list(df[["DN", "rewire"]].loc["ppi_change"])
b2 = np.add(y1, y2)
y3 = list(df[["DN", "rewire"]].loc["activ_change"])
b3 = np.add(b2, y3)
y4 = list(df[["DN", "rewire"]].loc["dbd_change"])

ax.bar(xs, y1, color=colors[0], label="∆ PDIs")
ax.bar(xs, y2, bottom=y1, color=colors[1], label="∆ PPIs")
ax.bar(xs, y3, bottom=b2, color=colors[2], label="∆ activity")
ax.bar(xs, y4, bottom=b3, color=colors[3], label="∆ DBD")

# add legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)

ax.set_ylabel("percent prevalence")
ax.set_xticklabels(["negative\nregulator", "rewirer"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/dn_stacked_bar.pdf", dpi="figure", bbox_inches="tight")


# In[40]:


# make stacked barchart situation of all assays (to compare)
tmp = pairs[pairs["dn_short"] != "NA"]
pdi_assessed = len(tmp[tmp["y1h_cat"] != "NA"])
ppi_assessed = len(tmp[tmp["y2h_cat"] != "NA"])
activ_assessed = len(tmp[tmp["m1h_cat"] != "NA"])

tot_assessed = pdi_assessed + ppi_assessed + activ_assessed

nc = pd.DataFrame.from_dict({"assessed": {"pdi": pdi_assessed/tot_assessed*100, 
                                    "ppi": ppi_assessed/tot_assessed*100,
                                    "activ": activ_assessed/tot_assessed*100}})
nc["assessed_cumsum"] = np.cumsum(nc["assessed"])
nc


# In[41]:


fig, ax = plt.subplots(figsize=(0.25, 1.5))

xs = ["assessed"]
y1 = list(nc[["assessed"]].loc["pdi"])
y2 = list(nc[["assessed"]].loc["ppi"])
b2 = np.add(y1, y2)
y3 = list(nc[["assessed"]].loc["activ"])

ax.bar(xs, y1, color=colors[0], label="PDIs")
ax.bar(xs, y2, bottom=y1, color=colors[1], label="PPIs")
ax.bar(xs, y3, bottom=b2, color=colors[2], label="activity")

# add legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)

ax.set_ylabel("percent prevalence")
ax.set_xticklabels(["total assessed"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/dn_stacked_bar.nc.pdf", dpi="figure", bbox_inches="tight")


# In[42]:


genes_w_dn = pairs[pairs["dn_short"] == "DN"][["gene_symbol", "family"]].drop_duplicates()
genes_w_rw = pairs[pairs["dn_short"] == "rewire"][["gene_symbol", "family"]].drop_duplicates()
genes_w_sim = pairs[pairs["dn_short"] == "similar"][["gene_symbol", "family"]].drop_duplicates()
tot_genes = pairs[["gene_symbol", "family"]].drop_duplicates()

tot_genes_per_f = tot_genes.groupby("family")["gene_symbol"].agg("count").reset_index()
dn_genes_per_f = genes_w_dn.groupby("family")["gene_symbol"].agg("count").reset_index()
rw_genes_per_f = genes_w_rw.groupby("family")["gene_symbol"].agg("count").reset_index()
sim_genes_per_f = genes_w_sim.groupby("family")["gene_symbol"].agg("count").reset_index()

family_cats = tot_genes_per_f.merge(dn_genes_per_f, 
                                    on="family", how="left").merge(rw_genes_per_f, 
                                                                   on="family", how="left").merge(sim_genes_per_f,
                                                                                                  on="family",
                                                                                                  how="left")
family_cats.fillna(0, inplace=True)
family_cats.columns = ["family", "tot", "dn", "rw", "sim"]

family_cats["tot_p"] = family_cats["tot"]/family_cats["tot"].sum(axis=0)*100
family_cats["dn_p"] = family_cats["dn"]/family_cats["dn"].sum(axis=0)*100
family_cats["rw_p"] = family_cats["rw"]/family_cats["rw"].sum(axis=0)*100
family_cats["sim_p"] = family_cats["sim"]/family_cats["sim"].sum(axis=0)*100
family_cats.sort_values(by="tot", ascending=False).head(11)


# In[43]:


colors = met_brewer.met_brew(name="Renoir", n=11, brew_type="discrete")
sns.palplot(colors)


# In[44]:


fig, ax = plt.subplots(figsize=(3, 1.75))

xs = ["total", "negative regulator", "rewirer", "similar"]

b = np.zeros(4)
c = 0
for i, row in family_cats.sort_values(by="tot", ascending=False).head(11).iterrows():
    y = list(row[["tot_p", "dn_p", "rw_p", "sim_p"]])
    ax.bar(xs, y, bottom=b, label=row.family, color=colors[c])
    b = np.add(b, y)
    c += 1

ax.bar(xs, np.subtract([100, 100, 100, 100], b), bottom=b, label="smaller families", color="lightgray")
ax.set_ylabel("% of categorized isoforms")

# add legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/dn_families_stacked_bar.pdf", dpi="figure", bbox_inches="tight")


# ## 5. join data with Joung et al. o/ex

# In[45]:


tfs = load_annotated_TFiso1_collection()


# In[46]:


def pd_translate(row):
    s = Seq(row["ORF sequence"])
    aa = s.translate()
    return str(aa)

joung_orf["seq_aa"] = joung_orf.apply(pd_translate, axis=1)


# In[47]:


tf_id_map = {}
for tf in tfs:
    isos = tfs[tf]
    for i, iso in enumerate(isos.isoforms):
        sub_dict = {}
        try: 
            iso_clone_acc = iso.clone_acc
        except AttributeError:
            continue

        try:
            iso_seq_aa = iso.aa_seq_GENCODE
        except AttributeError:
            iso_seq_aa = iso.aa_seq
        iso_ensts = iso.ensembl_transcript_ids

        
        # first try to match based on aa seq
        joung_sub = joung_orf[joung_orf["seq_aa"] == iso_seq_aa]
        
        if len(joung_sub) > 0:
            sub_dict["match_type"] = "seq_aa"
            sub_dict["joung_id"] = joung_sub["Name"].iloc[0]
        
        # if not found, then try ensts
        if len(joung_sub) == 0:
            if iso_ensts is None:
                continue
            
            for iso_enst in iso_ensts:
                joung_sub = joung_orf[joung_orf["RefSeq and Gencode ID"].str.contains(iso_enst)]
                if len(joung_sub) > 0:
                    continue
            
            if len(joung_sub) > 0:
                sub_dict["match_type"] = "enst"
                sub_dict["joung_id"] = joung_sub["Name"].iloc[0]
            else:
                continue
        
        sub_dict["enst"] = iso_ensts
        sub_dict["seq_aa"] = iso_seq_aa
        tf_id_map[iso_clone_acc] = sub_dict


# In[48]:


tf_id_map_df = pd.DataFrame.from_dict(tf_id_map, orient="index").reset_index()
print(len(tf_id_map_df))
tf_id_map_df.sample(5)


# In[49]:


joung_orf = joung_orf.merge(tf_id_map_df, left_on="Name", right_on="joung_id", how="left", suffixes=("_joung",
                                                                                                     "_tf1p0"))
joung_orf.sample(5)


# In[50]:


joung_data = joung_orf.merge(joung_data, on="Name", how="left")


# In[51]:


dn_ref = pairs[["gene_symbol", "family", "clone_acc_ref", "is_ref_novel_isoform", "is_MANE_select_isoform_cloned",
             "dn_short"]].drop_duplicates()
dn_ref.columns = ["gene_name", "family", "tf1p0_id", "is_novel", "is_MANE_select", "dn_cat"]
dn_ref["dn_cat"] = "ref"
dn_ref["iso_status"] = "ref"


# In[52]:


dn_alt = pairs[["gene_symbol", "family", "clone_acc_alt", "is_alt_novel_isoform", "is_MANE_select_isoform_cloned",
             "dn_short"]].drop_duplicates()
dn_alt.columns = ["gene_name", "family", "tf1p0_id", "is_novel", "is_MANE_select", "dn_cat"]
dn_alt["is_MANE_select"] = False # assuming none of the alts are the MANE select
dn_alt["iso_status"] = "alt"


# In[53]:


dn_cats = dn_ref.append(dn_alt).drop_duplicates()


# In[54]:


dn_cats = dn_cats.merge(joung_data, left_on="tf1p0_id", right_on="index", how="left")


# In[55]:


dn_cats.iso_status.value_counts()


# In[56]:


dn_cats[~pd.isnull(dn_cats["Name"])].iso_status.value_counts()


# In[57]:


refs_inc = len(dn_cats[(~pd.isnull(dn_cats["Name"])) & (dn_cats["iso_status"] == "ref")])
refs_tf1p0 = len(dn_cats[dn_cats["iso_status"] == "ref"])
print("%% of our ref seqs included in joung: %s" % (refs_inc/refs_tf1p0*100))


# In[58]:


alts_inc = len(dn_cats[(~pd.isnull(dn_cats["Name"])) & (dn_cats["iso_status"] == "alt")])
alts_tf1p0 = len(dn_cats[dn_cats["iso_status"] == "alt"])
print("%% of our alt seqs included in joung: %s" % (alts_inc/alts_tf1p0*100))


# In[59]:


dn_cats["orf_len"] = dn_cats["seq_aa_joung"].str.len()


# In[60]:


joung_down_tf1p0_map = joung_down_map.merge(dn_cats[["TF ORF", "tf1p0_id", "iso_status", "dn_cat", "orf_len"]],
                                           left_on="TF", right_on="TF ORF")
print(len(joung_down_tf1p0_map))
print(len(joung_down_tf1p0_map["TF ORF"].unique()))


# In[61]:


joung_down_tf1p0_map.fillna("NA", inplace=True)


# In[62]:


joung_tf1p0_cnts = joung_down_tf1p0_map.groupby(["TF", "tf1p0_id", "iso_status", 
                                                 "dn_cat", "orf_len"])["TF ORF"].agg("count").reset_index()
joung_tf1p0_cnts.columns = ["TF", "tf1p0_id", "iso_status", "dn_cat", "orf_len", "tot_cell_cnt"]


# In[63]:


dn_cats_nonan = dn_cats[~pd.isnull(dn_cats["Diffusion P-value"])]
len(dn_cats_nonan)


# In[64]:


dn_cats_nonan["neglog_diff_pval"] = -np.log10(dn_cats_nonan["Diffusion P-value"])
dn_cats_nonan.fillna("NA", inplace=True)


# In[65]:


dn_cats_nonan_ref = dn_cats_nonan[dn_cats_nonan["iso_status"] == "ref"]
dn_cats_nonan_alt = dn_cats_nonan[dn_cats_nonan["iso_status"] == "alt"]
dn_cats_nonan_diff = dn_cats_nonan_ref.merge(dn_cats_nonan_alt, on=["gene_name", "family", "RefSeq Gene Name"],
                                             how="left", suffixes=("_ref", "_alt"))
dn_cats_nonan_diff["diff_pval_diff"] = dn_cats_nonan_diff["Diffusion P-value_ref"] - dn_cats_nonan_diff["Diffusion P-value_alt"]
dn_cats_nonan_diff["diff_diff_diff"] = dn_cats_nonan_diff["Diffusion difference_ref"] - dn_cats_nonan_diff["Diffusion difference_alt"]

dn_cats_nonan_diff["abs_ddd"] = np.abs(dn_cats_nonan_diff["diff_diff_diff"])


# In[66]:


dn_cats_nonan[dn_cats_nonan["gene_name"] == "CREB1"]


# In[67]:


dn_cats_nonan[dn_cats_nonan["dn_cat"] != "NA"].sort_values(by="neglog_diff_pval", ascending=False).head(10)


# In[68]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(data=dn_cats_nonan[dn_cats_nonan["dn_cat"].isin(["ref"])], 
                x="Diffusion difference", y="neglog_diff_pval", 
                color="white", linewidth=0.5, edgecolor="black", alpha=0.8, zorder=10,
                **{"s": 9})

sns.scatterplot(data=dn_cats_nonan[dn_cats_nonan["dn_cat"].isin(["similar", "rewire", "DN"])], 
                x="Diffusion difference", y="neglog_diff_pval", 
                hue="dn_cat", palette=dn_pal, linewidth=0.25, edgecolor="black", alpha=0.8, zorder=10,
                **{"s": 9}, ax=ax)

for annot_clone, ha, va, offset, relpos, cs in zip(["GRHL3|3/7|08G09", "HMBOX1|3/5|03E06", 
                                                    "KLF7|3/8|10B10", "GRHL3|4/7|08F09", "PBX1|2/2|02C05",
                                                    "DLX1|2/2|07E09"],
                                                    ["center", "right", "right", "right", "right", "right"],
                                                    ["top", "bottom", "top", "center", "bottom", "bottom"],
                                                    [(0, -7), (3, 8), (-8, 7), (-7, 0), (-20, -15), (-8, -8)],
                                                    [(0.5, 1), (1, 0.5), (1, 0.5), (1, 0.5), (1, 0.5), (1, 0.5)],
                                                    ["arc3,rad=0", "arc3,rad=-0.3", "arc3,rad=-0.3", "arc3,rad=0",
                                                     "arc3,rad=0.3", "arc3,rad=-0.3"]):
    row = dn_cats_nonan[dn_cats_nonan["tf1p0_id"] == annot_clone].iloc[0]
    if row["dn_cat"] == "ref":
        color = "black"
    else:
        color = dn_pal[row["dn_cat"]]
    print("annot clone: %s | ha: %s | va: %s" % (annot_clone, ha, va))
    
    shorter_id = annot_clone.split("|")[0] + "-" + annot_clone.split("|")[1].split("/")[0]
    ax.annotate(shorter_id, xy=(row["Diffusion difference"], row["neglog_diff_pval"]), xytext=offset,
                color=color, ha=ha, va=va,
                textcoords='offset points', bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'),
                arrowprops=dict(arrowstyle="-", color=color, relpos=relpos, connectionstyle=cs))

ax.set_xlabel("over-expression effect size")
ax.set_ylabel("-log10(over-expression p-value)")
ax.set_title("effect of TF isoforms on differentiation\n(Joung et al.)")

ax.set_xlim((-0.04, 0.01))
ax.set_ylim((-0.01, 7.2))
ax.axhline(y=-np.log10(0.05), linestyle="dashed", color="black", linewidth=0.5)
ax.axvline(x=0, linestyle="dashed", color="black", linewidth=0.5)

ax.get_legend().remove()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/Joung_Volcano.pdf", dpi="figure", bbox_inches="tight")


# In[69]:


joung_cells["Name"] = joung_cells["TF"].str.strip().str.split("-", expand=True)[0]
print(len(joung_cells))

# # filter out anything with score < 0.2
joung_cells = joung_cells[joung_cells["prediction.score.max"] > 0.2]
print(len(joung_cells))


# In[70]:


joung_cells_grp = joung_cells.groupby(["Name", "TF", "predicted.id"])["batch"].agg("count").reset_index()


# In[71]:


joung_tf1p0_cnts["cell_cnt_qcut"] = pd.qcut(joung_tf1p0_cnts["tot_cell_cnt"], q=4, labels=[1, 2, 3, 4])
dn_cats_nonan = dn_cats_nonan.merge(joung_tf1p0_cnts[["TF", "tot_cell_cnt", "cell_cnt_qcut"]], 
                                    left_on="TF ORF", right_on="TF")


# In[72]:


tot_cell_cnt = dn_cats_nonan[["Name", "TF", "tot_cell_cnt"]].drop_duplicates()
diff_cell_cnt = joung_cells.groupby(["Name", "TF"])["batch"].agg("count").reset_index()
cell_cnt = tot_cell_cnt.merge(diff_cell_cnt, on=["Name", "TF"], how="left")
cell_cnt.fillna(0, inplace=True)
cell_cnt.columns = ["Name", "TF", "tot_cell_cnt", "diff_cell_cnt"]

orf_enr = cell_cnt.merge(joung_cells_grp, on=["Name", "TF"], how="left")
orf_enr["batch"].fillna(0, inplace=True)

orf_enr.columns = ["Name", "TF", "tot_cell_cnt", "diff_cell_cnt", "predicted.id", "id_cell_cnt"]
orf_enr["perc_cells_of_diff_tf"] = orf_enr["id_cell_cnt"]/orf_enr["diff_cell_cnt"]
orf_enr["perc_cells_of_tot_tf"] = orf_enr["id_cell_cnt"]/orf_enr["tot_cell_cnt"]


# In[73]:


orf_enr_dn = orf_enr.merge(dn_cats_nonan[["gene_name", "Name", "tf1p0_id",
                                          "dn_cat"]], on="Name").drop_duplicates(subset=["tf1p0_id",
                                                                                         "predicted.id",
                                                                                         "dn_cat"])


# In[74]:


has_alt = list(orf_enr_dn[orf_enr_dn["dn_cat"] != "ref"]["gene_name"].unique())
orf_enr_dn_filt = orf_enr_dn[orf_enr_dn["gene_name"].isin(has_alt)]

has_ref = list(orf_enr_dn_filt[orf_enr_dn_filt["dn_cat"] == "ref"]["gene_name"].unique())
orf_enr_dn_filt = orf_enr_dn_filt[orf_enr_dn_filt["gene_name"].isin(has_ref)]
len(orf_enr_dn_filt)


# In[75]:


orf_enr_dn_filt["dn_cat_s"] = pd.Categorical(orf_enr_dn_filt["dn_cat"], ["ref", "similar", "rewire", 
                                                                         "DN", "NA", "likely"])
orf_enr_dn_filt = orf_enr_dn_filt.sort_values(by=["gene_name", "dn_cat_s"])
orf_enr_dn_filt[orf_enr_dn_filt["gene_name"] == "PBX1"]


# In[76]:


cell_cnt["undiff_cell_cnt"] = cell_cnt["tot_cell_cnt"] - cell_cnt["diff_cell_cnt"]
len(cell_cnt)


# In[77]:


cell_cnt["diff_cell_perc"] = (cell_cnt["diff_cell_cnt"]/cell_cnt["tot_cell_cnt"])*100
cell_cnt["undiff_cell_perc"] = (cell_cnt["undiff_cell_cnt"]/cell_cnt["tot_cell_cnt"])*100
cell_cnt[cell_cnt["TF"].str.contains("GRHL3")]


# In[78]:


tmp = orf_enr_dn_filt[orf_enr_dn_filt["tf1p0_id"].str.contains("GRHL3")].pivot(index="tf1p0_id", 
                                                                               columns="predicted.id", 
                                                                               values="perc_cells_of_tot_tf")
tmp.drop(np.nan, axis=1, inplace=True)
tmp = tmp.loc[["GRHL3|3/7|08G09", "GRHL3|4/7|08F09", "GRHL3|1/7|08E10", "GRHL3|2/7|08A10"]]
tmp.fillna(0, inplace=True)

idx = pd.DataFrame(tmp.index)
idx = idx.merge(orf_enr_dn_filt[["tf1p0_id", "dn_cat"]], on="tf1p0_id").drop_duplicates().set_index("tf1p0_id")
idx["Isoform type"] = idx.dn_cat.map(dn_pal)

g = sns.clustermap(tmp, cmap="Greys", row_cluster=False, row_colors=idx["Isoform type"],
                   figsize=(4, 1), yticklabels=True, cbar_pos=(0, 1, 0.05, 0.2), linewidth=0.5, linecolor="grey")

g.savefig("../../figures/fig7/Joung_GRHL3_hm.pdf", bbox_inches="tight", dpi="figure")


# ## 6. plot expression profiles of isoform categories

# In[79]:


# use same downsample as prev figs
means_gtex_downsample = df_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean()


# In[80]:


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


# In[81]:


# calculate expression ratios - gtex
per_gene_gtex = ((2 ** df_gtex - 1)
                .groupby(genes_gtex)
                .transform('sum'))
f_gtex = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1)
        .mean())
f_gtex = f_gtex * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex = f_gtex * 100


# In[82]:


# calculate expression ratios -gtex downsampled
per_gene_gtex_ds = ((2 ** df_gtex.loc[:,metadata_gtex_dummy.index] - 1)
                   .groupby(genes_gtex)
                   .transform('sum'))

f_gtex_downsample = (((2 ** df_gtex - 1) / per_gene_gtex)
        .groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1)
        .mean())
f_gtex_downsample = f_gtex_downsample * (per_gene_gtex.groupby(df_gtex.columns.map(metadata_gtex_dummy['body_site']), axis=1).mean() >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])  # only count fractions if gene TPM is >= 1

f_gtex_downsample = f_gtex_downsample * 100


# In[83]:


# calculate gene-level tissue specificities
gene_dev_nonan_taus, gene_dev_nan_taus, gene_dev_array_max = calculate_tau(per_gene_dev.drop_duplicates())
gene_gtex_nonan_taus, gene_gtex_nan_taus, gene_gtex_array_max = calculate_tau(per_gene_gtex.drop_duplicates())
gene_gtex_ds_nonan_taus, gene_gtex_ds_nan_taus, gene_gtex_ds_array_max = calculate_tau(per_gene_gtex_ds.drop_duplicates())


# In[84]:


gene_taus = pd.DataFrame()
gene_taus["UID"] = per_gene_dev.drop_duplicates().index
gene_taus["dev_tau"] = gene_dev_nan_taus
gene_taus["gtex_tau"] = gene_gtex_nan_taus
gene_taus["gtex_ds_tau"] = gene_gtex_ds_nan_taus
gene_taus["gene_name"] = gene_taus["UID"].str.split("|", expand=True)[0]
gene_taus.sample(5)


# In[85]:


# join w pairs table
dev_ratios = f_dev.reset_index()
dev_ratios["clone_acc"] = dev_ratios["UID"].str.split(" ", expand=True)[0].astype(str)
dev_ratios = dev_ratios[dev_ratios["clone_acc"] != "noclone"]
len(dev_ratios)


# In[86]:


gtex_ratios = f_gtex.reset_index()
gtex_ratios["clone_acc"] = gtex_ratios["UID"].str.split(" ", expand=True)[0].astype(str)
gtex_ratios = gtex_ratios[gtex_ratios["clone_acc"] != "noclone"]
len(gtex_ratios)


# In[87]:


gtex_ds_ratios = f_gtex_downsample.reset_index()
gtex_ds_ratios["clone_acc"] = gtex_ds_ratios["UID"].str.split(" ", expand=True)[0].astype(str)
gtex_ds_ratios = gtex_ds_ratios[gtex_ds_ratios["clone_acc"] != "noclone"]
len(gtex_ds_ratios)


# In[88]:


dev_ratios = dev_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id")
gtex_ratios = gtex_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id")
gtex_ds_ratios = gtex_ds_ratios.merge(dn_cats, left_on="clone_acc", right_on="tf1p0_id")
print(len(dev_ratios))
print(len(gtex_ratios))
print(len(gtex_ds_ratios))


# In[89]:


dn_cats = dn_cats.merge(gene_taus, on="gene_name")
print(len(dn_cats))
dn_cats.head()


# In[90]:


ref_expr = dn_cats.groupby(["gene_name", "family", "dn_cat", "dev_tau",
                            "gtex_tau", "gtex_ds_tau"])["tf1p0_id"].agg("count").reset_index()
ref_expr = ref_expr.pivot(index="gene_name",
                          columns="dn_cat", values="tf1p0_id")
ref_expr.fillna(0, inplace=True)


# In[91]:


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


# In[92]:


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
ax.set_ylabel("gene-level tissue specificity (tau)", position=(0, 0.38))

# manually set left axis so it stops at 1.0
ax.set_ylim((0.45, 1.2))
ax.spines['left'].set_visible(False)
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
axes_to_data = ax.transAxes + ax.transData.inverted()
left_spine_in_data_coords = axes_to_data.transform((0, 0))
ax.plot([left_spine_in_data_coords[0], left_spine_in_data_coords[0]], [0.45, 1],
         color=ax.spines['bottom'].get_edgecolor(), linewidth=ax.spines['bottom'].get_linewidth())


fig.savefig("../../figures/fig7/DN_DevTau_Gene_Boxplot.pdf", dpi="figure", bbox_inches="tight")


# In[93]:


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
ax.set_ylabel("gene-level tissue specificity (tau)", position=(0, 0.38))

# manually set left axis so it stops at 1.0
ax.set_ylim((0.5, 1.2))
ax.spines['left'].set_visible(False)
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
axes_to_data = ax.transAxes + ax.transData.inverted()
left_spine_in_data_coords = axes_to_data.transform((0, 0))
ax.plot([left_spine_in_data_coords[0], left_spine_in_data_coords[0]], [0.5, 1],
         color=ax.spines['bottom'].get_edgecolor(), linewidth=ax.spines['bottom'].get_linewidth())

fig.savefig("../../figures/fig7/DN_GTExDsTau_Gene_Boxplot.pdf", dpi="figure", bbox_inches="tight")


# ## 7. load BRCA data

# In[94]:


brca_cnts_f = "../../data/processed/Nathans_analysis/Breast_cancer/isoCounts.BreastCancer.txt"
brca_tx_f = "../../data/processed/Nathans_analysis/Breast_cancer/transcript.BreastCancer.txt"
pam50_f = "../../data/processed/Nathans_analysis/Breast_cancer/groups.PAM50.txt"
bulk_f = "../../data/processed/Nathans_analysis/Breast_cancer/groups.BreastCancer_ratios.txt"


# In[95]:


skiprows=list(range(96320, 96387))+list(range(99680,99687))
brca = pd.read_table(brca_tx_f, sep="\t", skiprows=skiprows)
brca.shape


# In[96]:


pam50_samps = pd.read_table(pam50_f, header=None)
pam50_samps.columns = ["file", "samp_type_id", "samp_type"]
pam50_samps["tcga_id"] = pam50_samps["file"].str.split(".", expand=True)[0]
pam50_samps.samp_type.value_counts()


# In[97]:


bulk_samps = pd.read_table(bulk_f, header=None)
bulk_samps.columns = ["file", "samp_type_id", "samp_type"]
bulk_samps["tcga_id"] = bulk_samps["file"].str.split(".", expand=True)[0]
bulk_samps.samp_type.value_counts()


# In[98]:


# map brca sample types
brca_samps = list(bulk_samps[bulk_samps["samp_type"] != "controls"]["tcga_id"])
print("# breast cancer samples: %s" % len(brca_samps))


# same ctrls in both
ctrl_samps = list(pam50_samps[pam50_samps["samp_type"] == "controls"]["tcga_id"])
print("# control samples: %s" % len(ctrl_samps))

luma_samps = list(pam50_samps[pam50_samps["samp_type"] == "Luminal A"]["tcga_id"])
print("# Luminal A samples: %s" % len(luma_samps))

lumb_samps = list(pam50_samps[pam50_samps["samp_type"] == "Luminal B"]["tcga_id"])
print("# Luminal B samples: %s" % len(lumb_samps))

tn_samps = list(pam50_samps[pam50_samps["samp_type"] == "Basal-like"]["tcga_id"])
print("# Basal-like samples: %s" % len(tn_samps))

her2_samps = list(pam50_samps[pam50_samps["samp_type"] == "HER2-enriched"]["tcga_id"])
print("# HER2-enriched samples: %s" % len(her2_samps))

norm_samps = list(pam50_samps[pam50_samps["samp_type"] == "Normal-like"]["tcga_id"])
print("# Normal-like samples: %s" % len(norm_samps))


# In[99]:


# one brca samp is weirdly missing, remove
brca_samps = [x for x in brca_samps if x in brca.columns]
len(brca_samps)


# In[100]:


## patient is the 3rd value in the barcode
## source: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
bulk_samps["patient_id"] = bulk_samps["tcga_id"].str.split("-", expand=True)[2]
pam50_samps["patient_id"] = pam50_samps["tcga_id"].str.split("-", expand=True)[2]


# In[101]:


tcga_samps = bulk_samps.merge(pam50_samps, on=["tcga_id", "patient_id"], how="outer",
                              suffixes=("_brca", "_pam50"))
print(len(tcga_samps))


# In[102]:


tcga_ctrls = tcga_samps[(tcga_samps["samp_type_brca"] == "controls") | (tcga_samps["samp_type_pam50"] == "controls")]
len(tcga_ctrls)


# In[103]:


tcga_tumors = tcga_samps[(tcga_samps["samp_type_brca"] != "controls") | (tcga_samps["samp_type_pam50"] != "controls")]
len(tcga_tumors)


# In[104]:


tcga_paired = tcga_ctrls.merge(tcga_tumors, on=["patient_id"], suffixes=("_ctrl", "_tumor"))
print(len(tcga_paired))


# ## 8. aggregate TF iso expression across transcripts + calculate isoform ratios/med expr

# In[105]:


tf_id_map = pd.DataFrame()
gene_names = []
iso_ids = []
clone_accs = []
enst_ids = []

for tf in tfs:
    db = tfs[tf]
    for iso in db.isoforms:
        enst_id = iso.ensembl_transcript_ids
        try:
            clone_acc = iso.clone_acc
        except:
            clone_acc = "none"
        
        if enst_id is None:
            gene_names.append(tf)
            iso_ids.append(iso.name)
            clone_accs.append(clone_acc)
            enst_ids.append("none")
        else:
            for enst_id in iso.ensembl_transcript_ids:
                gene_names.append(tf)
                iso_ids.append(iso.name)
                clone_accs.append(clone_acc)
                enst_ids.append(enst_id)

tf_id_map["gene_name"] = gene_names
tf_id_map["iso_id"] = iso_ids
tf_id_map["clone_acc"] = clone_accs
tf_id_map["enst_id"] = enst_ids
print(len(tf_id_map))
tf_id_map.sample(5)


# In[106]:


def merge_id(row):
    if row.enst_id == "none":
        return row.clone_acc
    else:
        return row.enst_id
    
tf_id_map["merge_id"] = tf_id_map.apply(merge_id, axis=1)


# In[107]:


dd = tf_id_map[["iso_id", "gene_name"]].drop_duplicates()
print(len(dd))
gene_dict = {row.iso_id : row.gene_name for i, row in dd.iterrows()}


# In[108]:


brca_cols = [x for x in brca.columns if x != "UID"]
len(brca_cols)


# In[109]:


brca = brca.merge(tf_id_map, left_on="UID", right_on="merge_id")
len(brca)


# In[110]:


brca_isos = brca.groupby("iso_id")[brca_cols].agg("sum").reset_index()
len(brca_isos)


# In[111]:


# calculate isoform ratios, set anything w gene-level exp <= 1 to nan
brca_genes = pd.Series(index=brca_isos.iso_id, data=brca_isos.iso_id.map(gene_dict).values)

brca_idx = brca_isos.set_index("iso_id", inplace=False)
brca_idx = brca_idx[brca_cols]
brca_gene_sum = brca_idx.groupby(brca_genes).transform('sum')

f_brca = brca_idx/brca_gene_sum
f_brca_nan = f_brca * (brca_gene_sum >= 1).applymap(lambda x: {False: np.nan, True: 1}[x])


# In[112]:


tcga_paired_ctrls = list(tcga_paired["tcga_id_ctrl"].unique())
tcga_paired_tumors = list(tcga_paired["tcga_id_tumor"].unique())


# In[113]:


brca_isos["med_brca_tpm"] = brca_isos[brca_samps].median(axis=1)
brca_isos["med_ctrl_tpm"] = brca_isos[ctrl_samps].median(axis=1)
brca_isos["med_luma_tpm"] = brca_isos[luma_samps].median(axis=1)
brca_isos["med_lumb_tpm"] = brca_isos[lumb_samps].median(axis=1)
brca_isos["med_tn_tpm"] = brca_isos[tn_samps].median(axis=1)
brca_isos["med_her2_tpm"] = brca_isos[her2_samps].median(axis=1)
brca_isos["med_norm_tpm"] = brca_isos[norm_samps].median(axis=1)
brca_isos["med_paired-brca_tpm"] = brca_isos[tcga_paired_tumors].median(axis=1)
brca_isos["med_paired-ctrls_tpm"] = brca_isos[tcga_paired_ctrls].median(axis=1)


# In[114]:


f_brca_nan["med_brca_rationan"] = f_brca_nan[brca_samps].median(axis=1)
f_brca_nan["med_ctrl_rationan"] = f_brca_nan[ctrl_samps].median(axis=1)
f_brca_nan["med_luma_rationan"] = f_brca_nan[luma_samps].median(axis=1)
f_brca_nan["med_lumb_rationan"] = f_brca_nan[lumb_samps].median(axis=1)
f_brca_nan["med_tn_rationan"] = f_brca_nan[tn_samps].median(axis=1)
f_brca_nan["med_her2_rationan"] = f_brca_nan[her2_samps].median(axis=1)
f_brca_nan["med_norm_rationan"] = f_brca_nan[norm_samps].median(axis=1)
f_brca_nan["med_paired-brca_rationan"] = f_brca_nan[tcga_paired_tumors].median(axis=1)
f_brca_nan["med_paired-ctrls_rationan"] = f_brca_nan[tcga_paired_ctrls].median(axis=1)


# In[115]:


f_brca["med_brca_ratio"] = f_brca[brca_samps].median(axis=1)
f_brca["med_ctrl_ratio"] = f_brca[ctrl_samps].median(axis=1)
f_brca["med_luma_ratio"] = f_brca[luma_samps].median(axis=1)
f_brca["med_lumb_ratio"] = f_brca[lumb_samps].median(axis=1)
f_brca["med_tn_ratio"] = f_brca[tn_samps].median(axis=1)
f_brca["med_her2_ratio"] = f_brca[her2_samps].median(axis=1)
f_brca["med_norm_ratio"] = f_brca[norm_samps].median(axis=1)
f_brca["med_paired-brca_ratio"] = f_brca[tcga_paired_tumors].median(axis=1)
f_brca["med_paired-ctrls_ratio"] = f_brca[tcga_paired_ctrls].median(axis=1)


# ## 9. calculate expr/ratio change across paired samples

# In[116]:


paired_ctrl_samps = list(tcga_paired["tcga_id_ctrl"])
print(len(paired_ctrl_samps))
paired_tumor_samps = list(tcga_paired["tcga_id_tumor"])
print(len(paired_tumor_samps))


# In[117]:


## calculate p-value using wilcoxon
def paired_pval(row, ctrl_cols, tumor_cols):
    x = row[ctrl_cols]
    x = [x for x in x if not pd.isnull(x)]
    y = row[tumor_cols]
    y = [y for y in y if not pd.isnull(y)]
    
    try:
        stat, p = wilcoxon(x, y)
        return p
    except:
        return np.nan
    
## calculate p-value using wilcoxon
def paired_stat(row, ctrl_cols, tumor_cols):
    x = row[ctrl_cols]
    y = row[tumor_cols]
    
    try:
        stat, p = wilcoxon(x, y)
        return stat
    except:
        return np.nan

f_brca["wilcox_pval"] = f_brca.apply(paired_pval, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1)
f_brca["wilcox_stat"] = f_brca.apply(paired_stat, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1)
print(len(f_brca))

f_brca_filt = f_brca[~pd.isnull(f_brca["wilcox_pval"])]
print(len(f_brca_filt))

f_brca_filt["wilcox_padj"] = smt.multipletests(list(f_brca_filt["wilcox_pval"]), alpha=0.05, method="fdr_bh")[1]

f_brca_nan["wilcox_pval"] = f_brca_nan.apply(paired_pval, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1)
f_brca_nan["wilcox_stat"] = f_brca_nan.apply(paired_stat, ctrl_cols=paired_ctrl_samps, tumor_cols=paired_tumor_samps, axis=1)
print(len(f_brca_nan))

f_brca_nan_filt = f_brca_nan[~pd.isnull(f_brca_nan["wilcox_pval"])]
print(len(f_brca_nan_filt))

f_brca_nan_filt["wilcox_padj"] = smt.multipletests(list(f_brca_nan_filt["wilcox_pval"]), alpha=0.05, method="fdr_bh")[1]


# In[118]:


for i, row in tcga_paired.iterrows():
    f_brca_filt["paired-diff_%s_ratio" % (i+1)] = f_brca_filt[row.tcga_id_tumor]-f_brca_filt[row.tcga_id_ctrl]
    f_brca_nan_filt["paired-diff_%s_rationan" % (i+1)] = f_brca_nan_filt[row.tcga_id_tumor].fillna(0)-f_brca_nan[row.tcga_id_ctrl].fillna(0)


# In[119]:


paired_ratio_cols = [x for x in f_brca_filt.columns if "paired-diff_" in x]
paired_rationan_cols = [x for x in f_brca_nan_filt.columns if "paired-diff_" in x]


# In[120]:


f_brca_filt["med_paired-diff_ratio"] = f_brca_filt[paired_ratio_cols].median(axis=1)
f_brca_nan_filt["med_paired-diff_rationan"] = f_brca_nan_filt[paired_rationan_cols].median(axis=1)


# ## 10. merge BRCA data w/ DN cats

# In[121]:


f_brca_filt = f_brca_filt.merge(tf_id_map, on="iso_id")
f_brca_nan_filt = f_brca_nan_filt.merge(tf_id_map, on="iso_id")


# In[122]:


f_brca_nan_med_cols = ["clone_acc"] + [x for x in f_brca_nan_filt.columns if "med_" in x] + [x for x in f_brca_nan_filt.columns if "wilcox" in x]


# In[123]:


f_brca_med_cols = ["clone_acc"] + [x for x in f_brca_filt.columns if "med_" in x] + [x for x in f_brca_filt.columns if "wilcox" in x]


# In[124]:


f_brca_filt[f_brca_med_cols]


# In[125]:


dn_data_exp = dn_cats.merge(f_brca_nan_filt[f_brca_nan_med_cols], left_on="tf1p0_id", right_on="clone_acc")
dn_data_exp = dn_data_exp.merge(f_brca_filt[f_brca_med_cols], left_on="tf1p0_id", right_on="clone_acc", 
                                suffixes=("_nan", ""))
print(len(dn_data_exp))


# In[126]:


dn_data_exp.drop_duplicates(subset="tf1p0_id", inplace=True)
print(len(dn_data_exp))


# In[127]:


dn_data_exp["neglog_padj"] = -np.log10(dn_data_exp["wilcox_padj"])
dn_data_exp["neglog_padj_nan"] = -np.log10(dn_data_exp["wilcox_padj_nan"])


# In[128]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(data=dn_data_exp[dn_data_exp["dn_cat"].isin(["ref"])], 
                     x="med_paired-diff_ratio", y="neglog_padj", 
                     color="white", linewidth=0.5, edgecolor="black", alpha=0.8, zorder=10,
                     **{"s": 8})

sns.scatterplot(data=dn_data_exp[dn_data_exp["dn_cat"].isin(["similar", "rewire", "DN"])], 
                x="med_paired-diff_ratio", y="neglog_padj", 
                hue="dn_cat", palette=dn_pal, linewidth=0.25, edgecolor="black", alpha=0.8, zorder=10,
                **{"s": 8}, ax=ax)

for annot_clone, ha, va, offset, relpos, cs in zip(["STAT1|2/7|03H01", "NFIA|2/5|02G09", 
                                                    "TFDP2|1/4|03C12", "ZBTB25|4/5|04D09", "STAT1|1/7|01B05",
                                                    "ZNF451|8/8|06F02", "CREB1|2/2|01F12", "CREB1|1/2|02E01"],
                                                    ["center", "center", "left", "center", "center", "left",
                                                     "left", "right"],
                                                    ["top", "bottom", "center", "bottom", "top", "bottom",
                                                     "bottom", "center"],
                                                    [(7, -7), (0, 7), (7, 0), (-3, 7), (-5, -7), (3, 5),
                                                     (5, 5), (-7, 0)],
                                                    [(0.5, 1), (0.5, 1), (0, 0.5), (0.5, 0), (0.5, 1), (0, 0.5),
                                                     (0, 0.5), (1, 0.5)],
                                                    ["arc3,rad=0.3", "arc3,rad=0", "arc3,rad=0", "arc3,rad=0.3",
                                                     "arc3,rad=0.3", "arc3,rad=0.3", "arc3,rad=0.3", "arc3,rad=0"]):
    row = dn_data_exp[dn_data_exp["tf1p0_id"] == annot_clone].iloc[0]
    if row["dn_cat"] == "ref":
        color = "black"
    else:
        color = dn_pal[row["dn_cat"]]
    print("annot clone: %s | ha: %s | va: %s" % (annot_clone, ha, va))
    
    shorter_id = annot_clone.split("|")[0] + "-" + annot_clone.split("|")[1].split("/")[0]
    ax.annotate(shorter_id, xy=(row["med_paired-diff_ratio"], row["neglog_padj"]), xytext=offset,
                color=color, ha=ha, va=va,
                textcoords='offset points', bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'),
                arrowprops=dict(arrowstyle="-", color=color, relpos=relpos, connectionstyle=cs))

ax.set_xlabel("median isoform difference (tumor - normal)")
ax.set_ylabel("-log10(Wilcoxon adjusted p-value)")
ax.set_title("TF isoforms in breast cancer\n")

#ax.set_xlim((-0.25, 0.12))
#ax.set_ylim((-0.25, 14))
ax.axhline(y=-np.log10(0.05), linestyle="dashed", color="black", linewidth=0.5)
ax.axvline(x=0, linestyle="dashed", color="black", linewidth=0.5)

ax.get_legend().remove()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/BRCA_Volcano.pdf", dpi="figure", bbox_inches="tight")


# In[129]:


## are DNs enriched within those alt. isos that significantly change in BRCA?


# In[130]:


tots = pd.DataFrame(dn_data_exp.dn_cat.value_counts())
sig = pd.DataFrame(dn_data_exp[dn_data_exp["wilcox_padj"] < 0.05].dn_cat.value_counts())
brca_st = tots.join(sig, lsuffix="_tot", rsuffix="_sig")
brca_st = brca_st.loc[["DN", "rewire", "similar", "NA"]]
brca_st = brca_st/brca_st.sum(axis=0)
brca_st


# In[194]:


tots.sum()


# In[196]:


sig


# In[131]:


fe = np.zeros((2, 2))

alts = dn_data_exp[dn_data_exp["dn_cat"].isin(["DN", "rewire", "similar", "NA"])]

fe[0, 0] = len(alts[(alts["dn_cat"] == "DN") & 
                    (alts["wilcox_padj"] < 0.05)].tf1p0_id.unique())
fe[1, 0] = len(alts[(alts["dn_cat"] != "DN") & 
                    (alts["wilcox_padj"] < 0.05)].tf1p0_id.unique())
fe[0, 1] = len(alts[(alts["dn_cat"] == "DN") & 
                    (alts["wilcox_padj"] >= 0.05)].tf1p0_id.unique())
fe[1, 1] = len(alts[(alts["dn_cat"] != "DN") & 
                    (alts["wilcox_padj"] >= 0.05)].tf1p0_id.unique())
fe


# In[132]:


print(fisher_exact(fe))
p = fisher_exact(fe)[1]


# In[181]:


fig, ax = plt.subplots(figsize=(0.5, 1.5))

xs = ["total", "sig. in BRCA"]
y1 = list(brca_st[["dn_cat_tot", "dn_cat_sig"]].loc["NA"])
y2 = list(brca_st[["dn_cat_tot", "dn_cat_sig"]].loc["similar"])
b2 = np.add(y1, y2)
y3 = list(brca_st[["dn_cat_tot", "dn_cat_sig"]].loc["rewire"])
b3 = np.add(b2, y3)
y4 = list(brca_st[["dn_cat_tot", "dn_cat_sig"]].loc["DN"])

ax.bar(xs, y1, color=dn_pal["NA"], label="NA", edgecolor="black", linewidth=0.5)
ax.bar(xs, y2, bottom=y1, color=dn_pal["similar"], label="similar", edgecolor="black", linewidth=0.5)
ax.bar(xs, y3, bottom=b2, color=dn_pal["rewire"], label="rewire", edgecolor="black", linewidth=0.5)
ax.bar(xs, y4, bottom=b3, color=dn_pal["DN"], label="neg. reg.", edgecolor="black", linewidth=0.5)

# annotate pval
annotate_pval(ax, 0, 1, 1.025, 0, 1.025, p, fontsize-1)

# add legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)

ax.set_ylabel("% Alternative isoforms")
ax.set_xticklabels(xs, rotation=30, ha="right", va="top")

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/dn_stacked_bar_brca.pdf", dpi="figure", bbox_inches="tight")


# ### CREB1 vignette

# In[134]:


brca_isos = brca_isos.merge(tf_id_map[["iso_id", "gene_name"]], on="iso_id").drop_duplicates()
len(brca_isos)


# In[135]:


brca_isos_paired = brca_isos[["gene_name", "iso_id"] + tcga_paired_ctrls + tcga_paired_tumors]
new_ctrl_cols = ["normal - %s" % (i+1) for i, x in enumerate(tcga_paired_ctrls)]
new_tumor_cols = ["tumor - %s" % (i+1) for i, x in enumerate(tcga_paired_tumors)]
brca_isos_paired.columns = ["gene_name", "iso_id"] + new_ctrl_cols + new_tumor_cols


# In[136]:


def brca_expression_plot(gene_name, figsize, ylim, df, cols, fig_suffix, ctrls_line, tumor_line):
    df_sub = df[df["gene_name"] == gene_name]
    df_sub.set_index("iso_id", inplace=True)
    df_sub = df_sub[cols].drop_duplicates()
    #print(df_sub.head())
    n_isos = len(df_sub)
    palette = sns.color_palette("husl", as_cmap=False, n_colors=n_isos)
    fig, axes = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(figsize)
    ### bar chart ###
    (df_sub
          .T
          .plot.bar(ax=axes[0],
                    legend=False,
                    width=0.7,
                    color=list(palette)))
    ### percentages ###
    (df_sub.div(df_sub.sum(axis=0))
              .T.plot.bar(ax=axes[1], 
                          stacked=True,
                          legend=False,
                          color=list(palette)))
    axes[0].set_yscale("symlog")
    axes[0].set_ylabel('tpm')
    #axes[0].set_ylim(ylim)
    axes[1].set_ylabel('percent')
    axes[1].set_yticklabels(['{:.0%}'.format(t) for t in axes[1].get_yticks()])
    axes[1].legend(loc='lower left', bbox_to_anchor=(1, 0))
    axes[0].axhline(y=1, color='black', linewidth=0.5, linestyle="dashed")
    
    # add medians
    axes[1].plot(ctrls_line[0], ctrls_line[1], color="black", linewidth=0.5, linestyle="dashed")
    axes[1].plot(tumor_line[0], tumor_line[1], color="black", linewidth=0.5, linestyle="dashed")
    
    axes[0].spines['right'].set_visible(False)
    axes[0].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    
    plt.subplots_adjust(hspace=0.15)
    plt.savefig('../../figures/fig7/brca_' + gene_name + '_' + fig_suffix + '.pdf',
                bbox_inches='tight')


# In[137]:


creb1 = dn_data_exp[dn_data_exp["gene_name"] == "CREB1"][["gene_name", "tf1p0_id", "dn_cat", "med_paired-brca_ratio",
                                                          "med_paired-ctrls_ratio", "med_paired-diff_ratio"]]
creb1


# In[138]:


cols = new_ctrl_cols[0:35] + new_tumor_cols[0:35]
brca_expression_plot("CREB1", (9, 3), (0, 6), brca_isos_paired, cols, "paired",
                     ([0, 34.5], [creb1[creb1["dn_cat"] == "ref"]["med_paired-ctrls_ratio"].iloc[0], 
                                  creb1[creb1["dn_cat"] == "ref"]["med_paired-ctrls_ratio"].iloc[0]]), 
                     ([34.5, 70], [creb1[creb1["dn_cat"] == "ref"]["med_paired-brca_ratio"].iloc[0], 
                                   creb1[creb1["dn_cat"] == "ref"]["med_paired-brca_ratio"].iloc[0]]))


# In[139]:


f_brca_paired = f_brca_filt[["gene_name", "iso_id"] + tcga_paired_ctrls + tcga_paired_tumors]
new_ctrl_cols = ["normal - %s" % (i+1) for i, x in enumerate(tcga_paired_ctrls)]
new_tumor_cols = ["tumor - %s" % (i+1) for i, x in enumerate(tcga_paired_tumors)]
f_brca_paired.columns = ["gene_name", "iso_id"] + new_ctrl_cols + new_tumor_cols


# In[140]:


f_brca_paired_melt = pd.melt(f_brca_paired, id_vars=["gene_name", "iso_id"])
f_brca_paired_melt["samp"] = f_brca_paired_melt["variable"].str.split(" ", expand=True)[0]


# In[141]:


dn_data_exp[["gene_name", "tf1p0_id"]]


# In[142]:


dn_data_exp = dn_data_exp.merge(tf_id_map[["gene_name", "clone_acc", "iso_id", "merge_id"]], 
                                left_on=["gene_name", "tf1p0_id"],
                                right_on=["gene_name", "clone_acc"]).drop_duplicates(subset="tf1p0_id")
print(len(dn_data_exp))


# In[180]:


tmp = f_brca_paired_melt[f_brca_paired_melt["gene_name"] == "CREB1"]

fig = plt.figure(figsize=(1, 1.5))

ax = sns.boxplot(data=tmp, x="iso_id", y="value", hue="samp", fliersize=0,
                 palette={"normal": "gray", "tumor": sns.color_palette("Set2")[3]},
                 order=["CREB1-2", "CREB1-1"])
mimic_r_boxplot(ax)

# sns.swarmplot(data=tmp, x="iso_id", y="value", hue="samp",
#               palette={"normal": "gray", "tumor": sns.color_palette("Set2")[3]}, ax=ax,
#               size=1, edgecolor="black", linewidth=0.5, alpha=0.5, split=True,
#               order=["CREB1-2", "CREB1-1"])

# annotate w p-vals
ys = [0.52, 0.96]
for i, iso in enumerate(tmp.iso_id.unique()):
    print(iso)
    padj = dn_data_exp[dn_data_exp["iso_id"]==iso]["wilcox_padj"].iloc[0]
    print(padj)
    annotate_pval(ax, i-0.2, i+0.2, ys[i], 0, ys[i], padj, fontsize-1)

ax.set_xlabel("")
ax.set_xticklabels(["CREB1-2 (ref)", "CREB1-1 (alt)"], rotation=30, va="top", ha="right")
ax.set_ylabel("isoform ratio")
ax.set_title("CREB1 isoforms\nin breast cancer\n\n")
ax.get_legend().remove()
ax.set_ylim((-0.05, 1))

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig("../../figures/fig7/BRCA_CREB1_boxplot.pdf", dpi="figure", bbox_inches="tight")


# In[144]:


from plotting import (y2h_ppi_per_tf_gene_plot,
                      y1h_pdi_per_tf_gene_plot,
                      m1h_activation_per_tf_gene_plot)

from data_loading import (load_isoform_and_paralog_y2h_data,
                          load_m1h_activation_data,
                          load_y1h_pdi_data)


# In[145]:


y2h = load_isoform_and_paralog_y2h_data()
m1h = load_m1h_activation_data(add_missing_data=True)
y1h = load_y1h_pdi_data()


# In[177]:


gene_name = "CREB1"


# In[178]:


fig, ax = plt.subplots(figsize=(3, 0.6))

tfs[gene_name].exon_diagram(ax=ax, )
fig.savefig("../../figures/fig7/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[179]:


fig, ax = plt.subplots(figsize=(3, 0.6))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig7/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[182]:


fig, ax = plt.subplots(1, 1, figsize=(1, 0.7))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, iso_order=["CREB1-2", "CREB1-1"], xlim=(0, 8.2))
plt.savefig('../../figures/fig7/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[181]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(1.25, 1))
y1h_pdi_per_tf_gene_plot(tf.name, ax=ax, data=y1h, iso_order=["CREB1-2", "CREB1-1"])
plt.savefig('../../figures/fig7/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[174]:


def developmental_tissue_expression_plot(gene_name, palette_name, figsize, ylim, means, cols, fig_suffix, shorten_x=False):
    n_isos = len(means.loc[genes == gene_name])
    palette = sns.color_palette(palette_name, as_cmap=False, n_colors=n_isos)
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
    
    if shorten_x:
        xticks = list(axes[1].get_xticklabels())
        xticks_short = [x.get_text().split("-")[0].strip() for x in xticks]
        axes[1].set_xticklabels(xticks_short, rotation=90, va="top", ha="center")
        
    axes[0].set_ylabel('log2(tpm + 1)\n')
    axes[0].set_ylim(ylim)
    axes[1].set_ylabel('percent')
    axes[1].set_yticklabels(['{:.0%}'.format(t) for t in axes[1].get_yticks()])
    axes[1].legend(loc='lower left', bbox_to_anchor=(1, 0), frameon=False)
    axes[0].axhline(y=1, color='black', linewidth=0.5, linestyle="dashed")
    
    for spine in ['right', 'top']:
        axes[0].spines[spine].set_visible(False)
        axes[1].spines[spine].set_visible(False)
    
    plt.subplots_adjust(hspace=0.25)
    plt.savefig('../../figures/fig7/expression_' + gene_name + '_' + fig_suffix + '.pdf',
                bbox_inches='tight')


# In[153]:


if not (genes_gtex == genes_dev).all():
        raise UserWarning()
genes = genes_gtex


# In[154]:


heart_cols = [x for x in means_dev.columns if "heart" in x]
brain_cols = [x for x in means_dev.columns if "brain" in x]
liver_cols = [x for x in means_dev.columns if "liver" in x]
developmental_tissue_expression_plot("CREB1", "husl", (6, 1.75), (0, 6), means_dev, 
                                     heart_cols + brain_cols + liver_cols, 
                                     "means_dev_heart_brain_liver")


# In[155]:


developmental_tissue_expression_plot("CREB1", "husl", (5, 1.75), (0, 5), means_gtex, 
                                     means_gtex.columns, 
                                     "means_gtex_all")


# In[177]:


developmental_tissue_expression_plot("CREB1", "husl", (4.2, 1.5), (0, 5), means_gtex, 
                                     means_gtex.columns, 
                                     "means_gtex_all_short", shorten_x=True)


# ### other vignettes

# In[160]:


gene_name = "HSF2"


# In[161]:


fig, ax = plt.subplots(figsize=(5, 0.6))

tfs[gene_name].exon_diagram(ax=ax, )
fig.savefig("../../figures/fig7/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[162]:


fig, ax = plt.subplots(figsize=(6, 0.6))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig7/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[163]:


fig, ax = plt.subplots(1, 1, figsize=(2, 0.6))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(0, 2))
plt.savefig('../../figures/fig7/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[164]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(1.25, 0.6))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../../figures/fig7/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[165]:


gene_name = "NFYA"


# In[166]:


fig, ax = plt.subplots(figsize=(4, 0.8))

tfs[gene_name].exon_diagram(ax=ax, )
fig.savefig("../../figures/fig7/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[167]:


fig, ax = plt.subplots(figsize=(6, 0.6))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig7/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[168]:


fig, ax = plt.subplots(1, 1, figsize=(1, 0.8))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(0, 5))
plt.savefig('../../figures/fig7/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[169]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(3.5, 1.5))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../../figures/fig7/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[170]:


pairs[pairs["gene_symbol"] == "HSF2"][["dn_cat", "clone_acc_alt", "activation_fold_change_log2"]]


# In[171]:


gene_name = "PPARG"


# In[172]:


fig, ax = plt.subplots(figsize=(4, 2))

tfs[gene_name].exon_diagram(ax=ax, )
fig.savefig("../../figures/fig7/{}_exon_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[173]:


fig, ax = plt.subplots(figsize=(4, 2))

tfs[gene_name].protein_diagram(only_cloned_isoforms=True, draw_legend=False, ax=ax)
fig.savefig("../../figures/fig7/{}_protein_diagram.pdf".format(gene_name), bbox_inches="tight", dpi="figure")


# In[174]:


fig, ax = plt.subplots(1, 1, figsize=(1, 1))

df = m1h_activation_per_tf_gene_plot(gene_name, data=m1h, ax=ax, xlim=(-1.1, 3))
plt.savefig('../../figures/fig7/{}_m1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[175]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(3.25, 1))
y2h_ppi_per_tf_gene_plot(tf.name, ax=ax, data=y2h)
plt.savefig('../../figures/fig7/{}_y2h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[176]:


tf = tfs[gene_name]
fig, ax = plt.subplots(1, 1, figsize=(3.25, 1))
y1h_pdi_per_tf_gene_plot(tf.name, ax=ax, data=y1h)
plt.savefig('../../figures/fig7/{}_y1h-profile.pdf'.format(gene_name), bbox_inches='tight')


# In[232]:


pairs[pairs["gene_symbol"] == "PPARG"]


# In[ ]:




