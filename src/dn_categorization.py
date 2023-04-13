
# coding: utf-8

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

from scipy.stats import fisher_exact

import plotting
from plotting import PAPER_PRESET, PAPER_FONTSIZE

get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'svg'")
mpl.rcParams['figure.autolayout'] = False


# In[3]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# ## variables

# In[4]:


# note to kaia- make sure these files are in the 'data' directory


# In[5]:


dn_f = "../data/processed/TF-iso_ref-vs-alt.bug_fix.tsv"


# ## 1. import data

# In[6]:


dn = pd.read_table(dn_f)
dn


# ## 2. categorize M1H data

# In[7]:


def m1h_cat(row):
    
    # ref is activator
    if row.activation_ref >= 1:
        if row.activation_fold_change_log2 <= -1 and row.activation_alt <= 1 and row.activation_alt >= -1:
            return "activation loss"
        elif not pd.isnull(row.activation_fold_change_log2):
            return "rewire"
        else:
            return "NA"
    
    # ref is repressor
    elif row.activation_ref <= -1:
        if row.activation_fold_change_log2 >= 1 and row.activation_alt <= 1 and row.activation_alt >= -1:
            return "repression loss"
        elif not pd.isnull(row.activation_fold_change_log2):
            return "rewire"
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
        elif not pd.isnull(row.activation_fold_change_log2):
            return "rewire"
        else:
            return "NA"
        
dn["m1h_cat"] = dn.apply(m1h_cat, axis=1)
dn.m1h_cat.value_counts()


# ## 3. categorize Y1H data

# In[8]:


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
    
dn["y1h_cat"] = dn.apply(y1h_cat, axis=1)
dn.y1h_cat.value_counts()


# ## 4. categorize Y2H data

# In[9]:


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
        return "PPI rewire"
    
dn["y2h_cat"] = dn.apply(y2h_cat, axis=1)
dn.y2h_cat.value_counts()


# ## 5. categorize negative regulators (dominant negatives)

# In[10]:


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
                if row.dbd_pct_lost >= 50:
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
                if row.dbd_pct_lost >= 50:
                    n.append("DBD loss")
                
                if len(n) > 0:
                    s = ",".join(n)
                    s = "DN (%s)" % s
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
            if row.dbd_pct_lost >= 50:
                n.append("DBD loss")

            if len(n) > 0:
                s = ",".join(n)
                s = "DN (%s)" % s
            else:
                s = "rewire"

            return s
            
dn["dn_cat"] = dn.apply(dn_cat, axis=1)
dn.dn_cat.value_counts()


# In[11]:


dn[dn["gene_symbol"] == "GRHL3"]


# ## 6. extract more DN info for plotting

# In[12]:


dn["dn_short"] = dn["dn_cat"].str.split(" ", expand=True)[0]
dn.dn_short.value_counts()


# In[13]:


def mech_bool(row, mech_col):
    if "DN" in row.dn_cat:
        if mech_col in row.dn_cat:
            return True
        else:
            return False
    else:
        return np.nan
    
dn["dn_ppi"] = dn.apply(mech_bool, mech_col="PPIs", axis=1)
dn["dn_pdi"] = dn.apply(mech_bool, mech_col="PDIs", axis=1)
dn["dn_activ"] = dn.apply(mech_bool, mech_col="activ", axis=1)
dn["dn_dbd"] = dn.apply(mech_bool, mech_col="DBD loss", axis=1)
dn[dn["dn_short"] == "DN"].sample(5)


# ## 7. now make a few plots

# In[14]:


fig = plt.figure(figsize=(1.5, 1.75))

ax = sns.countplot(data=dn, x="dn_short", palette=sns.color_palette("Set2"),
                   order=["DN", "rewire", "NA"])
ax.set_xticklabels(["putative DN", "putative re-wirer", "NA"], ha="right", va="top", rotation=30)
ax.set_xlabel("")
ax.set_ylabel("count of alternative TF isoforms")

fig.savefig("../figures/DN_countplot.pdf", dpi="figure", bbox_inches="tight")


# In[15]:


from upsetplot import plot


# In[16]:


ppis = list(set(list(dn[dn["dn_ppi"] == True]["clone_acc_alt"])))
pdis = list(set(list(dn[dn["dn_pdi"] == True]["clone_acc_alt"])))
activ = list(set(list(dn[dn["dn_activ"] == True]["clone_acc_alt"])))
dbd = list(set(list(dn[dn["dn_dbd"] == True]["clone_acc_alt"])))

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

fig.savefig("../figures/DN_negreg_upset.pdf", dpi="figure", bbox_inches="tight")


# In[17]:


rw = dn[dn["dn_cat"] == "rewire"]
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

fig.savefig("../figures/DN_rewire_upset.pdf", dpi="figure", bbox_inches="tight")


# In[18]:


ppis = list(set(list(dn[dn["y2h_cat"] != "NA"]["clone_acc_alt"])))
pdis = list(set(list(dn[dn["y1h_cat"] != "NA"]["clone_acc_alt"])))
activ = list(set(list(dn[dn["m1h_cat"] != "NA"]["clone_acc_alt"])))

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

fig.savefig("../figures/DN_pairs_assessed_upset.pdf", dpi="figure", bbox_inches="tight")


# In[19]:


y = np.array([len(dn[dn["dn_short"] == "rewire"]), len(dn[dn["dn_short"] == "DN"]),
              len(dn[dn["dn_short"] == "NA"]), len(dn[dn["dn_short"] == "likely"])])
labels = ["rewire", "negative regulator", "NA", "likely non-functional"]
colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], "lightgray", "darkgray"]

fig = plt.figure(figsize=(1.75, 1.75))
ws, ls, ns = plt.pie(y, labels=labels, colors=colors, autopct='%1.1f%%', startangle=-45, explode=(0, 0.1, 0, 0))
for w in ws:
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
ns[3].set_color("white")

fig.savefig("../figures/dn_pie.incl_NA.pdf", dpi="figure", bbox_inches="tight")


# In[20]:


ys = np.array([len(dn[dn["dn_short"] == "rewire"]), len(dn[dn["dn_short"] == "DN"]),
              len(dn[dn["dn_short"] == "likely"])])
labels = ["rewirer", "negative regulator", "likely non-functional"]
colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1], "darkgray"]

fig, ax = plt.subplots(figsize=(2.0, 2.0), subplot_kw=dict(aspect="equal"))
ws, ls, ns = ax.pie(ys, colors=colors, labels=labels, autopct='%1.1f%%', startangle=90, explode=(0.05, 0.05, 0.1))

for n, w in zip(ns, ws):
    w.set_linewidth(0.5)
    w.set_edgecolor("black")
    n.set_fontweight("bold")
ns[2].set_text("")

fig.savefig("../figures/dn_pie.no_NA.pdf", dpi="figure", bbox_inches="tight")


# In[21]:


outer_ys = np.array([len(dn[(dn["dn_short"] == "rewire")]), 
                     len(dn[(dn["dn_short"] == "DN")])])
outer_labels = ["rewirer\n(%s%%)" % round((outer_ys[0]/np.sum(outer_ys)*100),1), 
                "negative regulator\n(%s%%)" % round((outer_ys[1]/np.sum(outer_ys)*100),1)]
outer_colors = [sns.color_palette("Set2")[2], sns.color_palette("Set2")[1]]

inner_ys = np.array([len(dn[(dn["dn_short"] == "rewire") & (dn["is_alt_novel_isoform"])]), 
               len(dn[(dn["dn_short"] == "rewire") & (~dn["is_alt_novel_isoform"])]), 
               len(dn[(dn["dn_short"] == "DN") & (dn["is_alt_novel_isoform"])]),
               len(dn[(dn["dn_short"] == "DN") & (~dn["is_alt_novel_isoform"])])])
inner_colors = [sns.light_palette(sns.color_palette("Set2")[2])[3], 
                sns.light_palette(sns.color_palette("Set2")[2])[1], 
                sns.light_palette(sns.color_palette("Set2")[1])[3], 
                sns.light_palette(sns.color_palette("Set2")[1])[1]]
hatches = ['++', '', '++', '']


fig, ax = plt.subplots(figsize=(2.0, 2.0), subplot_kw=dict(aspect="equal"))

o_ws, o_ls = ax.pie(outer_ys, colors=outer_colors, labels=outer_labels,
                          startangle=90, radius=1, wedgeprops=dict(width=0.3, edgecolor='w'))
i_ws, i_ls, i_ns = ax.pie(inner_ys, colors=inner_colors, autopct='%1.0f%%', 
                          startangle=90, radius=0.7, 
                          wedgeprops=dict(width=0.4, edgecolor='w'),
                          textprops={'size': 'smaller'}, pctdistance=0.7)

for i, w in enumerate(i_ws):
    w.set(hatch=hatches[i])
    
ax.set_title("alternative TF isoform categories\n(%s reference-alternative pairs)" % (np.sum(outer_ys)))

fig.savefig("../figures/dn_pie.novel_nested.pdf", dpi="figure", bbox_inches="tight")


# In[22]:


# create df for stacked bar chart
delta_pdis = dn[~dn["y1h_cat"].isin(["NA", "no PDI change"])]
pdis_vc = pd.DataFrame(delta_pdis.dn_short.value_counts()).reset_index()

delta_ppis = dn[~dn["y2h_cat"].isin(["NA", "no PPI change (important PPIs)"])]
ppis_vc = pd.DataFrame(delta_ppis.dn_short.value_counts()).reset_index()

delta_activ = dn[~dn["m1h_cat"].isin(["NA"])]
activ_vc = pd.DataFrame(delta_activ.dn_short.value_counts()).reset_index()

mrg = pdis_vc.merge(ppis_vc, on="index", how="outer").merge(activ_vc, on="index", how="outer")
mrg.fillna(0, inplace=True)
mrg.columns = ["index", "PDIs", "PPIs", "activity"]

to_plot = pd.melt(mrg, id_vars="index")
to_plot.sample(5)


# In[23]:


pal = {"likely": "darkgray", "NA": "lightgray", "rewire": sns.color_palette("Set2")[2],
       "DN": sns.color_palette("Set2")[1]}

fig = plt.figure(figsize=(2.8, 1.5))

ax = sns.barplot(data=to_plot[~to_plot["index"].isin(["likely", "NA"])], 
                 x="variable", hue="index", y="value", palette=pal)
plt.legend(loc=2, bbox_to_anchor=(1.01, 1))

ax.set_xlabel("")
ax.set_xticklabels(["∆ PDIs", "∆ PPIs", "∆ activity"])
ax.set_ylabel("count of ref/alt pairs with data")
ax.set_ylim((0, 100))

sub = to_plot[~to_plot["index"].isin(["likely", "NA"])]
for i, cat in enumerate(["PDIs", "PPIs", "activity"]):
    cat_sub = sub[(sub["variable"] == cat)]
    n_dn = cat_sub[cat_sub["index"] == "DN"].value.iloc[0]
    n_rw = cat_sub[cat_sub["index"] == "rewire"].value.iloc[0]
    p_dn = round(n_dn/np.sum(cat_sub["value"])*100, 1)
    p_rw = round(n_rw/np.sum(cat_sub["value"])*100, 1)
    
    ax.text(i-0.2, n_dn, "%s%%" % p_dn, va="bottom", ha="center", color=pal["DN"], fontsize=5.5)
    ax.text(i+0.2, n_rw, "%s%%" % p_rw, va="bottom", ha="center", color=pal["rewire"], fontsize=5.5)


# In[24]:


# make stacked barchart situation
tmp = dn[dn["dn_short"] == "DN"]
dn_pdi_change = len(tmp[tmp["dn_pdi"] == True])
dn_ppi_change = len(tmp[tmp["dn_ppi"] == True])
dn_activ_change = len(tmp[tmp["dn_activ"] == True])
dn_dbd_change = len(tmp[tmp["dn_dbd"] == True])
tot_dn = dn_pdi_change + dn_ppi_change + dn_activ_change + dn_dbd_change

tmp = dn[dn["dn_short"] == "rewire"]
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


# In[25]:


colors = met_brewer.met_brew(name="Hokusai3", n=4, brew_type="discrete")
sns.palplot(colors)


# In[26]:


fig, ax = plt.subplots(figsize=(2.1, 1.5))

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
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1))

ax.set_ylabel("percent prevalence")
ax.set_xticklabels(["negative regulator", "rewirer"])

fig.savefig("../figures/dn_stacked_bar.pdf", dpi="figure", bbox_inches="tight")


# In[27]:


dn.to_csv("../data/processed/TF-iso_ref-vs-alt.DN_cat.tsv", sep="\t", index=False)

