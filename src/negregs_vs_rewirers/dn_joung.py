# coding: utf-8

# In[1]:


import warnings

warnings.filterwarnings("ignore")


# In[2]:


import matplotlib as mpl
import matplotlib.pyplot as plt
import met_brewer
import pandas as pd
import numpy as np
import seaborn as sns
import sys
import upsetplot

import statsmodels.api as sm
import statsmodels.formula.api as smf

from Bio.Seq import Seq
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr

import plotting
from plotting import PAPER_PRESET, PAPER_FONTSIZE, nice_boxplot, mimic_r_boxplot


get_ipython().run_line_magic("matplotlib", "inline")
get_ipython().run_line_magic("config", "InlineBackend.figure_format = 'svg'")
mpl.rcParams["figure.autolayout"] = False


# In[3]:


from data_loading import (
    load_annotated_TFiso1_collection,
    load_valid_isoform_clones,
    load_developmental_tissue_expression_gencode,
)


# In[4]:


sns.set(**PAPER_PRESET)
fontsize = PAPER_FONTSIZE


# ## functions

# In[5]:


# base code
import numpy as np
import seaborn as sns
import statsmodels
from statsmodels.tools.tools import maybe_unwrap_results
from statsmodels.graphics.gofplots import ProbPlot
from statsmodels.stats.outliers_influence import variance_inflation_factor
import matplotlib.pyplot as plt
from typing import Type

style_talk = "seaborn-talk"  # refer to plt.style.available


class Linear_Reg_Diagnostic:
    """
    Diagnostic plots to identify potential problems in a linear regression fit.
    Mainly,
        a. non-linearity of data
        b. Correlation of error terms
        c. non-constant variance
        d. outliers
        e. high-leverage points
        f. collinearity

    Author:
        Prajwal Kafle (p33ajkafle@gmail.com, where 3 = r)
        Does not come with any sort of warranty.
        Please test the code one your end before using.
    """

    def __init__(
        self,
        results: Type[statsmodels.regression.linear_model.RegressionResultsWrapper],
    ) -> None:
        """
        For a linear regression model, generates following diagnostic plots:

        a. residual
        b. qq
        c. scale location and
        d. leverage

        and a table

        e. vif

        Args:
            results (Type[statsmodels.regression.linear_model.RegressionResultsWrapper]):
                must be instance of statsmodels.regression.linear_model object

        Raises:
            TypeError: if instance does not belong to above object

        Example:
        >>> import numpy as np
        >>> import pandas as pd
        >>> import statsmodels.formula.api as smf
        >>> x = np.linspace(-np.pi, np.pi, 100)
        >>> y = 3*x + 8 + np.random.normal(0,1, 100)
        >>> df = pd.DataFrame({'x':x, 'y':y})
        >>> res = smf.ols(formula= "y ~ x", data=df).fit()
        >>> cls = Linear_Reg_Diagnostic(res)
        >>> cls(plot_context="seaborn-paper")

        In case you do not need all plots you can also independently make an individual plot/table
        in following ways

        >>> cls = Linear_Reg_Diagnostic(res)
        >>> cls.residual_plot()
        >>> cls.qq_plot()
        >>> cls.scale_location_plot()
        >>> cls.leverage_plot()
        >>> cls.vif_table()
        """

        if (
            isinstance(
                results, statsmodels.regression.linear_model.RegressionResultsWrapper
            )
            is False
        ):
            raise TypeError(
                "result must be instance of statsmodels.regression.linear_model.RegressionResultsWrapper object"
            )

        self.results = maybe_unwrap_results(results)

        self.y_true = self.results.model.endog
        self.y_predict = self.results.fittedvalues
        self.xvar = self.results.model.exog
        self.xvar_names = self.results.model.exog_names

        self.residual = np.array(self.results.resid)
        influence = self.results.get_influence()
        self.residual_norm = influence.resid_studentized_internal
        self.leverage = influence.hat_matrix_diag
        self.cooks_distance = influence.cooks_distance[0]
        self.nparams = len(self.results.params)

    def __call__(self, plot_context="seaborn-paper"):
        # print(plt.style.available)
        with plt.style.context(plot_context):
            fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
            self.residual_plot(ax=ax[0, 0])
            self.qq_plot(ax=ax[0, 1])
            self.scale_location_plot(ax=ax[1, 0])
            self.leverage_plot(ax=ax[1, 1])
            plt.show()

        self.vif_table()
        return fig, ax

    def residual_plot(self, ax=None):
        """
        Residual vs Fitted Plot

        Graphical tool to identify non-linearity.
        (Roughly) Horizontal red line is an indicator that the residual has a linear pattern
        """
        if ax is None:
            fig, ax = plt.subplots()

        sns.residplot(
            x=self.y_predict,
            y=self.residual,
            lowess=True,
            scatter_kws={"alpha": 0.5},
            line_kws={"color": "red", "lw": 1, "alpha": 0.8},
            ax=ax,
        )

        # annotations
        residual_abs = np.abs(self.residual)
        abs_resid = np.flip(np.sort(residual_abs))
        abs_resid_top_3 = abs_resid[:3]
        for i, _ in enumerate(abs_resid_top_3):
            ax.annotate(i, xy=(self.y_predict[i], self.residual[i]), color="C3")

        ax.set_title("Residuals vs Fitted", fontweight="bold")
        ax.set_xlabel("Fitted values")
        ax.set_ylabel("Residuals")
        return ax

    def qq_plot(self, ax=None):
        """
        Standarized Residual vs Theoretical Quantile plot

        Used to visually check if residuals are normally distributed.
        Points spread along the diagonal line will suggest so.
        """
        if ax is None:
            fig, ax = plt.subplots()

        QQ = ProbPlot(self.residual_norm)
        QQ.qqplot(line="45", alpha=0.5, lw=1, ax=ax)

        # annotations
        abs_norm_resid = np.flip(np.argsort(np.abs(self.residual_norm)), 0)
        abs_norm_resid_top_3 = abs_norm_resid[:3]
        for r, i in enumerate(abs_norm_resid_top_3):
            ax.annotate(
                i,
                xy=(np.flip(QQ.theoretical_quantiles, 0)[r], self.residual_norm[i]),
                ha="right",
                color="C3",
            )

        ax.set_title("Normal Q-Q", fontweight="bold")
        ax.set_xlabel("Theoretical Quantiles")
        ax.set_ylabel("Standardized Residuals")
        return ax

    def scale_location_plot(self, ax=None):
        """
        Sqrt(Standarized Residual) vs Fitted values plot

        Used to check homoscedasticity of the residuals.
        Horizontal line will suggest so.
        """
        if ax is None:
            fig, ax = plt.subplots()

        residual_norm_abs_sqrt = np.sqrt(np.abs(self.residual_norm))

        ax.scatter(self.y_predict, residual_norm_abs_sqrt, alpha=0.5)
        sns.regplot(
            x=self.y_predict,
            y=residual_norm_abs_sqrt,
            scatter=False,
            ci=False,
            lowess=True,
            line_kws={"color": "red", "lw": 1, "alpha": 0.8},
            ax=ax,
        )

        # annotations
        abs_sq_norm_resid = np.flip(np.argsort(residual_norm_abs_sqrt), 0)
        abs_sq_norm_resid_top_3 = abs_sq_norm_resid[:3]
        for i in abs_sq_norm_resid_top_3:
            ax.annotate(
                i, xy=(self.y_predict[i], residual_norm_abs_sqrt[i]), color="C3"
            )
        ax.set_title("Scale-Location", fontweight="bold")
        ax.set_xlabel("Fitted values")
        ax.set_ylabel(r"$\sqrt{|\mathrm{Standardized\ Residuals}|}$")
        return ax

    def leverage_plot(self, ax=None):
        """
        Residual vs Leverage plot

        Points falling outside Cook's distance curves are considered observation that can sway the fit
        aka are influential.
        Good to have none outside the curves.
        """
        if ax is None:
            fig, ax = plt.subplots()

        ax.scatter(self.leverage, self.residual_norm, alpha=0.5)

        sns.regplot(
            x=self.leverage,
            y=self.residual_norm,
            scatter=False,
            ci=False,
            lowess=True,
            line_kws={"color": "red", "lw": 1, "alpha": 0.8},
            ax=ax,
        )

        # annotations
        leverage_top_3 = np.flip(np.argsort(self.cooks_distance), 0)[:3]
        for i in leverage_top_3:
            ax.annotate(i, xy=(self.leverage[i], self.residual_norm[i]), color="C3")

        xtemp, ytemp = self.__cooks_dist_line(0.5)  # 0.5 line
        ax.plot(xtemp, ytemp, label="Cook's distance", lw=1, ls="--", color="red")
        xtemp, ytemp = self.__cooks_dist_line(1)  # 1 line
        ax.plot(xtemp, ytemp, lw=1, ls="--", color="red")

        ax.set_xlim(0, max(self.leverage) + 0.01)
        ax.set_title("Residuals vs Leverage", fontweight="bold")
        ax.set_xlabel("Leverage")
        ax.set_ylabel("Standardized Residuals")
        ax.legend(loc="upper right")
        return ax

    def vif_table(self):
        """
        VIF table

        VIF, the variance inflation factor, is a measure of multicollinearity.
        VIF > 5 for a variable indicates that it is highly collinear with the
        other input variables.
        """
        vif_df = pd.DataFrame()
        vif_df["Features"] = self.xvar_names
        vif_df["VIF Factor"] = [
            variance_inflation_factor(self.xvar, i) for i in range(self.xvar.shape[1])
        ]

        print(vif_df.sort_values("VIF Factor").round(2))

    def __cooks_dist_line(self, factor):
        """
        Helper function for plotting Cook's distance curves
        """
        p = self.nparams
        formula = lambda x: np.sqrt((factor * p * (1 - x)) / x)
        x = np.linspace(0.001, max(self.leverage), 50)
        y = formula(x)
        return x, y


# In[6]:


pal = {
    "ref": sns.color_palette("Set2")[0],
    "rewire": sns.color_palette("Set2")[2],
    "DN": sns.color_palette("Set2")[1],
    "NA": "lightgray",
    "likely": "darkgray",
}


# ## variables

# In[7]:


dn_f = "../data/processed/TF-iso_ref-vs-alt.DN_cat.tsv"


# In[8]:


joung_orf_f = "../data/external/joung_files/Joung_ORF_lib.txt"
joung_data_f = "../data/external/joung_files/Joung_ORF_scores.txt"
joung_cells_f = "../data/external/joung_files/Figure3B_celltype_mapping.csv"


# In[9]:


joung_down_map_batch_f = "../data/external/joung_files/subsample_mapping_batch.txt"
joung_down_map_TF_f = "../data/external/joung_files/subsample_mapping_TF.txt"
joung_down_map_louvain_f = "../data/external/joung_files/subsample_mapping_louvain.txt"


# ## 1. import data

# In[10]:


dn = pd.read_table(dn_f)


# In[11]:


joung_orf = pd.read_table(joung_orf_f)
joung_orf["Name"] = joung_orf["Name"].str.strip()


# In[12]:


joung_data = pd.read_table(joung_data_f)
joung_data["Name"] = joung_data["TF ORF"].str.split("-", expand=True)[0].str.strip()


# In[13]:


joung_cells = pd.read_table(joung_cells_f, sep=",")


# In[14]:


joung_down_map_batch = pd.read_table(joung_down_map_batch_f, index_col=0)
print(len(joung_down_map_batch))
joung_down_map_TF = pd.read_table(joung_down_map_TF_f, index_col=0)
print(len(joung_down_map_TF))
joung_down_map_louvain = pd.read_table(joung_down_map_louvain_f, index_col=0)
print(len(joung_down_map_louvain))

joung_down_map = joung_down_map_batch.join(joung_down_map_TF).join(
    joung_down_map_louvain
)
print(len(joung_down_map))


# ## 2. map our clones to their ORFs

# In[15]:


tfs = load_annotated_TFiso1_collection()


# In[16]:


def pd_translate(row):
    s = Seq(row["ORF sequence"])
    aa = s.translate()
    return str(aa)


joung_orf["seq_aa"] = joung_orf.apply(pd_translate, axis=1)


# In[17]:


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
                joung_sub = joung_orf[
                    joung_orf["RefSeq and Gencode ID"].str.contains(iso_enst)
                ]
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


# In[18]:


tf_id_map_df = pd.DataFrame.from_dict(tf_id_map, orient="index").reset_index()
print(len(tf_id_map_df))
tf_id_map_df.sample(5)


# In[19]:


joung_orf = joung_orf.merge(
    tf_id_map_df,
    left_on="Name",
    right_on="joung_id",
    how="left",
    suffixes=("_joung", "_tf1p0"),
)
joung_orf.sample(5)


# ## 3. intersect our DN categorizations with their scores

# In[20]:


joung_data = joung_orf.merge(joung_data, on="Name", how="left")


# In[21]:


dn_ref = dn[
    [
        "gene_symbol",
        "family",
        "clone_acc_ref",
        "is_ref_novel_isoform",
        "is_MANE_select_isoform_cloned",
        "dn_short",
    ]
].drop_duplicates()
dn_ref.columns = [
    "gene_name",
    "family",
    "tf1p0_id",
    "is_novel",
    "is_MANE_select",
    "dn_cat",
]
dn_ref["dn_cat"] = "ref"
dn_ref["iso_status"] = "ref"


# In[22]:


dn_alt = dn[
    [
        "gene_symbol",
        "family",
        "clone_acc_alt",
        "is_alt_novel_isoform",
        "is_MANE_select_isoform_cloned",
        "dn_short",
    ]
].drop_duplicates()
dn_alt.columns = [
    "gene_name",
    "family",
    "tf1p0_id",
    "is_novel",
    "is_MANE_select",
    "dn_cat",
]
dn_alt["is_MANE_select"] = False  # assuming none of the alts are the MANE select
dn_alt["iso_status"] = "alt"


# In[23]:


dn_cats = dn_ref.append(dn_alt).drop_duplicates()


# In[24]:


dn_cats = dn_cats.merge(joung_data, left_on="tf1p0_id", right_on="index", how="left")
dn_cats.sample(5)


# ## 4. count overlap between collections

# In[25]:


dn_cats.iso_status.value_counts()


# In[26]:


dn_cats[~pd.isnull(dn_cats["Name"])].iso_status.value_counts()


# In[27]:


refs_inc = len(
    dn_cats[(~pd.isnull(dn_cats["Name"])) & (dn_cats["iso_status"] == "ref")]
)
refs_tf1p0 = len(dn_cats[dn_cats["iso_status"] == "ref"])
print("%% of our ref seqs included in joung: %s" % (refs_inc / refs_tf1p0 * 100))


# In[28]:


alts_inc = len(
    dn_cats[(~pd.isnull(dn_cats["Name"])) & (dn_cats["iso_status"] == "alt")]
)
alts_tf1p0 = len(dn_cats[dn_cats["iso_status"] == "alt"])
print("%% of our alt seqs included in joung: %s" % (alts_inc / alts_tf1p0 * 100))


# In[29]:


missing_ref = dn_cats[(dn_cats["iso_status"] == "ref") & (pd.isnull(dn_cats["Name"]))]
missing_ref.is_novel.value_counts()


# In[30]:


missing_ref[~missing_ref["is_novel"]].is_MANE_select.value_counts(dropna=False)


# In[31]:


missing_ref[(~missing_ref["is_novel"]) & (missing_ref["is_MANE_select"] == True)]


# In[32]:


# a lot more transcripts are missing here since the seq_aa includes muts
# asked luke to see if he can add the ref seq_aa as an attribute to merge on


# In[33]:


missing_alt = dn_cats[(dn_cats["iso_status"] == "alt") & (pd.isnull(dn_cats["Name"]))]
missing_alt.is_novel.value_counts()


# In[34]:


dn_cats["orf_len"] = dn_cats["seq_aa_joung"].str.len()


# ## 5. join DN categories with cell counts

# In[35]:


joung_down_tf1p0_map = joung_down_map.merge(
    dn_cats[["TF ORF", "tf1p0_id", "iso_status", "dn_cat", "orf_len"]],
    left_on="TF",
    right_on="TF ORF",
)
print(len(joung_down_tf1p0_map))
print(len(joung_down_tf1p0_map["TF ORF"].unique()))


# In[36]:


joung_down_tf1p0_map.fillna("NA", inplace=True)


# In[37]:


joung_tf1p0_cnts = (
    joung_down_tf1p0_map.groupby(["TF", "tf1p0_id", "iso_status", "dn_cat", "orf_len"])[
        "TF ORF"
    ]
    .agg("count")
    .reset_index()
)
joung_tf1p0_cnts.columns = [
    "TF",
    "tf1p0_id",
    "iso_status",
    "dn_cat",
    "orf_len",
    "tot_cell_cnt",
]
joung_tf1p0_cnts.head()


# ## 6. create dataframe for linear model

# In[38]:


lm = dn_cats[
    [
        "gene_name",
        "family",
        "tf1p0_id",
        "is_novel",
        "dn_cat",
        "iso_status",
        "Diffusion P-value",
        "Diffusion difference",
        "orf_len",
    ]
]
lm.columns = [
    "gene_name",
    "family",
    "tf1p0_id",
    "is_novel",
    "dn_cat",
    "iso_status",
    "diff_pval",
    "diff_diff",
    "orf_len",
]
lm["dn_cat"].replace("ref", "*ref", inplace=True)
lm["iso_status"].replace("ref", "*ref", inplace=True)
lm["dn_cat"].fillna("NA", inplace=True)
lm = lm.merge(joung_tf1p0_cnts[["tf1p0_id", "tot_cell_cnt"]], on=["tf1p0_id"])
lm = lm[lm["dn_cat"] != "likely"]
print(len(lm))
lm.sample(5)


# In[39]:


lm.dn_cat.value_counts()


# In[40]:


# log transform counts and p-values
lm["neglog_diff_pval"] = -np.log10(lm["diff_pval"])
lm["log_count"] = np.log10(lm["tot_cell_cnt"])


# ## 7. run regression to control for cell count variable

# In[41]:


mod = smf.ols(formula="neglog_diff_pval ~ log_count + dn_cat", data=lm)


# In[42]:


res = mod.fit()


# In[43]:


print(res.summary())


# In[44]:


cls = Linear_Reg_Diagnostic(res)


# In[45]:


fig, ax = cls()
fig.savefig("../figures/Joung_Model_QC.pdf", dpi="figure", bbox_inches="tight")


# ## 8. plot ORFeome p-values across categories

# In[46]:


dn_cats_nonan = dn_cats[~pd.isnull(dn_cats["Diffusion P-value"])]
len(dn_cats_nonan)


# In[47]:


dn_cats_nonan["neglog_diff_pval"] = -np.log10(dn_cats_nonan["Diffusion P-value"])
dn_cats_nonan.fillna("NA", inplace=True)


# In[48]:


nice_boxplot(
    dn_cats_nonan,
    "neglog_diff_pval",
    "dn_cat",
    pal,
    ["ref", "rewire", "DN", "NA"],
    [6.9, 7.6, 8.3, 5.2],
    -0.02,
    "",
    ["reference", "rewirer", "negative regulator", "NA"],
    "-log10(Diffusion p-value)",
    False,
    (-0.75, 9.2),
    "effect of TF isoforms on differentiation\n(Joung et al.)",
    "../figures/Joung_DiffP_Boxplot.pdf",
)


# In[49]:


nice_boxplot(
    dn_cats_nonan,
    "RNA Velocity difference",
    "dn_cat",
    pal,
    ["ref", "rewire", "DN", "NA"],
    [0, 0, 0, 0],
    -0.1,
    "",
    ["reference", "rewirer", "negative regulator", "NA"],
    "RNA Velocity difference",
    False,
    (-0.05, 0.02),
    "effect of TF isoforms on differentiation\n(Joung et al.)",
    "../figures/Joung_DiffDiff_Boxplot.pdf",
)


# In[50]:


sns.lmplot(
    data=dn_cats_nonan,
    x="neglog_diff_pval",
    y="Diffusion difference",
    hue="dn_cat",
    palette=pal,
    col="dn_cat",
    height=2,
    scatter_kws={"s": 7},
)


# ## 9. plot p-values across cell count quartiles

# In[51]:


joung_tf1p0_cnts["cell_cnt_qcut"] = pd.qcut(
    joung_tf1p0_cnts["tot_cell_cnt"], q=4, labels=[1, 2, 3, 4]
)


# In[52]:


dn_cats_nonan = dn_cats_nonan.merge(
    joung_tf1p0_cnts[["TF", "tot_cell_cnt", "cell_cnt_qcut"]],
    left_on="TF ORF",
    right_on="TF",
)


# In[53]:


fig = plt.figure(figsize=(5, 3))

ax = sns.boxplot(
    data=dn_cats_nonan,
    x="cell_cnt_qcut",
    y="neglog_diff_pval",
    hue="dn_cat",
    palette=pal,
    hue_order=["ref", "rewire", "DN", "NA"],
    fliersize=0,
)

sns.swarmplot(
    data=dn_cats_nonan,
    x="cell_cnt_qcut",
    y="neglog_diff_pval",
    hue="dn_cat",
    palette=pal,
    hue_order=["ref", "rewire", "DN", "NA"],
    ax=ax,
    split=True,
    size=4,
    edgecolor="black",
    linewidth=0.5,
    alpha=0.5,
)

mimic_r_boxplot(ax)

ax.set_xlabel("Subsampled cell count quartile")
ax.set_ylabel("Joung et al. over-expression -log10(Diffusion p-value)")

plt.legend(loc=2, bbox_to_anchor=(1.01, 1), facecolor="white")

fig.savefig(
    "../figures/Joung_DiffP_Boxplot_QCut.pdf", dpi="figure", bbox_inches="tight"
)


# In[54]:


fig = plt.figure(figsize=(5, 3))

ax = sns.boxplot(
    data=dn_cats_nonan,
    x="cell_cnt_qcut",
    y="Diffusion difference",
    hue="dn_cat",
    palette=pal,
    hue_order=["ref", "rewire", "DN", "NA"],
    fliersize=0,
)

sns.swarmplot(
    data=dn_cats_nonan,
    x="cell_cnt_qcut",
    y="Diffusion difference",
    hue="dn_cat",
    palette=pal,
    hue_order=["ref", "rewire", "DN", "NA"],
    ax=ax,
    split=True,
    size=4,
    edgecolor="black",
    linewidth=0.5,
    alpha=0.5,
)

mimic_r_boxplot(ax)

ax.set_xlabel("Subsampled cell count quartile")
ax.set_ylabel("Joung et al. over-expression Diffusion difference")

plt.legend(loc=2, bbox_to_anchor=(1.01, 1), facecolor="white")

fig.savefig(
    "../figures/Joung_DiffDiff_Boxplot_QCut.pdf", dpi="figure", bbox_inches="tight"
)


# ## 10. scatter for ref v alt

# In[55]:


dn_cats_nonan_ref = dn_cats_nonan[dn_cats_nonan["iso_status"] == "ref"]
dn_cats_nonan_alt = dn_cats_nonan[dn_cats_nonan["iso_status"] == "alt"]
dn_cats_nonan_diff = dn_cats_nonan_ref.merge(
    dn_cats_nonan_alt,
    on=["gene_name", "family", "RefSeq Gene Name"],
    how="left",
    suffixes=("_ref", "_alt"),
)
dn_cats_nonan_diff["diff_pval_diff"] = (
    dn_cats_nonan_diff["Diffusion P-value_ref"]
    - dn_cats_nonan_diff["Diffusion P-value_alt"]
)
dn_cats_nonan_diff["diff_diff_diff"] = (
    dn_cats_nonan_diff["Diffusion difference_ref"]
    - dn_cats_nonan_diff["Diffusion difference_alt"]
)

dn_cats_nonan_diff["abs_ddd"] = np.abs(dn_cats_nonan_diff["diff_diff_diff"])


# In[56]:


fig = plt.figure(figsize=(2.4, 2.2))

ax = sns.scatterplot(
    data=dn_cats_nonan_diff[dn_cats_nonan_diff["dn_cat_alt"] != "likely"],
    x="neglog_diff_pval_ref",
    y="neglog_diff_pval_alt",
    hue="dn_cat_alt",
    palette=pal,
    linewidth=0.5,
    edgecolor="black",
    alpha=0.8,
    zorder=10,
)
ax.set_xlabel("reference isoform -log10(p-value)")
ax.set_ylabel("alternative isoform -log10(p-value)")
ax.set_title("effect of TF isoforms on differentiation\n(Joung et al.)")
# ax.set_xlim((0, 7.2))
# ax.set_ylim((0, 7.2))
ax.plot([0, 7.2], [0, 7.2], linestyle="dashed", color="black", linewidth=1, zorder=1)

plt.legend(loc=2, bbox_to_anchor=(1.01, 1))

fig.savefig("../figures/Joung_Scatter.inc_NA.pdf", dpi="figure", bbox_inches="tight")


# In[57]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(
    data=dn_cats_nonan_diff[dn_cats_nonan_diff["dn_cat_alt"].isin(["rewire", "DN"])],
    x="neglog_diff_pval_ref",
    y="neglog_diff_pval_alt",
    hue="dn_cat_alt",
    palette=pal,
    linewidth=0.5,
    edgecolor="black",
    alpha=0.8,
    zorder=10,
)
ax.set_xlabel("reference isoform -log10(p-value)")
ax.set_ylabel("alternative isoform -log10(p-value)")
ax.set_title("effect of TF isoforms on differentiation\n(Joung et al.)")

ax.set_xlim((0, 7))
ax.set_ylim((0, 7))
ax.plot([0, 7], [0, 7], linestyle="dashed", color="black", linewidth=1, zorder=1)

ax.get_legend().remove()

fig.savefig("../figures/Joung_Scatter.pdf", dpi="figure", bbox_inches="tight")


# In[58]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(
    data=dn_cats_nonan_diff[dn_cats_nonan_diff["dn_cat_alt"].isin(["rewire", "DN"])],
    x="Diffusion difference_ref",
    y="Diffusion difference_alt",
    hue="dn_cat_alt",
    palette=pal,
    linewidth=0.5,
    edgecolor="black",
    alpha=0.8,
    zorder=10,
)
ax.set_xlabel("reference isoform effect size")
ax.set_ylabel("alternative isoform effect size")
ax.set_title("effect of TF isoforms on differentiation\n(Joung et al.)")

ax.set_xlim((-0.04, 0.01))
ax.set_ylim((-0.03, 0.01))
ax.plot(
    [-0.0025, -0.0025],
    [-0.0025, 0.0075],
    linestyle="dashed",
    color="black",
    linewidth=0.5,
)
ax.plot(
    [-0.0025, 0.0075],
    [-0.0025, -0.0025],
    linestyle="dashed",
    color="black",
    linewidth=0.5,
)
ax.plot(
    [-0.0025, 0.0075],
    [0.0075, 0.0075],
    linestyle="dashed",
    color="black",
    linewidth=0.5,
)
ax.plot(
    [0.0075, 0.0075],
    [-0.0025, 0.0075],
    linestyle="dashed",
    color="black",
    linewidth=0.5,
)
# ax.axvline(x=0, linestyle="dashed", color="black", linewidth=0.5)

ax.get_legend().remove()

fig.savefig("../figures/Joung_Scatter_EffSize.pdf", dpi="figure", bbox_inches="tight")


# In[59]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(
    data=dn_cats_nonan_diff[dn_cats_nonan_diff["dn_cat_alt"].isin(["rewire", "DN"])],
    x="Diffusion difference_ref",
    y="Diffusion difference_alt",
    hue="dn_cat_alt",
    palette=pal,
    linewidth=0.5,
    edgecolor="black",
    alpha=0.8,
    zorder=10,
)
ax.set_xlabel("reference isoform effect size")
ax.set_ylabel("alternative isoform effect size")
ax.set_title("effect of TF isoforms on differentiation\n(Joung et al.)")

ax.set_xlim((-0.0025, 0.0075))
ax.set_ylim((-0.0025, 0.0075))
ax.plot(
    [-0.0025, 0.0075],
    [-0.0025, 0.0075],
    linestyle="dashed",
    color="black",
    linewidth=1,
    zorder=1,
)

ax.get_legend().remove()

fig.savefig(
    "../figures/Joung_Scatter_EffSize_Zoom.pdf", dpi="figure", bbox_inches="tight"
)


# In[60]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(
    data=dn_cats_nonan_diff[
        (dn_cats_nonan_diff["dn_cat_alt"].isin(["rewire", "DN"]))
        & (dn_cats_nonan_diff["Diffusion P-value_alt"] < 0.1)
    ],
    x="Diffusion difference_ref",
    y="Diffusion difference_alt",
    hue="dn_cat_alt",
    palette=pal,
    linewidth=0.5,
    edgecolor="black",
    alpha=0.8,
    zorder=10,
)
ax.set_xlabel("reference isoform effect size")
ax.set_ylabel("alternative isoform effect size")
ax.set_title("effect of TF isoforms on differentiation\n(Joung et al.)")

ax.set_xlim((-0.04, 0.01))
ax.set_ylim((-0.03, 0.01))
ax.axhline(y=0, linestyle="dashed", color="black", linewidth=0.5)
ax.axvline(x=0, linestyle="dashed", color="black", linewidth=0.5)

ax.get_legend().remove()

# fig.savefig("../figures/Joung_Scatter.pdf", dpi="figure", bbox_inches="tight")


# In[61]:


tmp = dn_cats_nonan_diff[dn_cats_nonan_diff["dn_cat_alt"].isin(["rewire", "DN"])]
tmp.sort_values(by="Diffusion difference_alt", ascending=True)[
    [
        "tf1p0_id_ref",
        "tf1p0_id_alt",
        "dn_cat_alt",
        "Name_ref",
        "Name_alt",
        "Diffusion difference_ref",
        "Diffusion difference_alt",
        "neglog_diff_pval_ref",
        "neglog_diff_pval_alt",
        "tot_cell_cnt_ref",
        "tot_cell_cnt_alt",
    ]
].head(20)


# In[62]:


fig = plt.figure(figsize=(2, 2.2))

ax = sns.scatterplot(
    data=dn_cats_nonan[dn_cats_nonan["dn_cat"].isin(["ref", "rewire", "DN"])],
    x="Diffusion difference",
    y="neglog_diff_pval",
    hue="dn_cat",
    palette=pal,
    linewidth=0.25,
    edgecolor="black",
    alpha=0.8,
    zorder=10,
    **{"s": 8}
)

ax.set_xlabel("over-expression effect size")
ax.set_ylabel("-log10(over-expression p-value)")
ax.set_title("effect of TF isoforms on differentiation\n(Joung et al.)")

ax.set_xlim((-0.04, 0.01))
ax.set_ylim((-0.01, 7.2))
ax.axhline(y=-np.log10(0.05), linestyle="dashed", color="black", linewidth=0.5)
ax.axvline(x=0, linestyle="dashed", color="black", linewidth=0.5)

ax.get_legend().remove()

fig.savefig("../figures/Joung_Volcano.pdf", dpi="figure", bbox_inches="tight")


# In[63]:


tmp = dn_cats_nonan[dn_cats_nonan["dn_cat"] == "DN"]
# tmp = tmp[tmp["neglog_diff_pval"] > 3]
tmp.sort_values(by="neglog_diff_pval", ascending=False)[
    ["tf1p0_id", "dn_cat", "Diffusion difference", "neglog_diff_pval", "tot_cell_cnt"]
].head()


# ## 11. enrichment heatmaps

# In[64]:


joung_cells["Name"] = joung_cells["TF"].str.strip().str.split("-", expand=True)[0]
print(len(joung_cells))

# # filter out anything with score < 0.2
joung_cells = joung_cells[joung_cells["prediction.score.max"] > 0.2]
print(len(joung_cells))


# In[65]:


joung_cells_grp = (
    joung_cells.groupby(["Name", "TF", "predicted.id"])["batch"]
    .agg("count")
    .reset_index()
)


# In[66]:


tot_cell_cnt = dn_cats_nonan[["Name", "TF", "tot_cell_cnt"]].drop_duplicates()
diff_cell_cnt = joung_cells.groupby(["Name", "TF"])["batch"].agg("count").reset_index()
cell_cnt = tot_cell_cnt.merge(diff_cell_cnt, on=["Name", "TF"], how="left")
cell_cnt.fillna(0, inplace=True)
cell_cnt.columns = ["Name", "TF", "tot_cell_cnt", "diff_cell_cnt"]

orf_enr = cell_cnt.merge(joung_cells_grp, on=["Name", "TF"], how="left")
orf_enr["batch"].fillna(0, inplace=True)

orf_enr.columns = [
    "Name",
    "TF",
    "tot_cell_cnt",
    "diff_cell_cnt",
    "predicted.id",
    "id_cell_cnt",
]
orf_enr["perc_cells_of_diff_tf"] = orf_enr["id_cell_cnt"] / orf_enr["diff_cell_cnt"]
orf_enr["perc_cells_of_tot_tf"] = orf_enr["id_cell_cnt"] / orf_enr["tot_cell_cnt"]


# In[67]:


orf_enr_dn = orf_enr.merge(
    dn_cats_nonan[["gene_name", "Name", "tf1p0_id", "dn_cat"]], on="Name"
).drop_duplicates(subset=["tf1p0_id", "predicted.id", "dn_cat"])


# In[68]:


has_alt = list(orf_enr_dn[orf_enr_dn["dn_cat"] != "ref"]["gene_name"].unique())
orf_enr_dn_filt = orf_enr_dn[orf_enr_dn["gene_name"].isin(has_alt)]

has_ref = list(
    orf_enr_dn_filt[orf_enr_dn_filt["dn_cat"] == "ref"]["gene_name"].unique()
)
orf_enr_dn_filt = orf_enr_dn_filt[orf_enr_dn_filt["gene_name"].isin(has_ref)]
len(orf_enr_dn_filt)


# In[69]:


orf_enr_dn_filt["dn_cat_s"] = pd.Categorical(
    orf_enr_dn_filt["dn_cat"], ["ref", "rewire", "DN", "NA", "likely"]
)
orf_enr_dn_filt = orf_enr_dn_filt.sort_values(by=["gene_name", "dn_cat_s"])
orf_enr_dn_filt[orf_enr_dn_filt["gene_name"] == "PBX1"]


# In[70]:


cell_cnt["undiff_cell_cnt"] = cell_cnt["tot_cell_cnt"] - cell_cnt["diff_cell_cnt"]
len(cell_cnt)


# In[71]:


cell_cnt["diff_cell_perc"] = (
    cell_cnt["diff_cell_cnt"] / cell_cnt["tot_cell_cnt"]
) * 100
cell_cnt["undiff_cell_perc"] = (
    cell_cnt["undiff_cell_cnt"] / cell_cnt["tot_cell_cnt"]
) * 100
cell_cnt[cell_cnt["TF"].str.contains("GRHL3")]


# In[72]:


tmp = orf_enr_dn_filt[orf_enr_dn_filt["tf1p0_id"].str.contains("KLF7")].pivot(
    index="tf1p0_id", columns="predicted.id", values="perc_cells_of_tot_tf"
)
tmp.fillna(0, inplace=True)

idx = pd.DataFrame(tmp.index)
idx = (
    idx.merge(orf_enr_dn_filt[["tf1p0_id", "dn_cat"]], on="tf1p0_id")
    .drop_duplicates()
    .set_index("tf1p0_id")
)
idx["Isoform type"] = idx.dn_cat.map(pal)

g = sns.clustermap(
    tmp,
    cmap="Greys",
    row_cluster=False,
    row_colors=idx["Isoform type"],
    figsize=(4, 2.5),
    yticklabels=True,
    cbar_pos=(0, 1, 0.05, 0.2),
)
g.savefig("../figures/Joung_KLF7_hm.pdf", bbox_inches="tight", dpi="figure")


# In[73]:


tmp = orf_enr_dn_filt[orf_enr_dn_filt["tf1p0_id"].str.contains("PKNOX1")].pivot(
    index="tf1p0_id", columns="predicted.id", values="perc_cells_of_tot_tf"
)
tmp.fillna(0, inplace=True)

idx = pd.DataFrame(tmp.index)
idx = (
    idx.merge(orf_enr_dn_filt[["tf1p0_id", "dn_cat"]], on="tf1p0_id")
    .drop_duplicates()
    .set_index("tf1p0_id")
)
idx["Isoform type"] = idx.dn_cat.map(pal)

g = sns.clustermap(
    tmp,
    cmap="Greys",
    row_cluster=False,
    row_colors=idx["Isoform type"],
    figsize=(2, 0.75),
    yticklabels=True,
    cbar_pos=(0, 1, 0.05, 0.2),
    annot=True,
)
g.savefig("../figures/Joung_PKNOX1_hm.pdf", bbox_inches="tight", dpi="figure")


# In[74]:


cell_cnt[cell_cnt["TF"].str.contains("PKNOX1")]


# In[75]:


tmp = orf_enr_dn_filt[orf_enr_dn_filt["tf1p0_id"].str.contains("PBX1")].pivot(
    index="tf1p0_id", columns="predicted.id", values="perc_cells_of_tot_tf"
)
tmp.fillna(0, inplace=True)

idx = pd.DataFrame(tmp.index)
idx = (
    idx.merge(orf_enr_dn_filt[["tf1p0_id", "dn_cat"]], on="tf1p0_id")
    .drop_duplicates()
    .set_index("tf1p0_id")
)
idx["Isoform type"] = idx.dn_cat.map(pal)

g = sns.clustermap(
    tmp,
    cmap="Greys",
    row_cluster=False,
    row_colors=idx["Isoform type"],
    figsize=(4, 2.5),
    yticklabels=True,
    cbar_pos=(0, 1, 0.05, 0.2),
)
# g.savefig("../figures/Joung_HOXA1_hm.pdf", bbox_inches="tight", dpi="figure")


# In[76]:


tmp = orf_enr_dn_filt[orf_enr_dn_filt["tf1p0_id"].str.contains("GRHL3")].pivot(
    index="tf1p0_id", columns="predicted.id", values="perc_cells_of_tot_tf"
)
tmp.drop(np.nan, axis=1, inplace=True)
tmp = tmp.loc[
    ["GRHL3|3/7|08G09", "GRHL3|4/7|08F09", "GRHL3|1/7|08E10", "GRHL3|2/7|08A10"]
]
tmp.fillna(0, inplace=True)

idx = pd.DataFrame(tmp.index)
idx = (
    idx.merge(orf_enr_dn_filt[["tf1p0_id", "dn_cat"]], on="tf1p0_id")
    .drop_duplicates()
    .set_index("tf1p0_id")
)
idx["Isoform type"] = idx.dn_cat.map(pal)

g = sns.clustermap(
    tmp,
    cmap="Greys",
    row_cluster=False,
    row_colors=idx["Isoform type"],
    figsize=(4, 1),
    yticklabels=True,
    cbar_pos=(0, 1, 0.05, 0.2),
)

g.savefig("../figures/Joung_GRHL3_hm.pdf", bbox_inches="tight", dpi="figure")


# ## 12. write files

# In[77]:


joung_orf.to_csv("../data/processed/Joung_ORF_annotated.tsv", sep="\t", index=False)


# In[78]:


dn_cats.to_csv("../data/processed/DN_cats_Joung.tsv", sep="\t", index=False)


# In[81]:


joung_orf[joung_orf["RefSeq Gene Name"] == "GRHL3"]


# In[118]:


joung_orf[joung_orf["RefSeq Gene Name"] == "TBX5"]


# In[122]:


tfs["TBX5"]["TBX5|1/3|08E01"].aa_seq


# In[123]:


len(tfs["TBX5"]["TBX5|1/3|08E01"].aa_seq)


# In[125]:


tfs["TBX5"]["TBX5|2/3|08C02"].aa_seq


# In[126]:


len(tfs["TBX5"]["TBX5|2/3|08C02"].aa_seq)


# In[128]:


tfs["TBX5"].isoforms
