import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import patches
import pandas as pd
from scipy import stats
import seaborn as sns

from data_loading import paralog_pair_ppi_table
from scipy.stats import mannwhitneyu


COLOR_PURPLE = (155 / 255, 97 / 255, 153 / 255)

# kaia's added code
PAPER_PRESET = {
    "style": "ticks",
    "font": "Helvetica",
    "context": "paper",
    "rc": {
        "font.size": 7,
        "pdf.fonttype": 42,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "axes.linewidth": 0.5,
        "legend.fontsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "axes.edgecolor": "black",
        "xtick.major.pad": 3.0,
        "ytick.major.pad": 3.0,
    },
}
PAPER_FONTSIZE = 7


def violinplot_reflected(*args, lb=0, ub=1, bw_const=0.1, **kwargs):
    """
    monkeypatch from https://github.com/mwaskom/seaborn/issues/525

    Also use a constant bandwidth

    """
    fit_kde_func = sns.categorical._ViolinPlotter.fit_kde

    def reflected_once_kde(self, x, bw):
        kde, bw_used = fit_kde_func(self, x, bw)

        kde_evaluate = kde.evaluate

        def truncated_kde_evaluate(x):
            val = np.where((x >= lb) & (x <= ub), kde_evaluate(x), 0)
            val += np.where((x >= lb) & (x <= ub), kde_evaluate(lb - x), 0)
            val += np.where((x >= lb) & (x <= ub), kde_evaluate(ub - (x - ub)), 0)
            return val

        kde.evaluate = truncated_kde_evaluate
        return kde, bw_used

    def constant_bw(kde):
        """Have to scale by STD because the returned value is multiplied
        separately by the STD of each group."""
        return bw_const / np.std(kde.dataset)

    sns.categorical._ViolinPlotter.fit_kde = reflected_once_kde
    retval = sns.violinplot(*args, bw=constant_bw, **kwargs)
    sns.categorical._ViolinPlotter.fit_kde = fit_kde_func  # change back
    return retval


def binary_profile_matrix(
    data,
    ax=None,
    box_size=0.7,
    fill_color="black",
    border_color="black",
    column_label_rotation=40,
    bait_colors=None,
    bait_annot=None,
):
    """Used for edgotyping: displays binary results with a grid of boxes

    Empty box for negative, filled box for positive and
    missing box for missing values.

    (Copied over from ccsblib)

    Args:
        data (pandas.DataFrame): boolean values
        ax (matplotlib.axes.Axes): axes to draw on
        box_size (float): area of the boxes between 0 and 1
        fill_color (str): color of filled sqaure
        border_color (str): color of outside of square
        column_label_roation: angle in degrees of top labels
        bait_colors: optional, added to allow for custom colors of squares in PDI figure
        bait_annot: optional, added to allow for custom annotation of squares in PDI figure
    Examples:
        Display results of some dummy edgotyping data:
        .. plot::
            :context: close-figs
            >>> import pandas as pd
            >>> df = pd.DataFrame(index=['Allele A', 'Allele B', 'Allele C'],
            ...                   columns=['Partner W', 'Partner X', 'Parner Y', 'Partner Z'],
            ...                   data=[[True] * 4,
            ...                         [True, False, True, False],
            ...                         [False, False, np.nan, False]])
            >>> binary_profile_matrix(df)
        You can switch the rows and columns by transposing the input DataFrame:
        .. plot::
            :context: close-figs
             >>> binary_profile_matrix(df.T,
             ...                            fill_color='grey',
             ...                            border_color='black',
             ...                            box_size=0.6,
             ...                            column_label_rotation=90)
    """
    if box_size > 1.0 or box_size < 0.0:
        raise ValueError("box_size must be between 0-1")
    if ax is None:
        ax = plt.gca()
    ax.set_aspect("equal")
    positives = [
        (i, j)
        for i in range(data.shape[1])
        for j in range(data.shape[0])
        if pd.notnull(data.iloc[j, i]) and data.iloc[j, i] == 1
    ]
    negatives = [
        (i, j)
        for i in range(data.shape[1])
        for j in range(data.shape[0])
        if pd.notnull(data.iloc[j, i]) and data.iloc[j, i] == 0
    ]
    for x, y in negatives:
        if bait_colors is None:
            neg_fill = False
            facecolor = None
        else:
            neg_fill = True
            facecolor = bait_colors.iloc[y, x]

        r = patches.Rectangle(
            (x - box_size / 2, y - box_size / 2),
            box_size,
            box_size,
            fill=neg_fill,
            facecolor=facecolor,
            edgecolor=border_color,
            linewidth=1,
        )
        ax.add_patch(r)

        if bait_annot is not None:
            annot = bait_annot.iloc[y, x]
            ax.text(
                x,
                y,
                "{:.2f}".format(annot),
                fontsize=6,
                color="black",
                ha="center",
                va="center",
            )

    for x, y in positives:
        if bait_colors is None:
            facecolor = fill_color
            linewidth = 1
        else:
            facecolor = bait_colors.iloc[y, x]
            linewidth = 3

        r = patches.Rectangle(
            (x - box_size / 2, y - box_size / 2),
            box_size,
            box_size,
            fill=True,
            facecolor=facecolor,
            edgecolor=border_color,
            linewidth=linewidth,
        )
        ax.add_patch(r)

        if bait_annot is not None:
            annot = bait_annot.iloc[y, x]
            if annot != "NA":
                ax.text(
                    x,
                    y,
                    "{:.2f}".format(annot),
                    fontsize=6,
                    color="black",
                    ha="center",
                    va="center",
                )

    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.tick_top()
    ax.set_xticks(range(data.shape[1]))
    ax.xaxis.set_tick_params(length=0)
    ax.xaxis.set_ticklabels(data.columns, rotation=column_label_rotation, ha="center")
    ax.set_yticks(range(data.shape[0]))
    ax.yaxis.set_tick_params(length=0)
    ax.set_yticklabels(data.index)
    ax.set_ylim((-0.5, data.shape[0] - 0.5))
    ax.set_xlim((-0.5, data.shape[1] - 0.5))
    ax.invert_yaxis()


def isoform_display_name(s):
    """Convert clone accession ID to display friendly format"""
    return s.split("|")[0] + "-" + s.split("|")[1].split("/")[0]


def strikethrough(s):
    return "".join([c + "\u0336" for c in s])


def y2h_ppi_per_tf_gene_plot(
    gene_name,
    data,
    ax=None,
    min_n_isoforms=1,
    min_n_partners=1,
    iso_order=None,
):
    tf = data.loc[
        (data["ad_gene_symbol"] == gene_name),
        ["ad_clone_acc", "db_gene_symbol", "Y2H_result"],
    ].copy()
    tf["ad_clone_acc"] = tf["ad_clone_acc"].apply(isoform_display_name)
    tf = tf.pivot(index="ad_clone_acc", columns="db_gene_symbol", values="Y2H_result")
    # remove partners with no positives
    tf = tf.loc[:, tf.any(axis=0)]
    if ax is None:
        ax = plt.gca()
    if tf.shape[0] < min_n_isoforms or tf.shape[1] < min_n_partners:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No PPI data available",
            ha="center",
            va="center",
            fontsize=30,
            fontweight="bold",
            color="grey",
        )
        return False
    if iso_order is None:
        tf = tf
    else:
        tf = tf.loc[iso_order, :]
    binary_profile_matrix(tf, ax=ax, column_label_rotation=90)
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )
    return True


def y2h_ppi_per_paralog_pair_plot(
    tf_gene_a,
    tf_gene_b,
    data,
    ax=None,
    min_n_isoforms=1,
    min_n_partners=1,
    iso_order=None,
):
    """

    TODO: gap between the two genes?

    Arguments:
        tf_gene_a {str} -- [description]
        tf_gene_b {str} -- [description]
        data {pandas.DataFrame} -- [description]
        ax {pandas.DataFrame} -- [description] (default: {None})
        min_n_isoforms {int} -- [description] (default: {1})
        min_n_partners {int} -- [description] (default: {1})

    """
    tf = paralog_pair_ppi_table(data, tf_gene_a, tf_gene_b)
    tf["ad_clone_acc"] = tf["ad_clone_acc"].apply(isoform_display_name)
    tf = tf.pivot(index="ad_clone_acc", columns="db_gene_symbol", values="Y2H_result")
    if ax is None:
        ax = plt.gca()
    if tf.shape[0] < min_n_isoforms or tf.shape[1] < min_n_partners:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No PPI data available",
            ha="center",
            va="center",
            fontsize=30,
            fontweight="bold",
            color="grey",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    if iso_order is not None:
        tf = tf.loc[iso_order, :]
    binary_profile_matrix(tf, ax=ax, column_label_rotation=90)
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )


def y1h_pdi_per_tf_gene_plot(
    gene_name,
    data,
    ax=None,
    min_n_isoforms=1,
    min_n_partners=1,
    iso_order=None,
    bait_colors=None,
    bait_annot=None,
):
    tf = (
        data.loc[data["gene_symbol"] == gene_name, data.columns[1:]]
        .copy()
        .set_index("clone_acc")
    )
    tf.index = tf.index.map(isoform_display_name)
    tf = tf.loc[:, tf.any(axis=0)]
    if ax is None:
        ax = plt.gca()
    if tf.shape[0] < min_n_isoforms or tf.shape[1] < min_n_partners:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No PDI data available",
            ha="center",
            va="center",
            fontsize=30,
            fontweight="bold",
            color="grey",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return False
    if iso_order is None:
        tf = tf
    else:
        tf = tf.loc[iso_order, :]
    binary_profile_matrix(
        tf,
        ax=ax,
        column_label_rotation=90,
        bait_colors=bait_colors,
        bait_annot=bait_annot,
    )
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )
    return True


def y1h_pdi_per_paralog_pair_plot(
    gene_name_a,
    gene_name_b,
    data,
    ax=None,
    min_n_isoforms=1,
    min_n_partners=1,
    iso_order=None,
    bait_colors=None,
    bait_annot=None,
):
    tf = (
        data.loc[data["gene_symbol"].isin([gene_name_a, gene_name_b]), data.columns[1:]]
        .copy()
        .set_index("clone_acc")
    )
    tf.index = tf.index.map(isoform_display_name)
    tf = tf.loc[:, tf.any(axis=0)]
    if ax is None:
        ax = plt.gca()
    if tf.shape[0] < min_n_isoforms or tf.shape[1] < min_n_partners:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No PDI data available",
            ha="center",
            va="center",
            fontsize=30,
            fontweight="bold",
            color="grey",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return False
    if iso_order is not None:
        tf = tf.loc[iso_order, :]
    binary_profile_matrix(
        tf,
        ax=ax,
        column_label_rotation=90,
        bait_colors=bait_colors,
        bait_annot=bait_annot,
    )
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )
    return True


def m1h_activation_per_tf_gene_plot(
    tf_gene_name, data, ax=None, xlim=None, iso_order=None
):
    if ax is None:
        ax = plt.gca()
    rep_columns = [c for c in data.columns if c.startswith("M1H_rep")]
    is_all_na = (
        data[rep_columns]
        .isnull()
        .groupby(data["gene_symbol"])
        .all()
        .all(axis=1)[tf_gene_name]
    )
    if tf_gene_name not in data["gene_symbol"].values or is_all_na:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No activation data available",
            ha="center",
            va="center",
            fontsize=30,
            fontweight="bold",
            color="grey",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return False

    clones = [
        isoform_display_name(acc)
        for acc in data.loc[data["gene_symbol"] == tf_gene_name, "clone_acc"].values
        for __ in range(len(rep_columns))
    ]
    values = data.loc[data["gene_symbol"] == tf_gene_name, rep_columns].values.flatten()
    n_reps = len(rep_columns)

    if iso_order is None:
        y = clones[::n_reps]
        width = (
            data.loc[data["gene_symbol"] == tf_gene_name, rep_columns]
            .mean(axis=1)
            .values
        )
    else:
        y = iso_order
        width = []
        for iso in iso_order:
            idxs = [i for i, x in enumerate(clones) if x == iso]
            width.append(np.mean(values[idxs]))

    ax.barh(
        y=y,
        width=width,
        edgecolor="black",
        color="slategrey",
        alpha=0.5,
        height=0.6,
    )

    # swarmplot to make points more visible
    df = pd.DataFrame()
    df["clone"] = clones
    df["value"] = values
    if iso_order is None:
        order = y
    else:
        order = iso_order

    sns.stripplot(
        data=df,
        x="value",
        y="clone",
        ax=ax,
        color="white",
        linewidth=1,
        edgecolor="black",
        size=4,
        order=order,
    )
    ax.set_ylabel("")

    ax.set_yticks(
        clones[::n_reps]
    )  # needed to avoid truncating clones with missing data
    ax.set_yticklabels(
        [
            c if not np.isnan(v) else strikethrough(c)
            for c, v in zip(clones[::n_reps], values[::n_reps])
        ]
    )

    if xlim is None:
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(min([-3, xmin]), max([3, xmax]))
    else:
        ax.set_xlim(xlim)
    ax.set_ylim(-0.5, len(clones[::n_reps]) - 0.5)
    ax.set_xlabel("log2(activation fold change)")
    ax.axvline(0, linestyle="-", color="black")
    ax.axvline(-1, linestyle="--", color="black")
    ax.axvline(1, linestyle="--", color="black")
    ax.invert_yaxis()
    for pos in ["top", "left", "right"]:
        ax.spines[pos].set_visible(False)
    ax.yaxis.set_tick_params(length=0)
    return True


def m1h_activation_per_paralog_pair_plot(
    gene_name_a, gene_name_b, data, ax=None, xlim=None, iso_order=None
):
    tf = (
        data.loc[data["gene_symbol"].isin([gene_name_a, gene_name_b]), data.columns[1:]]
        .copy()
        .set_index("clone_acc")
    )

    if ax is None:
        ax = plt.gca()
    rep_columns = [c for c in data.columns if c.startswith("M1H_rep")]

    is_all_na = (
        data[rep_columns]
        .isnull()
        .groupby(data["gene_symbol"])
        .all()
        .all(axis=1)[[gene_name_a, gene_name_b]]
        .sum()
        .astype(bool)
    )

    if (
        gene_name_a not in data["gene_symbol"].values
        and gene_name_b not in data["gene_symbol"].values
    ) or is_all_na:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No activation data available",
            ha="center",
            va="center",
            fontsize=30,
            fontweight="bold",
            color="grey",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return False

    clones_a = [
        isoform_display_name(acc)
        for acc in data.loc[data["gene_symbol"] == gene_name_a, "clone_acc"].values
        for __ in range(len(rep_columns))
    ]
    values_a = data.loc[
        data["gene_symbol"] == gene_name_a, rep_columns
    ].values.flatten()

    clones_b = [
        isoform_display_name(acc)
        for acc in data.loc[data["gene_symbol"] == gene_name_b, "clone_acc"].values
        for __ in range(len(rep_columns))
    ]
    values_b = data.loc[
        data["gene_symbol"] == gene_name_b, rep_columns
    ].values.flatten()

    n_reps = len(rep_columns)

    clones = clones_a + clones_b
    values = np.asarray(list(values_a) + list(values_b))

    if iso_order is None:
        y = clones[::n_reps]
        width = (
            data.loc[data["gene_symbol"].isin([gene_name_a, gene_name_b]), rep_columns]
            .mean(axis=1)
            .values
        )
    else:
        y = iso_order
        width = []
        for iso in iso_order:
            idxs = [i for i, x in enumerate(clones) if x == iso]
            width.append(np.mean(values[idxs]))

    ax.barh(
        y=y,
        width=width,
        edgecolor="black",
        color="slategrey",
        alpha=0.5,
        height=0.6,
    )

    # swarmplot to make points more visible
    df = pd.DataFrame()
    df["clone"] = clones
    df["value"] = values
    if iso_order is None:
        order = y
    else:
        order = iso_order

    sns.stripplot(
        data=df,
        x="value",
        y="clone",
        ax=ax,
        color="white",
        linewidth=1,
        edgecolor="black",
        size=4,
        order=order,
    )
    ax.set_ylabel("")

    ax.set_yticks(
        clones[::n_reps]
    )  # needed to avoid truncating clones with missing data
    ax.set_yticklabels(
        [
            c if not np.isnan(v) else strikethrough(c)
            for c, v in zip(clones[::n_reps], values[::n_reps])
        ]
    )

    if xlim is None:
        ax.set_xlim(-3, 12)
    else:
        ax.set_xlim(xlim)
    ax.set_ylim(-0.5, len(clones[::n_reps]) - 0.5)
    ax.set_xlabel("Log2 M1H readout")
    ax.axvline(0, linestyle="-", color="black")
    ax.axvline(-1, linestyle="--", color="black")
    ax.axvline(1, linestyle="--", color="black")
    ax.invert_yaxis()
    for pos in ["top", "left", "right"]:
        ax.spines[pos].set_visible(False)
    ax.yaxis.set_tick_params(length=0)
    return True


def validation_plot(
    positives=None,
    n_tested=None,
    data=None,
    selections=None,
    result_column="result",
    labels=None,
    colors=None,
    ax=None,
    bayes_errors=True,
    y_max=1.0,
    draw_numbers=True,
    xlabel_rotation=0,
    errorbar_capsize=0.9,
    errorbar_thickness=None,
    bar_spacing=0.2,
    fontsize=10,
):
    """Compare the validation rate of different cateogies.

    Missing values are not used in the denominator.

    See the :ref:`tutorial </validation_plots.ipynb>` for more information.

    Args:
        positives (list(int)): number tested positive in each category
        n_tested (list(int)): number successfully tested in each category
        data (pandas.DataFrame): results of validation experiment
        selections (list(pandas.Series)): boolean rows index for each category
        ax (matplotlib.axes.Axes): Axes to draw plot onto
        bayes_errors (bool): do Bayesian error bars, if false use standard error
                            on proportion
        result_column (str): column containing 0/1/nan for result of test
        labels (list(str)): name of each category
        colors (list(str)): color for each bar
        y_max (float): y axis upper limit
        draw_numbers (bool): flag to print the numbers on top of the bars
        xlabel_rotation (float): rotating the x axis labels
        errorbar_capsize (float): as fraction of the width of the bars
        errorbar_thickness (float): width of error bar lines
        bar_spacing (float): must be between 0 and 1

    Examples:
        There are two ways to call the function. Either give it the raw data:



        .. plot::
            :context: close-figs

            >>> from plotting import validation_plot
            >>> validation_plot(positives=[20, 1, 19],
            ...                 n_tested=[100, 100, 100],
            ...                 labels=['PRS', 'RRS', 'Y2H'],
            ...                 y_max=0.3)

        Or pass it the validation results as a DataFrame and a list of the rows
        for each category:

        .. plot::
            :context: close-figs

            >>> from ccsblib import huri
            >>> data = huri.load_validation_data()
            >>> exp = (data['assay'] == 'MAPPIT') & (data['standard_batch'] == 'Hvs01')
            >>> sources = ['lit_bm_2013_rand250', 'Hs01', 'RRS']
            >>> categories = [exp & (data['source'] == cat) for cat in sources]
            >>> validation_plot(data=data,
            ...                 selections=categories,
            ...                 labels=sources,
            ...                 y_max=0.25)

    """
    signature_a = positives is not None and n_tested is not None
    signature_b = data is not None and selections is not None
    if signature_a == signature_b:
        msg = """Must supply only one of both positives and n_tested or
                 both data and selections."""
        raise ValueError(msg)
    if signature_b:
        if (
            not data.loc[data[result_column].notnull(), result_column]
            .isin({0, 1})
            .all()
        ):
            raise ValueError("Only expect 0/1/missing in result column")
        positives = [
            (data.loc[rows, :][result_column] == 1).sum() for rows in selections
        ]
        n_tested = [
            data.loc[rows, :][result_column].notnull().sum() for rows in selections
        ]
    _validation_plot(
        positives=positives,
        n_tested=n_tested,
        colors=colors,
        labels=labels,
        ax=ax,
        bayes_errors=bayes_errors,
        y_max=y_max,
        draw_numbers=draw_numbers,
        xlabel_rotation=xlabel_rotation,
        errorbar_capsize=errorbar_capsize,
        errorbar_thickness=errorbar_thickness,
        bar_spacing=bar_spacing,
        fontsize=fontsize,
    )


def _validation_plot(
    positives,
    n_tested,
    colors=None,
    labels=None,
    ax=None,
    bayes_errors=True,
    y_max=1.0,
    draw_numbers=True,
    xlabel_rotation=0.0,
    errorbar_capsize=5,
    errorbar_thickness=None,
    bar_spacing=0.2,
    fontsize=10,
):
    if len(positives) != len(n_tested):
        raise ValueError("Lengths of positives and n_tested must be equal")
    if any([p > n for p, n in zip(positives, n_tested)]):
        raise ValueError("Number of positives must be <= number tested")
    if bar_spacing > 1.0 or bar_spacing < 0.0:
        msg = "bar_spacing={}\nbar_spacing must be between 0 and 1"
        msg = msg.format(bar_spacing)
        raise ValueError(msg)
    bar_width = 1.0 - bar_spacing
    if ax is None:
        ax = plt.gca()
    if labels is None:
        labels = [""] * len(positives)
    if colors is None:
        colors = [None] * len(positives)
    ax.set_yticks(np.arange(0.0, 1.0, 0.1), minor=False)
    ax.set_yticks(np.arange(0.05, 1.0, 0.1), minor=True)
    # ax.set_facecolor("0.96")
    ax.set_axisbelow(True)
    ax.grid(color="white", axis="y", which="both", zorder=5)
    pos = np.array(positives)
    tested = np.array(n_tested)
    neg = tested - pos
    fracs = pos / tested
    if bayes_errors:
        intv = stats.beta.interval(0.6827, pos + 1, neg + 1)
        errs = [fracs - intv[0], intv[1] - fracs]
        errs[0][pos == 0] = 0.0
        errs[1][neg == 0] = 0.0
    else:
        stdErrProp = np.sqrt((fracs * (1.0 - fracs)) / (pos + neg))
        errs = [stdErrProp, stdErrProp]
    for i in range(len(positives)):
        ax.bar(i, fracs[i], color=colors[i], label=labels[i], width=bar_width)
        if draw_numbers:
            c = "white"
            h = 0.02  # default height to draw numbers
            if fracs[i] < h:
                c = "black"
            if (errs[1][i] + fracs[i]) > h and (fracs[i] - errs[0][i]) < (h + 0.04):
                c = "black"
                h = fracs[i] + errs[1][i] + 0.02
            ax.text(
                i,
                h,
                "{}/{}".format(pos[i], pos[i] + neg[i]),
                color=c,
                ha="center",
                fontsize=fontsize,
            )
    bar_width_pixels = (
        ax.transData.transform((bar_width, 0)) - ax.transData.transform((0, 0))
    )[0]
    # 72 comes from definition of a point as 1 / 72 inches
    bar_width_points = (72.0 / ax.figure.dpi) * bar_width_pixels
    ax.errorbar(
        range(fracs.shape[0]),
        fracs,
        yerr=errs,
        color="black",
        fmt="none",
        # I don't understand why I needed the factor of 0.5 below
        capsize=bar_width_points * errorbar_capsize * 0.5,
        elinewidth=errorbar_thickness,
        capthick=errorbar_thickness,
    )
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=xlabel_rotation)
    ax.set_ylim((0.0, y_max))
    ax.set_ylabel("Fraction positive")


def validation_titration_plot(
    data,
    selections,
    threshold=None,
    xmin=None,
    xmax=None,
    ymax=None,
    score_column="score",
    labels=None,
    colors=None,
    line_styles=None,
    ax=None,
    threshold_label=None,
    threshold_color="grey",
    plot_kwargs=None,
):
    """
    - error bars

    """
    if ax is None:
        ax = plt.gca()
    if labels is None:
        labels = [""] * len(selections)
    if colors is None:
        colors = [None] * len(selections)
    if line_styles is None:
        line_styles = ["-"] * len(selections)
    if xmin is None:
        xmin = data[score_column].min()
    if xmax is None:
        xmax = data[score_column].max()
    n_points = 1000  # TODO: this is a bad way to do it, would be better to just get every point where there is a pair
    points = np.linspace(xmin, xmax, n_points)
    for selection, label, color, line_style in zip(
        selections, labels, colors, line_styles
    ):
        n = data.loc[selection, score_column].notnull().sum()
        pos = np.array([(data.loc[selection, score_column] > x).sum() for x in points])
        neg = n - pos
        fracs = pos / n
        ax.plot(
            points, fracs, label=label, color=color, linestyle=line_style
        )  # BUG , **plot_kwargs)
        intv = stats.beta.interval(0.6827, pos + 1, neg + 1)
        errs = [fracs - intv[0], intv[1] - fracs]
        errs[0][pos == 0] = 0.0
        errs[1][neg == 0] = 0.0
        ax.fill_between(
            points,
            fracs - errs[0],
            fracs + errs[1],
            color=color,
            alpha=0.2,
            linewidth=0,
        )
    ax.set_ylim(0, ymax)
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel("Fraction positive")
    ax.set_xlabel("Score threshold")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    if threshold is not None:
        ax.axvline(
            x=threshold,
            ymin=0,
            ymax=1,
            linestyle="--",
            color=threshold_color,
            linewidth=1,
        )
        if threshold_label is not None:
            ax.text(
                x=threshold + (xmax - xmin) * 0.02,
                y=ymax,
                s=threshold_label,
                color=threshold_color,
                verticalalignment="top",
                horizontalalignment="left",
                fontsize=8,
            )


def table_circle_size_plot(df, ax=None, scale=3000, fontsize=12):
    if ax is None:
        ax = plt.gca()
    scale_factor = scale / df.sum().sum()
    ax.scatter(
        x=[j for i in range(df.shape[0]) for j in range(df.shape[1])],
        y=[i for i in range(df.shape[0]) for j in range(df.shape[1])],
        s=df.values.flatten() * scale_factor,
        clip_on=False,
    )
    # write numbers
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            ax.text(
                s=df.values[i, j],
                x=j,
                y=i,
                ha="center",
                va="center",
                fontsize=fontsize,
                color="white",
            )
    ax.set_xticks(range(df.shape[1]))
    ax.set_yticks(range(df.shape[0]))
    ax.set_xlabel(df.columns.name)
    ax.set_ylabel(df.index.name)
    ax.set_xticklabels(df.columns.values, rotation=45, ha="right")
    ax.set_yticklabels(df.index.values)
    ax.set_xlim(-0.5, df.shape[1] - 0.5)
    ax.set_ylim(df.shape[0] - 0.5, -0.5)  # NOTE: flipping axis
    ax.set_aspect("equal")
    ax.xaxis.set_tick_params(length=0)
    ax.yaxis.set_tick_params(length=0)
    for loc in ["left", "right", "top", "bottom"]:
        ax.spines[loc].set_visible(False)


def mimic_r_boxplot(ax):
    for i, patch in enumerate(ax.artists):
        r, g, b, a = patch.get_facecolor()
        col = (r, g, b, 1)
        patch.set_facecolor((r, g, b, 0.5))
        patch.set_edgecolor((r, g, b, 1))

        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        line_order = ["lower", "upper", "whisker_1", "whisker_2", "med", "fliers"]
        for j in range(i * 6, i * 6 + 6):
            elem = line_order[j % 6]
            line = ax.lines[j]
            if "whisker" in elem:
                line.set_visible(False)
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)
            if "fliers" in elem:
                line.set_alpha(0.5)


def annotate_pval(ax, x1, x2, y, h, text_y, val, fontsize, text=None):
    from decimal import Decimal

    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], c="black", linewidth=0.5)
    if val < 0.0001:
        s = "{:.2e}".format(Decimal(val))
        # text = "**"
    elif val < 0.05:
        s = "%.4f" % val
        # text = "*"
    else:
        s = "%.2f" % val
        # text = "n.s."

    if text is not None:
        s = text

    ax.text(
        (x1 + x2) * 0.5,
        text_y,
        s,
        ha="center",
        va="bottom",
        color="black",
        size=fontsize,
    )


def nice_boxplot(
    df,
    ycat,
    xcat,
    pal,
    xorder,
    pys,
    ay,
    xlabel,
    xticklabels,
    ylabel,
    log_scale,
    ylim,
    title,
    figsize=(2,2.5),
):
    fig = plt.figure(figsize=figsize)

    ax = sns.boxplot(data=df, y=ycat, x=xcat, order=xorder, palette=pal, fliersize=0)
    mimic_r_boxplot(ax)

    sns.swarmplot(
        data=df,
        y=ycat,
        x=xcat,
        order=xorder,
        palette=pal,
        ax=ax,
        size=2,
        edgecolor="black",
        linewidth=0.5,
        alpha=0.5,
    )

    # calculate differences
    for comp, xs, y, d_y in zip(
        [
            (xorder[0], xorder[1]),
            (xorder[0], xorder[2]),
            (xorder[0], xorder[3]),
            (xorder[1], xorder[2]),
        ],
        [(0, 1), (0, 2), (0, 3), (1, 2)],
        pys,
        [0, 0, 0, 0],
    ):
        cat_a = comp[0]
        cat_b = comp[1]
        dist_a = list(df[(df[xcat] == cat_a)][ycat])
        dist_b = list(df[(df[xcat] == cat_b)][ycat])

        u, p = stats.mannwhitneyu(dist_a, dist_b, alternative="two-sided")
        print(p)

        annotate_pval(ax, xs[0], xs[1], y, 0, y - (y * d_y), p, PAPER_FONTSIZE - 1)
    return fig, ax


def nice_violinplot(
    df,
    ycat,
    xcat,
    pal,
    xorder,
    pys,
    ay,
    xlabel,
    xticklabels,
    ylabel,
    log_scale,
    ylim,
    title,
    figf,
):
    fig = plt.figure(figsize=(2, 2))

    ax = sns.violinplot(
        data=df,
        y=ycat,
        x=xcat,
        order=xorder,
        palette=pal,
        cut=0,
        inner="quartiles",
        scale="width",
    )

    # edit quartile lines
    for l in ax.lines:
        l.set_linestyle("--")
        l.set_linewidth(0.6)
        l.set_color("black")
        l.set_alpha(0.5)
    for l in ax.lines[1::3]:
        l.set_linestyle("-")
        l.set_linewidth(1.0)
        l.set_color("black")
        l.set_alpha(1)

    # calculate differences
    for comp, xs, y, d_y in zip(
        [
            (xorder[0], xorder[1]),
            (xorder[0], xorder[2]),
            (xorder[0], xorder[3]),
            (xorder[1], xorder[2]),
        ],
        [(0, 1), (0, 2), (0, 3), (1, 2)],
        pys,
        [0, 0, 0, 0],
    ):
        cat_a = comp[0]
        cat_b = comp[1]
        dist_a = list(df[(df[xcat] == cat_a)][ycat])
        dist_b = list(df[(df[xcat] == cat_b)][ycat])

        u, p = stats.mannwhitneyu(dist_a, dist_b, alternative="two-sided")
        print(p)

        annotate_pval(ax, xs[0], xs[1], y, 0, y - (y * d_y), p, PAPER_FONTSIZE)

    # add N to plot
    for i, label in enumerate(xorder):
        n = len(df[(df[xcat] == label) & (~pd.isnull(ycat))])
        print(n)
        ax.annotate(
            str(n),
            xy=(i, ay),
            xycoords="data",
            xytext=(0, 0),
            textcoords="offset pixels",
            ha="center",
            va="top",
            color=pal[label],
            size=PAPER_FONTSIZE,
        )

    ax.set_xlabel(xlabel)
    ax.set_xticklabels(xticklabels, ha="right", va="top", rotation=30)
    ax.set_ylabel(ylabel)
    if log_scale:
        ax.set_yscale("log")
    ax.set_ylim(ylim)

    ax.set_title(title)

    fig.savefig(figf, dpi="figure", bbox_inches="tight")


def checkerboard(
    data,
    protein_a_column=None,
    protein_b_column=None,
    detection_columns=None,
    sort=True,
    assay_labels=None,
    positive_color="yellow",
    negative_color="white",
    ax=None,
):
    """Plot yes/no detection for benchmark PPI set with different assays

    See Braun et al, Nature Methods, 2010 for examples.

    Args:
        data (pandas.DataFrame): PPI test results. No missing values.
        protein_a/b_column (str): columns with protein names
        detection_columns (list(str)): name of columns containing boolean results
        sort (bool): whether to sort pairs by number of assays detected and assay order
        assay_labels (list(str)): names of assays to print
        positive_color (str/RGB/RGBA or list(colors)): single color or list of colors for each different assay
        negative_color (str/RGB/RGBA): color to indicate undetected pairs
        ax (matplotlib.axes.Axes): Axes to draw plot onto

    Examples:
        Make a checkerboard of some dummy data:

        .. plot::
            :context: close-figs

            >>> import pandas as pd
            >>> from plotting import checkerboard
            >>> prs_results = pd.DataFrame(columns=['gene_a', 'gene_b', 'Y2H', 'MAPPIT', 'GPCA'],
            ...                            data=[['ABC1', 'ABC2', False, False, False],
            ...                                  ['EFG1', 'EFG2', False, False, False],
            ...                                  ['HIJ1', 'HIJ2', True, False, False],
            ...                                  ['KLM1', 'KLM2', False, False, False],
            ...                                  ['NOP1', 'NOP2', True, True, True],
            ...                                  ['QRS1', 'QRS2', True, False, True],
            ...                                  ['TUV1', 'TUV2', False, False, True],
            ...                                  ['XYZ1', 'XYZ2', False, False, False]])
            >>> checkerboard(data=prs_results,
            ...              protein_a_column='gene_a',
            ...              protein_b_column='gene_b',
            ...              detection_columns=['Y2H',
            ...                                 'MAPPIT',
            ...                                 'GPCA'])

    """
    df = data.copy()
    if ax is None:
        ax = plt.gca()
    if assay_labels is None:
        assay_labels = detection_columns
    if protein_a_column is None:
        protein_a_column = df.columns[0]
    if protein_b_column is None:
        protein_a_column = df.columns[1]
    if detection_columns is None:
        detection_columns = list(df.columns[df.dtypes == bool])
    elif isinstance(detection_columns, str):
        detection_columns = [detection_columns]
    df["total_positives"] = df[detection_columns].sum(axis=1)
    if sort:
        df = df.sort_values(by=["total_positives"] + detection_columns, ascending=False)
    if isinstance(positive_color, list) and len(positive_color) == len(
        detection_columns
    ):
        results = df[detection_columns].values.T
        for i, color in enumerate(positive_color):
            m = np.zeros(shape=results.shape, dtype=bool)
            m[i, :] = True
            ax.imshow(
                results * m,
                cmap=matplotlib.colors.ListedColormap(
                    [negative_color if i == 0 else (0, 0, 0, 0), color]
                ),
            )
    else:
        ax.imshow(
            df[detection_columns].values.T,
            cmap=matplotlib.colors.ListedColormap([negative_color, positive_color]),
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_yticks(range(len(detection_columns)))
    ax.set_yticklabels(assay_labels)
    ax.yaxis.set_tick_params(length=0)
    ax.set_xticks([])
    len_longest_name = df[protein_a_column].str.len().max()
    for i, (name_a, name_b) in enumerate(
        zip(df[protein_a_column].values, df[protein_b_column].values)
    ):
        ax.text(
            i,
            -0.6,
            name_a + " " * (len_longest_name - len(name_a) + 2) + name_b,
            rotation=90,
            fontfamily="monospace",
            va="bottom",
            ha="center",
        )
