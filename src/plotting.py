import os
import sys
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
import pandas as pd
from scipy import stats
import seaborn as sns

from data_loading import paralog_pair_ppi_table


COLOR_PURPLE = (155 / 255, 97 / 255, 153 / 255)

## kaia's added code
PAPER_PRESET = {
    "style": "ticks",
    "font": "Helvetica",
    "context": "paper",
    "rc": {
        "font.size": 7,
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


def violinplot_reflected(*args, **kwargs):
    """
    monkeypatch from https://github.com/mwaskom/seaborn/issues/525
    """
    fit_kde_func = sns.categorical._ViolinPlotter.fit_kde

    def reflected_once_kde(self, x, bw):
        lb = 0
        ub = 1

        kde, bw_used = fit_kde_func(self, x, bw)

        kde_evaluate = kde.evaluate

        def truncated_kde_evaluate(x):
            val = np.where((x >= lb) & (x <= ub), kde_evaluate(x), 0)
            val += np.where((x >= lb) & (x <= ub), kde_evaluate(lb - x), 0)
            val += np.where((x > lb) & (x <= ub), kde_evaluate(ub - (x - ub)), 0)
            return val

        kde.evaluate = truncated_kde_evaluate
        return kde, bw_used

    sns.categorical._ViolinPlotter.fit_kde = reflected_once_kde
    retval = sns.violinplot(*args, **kwargs)
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
        (data["category"] == "tf_isoform_ppis") & (data["ad_gene_symbol"] == gene_name),
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
        return
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
        data.loc[data["tf"] == gene_name, data.columns[1:]]
        .copy()
        .set_index("unique_acc")
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
        return
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


def m1h_activation_per_tf_gene_plot(tf_gene_name, data, ax=None, xlim=None):
    if ax is None:
        ax = plt.gca()
    rep_columns = [c for c in data.columns if c.startswith("M1H_rep")]
    is_all_na = (
        data[rep_columns].isnull().groupby(data["gene"]).all().all(axis=1)[tf_gene_name]
    )
    if tf_gene_name not in data["gene"].values or is_all_na:
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
        return
    clones = [
        isoform_display_name(acc)
        for acc in data.loc[data["gene"] == tf_gene_name, "clone_acc"].values
        for __ in range(len(rep_columns))
    ]
    values = data.loc[data["gene"] == tf_gene_name, rep_columns].values.flatten()
    n_reps = len(rep_columns)
    ax.barh(
        y=clones[::n_reps],
        width=data.loc[data["gene"] == tf_gene_name, rep_columns].mean(axis=1).values,
        edgecolor="black",
        color="slategrey",
        alpha=0.5,
        height=0.6,
    )

    # swarmplot to make points more visible
    df = pd.DataFrame()
    df["clone"] = clones
    df["value"] = values
    sns.stripplot(
        data=df,
        x="value",
        y="clone",
        ax=ax,
        color="white",
        linewidth=1,
        edgecolor="black",
        size=4,
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

    if xlim == None:
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


def annotate_pval(ax, x1, x2, y, h, text_y, val, fontsize):
    from decimal import Decimal

    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c="black", linewidth=0.5)
    if val < 0.0001:
        text = "{:.2e}".format(Decimal(val))
        # text = "**"
    elif val < 0.05:
        text = "%.4f" % val
        # text = "*"
    else:
        text = "%.4f" % val
        # text = "n.s."
    ax.text(
        (x1 + x2) * 0.5,
        text_y,
        text,
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
    figf,
):
    fig = plt.figure(figsize=(2, 2.5))

    ax = sns.boxplot(data=df, y=ycat, x=xcat, order=xorder, palette=pal, fliersize=0)

    sns.swarmplot(
        data=df,
        y=ycat,
        x=xcat,
        order=xorder,
        palette=pal,
        ax=ax,
        size=4,
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

        u, p = mannwhitneyu(dist_a, dist_b, alternative="two-sided")
        print(p)

        annotate_pval(ax, xs[0], xs[1], y, 0, y - (y * d_y), p, PAPER_FONTSIZE)


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

        u, p = mannwhitneyu(dist_a, dist_b, alternative="two-sided")
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
