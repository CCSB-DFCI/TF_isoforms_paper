import os
import sys
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
import pandas as pd
from scipy import stats

sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)), "../.."))
from isoform_pairwise_metrics import paralog_pair_ppi_table


COLOR_PURPLE = (155 / 255, 97 / 255, 153 / 255)


def binary_profile_matrix(
    data,
    ax=None,
    box_size=0.7,
    fill_color="black",
    border_color="black",
    column_label_rotation=40,
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
        r = patches.Rectangle(
            (x - box_size / 2, y - box_size / 2),
            box_size,
            box_size,
            fill=False,
            edgecolor=border_color,
        )
        ax.add_patch(r)
    for x, y in positives:
        r = patches.Rectangle(
            (x - box_size / 2, y - box_size / 2),
            box_size,
            box_size,
            fill=True,
            facecolor=fill_color,
            edgecolor=border_color,
        )
        ax.add_patch(r)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.xaxis.tick_top()
    ax.set_xticks(range(data.shape[1]))
    ax.xaxis.set_tick_params(length=0)
    ax.xaxis.set_ticklabels(data.columns, rotation=column_label_rotation, ha="left")
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
    gene_name, data, ax=None, min_n_isoforms=1, min_n_partners=1
):
    tf = data.loc[
        (data["category"] == "tf_isoform_ppis") & (data["ad_gene_symbol"] == gene_name),
        ["ad_clone_acc", "db_gene_symbol", "Y2H_result"],
    ].copy()
    # tf['score'] = tf['score'].map({'1': True,
    #                               '0': False,
    #                               'AA': np.nan,
    #                               'NC': np.nan})
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
    binary_profile_matrix(tf, ax=ax, column_label_rotation=90)
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )


def y2h_ppi_per_paralog_pair_plot(
    tf_gene_a, tf_gene_b, data, ax=None, min_n_isoforms=1, min_n_partners=1
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
    tf["score"] = tf["score"].map({"1": True, "0": False, "AA": np.nan, "NC": np.nan})
    tf["ad_clone_acc"] = tf["ad_clone_acc"].apply(isoform_display_name)
    tf = tf.pivot(index="ad_clone_acc", columns="db_gene_symbol", values="score")
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
    binary_profile_matrix(tf, ax=ax, column_label_rotation=90)
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )


def y1h_pdi_per_tf_gene_plot(
    gene_name, data, ax=None, min_n_isoforms=1, min_n_partners=1
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
    binary_profile_matrix(tf, ax=ax, column_label_rotation=90)
    ax.set_yticklabels(
        [
            strikethrough(name) if all_na else name
            for name, all_na in tf.isnull().all(axis=1).items()
        ]
    )


def m1h_activation_per_tf_gene_plot(tf_gene_name, data, ax=None):
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
        color="grey",
        alpha=0.5,
        height=0.6,
    )
    ax.scatter(y=clones, x=values, alpha=0.5, color="C0")
    ax.set_yticks(
        clones[::n_reps]
    )  # needed to avoid truncating clones with missing data
    ax.set_yticklabels(
        [
            c if not np.isnan(v) else strikethrough(c)
            for c, v in zip(clones[::n_reps], values[::n_reps])
        ]
    )
    ax.set_xlim(-3, 12)
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
