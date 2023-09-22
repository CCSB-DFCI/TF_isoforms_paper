from itertools import combinations

import numpy as np
import pandas as pd


def paralog_pair_ppi_table(data, tf_gene_a, tf_gene_b):
    gene_a_partners = data.loc[
        (data["ad_gene_symbol"] == tf_gene_a)
        & (data["category"] == "tf_isoform_ppis")
        & (data["Y2H_result"] == True),
        "db_gene_symbol",
    ].unique()
    gene_b_partners = data.loc[
        (data["ad_gene_symbol"] == tf_gene_b)
        & (data["category"] == "tf_isoform_ppis")
        & (data["Y2H_result"] == True),
        "db_gene_symbol",
    ].unique()
    partners = np.concatenate([gene_a_partners, gene_b_partners])
    tf = data.loc[
        ((data["ad_gene_symbol"] == tf_gene_a) | (data["ad_gene_symbol"] == tf_gene_b))
        & data["category"].isin(
            [
                "tf_isoform_ppis",
                "tf_paralog_ppis",
                "non_paralog_control",
                "paralog_with_PDI",
            ]
        )
        & data["db_gene_symbol"].isin(partners),
        ["ad_clone_acc", "db_gene_symbol", "Y2H_result"],
    ].copy()
    return tf


def isoforms_of_paralogs_pairs(paralog_pairs, isoforms):
    """Given a list of paralog gene pairs and a list of isoforms per gene,
    return all pairwise isoform combinations between paralogous genes.

    Args:
        pairs ([type]): [description]
        isoforms (bool): [description]

    Returns:
        [type]: [description]
    """
    pairs = pd.merge(
        paralog_pairs.loc[:, ["tf_gene_a", "tf_gene_b", "is_paralog_pair"]],
        isoforms.loc[:, ["gene", "clone_acc", "aa_seq"]],
        how="inner",
        left_on="tf_gene_a",
        right_on="gene",
    ).rename(columns={"aa_seq": "aa_seq_a", "clone_acc": "clone_acc_a"})
    pairs = pd.merge(
        pairs,
        isoforms.loc[:, ["gene", "clone_acc", "aa_seq"]],
        how="inner",
        left_on="tf_gene_b",
        right_on="gene",
    ).rename(columns={"aa_seq": "aa_seq_b", "clone_acc": "clone_acc_b"})
    pairs = pairs.drop(columns=["gene_x", "gene_y"])
    return pairs


def _pairs_comparison_table(pairs, y2h, y1h, m1h):
    if y2h is not None:
        pairs["ppi_n_tested"] = pairs.apply(
            ppi_metric, data=y2h, function=number_tested_partners, axis=1
        )
        pairs["ppi_n_shared"] = pairs.apply(
            ppi_metric, data=y2h, function=number_shared_partners, axis=1
        )
        pairs["ppi_n_min"] = pairs.apply(
            ppi_metric, data=y2h, function=number_min_partners, axis=1
        )
        pairs["ppi_n_min_diff"] = pairs.apply(
            ppi_metric, data=y2h, function=min_difference, axis=1
        )
        pairs["ppi_jaccard"] = pairs.apply(
            ppi_metric, data=y2h, function=jaccard_index, axis=1
        )
        pairs["ppi_simpson"] = pairs.apply(
            ppi_metric, data=y2h, function=simpsons_index, axis=1
        )
        pairs["ppi_n_diff"] = pairs["ppi_n_tested"] - pairs["ppi_n_shared"]
        pairs["ppi_delta_n"] = pairs.apply(
            ppi_metric, data=y2h, function=size_b_minus_size_a, axis=1
        )

    if y1h is not None:
        pairs["pdi_n_tested"] = pairs.apply(
            pdi_metric, data=y1h, function=number_tested_partners, axis=1
        )
        pairs["pdi_n_shared"] = pairs.apply(
            pdi_metric, data=y1h, function=number_shared_partners, axis=1
        )
        pairs["pdi_n_min"] = pairs.apply(
            pdi_metric, data=y1h, function=number_min_partners, axis=1
        )
        pairs["pdi_n_min_diff"] = pairs.apply(
            pdi_metric, data=y1h, function=min_difference, axis=1
        )
        pairs["pdi_jaccard"] = pairs.apply(
            pdi_metric, data=y1h, function=jaccard_index, axis=1
        )
        pairs["pdi_simpson"] = pairs.apply(
            pdi_metric, data=y1h, function=simpsons_index, axis=1
        )
        pairs["pdi_n_diff"] = pairs["pdi_n_tested"] - pairs["pdi_n_shared"]
        pairs["pdi_delta_n"] = pairs.apply(
            pdi_metric, data=y1h, function=size_b_minus_size_a, axis=1
        )
    if m1h is not None:
        pairs["m1h_min"] = pairs.apply(m1h_min, data=m1h, axis=1)
        pairs["m1h_max"] = pairs.apply(m1h_max, data=m1h, axis=1)
        pairs["activation_fold_change"] = pairs.apply(fold_change_m1h, data=m1h, axis=1)
        pairs["activation_abs_fold_change"] = pairs["activation_fold_change"].abs()

    aa_ident = load_seq_comparison_data()
    pairs["aa_seq_pct_id"] = pairs.apply(
        lambda x: "_".join(sorted([x["clone_acc_a"], x["clone_acc_b"]])), axis=1
    ).map(aa_ident)
    # TMP DEBUG
    # if pairs["aa_seq_pct_id"].isnull().any():
    #    raise UserWarning("Problem with sequence similarity data")
    pairs = pairs.dropna(subset=["aa_seq_pct_id"])  # HACK
    pairs = pairs.set_index("pair")
    return pairs


def pairs_of_paralogs_and_isoforms_comparison_table(
    isoforms,
    paralog_pairs,
    y2h,
    y1h,
    m1h,
    restrict_isoforms_to_those_with_paralogs=False,
):
    """

    Restricted to isoforms pairs in the paralogs dataset.

    Args:
        isoforms (): [description]
        paralog_pairs ([type]): [description]
        y2h ([type]): [description]
        y1h ([type]): [description]
        m1h ([type]): [description]

    Returns:
        [type]: [description]

    """
    pairs = isoforms_of_paralogs_pairs(paralog_pairs, isoforms)
    pairs["category"] = pairs["is_paralog_pair"].map(
        {True: "paralogs", False: "non-paralog-control"}
    )
    pairs = pairs.drop(columns=["aa_seq_a", "aa_seq_b", "is_paralog_pair"])
    iso_pairs = []
    for tf_gene in isoforms["gene"].unique():
        tf_iso = isoforms.loc[isoforms["gene"] == tf_gene, "clone_acc"].values
        for iso_a, iso_b in combinations(tf_iso, 2):
            iso_pairs.append((tf_gene, tf_gene, iso_a, iso_b))
    iso_pairs = pd.DataFrame(
        data=iso_pairs, columns=["tf_gene_a", "tf_gene_b", "clone_acc_a", "clone_acc_b"]
    )
    iso_pairs["category"] = "isoforms"
    if restrict_isoforms_to_those_with_paralogs:
        iso_pairs = iso_pairs.loc[
            iso_pairs["tf_gene_a"].isin(pairs["tf_gene_a"])
            | iso_pairs["tf_gene_a"].isin(pairs["tf_gene_b"]),
            :,
        ]
    pairs = pd.concat([pairs, iso_pairs])
    pairs["pair"] = pairs.apply(
        lambda x: "_".join(sorted([x["clone_acc_a"], x["clone_acc_b"]])), axis=1
    )
    return _pairs_comparison_table(pairs, y2h, y1h, m1h)


def tissue_expression_similarity_per_gene():
    pass


def ppi_metric(row, data, function, suffixes=("_a", "_b")):
    ad_a = row["clone_acc" + suffixes[0]]
    ad_b = row["clone_acc" + suffixes[1]]
    pair = data.loc[data["ad_clone_acc"].isin([ad_a, ad_b]), :].pivot(
        values="Y2H_result", index="db_gene_symbol", columns="ad_clone_acc"
    )
    if ad_a not in pair.columns or ad_b not in pair.columns:
        return np.nan
    # remove any partner with AA / NC / NS / NaN in either
    pair = pair.loc[pair.notnull().all(axis=1), :].astype(int).astype(bool)
    # remove partners that tested negative in both
    pair = pair.loc[pair.any(axis=1), :]
    if pair.shape[0] > 0:
        return function(
            set(pair.index[pair[ad_a]].values), set(pair.index[pair[ad_b]].values)
        )
    else:
        return np.nan


def pdi_metric(row, data, function):
    df = data.loc[
        (data["unique_acc"] == row["clone_acc_a"])
        | (data["unique_acc"] == row["clone_acc_b"]),
        data.columns[2:],
    ].copy()
    if df.shape[0] < 2:
        return np.nan
    df = df.loc[:, df.any(axis=0)]
    if df.shape[1] == 0:
        return np.nan
    # kaia edited these 2 lines to drop any baits with NA
    # my version of pandas wasn't auto-filtering these out, think it got fixed in later versions
    a = set(df.columns[df.iloc[0].fillna(False)])
    b = set(df.columns[df.iloc[1].fillna(False)])
    return function(a, b)


def jaccard_index(a, b):
    if len(a) == 0 and len(b) == 0:
        return np.nan
    return len(a.intersection(b)) / float(len(a.union(b)))


def simpsons_index(a, b):
    min_size = min(len(a), len(b))
    if min_size == 0:
        return np.nan
    else:
        return len(a.intersection(b)) / float(min_size)


def number_tested_partners(a, b):
    """Comes up with nan when it should be 0?"""
    return len(a.union(b))


def number_shared_partners(a, b):
    return len(a.intersection(b))


def number_min_partners(a, b):
    return min(len(a), len(b))


def min_difference(a, b):
    return min(len(a.difference(b)), len(b.difference(a)))


def size_b_minus_size_a(a, b):
    return len(b) - len(a)


def fold_change_m1h(row, data):
    if (
        row["clone_acc_a"] not in data["clone_acc"].values
        or row["clone_acc_b"] not in data["clone_acc"].values
    ):
        return np.nan
    a = (
        data.loc[
            data["clone_acc"] == row["clone_acc_a"],
            [c for c in data.columns if c.startswith("M1H_rep")],
        ]
        .mean(axis=1)
        .values[0]
    )
    b = (
        data.loc[
            data["clone_acc"] == row["clone_acc_b"],
            [c for c in data.columns if c.startswith("M1H_rep")],
        ]
        .mean(axis=1)
        .values[0]
    )
    return b - a


def m1h_min(row, data):
    if (
        row["clone_acc_a"] not in data["clone_acc"].values
        or row["clone_acc_b"] not in data["clone_acc"].values
    ):
        return np.nan
    a = (
        data.loc[
            data["clone_acc"] == row["clone_acc_a"],
            [c for c in data.columns if c.startswith("M1H_rep")],
        ]
        .mean(axis=1)
        .values[0]
    )
    b = (
        data.loc[
            data["clone_acc"] == row["clone_acc_b"],
            [c for c in data.columns if c.startswith("M1H_rep")],
        ]
        .mean(axis=1)
        .values[0]
    )
    return min([a, b])


def m1h_max(row, data):
    if (
        row["clone_acc_a"] not in data["clone_acc"].values
        or row["clone_acc_b"] not in data["clone_acc"].values
    ):
        return np.nan
    a = (
        data.loc[
            data["clone_acc"] == row["clone_acc_a"],
            [c for c in data.columns if c.startswith("M1H_rep")],
        ]
        .mean(axis=1)
        .values[0]
    )
    b = (
        data.loc[
            data["clone_acc"] == row["clone_acc_b"],
            [c for c in data.columns if c.startswith("M1H_rep")],
        ]
        .mean(axis=1)
        .values[0]
    )
    return max([a, b])
