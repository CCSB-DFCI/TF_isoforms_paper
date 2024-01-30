from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import tqdm

from data_loading import (
    load_dbd_accessions,
    load_annotated_TFiso1_collection,
    load_y2h_isoform_data,
    load_y1h_pdi_data,
    load_m1h_activation_data,
    DIMERIZING_TF_FAMILIES,
    load_tf_families,
    load_human_tf_db,
    load_ppi_partner_categories,
)
from .utils import cache_with_pickle, DATA_DIR


@cache_with_pickle
def load_seq_id_between_cloned_isoforms():
    """ """
    tfs = load_annotated_TFiso1_collection()
    tfs = {k: v for k, v in tfs.items() if len(v.cloned_isoforms) >= 2}

    def pairwise_global_aa_sequence_similarity(aa_seq_a, aa_seq_b):
        alignment = pairwise2.align.globalds(
            aa_seq_a, aa_seq_b, blosum62, -10, -0.5, penalize_end_gaps=False
        )
        alignment = pairwise2.format_alignment(*alignment[0]).split("\n")[1]
        return alignment.count("|") / len(alignment) * 100

    aa_id = []
    tfs = {k: tfs[k] for k in sorted(tfs.keys())}
    for tf in tqdm.tqdm(list(tfs.values())):
        for iso_a, iso_b in combinations(tf.cloned_isoforms, 2):
            aa_id.append(
                (
                    tf.name,
                    iso_a.clone_acc,
                    iso_b.clone_acc,
                    pairwise_global_aa_sequence_similarity(iso_a.aa_seq, iso_b.aa_seq),
                )
            )
    aa_id = pd.DataFrame(
        aa_id,
        columns=[
            "gene_symbol",
            "clone_acc_a",
            "clone_acc_b",
            "aa_seq_pct_identity",
        ],
    )
    aa_id = pd.concat(
        [
            aa_id,
            aa_id.copy().rename(
                columns={"clone_acc_a": "clone_acc_b", "clone_acc_b": "clone_acc_a"}
            ),
        ]
    )
    return aa_id


def load_ref_vs_alt_isoforms_table():
    fpath = Path(__file__).resolve().parents[2] / "output" / "TF-iso_ref-vs-alt.tsv"
    if not fpath.exists():
        print("generating reference vs alternative table")
        _write_TF_iso_ref_vs_alt_table(fpath)
    return pd.read_csv(fpath, sep="\t")


def _pairs_of_ref_vs_alt_isoforms_comparison_table(tfs):
    iso_pairs = []
    for tf in tfs.values():
        ref = tf.cloned_reference_isoform
        for alt in tf.cloned_isoforms:
            if alt.name == ref.name:
                continue
            iso_pairs.append(
                (
                    tf.name,
                    tf.ensembl_gene_id,
                    tf.tf_family,
                    tf.tf_family in DIMERIZING_TF_FAMILIES,
                    ref.clone_acc,
                    alt.clone_acc,
                    "|".join(ref.ensembl_transcript_ids)
                    if ref.ensembl_transcript_ids is not None
                    else np.nan,
                    "|".join(alt.ensembl_transcript_ids)
                    if alt.ensembl_transcript_ids is not None
                    else np.nan,
                    ref.is_novel_isoform(),
                    alt.is_novel_isoform(),
                    tf.cloned_MANE_select_isoform,
                    len(ref.aa_seq),
                    len(alt.aa_seq),
                    len(ref.exons),
                    len(alt.exons),
                    tf.splicing_categories(ref.name, alt.name)[
                        "alternative N-terminal"
                    ],
                    tf.splicing_categories(ref.name, alt.name)[
                        "alternative C-terminal"
                    ],
                    tf.splicing_categories(ref.name, alt.name)[
                        "alternative internal exon"
                    ],
                    tf.splicing_categories(ref.name, alt.name)[
                        "alternative 5' splice site"
                    ],
                    tf.splicing_categories(ref.name, alt.name)[
                        "alternative 3' splice site"
                    ],
                    tf.splicing_categories(ref.name, alt.name)["exon skipping"],
                    tf.splicing_categories(ref.name, alt.name)[
                        "mutually exclusive exons"
                    ],
                    tf.splicing_categories(ref.name, alt.name)["intron retention"],
                )
            )
    iso_pairs = pd.DataFrame(
        data=iso_pairs,
        columns=[
            "gene_symbol",
            "Ensembl_gene_ID",
            "family",
            "is_dimerizing_TF_family",
            "clone_acc_ref",
            "clone_acc_alt",
            "Ensembl_transcript_IDs_ref",
            "Ensembl_transcript_IDs_alt",
            "is_ref_novel_isoform",
            "is_alt_novel_isoform",
            "is_MANE_select_isoform_cloned",
            "n_aa_ref",
            "n_aa_alt",
            "n_exons_ref",
            "n_exons_alt",
            "is_alternative_N_terminal",
            "is_alternative_C_terminal",
            "is_alternative_internal_exon",
            "is_alternative_5_prime_donor",
            "is_alternative_3_prime_acceptor",
            "is_exon_skipping",
            "is_mutually_exclusive_exons",
            "is_intron_retention",
        ],
    )
    return iso_pairs


def _write_TF_iso_ref_vs_alt_table(outpath):
    tfs = load_annotated_TFiso1_collection()
    df = _pairs_of_ref_vs_alt_isoforms_comparison_table(tfs)
    m1h = load_m1h_activation_data()
    m1h["mean"] = m1h[["M1H_rep1", "M1H_rep2", "M1H_rep3"]].mean(axis=1)
    m1h["abs_mean"] = m1h["mean"].abs()
    y1h = load_y1h_pdi_data()
    if y1h["clone_acc"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    y1h = y1h.set_index("clone_acc")

    def load_dbd_affected():
        df = pd.concat(
            [
                g.aa_feature_disruption(g.cloned_reference_isoform.name)
                for g in tfs.values()
            ]
        )
        df["is_DBD"] = df["accession"].isin(load_dbd_accessions())
        df_new = (
            df.loc[df["is_DBD"], :]
            .groupby(["gene_symbol", "ref_iso", "alt_iso"])[["deletion", "frameshift"]]
            .sum()
            .sum(axis=1)
            / df.loc[df["is_DBD"], :]
            .groupby(["gene_symbol", "ref_iso", "alt_iso"])["length"]
            .sum()
        ).to_frame(name="dbd_fraction")
        df_new["dbd_insertion_n_aa"] = (
            df.loc[df["is_DBD"], :]
            .groupby(["gene_symbol", "ref_iso", "alt_iso"])["insertion"]
            .sum()
        )
        df = df_new.reset_index()
        df["dbd_pct_lost"] = df["dbd_fraction"] * 100.0
        df = df.drop(columns=["dbd_fraction"])
        return df

    dbd = load_dbd_affected()
    dbd["clone_acc_ref"] = dbd["ref_iso"].map(
        {iso.name: iso.clone_acc for tf in tfs.values() for iso in tf.cloned_isoforms}
    )
    dbd["clone_acc_alt"] = dbd["alt_iso"].map(
        {iso.name: iso.clone_acc for tf in tfs.values() for iso in tf.cloned_isoforms}
    )
    dbd = dbd.drop(columns=["gene_symbol", "ref_iso", "alt_iso"])
    df = pd.merge(df, dbd, how="left", on=["clone_acc_ref", "clone_acc_alt"])
    df["dbd_affected"] = df["dbd_pct_lost"] > 0
    df.loc[df["dbd_pct_lost"].isnull(), "dbd_affected"] = np.nan

    aa_ident = load_seq_id_between_cloned_isoforms()
    df = pd.merge(
        df,
        aa_ident.loc[:, ["clone_acc_a", "clone_acc_b", "aa_seq_pct_identity"]],
        how="left",
        left_on=["clone_acc_ref", "clone_acc_alt"],
        right_on=["clone_acc_a", "clone_acc_b"],
    ).drop(columns=["clone_acc_a", "clone_acc_b"])
    if df["aa_seq_pct_identity"].isnull().any():
        raise UserWarning("Unexpected missing sequence similarity values")

    y2h = load_y2h_isoform_data()
    y2h_complete = load_y2h_isoform_data(require_at_least_one_ppi_per_isoform=False)
    _add_PPI_columns(df=df, y2h=y2h_complete)
    ppi_partner_cats = load_ppi_partner_categories()
    tfdb = load_human_tf_db()
    fam = load_tf_families()
    y2h["ad_tf_family"] = y2h["ad_gene_symbol"].map(fam)
    y2h["db_tf_family"] = y2h["db_gene_symbol"].map(fam)
    y2h["is_dimerizing_ppi"] = y2h["ad_tf_family"].isin(DIMERIZING_TF_FAMILIES) & (
        y2h["ad_tf_family"] == y2h["db_tf_family"]
    )
    y2h["is_tf_tf_ppi"] = y2h["db_gene_symbol"].isin(tfdb["HGNC symbol"].unique())

    # of reference dimer PPI, are all lost, some lost, none lost
    def ppi_pertubation(row, ppi):
        ref_clone_acc = row["clone_acc_ref"]
        alt_clone_acc = row["clone_acc_alt"]
        if (
            ref_clone_acc not in ppi["ad_clone_acc"].unique()
            or alt_clone_acc not in ppi["ad_clone_acc"].unique()
        ):
            return np.nan
        df = (
            ppi.loc[
                ppi["ad_clone_acc"].isin([ref_clone_acc, alt_clone_acc]),
                ["ad_clone_acc", "db_gene_symbol", "Y2H_result"],
            ]
            .pivot(values="Y2H_result", index="db_gene_symbol", columns="ad_clone_acc")
            .dropna()
        )
        df = df.loc[df.any(axis=1), :]
        if df.shape[0] == 0:
            return np.nan
        if df.all().all():
            return "retains all"
        elif not df[alt_clone_acc].any():
            return "loses all"
        elif df[alt_clone_acc].sum() > df[ref_clone_acc].sum():
            return "gains some"
        else:
            return "loses some"

    df["dimer_ppi"] = df.apply(
        ppi_pertubation, ppi=y2h.loc[y2h["is_dimerizing_ppi"], :], axis=1
    )
    df["other_than_dimer_ppi"] = df.apply(
        ppi_pertubation, ppi=y2h.loc[~y2h["is_dimerizing_ppi"], :], axis=1
    )
    df["tf_tf_ppi"] = df.apply(
        ppi_pertubation, ppi=y2h.loc[y2h["is_tf_tf_ppi"], :], axis=1
    )

    df["tf_cofactor_ppi"] = df.apply(
        ppi_pertubation,
        ppi=y2h.loc[
            y2h["db_gene_symbol"].isin(
                ppi_partner_cats.loc[
                    ppi_partner_cats["category"] == "cofactor", "gene_symbol_partner"
                ].unique()
            ),
            :,
        ],
        axis=1,
    )
    df["tf_signalling_ppi"] = df.apply(
        ppi_pertubation,
        ppi=y2h.loc[
            y2h["db_gene_symbol"].isin(
                ppi_partner_cats.loc[
                    ppi_partner_cats["category"] == "signaling",
                    "gene_symbol_partner",
                ].unique()
            ),
            :,
        ],
        axis=1,
    )

    _add_PDI_columns(df=df, y1h=y1h)
    _add_activation_columns(df=df, m1h=m1h)
    df.to_csv(outpath, sep="\t", index=False)


def _add_PPI_columns(df, y2h, suffixes=("_ref", "_alt")):
    n_ppi = y2h.loc[(y2h["Y2H_result"] == True), :].groupby("ad_clone_acc").size()
    df[f"n_positive_PPI{suffixes[0]}"] = df[f"clone_acc{suffixes[0]}"].map(n_ppi)
    df[f"n_positive_PPI{suffixes[1]}"] = df[f"clone_acc{suffixes[1]}"].map(n_ppi)
    # BUG MISSING 0's here!
    df.loc[
        df[f"n_positive_PPI{suffixes[0]}"].isnull()
        & df[f"clone_acc{suffixes[0]}"].isin(
            y2h.loc[(y2h["Y2H_result"] == False), "ad_clone_acc"].unique()
        ),
        f"n_positive_PPI{suffixes[0]}",
    ] = 0
    df.loc[
        df[f"n_positive_PPI{suffixes[1]}"].isnull()
        & df[f"clone_acc{suffixes[1]}"].isin(
            y2h.loc[(y2h["Y2H_result"] == False), "ad_clone_acc"].unique()
        ),
        f"n_positive_PPI{suffixes[1]}",
    ] = 0

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

    df[f"n_PPI_successfully_tested_in{suffixes[0]}_and{suffixes[1]}"] = df.apply(
        ppi_metric,
        data=y2h,
        suffixes=suffixes,
        function=number_tested_partners,
        axis=1,
    )
    df[f"n_positive_PPI{suffixes[0]}_filtered"] = df.apply(
        ppi_metric,
        data=y2h,
        suffixes=suffixes,
        function=lambda a, b: len(a),
        axis=1,
    )
    df[f"n_positive_PPI{suffixes[1]}_filtered"] = df.apply(
        ppi_metric,
        data=y2h,
        suffixes=suffixes,
        function=lambda a, b: len(b),
        axis=1,
    )
    df["n_shared_PPI"] = df.apply(
        ppi_metric,
        data=y2h,
        suffixes=suffixes,
        function=number_shared_partners,
        axis=1,
    )
    df["n_PPI_diff"] = (
        df[f"n_PPI_successfully_tested_in{suffixes[0]}_and{suffixes[1]}"]
        - df["n_shared_PPI"]
    )
    df["PPI_delta_n"] = (
        df[f"n_positive_PPI{suffixes[1]}_filtered"]
        - df[f"n_positive_PPI{suffixes[0]}_filtered"]
    )
    df["PPI_jaccard"] = df.apply(
        ppi_metric,
        data=y2h,
        suffixes=suffixes,
        function=jaccard_index,
        axis=1,
    )


def _add_PDI_columns(df, y1h, suffixes=("_ref", "_alt")):
    def pdi_metric(row, data, function, suffixes=("_a", "_b")):
        clone_acc_a = row["clone_acc" + suffixes[0]]
        clone_acc_b = row["clone_acc" + suffixes[1]]
        df = data.loc[
            (data.index == clone_acc_a) | (data.index == clone_acc_b),
            data.columns[1:],
        ].copy()
        if df.shape[0] < 2:
            return np.nan
        df = df.loc[
            [clone_acc_a, clone_acc_b], df.any(axis=0) & ~df.isnull().any(axis=0)
        ]
        if df.shape[1] == 0:
            return np.nan
        a = set(df.columns[df.iloc[0]])
        b = set(df.columns[df.iloc[1]])
        return function(a, b)

    n_pdi = y1h.drop(columns=["gene_symbol"]).sum(axis=1)
    df[f"n_positive_PDI{suffixes[0]}"] = df[f"clone_acc{suffixes[0]}"].map(n_pdi)
    df[f"n_positive_PDI{suffixes[1]}"] = df[f"clone_acc{suffixes[1]}"].map(n_pdi)
    df[f"n_PDI_successfully_tested_in{suffixes[0]}_and{suffixes[1]}"] = df.apply(
        pdi_metric,
        data=y1h,
        suffixes=suffixes,
        function=number_tested_partners,
        axis=1,
    )

    df[f"n_positive_PDI{suffixes[0]}_filtered"] = df.apply(
        pdi_metric,
        data=y1h,
        suffixes=suffixes,
        function=lambda a, b: len(a),
        axis=1,
    )
    df[f"n_positive_PDI{suffixes[1]}_filtered"] = df.apply(
        pdi_metric,
        data=y1h,
        suffixes=suffixes,
        function=lambda a, b: len(b),
        axis=1,
    )

    df["n_shared_PDI"] = df.apply(
        pdi_metric,
        data=y1h,
        suffixes=suffixes,
        function=number_shared_partners,
        axis=1,
    )
    df["PDI_jaccard"] = df.apply(
        pdi_metric, data=y1h, suffixes=suffixes, function=jaccard_index, axis=1
    )


def _add_activation_columns(df, m1h, suffixes=("_ref", "_alt")):
    if "gene_symbol" in df.columns:
        df["at_least_one_isoform_in_gene_abs_activation_gte_2fold"] = df[
            "gene_symbol"
        ].map(m1h.groupby("gene_symbol")["abs_mean"].max() >= 1)
    else:
        df[f"at_least_one_isoform_in_gene{suffixes[0]}_abs_activation_gte_2fold"] = df[
            f"gene_symbol{suffixes[0]}"
        ].map(m1h.groupby("gene_symbol")["abs_mean"].max() >= 1)
        df[f"at_least_one_isoform_in_gene{suffixes[1]}_abs_activation_gte_2fold"] = df[
            f"gene_symbol{suffixes[1]}"
        ].map(m1h.groupby("gene_symbol")["abs_mean"].max() >= 1)

    df[f"activation{suffixes[0]}"] = df[f"clone_acc{suffixes[0]}"].map(
        m1h.set_index("clone_acc")["mean"]
    )
    df[f"activation{suffixes[1]}"] = df[f"clone_acc{suffixes[1]}"].map(
        m1h.set_index("clone_acc")["mean"]
    )
    df["activation_fold_change_log2"] = (
        df[f"activation{suffixes[1]}"] - df[f"activation{suffixes[0]}"]
    )
    df["activation_abs_fold_change_log2"] = df["activation_fold_change_log2"].abs()


def tissue_expression_similarity_per_isoform():
    """

    work in progress

    """
    df, metadata, genes = load_gtex_remapped()
    means = df.groupby(df.columns.map(metadata["body_site"]), axis=1).mean()
    similarities = []
    for gene in genes.unique():
        vals = means.loc[genes == gene, :].copy()
        # require at least 1 TPM at gene level per tissue
        vals.loc[:, vals.sum(axis=0) < 1] = np.nan
        # require at least 1 TPM in one tissue per isoform
        vals.loc[(vals >= 1).sum(axis=1) == 0, :] = np.nan
        pcc_corrected = vals.T.corr()
        for iso_a, iso_b in combinations(genes.loc[genes == gene].index.values, 2):
            similarities.append([iso_a, iso_b, pcc_corrected.loc[iso_a, iso_b]])
    similarities = pd.DataFrame(
        data=similarities, columns=["iso_a", "iso_b", "PCC_corrected"]
    )
    clones = load_valid_isoform_clones()

    def get_pair_id(row):
        return "_".join(
            sorted(
                [
                    [
                        x
                        for x in row["iso_a"].split(" ")[0].split("_")
                        if x in clones["clone_acc"].unique()
                    ][0],
                    [
                        x
                        for x in row["iso_b"].split(" ")[0].split("_")
                        if x in clones["clone_acc"].unique()
                    ][0],
                ]
            )
        )

    cloned = ~similarities["iso_a"].str.startswith("noclone") & ~similarities[
        "iso_b"
    ].str.startswith("noclone")
    similarities = similarities.loc[cloned, :]
    similarities["pair"] = similarities.apply(get_pair_id, axis=1)
    if similarities["pair"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    if not pairs["pair"].isin(similarities["pair"].values).all():
        raise UserWarning("unexpected missing pairs")
    pairs["tissue_expression_pcc"] = pairs["pair"].map(
        similarities.set_index("pair")["PCC_corrected"]
    )


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


def load_condensate_data():
    """
    Reference vs alternative isoforms table with additional info
    about condensates, restricted to pairs that were tested in the
    condensates assay.
    """
    df = pd.read_excel(DATA_DIR / "internal/TFiso_LLPS_scores_20231215.xlsx")
    df["gene_symbol"] = df["isoform_acc"].map(lambda x: x.split("|")[0])
    df["condensates_observed_HEK"] = df["LLPS_HEK293"] != 0
    df["condensates_observed_U2OS"] = df["LLPS_U2OS"] != 0
    df["is_cloned_reference"] = df["Ref_isoform"].map(
        {"Reference": True, "Alternative": False}
    )
    df = df.rename(
        columns={
            "LLPS_HEK293": "HEK_Condensate",
            "LLPS_U2OS": "U2OS_Condensate",
            "Localization_U2OS": "localization_U2OS",
        }
    )
    df["HEK_Condensate"] = (
        df["HEK_Condensate"]
        .map(
            lambda x: {
                "cc/nm": "BOTH",
                "few nc": "NC",
                "few cc": "CC",
            }.get(x, x)
        )
        .str.strip()
        .str.upper()
    )
    df["U2OS_Condensate"] = (
        df["U2OS_Condensate"]
        .map(
            lambda x: {
                "cc/nm": "BOTH",
                "few nc": "NC",
                "few cc": "CC",
            }.get(x, x)
        )
        .str.strip()
        .str.upper()
    )
    rename_loc = {"mostly in nucleus": "both", "in nucleus": "nucleus"}
    for cl in ["HEK", "U2OS"]:
        df[f"localization_{cl}"] = (
            df[f"localization_{cl}"]
            .str.strip()
            .str.lower()
            .map(lambda x: rename_loc.get(x, x))
        )
    if df["isoform_acc"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    df = df.set_index("isoform_acc")

    pairs = load_ref_vs_alt_isoforms_table()
    pairs.loc[
        (pairs["n_positive_PPI_ref"] == 0) | (pairs["n_positive_PPI_alt"] == 0),
        "PPI_jaccard",
    ] = np.nan
    pairs["PPI_Jaccard_d"] = 1 - pairs["PPI_jaccard"]
    pairs["PDI_Jaccard_d"] = 1 - pairs["PDI_jaccard"]
    for x in ["ref", "alt"]:
        for var in [
            "condensates_observed_HEK",
            "condensates_observed_U2OS",
            "HEK_Condensate",
            "U2OS_Condensate",
            "localization_HEK",
            "localization_U2OS",
        ]:
            pairs[var + "_" + x] = pairs["clone_acc_" + x].map(df[var])

    # Only take pairs with condensate info
    pairs = pairs.loc[
        pairs["condensates_observed_HEK_ref"].notnull()
        & pairs["condensates_observed_HEK_alt"].notnull(),
        :,
    ]

    for cl in ["HEK", "U2OS"]:
        pairs.loc[:, f"condensate_cat_{cl}"] = "Unchanged"
        pairs.loc[
            pairs[f"{cl}_Condensate_ref"].notnull()
            & pairs[f"{cl}_Condensate_alt"].isnull(),
            f"condensate_cat_{cl}",
        ] = "LOC"
        pairs.loc[
            pairs[f"{cl}_Condensate_ref"].isnull()
            & pairs[f"{cl}_Condensate_alt"].notnull(),
            f"condensate_cat_{cl}",
        ] = "GOC"
        pairs.loc[
            (
                (pairs[f"{cl}_Condensate_ref"] == "NC")
                & (pairs[f"{cl}_Condensate_alt"] == "CC")
            )
            | (
                (pairs[f"{cl}_Condensate_ref"] == "CC")
                & (pairs[f"{cl}_Condensate_alt"] == "NC")
            ),
            f"condensate_cat_{cl}",
        ] = "Changed localization"
    pairs["condensate_cat_merged_HEK"] = pairs["condensate_cat_HEK"].map(
        {
            "Unchanged": "No difference",
            "LOC": "Difference",
            "GOC": "Difference",
            "Changed localization": "Difference",
        }
    )
    pairs["condensate_cat_merged_U2OS"] = pairs["condensate_cat_U2OS"].map(
        {
            "Unchanged": "No difference",
            "LOC": "Difference",
            "GOC": "Difference",
            "Changed localization": "Difference",
        }
    )

    for cl in ["HEK", "U2OS"]:
        pairs[f"condensate_cat_only_{cl}"] = pairs[f"condensate_cat_{cl}"].map(
            {
                "Unchanged": "No difference",
                "LOC": "Difference",
                "GOC": "Difference",
                "Changed localization": "No difference",
            }
        )

    for cl in ["HEK", "U2OS"]:
        var = f"localization_cat_{cl}"
        pairs[var] = np.nan
        pairs.loc[
            pairs[f"localization_{cl}_ref"] == pairs[f"localization_{cl}_alt"], var
        ] = "No difference"
        pairs.loc[
            pairs[f"localization_{cl}_ref"] != pairs[f"localization_{cl}_alt"], var
        ] = "Difference"
        if pairs[var].isnull().any():
            raise UserWarning("bug in code")

    for cl in ["HEK", "U2OS"]:
        pairs[f"condensate_cat_only_detailed_{cl}"] = pairs[f"condensate_cat_{cl}"].map(
            {
                "Unchanged": "Both form condensates",
                "LOC": "Alternative loses condensate",
                "GOC": "Alernative gains condensate",
                "Changed localization": "Both form condensates",
            }
        )
        pairs.loc[
            (pairs[f"condensates_observed_{cl}_ref"] == False)
            & (pairs[f"condensates_observed_{cl}_alt"] == False),
            f"condensate_cat_only_detailed_{cl}",
        ] = "Neither form condensates"
    for cl in ["HEK", "U2OS"]:
        pairs[f"condensate_or_loc_change_{cl}"] = (
            (
                pairs[f"condensates_observed_{cl}_ref"]
                == pairs[f"condensates_observed_{cl}_alt"]
            )
            & (pairs[f"localization_{cl}_ref"] == pairs[f"localization_{cl}_alt"])
        ).map({True: "No difference", False: "Difference"})
    for cl in ["HEK", "U2OS"]:
        c = f"combined_cat_{cl}"
        pairs[c] = np.nan
        diff_loc = pairs[f"localization_{cl}_ref"] != pairs[f"localization_{cl}_alt"]
        diff_cond = (
            pairs[f"condensates_observed_{cl}_ref"]
            != pairs[f"condensates_observed_{cl}_alt"]
        )
        pairs.loc[diff_loc & ~diff_cond, c] = "Difference in localization"
        pairs.loc[~diff_loc & diff_cond, c] = "Difference in condensate formation"
        pairs.loc[diff_loc & diff_cond, c] = "Difference in both"
        pairs.loc[
            ~diff_loc & ~diff_cond, c
        ] = "Same localization and condensate formation"
        if pairs[c].isnull().any():
            raise UserWarning("Bug in code")
    pairs["condensate_or_loc_change_both"] = pairs["condensate_or_loc_change_HEK"]
    pairs.loc[
        pairs["condensate_or_loc_change_HEK"] != pairs["condensate_or_loc_change_U2OS"],
        "condensate_or_loc_change_both",
    ] = np.nan

    def detailed_condensate_cat(row, cl):
        a = row[f"{cl}_Condensate_ref"]
        if pd.isnull(a):
            a = "None"
        b = row[f"{cl}_Condensate_alt"]
        if pd.isnull(b):
            b = "None"
        return "{} -> {}".format(a, b)

    for cl in ["HEK", "U2OS"]:
        pairs[f"condensate_cat_detailed_{cl}"] = pairs.apply(
            detailed_condensate_cat, cl=cl, axis=1
        )

    return pairs, df
