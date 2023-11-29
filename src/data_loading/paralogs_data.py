from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import tqdm

from .clones_and_assays_data import (
    load_full_y2h_data_including_controls,
    load_y1h_pdi_data,
    load_m1h_activation_data,
    load_valid_isoform_clones,
    load_paralog_pairs_tested_in_y2h,
)
from .load_annotated_TFs import load_annotated_TFiso1_collection
from .utils import cache_with_pickle, DATA_DIR
from .isoform_pairwise_metrics import (
    _add_PPI_columns,
    _add_PDI_columns,
    _add_activation_columns,
    load_seq_id_between_cloned_isoforms,
)
from .protein_data import load_tf_families


@cache_with_pickle
def load_seq_id_between_cloned_genes():
    """
    Using reference isoforms. Used to define paralog pairs.

    """
    tfs = load_annotated_TFiso1_collection(include_single_isoform_genes=True)

    def pairwise_global_aa_sequence_similarity(aa_seq_a, aa_seq_b):
        alignment = pairwise2.align.globalds(
            aa_seq_a, aa_seq_b, blosum62, -10, -0.5, penalize_end_gaps=False
        )
        alignment = pairwise2.format_alignment(*alignment[0]).split("\n")[1]
        return alignment.count("|") / len(alignment) * 100

    aa_id = []
    tfs = {k: tfs[k] for k in sorted(tfs.keys())}
    for tf_a, tf_b in tqdm.tqdm(list(combinations(tfs.values(), 2))):
        # NOTE: here we're using the reference isoform, which may not be cloned
        aa_id.append(
            (
                tf_a.name,
                tf_b.name,
                "|".join(tf_a.reference_isoform.ensembl_transcript_names),
                "|".join(tf_b.reference_isoform.ensembl_transcript_names),
                pairwise_global_aa_sequence_similarity(
                    tf_a.reference_isoform.aa_seq, tf_b.reference_isoform.aa_seq
                ),
            )
        )
    aa_id = pd.DataFrame(
        aa_id,
        columns=[
            "gene_symbol_a",
            "gene_symbol_b",
            "isoform_a",
            "isoform_b",
            "aa_seq_pct_identity",
        ],
    )
    return aa_id


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


def load_paralogs_vs_isoforms_comparison_table():
    """

    TODO:
    - add any other columns???

    """

    fpath = (
        Path(__file__).resolve().parents[2] / "output" / "TF-paralogs-vs-isoforms.tsv"
    )
    if not fpath.exists():
        print("generating isoforms vs paralogs table")
        _write_TF_iso_vs_paralogs_table(fpath)
    return pd.read_csv(fpath, sep="\t")


def load_ensembl_compara_paralog_pairs_for_cloned_tfs():
    df = pd.read_csv(DATA_DIR / "external/Ensembl_TF_paralogs.txt", sep="\t")
    df = df.dropna()

    tfs = load_annotated_TFiso1_collection(include_single_isoform_genes=True)
    ensembl_to_hgnc = {tf.ensembl_gene_id: tf.name for tf in tfs.values()}

    df["gene1"] = df["Gene stable ID"].map(ensembl_to_hgnc)
    df["gene2"] = df["Human paralogue gene stable ID"].map(ensembl_to_hgnc)
    df = df.dropna()
    df["gene_symbol_a"] = df[["gene1", "gene2"]].min(axis=1)
    df["gene_symbol_b"] = df[["gene1", "gene2"]].max(axis=1)
    df = df.loc[:, ["gene_symbol_a", "gene_symbol_b"]].drop_duplicates()
    return df


def load_all_within_family_gene_pairs_for_cloned_tfs():
    tfs = load_annotated_TFiso1_collection(include_single_isoform_genes=True)
    df = pd.DataFrame(
        [(tf.name, tf.tf_family) for tf in tfs.values()],
        columns=["gene_symbol", "family"],
    )
    df = df.loc[df["family"] != "Unknown", :]
    # dealing with:
    # Homeodomain; POU
    # Homeodomain; Paired box
    # C2H2 ZF; AT hook
    # C2H2 ZF; BED ZF
    df["family"] = df["family"].str.split("; ")
    df = df.explode(column="family")
    pairs = []
    for fam in df["family"].unique():
        for a, b in combinations(df.loc[df["family"] == fam, "gene_symbol"].values, 2):
            pairs.append(list(sorted([a, b])))
    df = pd.DataFrame(pairs, columns=["gene_symbol_a", "gene_symbol_b"])
    df = df.drop_duplicates()
    return df


def load_all_different_family_gene_pairs_for_cloned_tfs():
    """
    TODO: delete this, since we didn't end up using it...
    """
    clones = load_valid_isoform_clones(include_single_isoform_genes=True)
    cloned_genes = (
        clones.sort_values("gene_symbol")["gene_symbol"].drop_duplicates().values
    )
    df = pd.DataFrame(
        combinations(cloned_genes, 2), columns=["gene_symbol_a", "gene_symbol_b"]
    )
    within = load_all_within_family_gene_pairs_for_cloned_tfs()
    df = df.loc[
        ~(df["gene_symbol_a"] + "_" + df["gene_symbol_b"]).isin(
            within["gene_symbol_a"] + "_" + within["gene_symbol_b"]
        ),
        :,
    ]

    if df.isnull().any().any():
        raise UserWarning("unexpected missing values")
    if df.duplicated().any():
        raise UserWarning("unexpected duplicates")
    if (df["gene_symbol_a"] == df["gene_symbol_b"]).any():
        raise UserWarning("should be different genes")
    return df


def load_non_paralog_control():
    df = load_ensembl_compara_paralog_pairs_for_cloned_tfs()
    fam = load_tf_families()
    df["families_a"] = df["gene_symbol_a"].map(fam).apply(lambda x: set(x.split("; ")))
    df["families_b"] = df["gene_symbol_b"].map(fam).apply(lambda x: set(x.split("; ")))
    np.random.seed(3489734)
    genes_a = list(
        df[["gene_symbol_a", "families_a"]].itertuples(name=None, index=False)
    )
    genes_b = list(
        df[["gene_symbol_b", "families_b"]].itertuples(name=None, index=False)
    )
    rnd_rows = []
    while len(rnd_rows) < df.shape[0]:
        to_add = []
        np.random.shuffle(genes_b)
        for (gene_a, families_a), (gene_b, families_b) in zip(genes_a, genes_b):
            if len(families_a.intersection(families_b)) == 0:
                to_add.append((gene_a, gene_b, families_a, families_b))
        for gene_a, gene_b, families_a, families_b in to_add:
            # NOTE: this just removes the first occurance
            genes_a.remove((gene_a, families_a))
            genes_b.remove((gene_b, families_b))
        rnd_rows += to_add
        if len(to_add) == 0:
            break
    rnd_rows = pd.DataFrame(data=rnd_rows, columns=df.columns)
    a = rnd_rows[["gene_symbol_a", "gene_symbol_b"]].min(axis=1)
    b = rnd_rows[["gene_symbol_a", "gene_symbol_b"]].max(axis=1)
    rnd_rows["gene_symbol_a"] = a
    rnd_rows["gene_symbol_b"] = b
    if (rnd_rows["families_a"] == rnd_rows["families_b"]).any():
        raise UserWarning("something went wrong")
    if (rnd_rows["gene_symbol_a"] == rnd_rows["gene_symbol_b"]).any():
        raise UserWarning("something went wrong")
    return rnd_rows.loc[:, ["gene_symbol_a", "gene_symbol_b"]]


def _write_TF_iso_vs_paralogs_table(fpath):
    df = load_ensembl_compara_paralog_pairs_for_cloned_tfs()
    df["is_paralog_pair"] = True
    ctrl = load_non_paralog_control()
    ctrl["is_paralog_pair"] = False
    df = pd.concat([df, ctrl])
    pairs = load_paralog_pairs_tested_in_y2h()
    pairs["is_tested_in_Y2H"] = True
    # NOTE: not including all pairs tested as paralogs in Y2H
    # since they don't meet the ensembl compara definition
    df = pd.merge(
        df,
        pairs,  # .loc[pairs["is_paralog_pair"], :],
        how="left",
        on=["gene_symbol_a", "gene_symbol_b", "is_paralog_pair"],
    )
    df = pd.concat([df, pairs.loc[~pairs["is_paralog_pair"], :]])
    df["is_tested_in_Y2H"] = df["is_tested_in_Y2H"].fillna(False)

    df = df.drop_duplicates(
        subset=["gene_symbol_a", "gene_symbol_b"]
    )  # only necessary because ZBTB48 ZIC3 in non-paralog control

    aa_id = load_seq_id_between_cloned_genes()
    paralog_gene_pairs = pd.merge(
        df,
        aa_id.loc[:, ["gene_symbol_a", "gene_symbol_b", "aa_seq_pct_identity"]],
        how="left",
        on=["gene_symbol_a", "gene_symbol_b"],
    )

    if not (
        paralog_gene_pairs["gene_symbol_a"] < paralog_gene_pairs["gene_symbol_b"]
    ).all():
        raise UserWarning("expected gene symbols to be sorted")
    if paralog_gene_pairs.isnull().any().any():
        raise UserWarning("unexpected missing values")

    y2h = load_full_y2h_data_including_controls()
    # need to include tested isoforms that showed no PDI hits
    y1h = load_y1h_pdi_data(add_missing_data=True, include_pY1H_data=False)
    if y1h["clone_acc"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    y1h = y1h.set_index("clone_acc")
    m1h = load_m1h_activation_data()
    m1h["mean"] = m1h[["M1H_rep1", "M1H_rep2", "M1H_rep3"]].mean(axis=1)
    m1h["abs_mean"] = m1h["mean"].abs()
    df = _pairs_of_paralogs_and_isoforms_comparison_table(
        paralog_pairs=paralog_gene_pairs,
        y2h=y2h,
        y1h=y1h,
        m1h=m1h,
    )
    if df.index.duplicated().any():
        raise UserWarning("unexpected duplicates")
    df.to_csv(fpath, index=False, sep="\t")


def _pairs_of_paralogs_and_isoforms_comparison_table(
    paralog_pairs,
    y2h,
    y1h,
    m1h,
):
    """

    Restricted to isoforms pairs in the paralogs dataset.

    Args:
        paralog_pairs ([type]): [description]
        y2h ([type]): [description]
        y1h ([type]): [description]
        m1h ([type]): [description]
        restrict_isoforms_to_those_with_paralogs (bool):


    Returns:
        [type]: [description]

    """
    if (paralog_pairs["gene_symbol_a"] == paralog_pairs["gene_symbol_b"]).any():
        raise ValueError("Gene incorrectly paired with itself as a paralog")

    tfs = load_annotated_TFiso1_collection(include_single_isoform_genes=True)
    pairs = paralog_pairs.copy()
    pairs["clone_acc_a"] = pairs["gene_symbol_a"].apply(
        lambda x: tfs[x].cloned_reference_isoform.clone_acc
    )
    pairs["clone_acc_b"] = pairs["gene_symbol_b"].apply(
        lambda x: tfs[x].cloned_reference_isoform.clone_acc
    )

    pairs["category"] = pairs["is_paralog_pair"].map(
        {True: "paralogs", False: "non-paralog-control"}
    )
    pairs = pairs.drop(columns=["is_paralog_pair"])

    iso_pairs = []
    for tf in tfs.values():
        ref_acc = tf.cloned_reference_isoform.clone_acc
        for alt_iso in tf.cloned_isoforms:
            if alt_iso.clone_acc == ref_acc:
                continue
            iso_pairs.append(
                (
                    tf.name,
                    tf.name,
                    ref_acc,
                    alt_iso.clone_acc,
                )
            )

    iso_pairs = pd.DataFrame(
        data=iso_pairs,
        columns=["gene_symbol_a", "gene_symbol_b", "clone_acc_a", "clone_acc_b"],
    )

    aa_ident = load_seq_id_between_cloned_isoforms()
    iso_pairs = pd.merge(
        iso_pairs,
        aa_ident.loc[:, ["clone_acc_a", "clone_acc_b", "aa_seq_pct_identity"]],
        how="left",
        on=["clone_acc_a", "clone_acc_b"],
    )
    if iso_pairs["aa_seq_pct_identity"].isnull().any():
        raise UserWarning("Unexpected missing sequence similarity values")

    iso_pairs["category"] = "isoforms"
    pairs = pd.concat([pairs, iso_pairs])
    pairs["pair"] = pairs.apply(
        lambda x: "_".join(sorted([x["clone_acc_a"], x["clone_acc_b"]])), axis=1
    )
    if pairs["aa_seq_pct_identity"].isnull().any():
        err_msg = "Problem with sequence similarity data."
        err_msg += f"\nMissing for: "
        err_msg += str(
            pairs.loc[
                pairs["aa_seq_pct_identity"].isnull(), ["clone_acc_a", "clone_acc_b"]
            ]
        )
        raise UserWarning(err_msg)

    family = {gene_symbol: tf.tf_family for gene_symbol, tf in tfs.items()}
    pairs["family_a"] = pairs["gene_symbol_a"].map(family)
    pairs["family_b"] = pairs["gene_symbol_b"].map(family)
    is_mane_cloned = {
        gene_symbol: tf.cloned_MANE_select_isoform for gene_symbol, tf in tfs.items()
    }
    pairs["is_MANE_select_isoform_cloned_a"] = pairs["gene_symbol_a"].map(
        is_mane_cloned
    )
    pairs["is_MANE_select_isoform_cloned_b"] = pairs["gene_symbol_b"].map(
        is_mane_cloned
    )
    pairs["is_MANE_select_isoform_cloned_both"] = (
        pairs["is_MANE_select_isoform_cloned_a"]
        & pairs["is_MANE_select_isoform_cloned_b"]
    )

    return _pairs_comparison_table(pairs, y2h, y1h, m1h)


def _pairs_comparison_table(pairs, y2h, y1h, m1h):
    _add_PPI_columns(pairs, y2h, suffixes=("_a", "_b"))
    _add_PDI_columns(pairs, y1h, suffixes=("_a", "_b"))
    _add_activation_columns(pairs, m1h, suffixes=("_a", "_b"))
    pairs = pairs.set_index("pair")
    # remove non-symmetric metrics, which might make sense for
    # ref vs alt isoforms but don't for gene a vs gene b paralogs
    pairs = pairs.drop(columns=["activation_fold_change_log2", "PPI_delta_n"])
    return pairs


def tissue_expression_similarity_per_gene():
    pass
