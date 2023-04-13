"""
We use GENCODE v30 for transcripts but we use the TF db HGNC symbols
which leads to 13 cases where the gene and transcript names aren't consistent
E.g. ZZZ3 has transcripts AC118549.1-202 and AC118549.1-204
"""

from pathlib import Path
import itertools
import warnings
from collections import defaultdict
import functools
import re
import pickle

import numpy as np
import pandas as pd
from Bio import SeqIO
import tqdm

from ccsblib import paros_connection  # TMP

import isolib


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
CACHE_DIR = Path(__file__).resolve().parents[1] / "cache"


def cache_with_pickle(func):
    """NOTE: only works for functions without arguments"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        cached_file = CACHE_DIR / (func.__name__ + ".pkl")
        # TODO: this only works if the underlying function takes kwargs
        use_cache = kwargs.get("use_cache", True)
        if cached_file.exists() and use_cache:
            with open(cached_file, "rb") as f:
                print("reading from cache")
                ret = pickle.load(f)
            return ret
        else:
            ret = func(*args, **kwargs)
            with open(cached_file, "wb") as f:
                pickle.dump(ret, f)
            return ret

    return wrapper


# TODO: this is a mess because there are circular function dependencies here
def load_valid_isoform_clones():
    """The subset of TF isoform clones that have passed the stringent
       annotation process (at-length AA match to GENCODE, or approved by
       Gloria/GENCODE team).

    Returns:
        set(str): clone accession IDs

    """
    df = pd.read_csv(DATA_DIR / "internal/isoform_clones.tsv", sep="\t")
    df["clone_name"] = df["clone_acc"].map(
        lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
    )
    y2h = pd.read_csv(
        DATA_DIR / "internal/Y2H-data_2022-03-08.tsv",
        sep="\t",
        na_values=[""],
        keep_default_na=False,
    )
    y1h = load_y1h_pdi_data(include_pY1H_data=False)
    m1h = load_m1h_activation_data()
    df["in_m1h"] = df["clone_acc"].isin(m1h["clone_acc"])
    df["in_y1h"] = df["clone_acc"].isin(y1h["unique_acc"])
    df["in_y2h"] = df["clone_name"].isin(
        y2h.loc[
            (y2h["category"] == "tf_isoform_ppis")
            & y2h["Y2H_result"].notnull(),  # y2h['score'].isin(['0', '1']),
            "ad_clone_name",
        ]
    )
    # dropping duplicates with identical AA seqs, keeping those with M1H data
    df = df.sort_values(
        ["gene", "in_m1h", "in_y2h", "in_y1h"], ascending=[True] + [False] * 3
    )
    df = df.loc[
        ~df["dup_idx"].duplicated() | df["dup_idx"].isnull(), ["gene", "clone_acc"]
    ].sort_values(["gene", "clone_acc"])

    nt_seq_file = DATA_DIR / "internal/j2_6k_unique_isoacc_and_nt_seqs.fa"
    nt = {r.id: str(r.seq) for r in SeqIO.parse(nt_seq_file, format="fasta")}
    df["cds"] = df["clone_acc"].map(nt)

    aa_seq_file = DATA_DIR / "internal/j2_6k_unique_isoacc_and_prot_seqs.fa"
    aa = {
        r.id.split("xxx")[1]: str(r.seq)
        for r in SeqIO.parse(aa_seq_file, format="fasta")
    }
    df["aa_seq"] = df["clone_acc"].map(aa)
    df["num_aa"] = df["aa_seq"].str.len()
    iso_annot = pd.read_csv(
        DATA_DIR / "internal/c_conso_annot_table_man_annot.tsv", sep="\t"
    )
    df["is_novel_isoform"] = df["clone_acc"].map(
        iso_annot.set_index("unique_acc")["gc_match"] == 0
    )
    exclude_both_strands_genes = ["FOXD4L3", "NANOG"]
    df = df.loc[~df["gene"].isin(exclude_both_strands_genes), :]
    df["clone_name"] = df["clone_acc"].map(
        lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
    )
    return df


def load_tf_isoform_y2h_screen_results():
    """There were two screens performed:

    The cloned TF isoforms as AD-fusions against DB-fusions of:
    (1) ORFeome 9.1
    (2) Subset of TFs and co-factors

    Returns:
        pandas.DataFrame: for each pair, was it found in the first and second screens

    """
    df = pd.read_csv(DATA_DIR / "internal/tf_isoform_y2h_screen.tsv", sep="\t")
    if not (df["in_orfeome_screen"] | df["in_focussed_screen"]).all():
        raise UserWarning("Something went wrong...")
    return df


def load_y2h_isoform_data(
    require_at_least_one_ppi_per_isoform=True,
    add_missing_data=False,
    filter_for_valid_clones=True,
    require_at_least_two_partners=False,
):
    """

    filtered dataset:
        remove partners with no positive tests per tf gene
        remove isoforms with no successful tests
        remove genes with < 2 isoforms

    Args:
        require_at_least_one_ppi_per_isoform (bool, optional): Defaults to True.
        add_missing_data (bool, optional): Defaults to False.
        filter_for_valid_clones (bool, optional): Defaults to True.

    Returns:
        pandas.DataFrame

    """
    y2h = load_isoform_and_paralog_y2h_data(add_missing_data, filter_for_valid_clones)
    ppi = y2h.loc[
        (y2h["category"] == "tf_isoform_ppis"),
        ["ad_clone_acc", "ad_gene_symbol", "db_gene_symbol", "Y2H_result"],
    ].copy()
    # at least one positive with an isoform per partner, for each TF gene
    ppi = ppi.loc[
        ppi.groupby(["ad_gene_symbol", "db_gene_symbol"])["Y2H_result"].transform(
            lambda row: (row == True).any()
        ),
        :,
    ]
    # at least one succussful test per isoform
    ppi = ppi.loc[
        ppi.groupby("ad_clone_acc")["Y2H_result"].transform(
            lambda x: (x.notnull().any())
        ),
        :,
    ]
    if require_at_least_one_ppi_per_isoform:
        ppi = ppi.loc[
            ppi.groupby("ad_clone_acc")["Y2H_result"].transform(
                lambda x: (x == True).any()
            ),
            :,
        ]
    # at least two isoforms per TF gene
    ppi = ppi.loc[
        ppi.groupby("ad_gene_symbol")["ad_clone_acc"].transform(
            lambda x: x.nunique() >= 2
        ),
        :,
    ]
    # successful tests with at least two isoforms of a TF gene, per partner
    ppi = ppi.loc[
        ppi.groupby(["ad_gene_symbol", "db_gene_symbol"])["Y2H_result"].transform(
            lambda x: x.notnull().sum() >= 2
        ),
        :,
    ]
    # require at least two successfully tested partners per TF gene
    if require_at_least_two_partners:
        ppi = ppi.loc[
            ppi["ad_gene_symbol"].map(
                (
                    ppi.loc[ppi["Y2H_result"].notnull()]
                    .groupby(["ad_gene_symbol", "db_gene_symbol"])
                    .size()
                    >= 2
                )
                .groupby("ad_gene_symbol")
                .sum()
                >= 2
            ),
            :,
        ]
    return ppi


def load_y2h_paralogs_additional_data():
    """Pairs tested in Y2H for the paralogs data, in addition to isoform pairs."""
    y2h = load_isoform_and_paralog_y2h_data()
    pairs = load_paralog_pairs()
    y2h_paralog = y2h.loc[
        y2h["category"].isin(
            ["tf_paralog_ppis", "paralog_with_PDI", "non_paralog_control"]
        ),
        :,
    ].copy()
    pair_map = defaultdict(set)
    for _i, row in pairs.iterrows():
        a, b = row["tf_gene_a"], row["tf_gene_b"]
        pair_map[a].add(b)
        pair_map[b].add(a)

    def find_matching_gene(row):
        matches = pair_map[row["ad_gene_symbol"]]
        matches = set(
            y2h.loc[
                (y2h["category"] == "tf_isoform_ppis")
                & y2h["ad_gene_symbol"].isin(matches)
                & (y2h["db_gene_symbol"] == row["db_gene_symbol"]),
                "ad_gene_symbol",
            ].unique()
        )
        if len(matches) == 0:
            return np.nan
        else:
            return "|".join(matches)

    y2h_paralog["paired_tf_gene"] = y2h_paralog.apply(find_matching_gene, axis=1)

    gte2iso = (
        y2h.loc[y2h["category"] == "tf_isoform_ppis", :]
        .groupby("ad_gene_symbol")["ad_clone_acc"]
        .nunique()
        >= 2
    )
    gte2iso = set(gte2iso.index[gte2iso])
    y2h_paralog["at_least_2_isoforms"] = y2h_paralog["ad_gene_symbol"].isin(
        gte2iso
    ) & y2h_paralog["paired_tf_gene"].apply(
        lambda x: any(g in gte2iso for g in x.split("|")) if pd.notnull(x) else False
    )

    gte2partner = (
        y2h.loc[y2h["category"] == "tf_isoform_ppis", :]
        .groupby("ad_gene_symbol")["db_gene_symbol"]
        .nunique()
        >= 2
    )
    gte2partner = set(gte2partner.index[gte2partner])
    y2h_paralog["at_least_2_partners"] = y2h_paralog["ad_gene_symbol"].isin(
        gte2partner
    ) & y2h_paralog["paired_tf_gene"].apply(
        lambda x: any(g in gte2partner for g in x.split("|"))
        if pd.notnull(x)
        else False
    )

    return y2h_paralog


def load_isoform_and_paralog_y2h_data(
    add_missing_data=False,
    add_partner_cateogories=False,
    remove_keratin_associated_proteins=True,
):
    """
    TODO: add_missing_data doesn't seem to do anything......

    remove_keratin_associated_proteins (bool): kertin associated proteins and
        the related late cornified envelope proteins are expressed in hair and
        skin, respectively. They bind many partners, and so their interactions
        are less likely to be relevant.

    """
    df = pd.read_csv(
        DATA_DIR / "internal/Y2H-data_2022-03-08.tsv",  # 'internal/y2h_data.tsv',
        sep="\t",
        na_values=[""],
        keep_default_na=False,
    )
    valid_clones = load_valid_isoform_clones()
    df["ad_clone_acc"] = df["ad_clone_name"].map(
        valid_clones.set_index("clone_name")["clone_acc"]
    )

    if remove_keratin_associated_proteins:
        df = df.loc[
            ~(
                df["db_gene_symbol"].str.startswith("KRTAP")
                | df["db_gene_symbol"].str.startswith("LCE")
            ),
            :,
        ]

    if add_missing_data:
        all_possible_ints = (
            pd.merge(
                df.loc[
                    df["category"] == "tf_isoform_ppis",
                    ["ad_gene_symbol", "db_gene_symbol", "category"],
                ].drop_duplicates(),
                valid_clones,
                left_on="ad_gene_symbol",
                right_on="gene",
            )
            .drop(columns="gene")
            .rename(columns={"clone_name": "ad_clone_name"})
        )
        df = pd.merge(
            df,
            all_possible_ints,
            on=["ad_gene_symbol", "db_gene_symbol", "ad_clone_name", "category"],
            how="outer",
        )
    if add_partner_cateogories:
        cat_info = load_ppi_partner_categories()
        cats = cat_info.groupby("category")["partner"].apply(set).to_dict()
        for cat, members in cats.items():
            df["is_partner_category_" + "_".join(cat.split())] = df[
                "db_gene_symbol"
            ].isin(members)
        cofac_type = cat_info.groupby("cofactor_type")["partner"].apply(set).to_dict()
        for subtype, members in cofac_type.items():
            df["is_cofactor_subtype_" + subtype] = df["db_gene_symbol"].isin(members)
    return df


def load_y1h_pdi_data(add_missing_data=False, include_pY1H_data=True):
    df = pd.read_csv(DATA_DIR / "internal/a2_juan_pdi_w_unique_isoacc.tsv", sep="\t")
    df = (
        pd.concat([df.loc[:, ["tf", "unique_acc"]], pd.get_dummies(df["bait"])], axis=1)
        .groupby(["tf", "unique_acc"])
        .sum()
        > 0
    ).reset_index()
    zeros = pd.read_csv(DATA_DIR / "internal/a2_juan_isoforms_wo_pdi.tsv", sep="\t")
    df = pd.concat([df, zeros], axis=0, sort=False).reset_index(drop=True).fillna(False)
    if add_missing_data:
        isoforms = load_valid_isoform_clones()
        df = pd.merge(
            df,
            isoforms.loc[:, ["gene", "clone_acc"]].rename(
                columns={"gene": "tf", "clone_acc": "unique_acc"}
            ),
            on=["tf", "unique_acc"],
            how="outer",
        )
    df[df.columns[2:]] = df[df.columns[2:]].astype("boolean")
    if include_pY1H_data:
        pY1H = load_additional_PDI_data_from_unpaired_cases_in_paired_Y1H_experiment()
        df = pd.merge(df, pY1H, on=["tf", "unique_acc"], how="outer")
    df = df.sort_values(["tf", "unique_acc"])
    return df


def load_additional_PDI_data_from_unpaired_cases_in_paired_Y1H_experiment():
    clones = load_valid_isoform_clones().set_index("clone_name")["clone_acc"]

    tested_genes = ["MAX", "PPARG", "RARG", "RXRG", "STAT1", "STAT3"]
    df = pd.concat(
        pd.read_excel(
            DATA_DIR / "internal/pY1H isoform single-TF results v2.xlsx",
            sheet_name=gene,
            usecols=[4, 5, 6, 7],
        )
        for gene in tested_genes
    )
    df["interaction"] = (
        df["interaction"]
        .fillna("False")
        .map(
            {
                "False": False,
                "moderate/strong": True,
                "very strong": True,
                "inconclusive": np.nan,
            }
        )
    )
    if (df.groupby("bait gene")["bait #"].nunique() != 1).any():
        raise UserWarning("unexpected multiple baits per gene")
    df["isoform"] = df["isoform"].map(lambda x: {"MAX": "MAX-1"}.get(x, x))
    df["tf"] = df["isoform"].apply(lambda x: x.split("-")[0])
    df["unique_acc"] = df["isoform"].map(clones)
    df = df.pivot(index=["tf", "unique_acc"], columns="bait gene", values="interaction")
    df = df.astype("boolean")
    df = df.reset_index()
    return df


def latest_y2h_data():
    qry = """      select b.large_plate_name,
            retest_pla, retest_pos,
            ad_clone_acc, ad_orf_id,
                db_gene_symbol, db_orf_id,
                    b.manual_score_growth as 3AT,
                            d.manual_score_growth as growth_score_LW_day1,
                    e.manual_score_growth as growth_score_LW_day2
        from tf_screen.final_retest as a
        inner join (select * from iscore.growth_result where large_plate_name like 'TFfnlr07%') as b
        on a.retest_pla = CONVERT(SUBSTRING(b.scoring_pla, -3), SIGNED)
        and a.retest_pos = b.scoring_pos
        
        LEFT join (select * from iscore.growth_result where large_plate_name like 'TFfnlr02%') as d
        on a.retest_pla = CONVERT(SUBSTRING(d.scoring_pla, -3), SIGNED)
        and a.retest_pos = d.scoring_pos
        
            LEFT join (select * from iscore.growth_result where large_plate_name like 'TFfnlr06%') as e
        on a.retest_pla = CONVERT(SUBSTRING(e.scoring_pla, -3), SIGNED)
        and a.retest_pos = e.scoring_pos

        
        where ad_clone_acc is not NULL;"""
    df = pd.read_sql(qry, paros_connection())
    df["plate_is_portrait_orientation"] = df["retest_pla"].map(
        df.loc[df["ad_clone_acc"] == "control-1", :]
        .drop_duplicates()
        .set_index("retest_pla")["retest_pos"]
        != "H01"
    )
    df = df.loc[~df["ad_clone_acc"].str.startswith("control-"), :]
    df["growth_score_LW_day2"] = df["growth_score_LW_day2"].fillna(
        df["growth_score_LW_day1"]
    )
    df = df.drop(columns=["growth_score_LW_day1"]).rename(
        columns={"growth_score_LW_day2": "LW"}
    )
    df = df.drop_duplicates()

    old = load_isoform_and_paralog_y2h_data()
    old = old.loc[~old["category"].isin(["rrs_paralogs", "lit_bm_paralogs"]), :]
    df = pd.merge(
        df,
        old.loc[:, ["ad_orf_id", "db_orf_id", "category"]],
        on=["ad_orf_id", "db_orf_id"],
        how="left",
    )
    wrong_cat = (
        df["category"].notnull()
        & ~df["category"].isin({"lit_bm_isoforms", "rrs_isoforms"})
        & ~df["plate_is_portrait_orientation"]
        & (df["retest_pos"].str.slice(0, 1) == "H")
    )

    wrong_cat = wrong_cat | (
        df["category"].notnull()
        & ~df["category"].isin({"lit_bm_isoforms", "rrs_isoforms"})
        & df["plate_is_portrait_orientation"]
        & df["retest_pos"].str.slice(1).isin({"10", "11", "12"})
    )
    wrong_cat = wrong_cat | (
        df["category"].notnull()
        & df["category"].isin({"lit_bm_isoforms", "rrs_isoforms"})
        & ~df["plate_is_portrait_orientation"]
        & (df["retest_pos"].str.slice(0, 1) != "H")
    )
    wrong_cat = wrong_cat | (
        df["category"].notnull()
        & df["category"].isin({"lit_bm_isoforms", "rrs_isoforms"})
        & df["plate_is_portrait_orientation"]
        & ~df["retest_pos"].str.slice(1).isin({"10", "11", "12"})
    )
    df = df.loc[~wrong_cat, :]
    if (df["category"].isnull() & (df["ad_orf_id"] > 0)).sum() > 0:
        raise UserWarning("missing categories")

    if df.duplicated(["retest_pla", "retest_pos"]).any():
        raise UserWarning("Unexpected duplicates")

    def find_empty_ad(row, media):
        if row["category"] not in {"lit_bm_isoforms", "rrs_isoforms"}:
            matches = df.loc[
                (df["ad_clone_acc"] == "empty-AD")
                & (df["retest_pla"] == row["retest_pla"])
                & (df["db_orf_id"] == row["db_orf_id"]),
                :,
            ]
            if matches.shape[0] > 0:
                vals = set(matches[media].values)
                if len(vals) > 1:
                    matches = df.loc[
                        (df["ad_clone_acc"] == "empty-AD")
                        & (df["retest_pla"] == row["retest_pla"])
                        & (df["retest_pos"].str.slice(1) == row["retest_pos"][1:])
                        & (df["retest_pos"] > row["retest_pos"])
                        & (df["db_orf_id"] == row["db_orf_id"]),
                        :,
                    ]
                    return matches.sort_values("retest_pos")[media].values[0]
                return list(vals)[0]
        else:
            matches = df.loc[
                (df["ad_clone_acc"] == "empty-AD")
                & (df["db_orf_id"] == row["db_orf_id"])
                & (
                    (
                        ~df["plate_is_portrait_orientation"]
                        & (df["retest_pos"].str.slice(0, 1) == "H")
                    )
                    | (
                        df["plate_is_portrait_orientation"]
                        & df["retest_pos"].str.slice(1).isin({"10", "11", "12"})
                    )
                ),
                :,
            ]
            if matches.shape[0] > 0:
                vals = set(matches[media].values)
                if len(vals) > 1:
                    print("error", row["retest_pla", "retest_pos"])
                return list(vals)[0]

    df["empty_AD_3AT"] = df.loc[df["ad_orf_id"] > 0].apply(
        find_empty_ad, args=["3AT"], axis=1
    )
    df["empty_AD_LW"] = df.loc[df["ad_orf_id"] > 0].apply(
        find_empty_ad, args=["LW"], axis=1
    )
    df = df.loc[df["ad_orf_id"] > 0, :]

    def score_pair(row):
        if any(
            row[spot] is None for spot in ["3AT", "LW", "empty_AD_3AT", "empty_AD_LW"]
        ):
            return np.nan
        if any(
            row[spot] == "NA" for spot in ["3AT", "LW", "empty_AD_3AT", "empty_AD_LW"]
        ):
            return np.nan
        score_3at, score_lw, score_empty_ad_3at, score_empty_ad_lw = (
            int(row[spot]) for spot in ["3AT", "LW", "empty_AD_3AT", "empty_AD_LW"]
        )
        if score_lw < 2 or score_empty_ad_lw < 2:
            return np.nan
        if score_empty_ad_3at == 4:
            return np.nan
        if score_3at <= 1:
            return False
        if (score_3at > score_empty_ad_3at) and (score_3at >= 2):
            return True
        if score_3at <= score_empty_ad_3at:
            return False
        raise UserWarning("Unexpected case:", row)

    df["Y2H_result"] = df.apply(score_pair, axis=1)

    df["ad_gene_symbol"] = df["ad_clone_acc"].apply(lambda x: x[: x.rfind("-")])
    df = df.loc[
        :,
        [
            "large_plate_name",
            "retest_pla",
            "retest_pos",
            "ad_gene_symbol",
            "ad_clone_acc",
            "ad_orf_id",
            "db_gene_symbol",
            "db_orf_id",
            "category",
            "3AT",
            "LW",
            "empty_AD_3AT",
            "empty_AD_LW",
            "Y2H_result",
        ],
    ]

    return df


def load_m1h_activation_data(add_missing_data=False):
    """ """
    df = pd.read_csv(DATA_DIR / "internal/a_m1h_final_table.tsv", sep="\t")
    df = df.rename(columns={"pos_acc": "clone_acc"})
    for column in df.columns:
        if column.startswith("M1H_rep"):
            df[column] = np.log2(df[column])
    if add_missing_data:
        isoforms = load_valid_isoform_clones()
        df = pd.merge(df, isoforms, on=["gene", "clone_acc"], how="outer")
    df = df.sort_values(["gene", "clone_acc"])
    return df


@cache_with_pickle
def load_gtex_gencode():
    """GTEx mapped to TFs in GENCODE v30, i.e. not including our novel isoforms"""
    df = pd.read_csv(
        DATA_DIR / "processed/expression_2022-09-01/transcript.GTEx-GC30_Isoforms.txt",
        sep="\t",
    )
    metadata = pd.read_csv(
        DATA_DIR / "processed/gtex_2022/GTEx_SRARunTable.txt", sep="\t"
    )
    if metadata["Run"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    metadata = metadata.set_index("Run")
    df = df.set_index("UID")
    df = (df + 1.0).apply(np.log2)

    def extract_ensembl_gene_id(s):
        return s.split("|")[1].split(".")[0]

    genes = pd.Series(
        index=df.index,
        data=df.index.map(extract_ensembl_gene_id).values,
    )
    tfs = load_annotated_gencode_tfs()
    ensembl_gene_id_to_hgnc_symbol = {
        tf.ensembl_gene_id: tf.name for tf in tfs.values()
    }
    tfs = {tf.ensembl_gene_id: tf for tf in tfs.values()}
    df = df.loc[genes.isin({tf.ensembl_gene_id for tf in tfs.values()}), :]
    if genes[df.index].nunique() != len(tfs):
        raise UserWarning("Unexpected missing genes")
    df = df.loc[genes != "PCGF6", :]
    df.index = df.index.map(lambda x: _convert_to_merged_protein_isoform_ids(x, tfs))
    df = df.loc[df.index != "fail", :]
    df = df.groupby("UID").sum()
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.orfs]):
        raise UserWarning("Something went wrong")
    genes = genes.loc[genes.isin({tf.ensembl_gene_id for tf in tfs.values()})]
    genes.index = genes.index.map(
        lambda x: _convert_to_merged_protein_isoform_ids(x, tfs)
    )
    genes = genes[~genes.index.duplicated(keep="first")]
    genes = genes.loc[genes.index.isin(df.index.values)].map(
        ensembl_gene_id_to_hgnc_symbol
    )
    if not df.columns.isin(metadata.index).all():
        raise UserWarning("Missing metadata")
    metadata = metadata.loc[metadata.index.isin(df.columns), :]
    return df, metadata, genes


def _convert_to_merged_protein_isoform_ids(s, tfs):
    """
    Map an ensembl transcript of a TF gene to the IDs when merging
    transcripts with identical protein sequences

    tfs: dict with keys of ensembl_gene_id
    """
    ensembl_trancript_name = s.split("|")[4]
    ensembl_gene_id = s.split("|")[1].split(".")[0]
    if ensembl_trancript_name not in tfs[ensembl_gene_id]:
        # print('failed for key: ' + ensembl_trancript_name)
        # these are not in GENCODE Basic set, i.e. not reliable transcripts
        return "fail"
    return "_".join(
        sorted(tfs[ensembl_gene_id][ensembl_trancript_name].ensembl_transcript_names)
    )


def _convert_to_joint_clone_and_ensembl_id(s, tfs):
    """
    Joint clone and ensembl transcript IDs:
    clone then ensembl, seperated by space
    underscores join multiple IDs, sorted,
    noclone and nomatch
    """
    if len(s.split("|")) == 3:  # novel clones
        return s + " nomatch"  # TODO: check duplicates
    elif len(s.split("|")) in [9, 10, 11]:
        ensembl_trancript_name = s.split("|")[4]
        gene = s.split("|")[5]
        if gene == "AC092072.1":  # HACK
            gene = "ZNF223"
        if ensembl_trancript_name not in tfs[gene]:
            # print('failed for key: ' + ensembl_trancript_name)
            # these are not in GENCODE Basic set, i.e. not reliable transcripts
            return "fail"
        iso = tfs[gene][ensembl_trancript_name]
        if hasattr(iso, "clone_acc"):
            clone_acc = iso.clone_acc
        else:
            clone_acc = "noclone"
        return (
            clone_acc
            + " "
            + "_".join(
                sorted(tfs[gene][ensembl_trancript_name].ensembl_transcript_names)
            )
        )
    else:
        raise UserWarning("Failed to map ID from " + s)


@cache_with_pickle
def load_gtex_remapped():
    """Load the GTEx data mapped to gencode v30 + our novel isoform clones

    Transcripts with identical amino acid sequences are merged.

    Requires a lot of processing because of problems with the mapping file
        - there are some clones that are no longer in our list
        - some cloned isoform + gencode isoform matched pairs are not merged
        - the identification is a mess

    """
    metadata = pd.read_csv(
        DATA_DIR / "processed/gtex_2022/GTEx_SRARunTable.txt", sep="\t"
    )
    df = pd.read_csv(
        DATA_DIR / "processed/gtex_2022/transcript.BreastCancer.txt", sep="\t"
    )
    if metadata["Run"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    metadata = metadata.set_index("Run")
    df = df.set_index("UID")
    df = (df + 1.0).apply(np.log2)

    def extract_gene_name(s):
        fields = s.split("|")
        if len(fields) == 3:
            return fields[0]
        elif len(fields) in {9, 10, 11}:
            return fields[5]
        else:
            raise UserWarning("Failed to extract gene name from " + s)

    genes = pd.Series(index=df.index, data=df.index.map(extract_gene_name).values)
    genes = genes.map(
        lambda x: "ZNF223" if x == "AC092072.1" else x
    )  # HACK for renamed gene
    clones = load_valid_isoform_clones()
    df = df.loc[genes.isin(clones["gene"].unique()), :]
    if genes[df.index].nunique() != clones["gene"].nunique():
        raise UserWarning("Unexpected missing genes")
    df = df.loc[genes != "PCGF6", :]
    tfs = load_annotated_6k_collection()
    df.index = df.index.map(lambda x: _convert_to_joint_clone_and_ensembl_id(x, tfs))
    df = df.loc[df.index != "fail", :]
    df = df.groupby("UID").sum()
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.orfs]):
        raise UserWarning("Something went wrong")

    def extract_gene_name_from_joint_id(s):
        if s.split(" ")[0] != "noclone":
            return s.split(" ")[0].split("|")[0]
        else:
            return s.split(" ")[1].split("_")[0][:-4]

    genes = pd.Series(
        index=df.index, data=df.index.map(extract_gene_name_from_joint_id).values
    )
    if not df.columns.isin(metadata.index).all():
        raise UserWarning("Missing metadata")
    metadata = metadata.loc[metadata.index.isin(df.columns), :]
    return df, metadata, genes


@cache_with_pickle
def load_developmental_tissue_expression_gencode():
    """
    Cardoso-Moreira et al. Nature


    """
    metadata = pd.read_csv(DATA_DIR / "processed/Cardoso-Moreira_et_al/metadata.txt")
    df = pd.read_csv(
        DATA_DIR
        / "processed/expression_2022-09-01/transcript.Cardoso-Moreira-et-al-Nature-2019-GC30_Isoforms.txt",
        sep="\t",
    )
    if metadata["Run"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    metadata = metadata.set_index("Run")
    df = df.set_index("UID")
    df = (df + 1.0).apply(np.log2)

    def extract_ensembl_gene_id(s):
        return s.split("|")[1].split(".")[0]

    genes = pd.Series(
        index=df.index,
        data=df.index.map(extract_ensembl_gene_id).values,
    )
    tfs = load_annotated_gencode_tfs()
    ensembl_gene_id_to_hgnc_symbol = {
        tf.ensembl_gene_id: tf.name for tf in tfs.values()
    }
    tfs = {tf.ensembl_gene_id: tf for tf in tfs.values()}
    df = df.loc[genes.isin({tf.ensembl_gene_id for tf in tfs.values()}), :]
    if genes[df.index].nunique() != len(tfs):
        raise UserWarning("Unexpected missing genes")
    df = df.loc[genes != "PCGF6", :]
    df.index = df.index.map(lambda x: _convert_to_merged_protein_isoform_ids(x, tfs))
    df = df.loc[df.index != "fail", :]
    df = df.groupby("UID").sum()
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.orfs]):
        raise UserWarning("Something went wrong")
    genes = genes.loc[genes.isin({tf.ensembl_gene_id for tf in tfs.values()})]
    genes.index = genes.index.map(
        lambda x: _convert_to_merged_protein_isoform_ids(x, tfs)
    )
    genes = genes[~genes.index.duplicated(keep="first")]
    genes = genes.loc[genes.index.isin(df.index.values)].map(
        ensembl_gene_id_to_hgnc_symbol
    )

    # the file has ERR2598062.fastq.gz instead of ERR2598062
    df.columns = df.columns.str.slice(0, -len(".fastq.gz"))

    if not df.columns.isin(metadata.index).all():
        raise UserWarning("Missing metadata")
    metadata = metadata.loc[metadata.index.isin(df.columns), :]
    return df, metadata, genes


@cache_with_pickle
def load_developmental_tissue_expression_remapped():
    """
    Cardoso-Moreira et al. Nature, remapped to include our novel isoforms


    """
    data_dir = DATA_DIR / "processed/Cardoso-Moreira_et_al"
    metadata = pd.read_csv(data_dir / "metadata.txt")
    df = pd.read_csv(data_dir / "transcript.Hs-Dev-Timecourse.txt", sep="\t")
    if metadata["Run"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    metadata = metadata.set_index("Run")
    df = df.set_index("UID")
    df = (df + 1.0).apply(np.log2)

    def extract_gene_name(s):
        fields = s.split("|")
        if len(fields) == 3:
            return fields[0]
        elif len(fields) in {9, 10, 11}:
            return fields[5]
        else:
            raise UserWarning("Failed to extract gene name from " + s)

    genes = pd.Series(index=df.index, data=df.index.map(extract_gene_name).values)
    genes = genes.map(
        lambda x: "ZNF223" if x == "AC092072.1" else x
    )  # HACK for renamed gene
    clones = load_valid_isoform_clones()
    df = df.loc[genes.isin(clones["gene"].unique()), :]
    if genes[df.index].nunique() != clones["gene"].nunique():
        raise UserWarning("Unexpected missing genes")
    df = df.loc[genes != "PCGF6", :]
    tfs = load_annotated_6k_collection()
    df.index = df.index.map(lambda x: _convert_to_joint_clone_and_ensembl_id(x, tfs))
    df = df.loc[df.index != "fail", :]
    df = df.groupby("UID").sum()
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.orfs]):
        raise UserWarning("Something went wrong")

    def extract_gene_name_from_joint_id(s):
        if s.split(" ")[0] != "noclone":
            return s.split(" ")[0].split("|")[0]
        else:
            return s.split(" ")[1].split("_")[0][:-4]

    genes = pd.Series(
        index=df.index, data=df.index.map(extract_gene_name_from_joint_id).values
    )
    return df, metadata, genes


def load_seq_comparison_data():
    """
    Pairwise sequence comparisons of AA.
    Needleman algorithm, global alignment.
    Note - I checked and there are no duplicate rows.
    """
    df = pd.read_table(
        DATA_DIR
        / "processed/a_2019-09-10_AA_seq_identity_Isoform_series_for_all_6Kpairs_unique_acc.txt"
    )
    df["pair"] = df.apply(lambda x: "_".join(sorted([x.iso1, x.iso2])), axis=1)
    df = df[["pair", "AAseq_identity%"]]
    df.columns = ["pair", "aa_seq_pct_id"]

    df_b = pd.read_csv(
        DATA_DIR / "processed/paralog_non_paralog_seq_id.tsv", sep="\t"
    ).drop_duplicates()
    df_b["pair"] = df_b.apply(
        lambda x: "_".join(sorted([x["clone_acc_a"], x["clone_acc_b"]])), axis=1
    )
    df_b = df_b.loc[:, ["pair", "aa_seq_pct_id"]]
    df = pd.concat([df, df_b])

    if df["pair"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    df = df.set_index("pair")
    if (df["aa_seq_pct_id"] < 0).any() or (df["aa_seq_pct_id"] > 100).any():
        raise UserWarning("Percent values outside 0-100")
    return df["aa_seq_pct_id"]


def load_paralog_pairs(filter_for_valid_clones=True):
    """Pairs of TF gene paralogs and non-paralogs that were tested in Y2H pairwise tests.

    TODO: change name, since it's paralog pairs tested in Y2H

    WARNING: the aa sequence identity is just whatever is first in the file.
    Need to settle on which one to use.

    Returns:
        pandas.DataFrame: one row for each pair

    """
    df = pd.read_csv(
        DATA_DIR / "internal/a_tf_iso_paralog_nonparalogs_tested.tsv", sep="\t"
    )
    df["is_paralog_pair"] = df["cat2"] == "paralog"
    aa = pd.read_csv(
        DATA_DIR
        / "processed/b_2018-11-30_AA_seq_identity_Paralog_comparisons_unique_acc.txt",
        sep="\t",
    )
    if not (df["tf1"] < df["tf2"]).all():
        raise UserWarning("Expected genes to be ordered")
    (
        aa[["gene1", "gene2"]].min(axis=1) + "_" + aa[["gene1", "gene2"]].max(axis=1)
    ).duplicated().any()
    aa["tf1"] = aa[["gene1", "gene2"]].min(axis=1)
    aa["tf2"] = aa[["gene1", "gene2"]].max(axis=1)
    df = (
        pd.merge(df, aa, how="left", on=["tf1", "tf2"])
        .loc[:, ["tf1", "tf2", "is_paralog_pair", "AAseq_identity%"]]
        .rename(
            columns={
                "tf1": "tf_gene_a",
                "tf2": "tf_gene_b",
                "AAseq_identity%": "pct_aa_seq_identity",
            }
        )
    )
    df = df.drop_duplicates()
    if filter_for_valid_clones:
        valid_clones = load_valid_isoform_clones()
        df = df.loc[
            df["tf_gene_a"].isin(valid_clones["gene"])
            & df["tf_gene_b"].isin(valid_clones["gene"]),
            :,
        ]
    if (df["tf_gene_a"] == df["tf_gene_b"]).any():
        raise ValueError("Same gene twice, should be two different genes")
    return df


def load_ppi_partner_categories():
    """Juan's manual classification of the PPI interaction partners.

     Note that a gene can be in multiple categories.

    Returns:
        pandas.DataFrame: gene and category

    """
    df = pd.read_excel(
        DATA_DIR / "internal/20191028- Uniprot functions for interactors.xlsx",
        sheet_name="Final",
    )
    if df["Function class"].isnull().any():
        raise UserWarning("Unexpected missing values")
    if df["partner"].duplicated().any():
        raise UserWarning("Unexpected duplicate entries")
    df = df.set_index("partner")
    cofac_type = df["Cofactor type?"].copy()
    df = (
        df["Function class"]
        .str.split(", ", expand=True)
        .stack()
        .reset_index()
        .loc[:, ["partner", 0]]
        .rename(columns={0: "category"})
    )
    df["category"] = df["category"].str.strip()
    cf_rows = df["category"] == "cofactor"
    df.loc[cf_rows, "cofactor_type"] = df.loc[cf_rows, "partner"].map(cofac_type)
    if not (df["partner"] == df["partner"].str.strip()).all():
        raise UserWarning("Possibly something wrong with gene names column")
    return df


def load_tf_families():
    """From the Lambert et al. review in Cell 2018

    Returns:
        pandas.Series: HGNC gene symbol to TF family

    """
    tf_fam = pd.read_csv(DATA_DIR / "external/Human_TF_DB_v_1.01.csv")
    tf_fam = tf_fam.loc[tf_fam["Is TF?"] == "Yes", ["HGNC symbol", "DBD"]]
    if tf_fam["HGNC symbol"].duplicated().any():
        raise UserWarning("Unexpected duplicates")
    tf_fam = tf_fam.set_index("HGNC symbol")["DBD"]
    return tf_fam


def load_human_tf_db():
    """From the Lambert et al. Cell 2018

    Returns:
        pandas.DataFrame

    """
    tf_db = pd.read_csv(DATA_DIR / "external/Human_TF_DB_v_1.01.csv")
    if tf_db["Is TF?"].isnull().any() or ~tf_db["Is TF?"].isin({"Yes", "No"}).all():
        raise UserWarning("Problem with column")
    # file contains all 2,765 proteins examined, of which 1,639 are classified as TFs
    tf_db["Is TF?"] = tf_db["Is TF?"].map({"Yes": True, "No": False})
    tf_db = tf_db.loc[tf_db["Is TF?"], :]
    return tf_db


def convert_old_id_to_new_id(old_id):
    if len(old_id.split("xxx")) != 3:
        raise UserWarning("Unrecognized old isoform clone ID")
    return old_id.split("xxx")[1]


def read_hmmer3_domtab(filepath):
    """Parser for HMMER 3 hmmscan with the --domtblout option.
    Args:
        filepath (str): location of file
    Returns:
        pandas.DataFrame: see the HMMER User's Guide for a description of the
                          columns.
    """
    columns = [
        "target name",
        "target accession",
        "tlen",
        "query name",
        "query accession",
        "qlen",
        "E-value",
        "score",
        "bias",
        "#",
        "of",
        "c-Evalue",
        "i-Evalue",
        "full_sequence_score",
        "full_sequence_bias",
        "hmm_coord_from",
        "hmm_coord_to",
        "ali_coord_from",
        "ali_coord_to",
        "env_coord_from",
        "env_coord_to",
        "acc",
        "description of target",
    ]
    ncols = len(columns)
    # can't just use pandas.read_csv because spaces are used both as the
    # delimeter and are present in the values of the last columns
    with open(filepath, "r") as f:
        data = []
        for line in f:
            if line.startswith("#"):
                continue
            s = line.split()
            data.append(s[: ncols - 1] + [" ".join(s[ncols - 1 :])])
    df = pd.DataFrame(data=data, columns=columns)
    int_cols = [
        "tlen",
        "qlen",
        "#",
        "of",
        "hmm_coord_from",
        "hmm_coord_to",
        "ali_coord_from",
        "ali_coord_to",
        "env_coord_from",
        "env_coord_to",
    ]
    float_cols = [
        "E-value",
        "score",
        "bias",
        "c-Evalue",
        "i-Evalue",
        "full_sequence_score",
        "full_sequence_bias",
        "acc",
    ]
    df[int_cols] = df[int_cols].astype(int)
    df[float_cols] = df[float_cols].astype(float)
    return df


def _load_pfam_domains(fpath, cutoff_seq=0.01, cutoff_dom=0.01):
    pfam = read_hmmer3_domtab(fpath)
    pfam["pfam_ac"] = pfam["target accession"].str.replace(r"\..*", "")
    pfam = pfam.loc[(pfam["E-value"] < cutoff_seq) & (pfam["c-Evalue"] < cutoff_dom), :]
    pfam = _remove_overlapping_domains(pfam)
    return pfam


def load_pfam_domains_6k():
    filtered_pfam_path = CACHE_DIR / "pfam_6K.tsv"
    if filtered_pfam_path.exists():
        return pd.read_csv(filtered_pfam_path, sep="\t")
    pfam = _load_pfam_domains(DATA_DIR / "processed/6K_hmmer_2020-04-16_domtabl.txt")
    pfam["query name"] = pfam["query name"].map(convert_old_id_to_new_id)
    pfam.to_csv(filtered_pfam_path, index=False, sep="\t")
    return pfam


def load_pfam_domains_gencode():
    filtered_pfam_path = CACHE_DIR / "pfam_gencode.v30.tsv"
    if filtered_pfam_path.exists():
        pfam = pd.read_csv(filtered_pfam_path, sep="\t")
        pfam["query name"] = pfam["query name"].str.slice(8)
        pfam["query name"] = pfam["query name"].apply(
            lambda x: "|".join(sorted(x.split("|")))
        )
        return pfam
    pfam = _load_pfam_domains(DATA_DIR / "processed/2020-04-23_GC30_6K_domtabl.txt")
    pfam.to_csv(filtered_pfam_path, index=False, sep="\t")
    pfam["query name"] = pfam["query name"].str.slice(8)  # remove 'GC_grp:_'
    pfam["query name"] = pfam["query name"].apply(
        lambda x: "|".join(sorted(x.split("|")))
    )
    return pfam


def _is_overlapping(dom_a, dom_b):
    if (
        dom_b["env_coord_from"] > dom_a["env_coord_to"]
        or dom_a["env_coord_from"] > dom_b["env_coord_to"]
    ):
        return False
    else:
        return True


def load_pfam_clans():
    """
    TODO: return table instead of dict

    """
    clans = pd.read_csv(DATA_DIR / "external/Pfam-A.clans.tsv", sep="\t", header=None)
    if clans.loc[:, 0].duplicated().any():
        raise UserWarning()
    clans = clans.iloc[:, [0, 1]].dropna().set_index(0)[1].to_dict()
    return clans


def _remove_overlapping_domains_from_same_clan(pfam_in):
    """
    Trying to follow the method in the pfam_scan pearl script from EBI.
    """
    pfam = pfam_in.copy()
    clans = load_pfam_clans()
    to_remove = set()
    pfam = pfam.sort_values(["query name", "E-value"])
    # This is a bit slow
    for iso_id in tqdm.tqdm(pfam["query name"].unique()):
        doms = pfam.loc[
            (pfam["query name"] == iso_id) & (pfam["pfam_ac"].isin(clans)), :
        ]
        for i, j in itertools.combinations(doms.index, 2):
            if i in to_remove or j in to_remove:
                continue
            dom_a = doms.loc[i, :]
            dom_b = doms.loc[j, :]
            if clans[dom_a["pfam_ac"]] == clans[dom_b["pfam_ac"]]:
                if _is_overlapping(dom_a, dom_b):
                    to_remove.add(j)
    pfam = pfam.drop(to_remove)
    return pfam


def _remove_overlapping_domains(pfam_in):
    """
    I got a lot of overlaps of domains not in the same clan, so trying
    removing all overlapping pfam domains.
    """
    pfam = pfam_in.copy()
    to_remove = set()
    pfam = pfam.sort_values(["query name", "E-value"])
    # This is a bit slow
    for iso_id in tqdm.tqdm(pfam["query name"].unique()):
        doms = pfam.loc[(pfam["query name"] == iso_id), :]
        for i, j in itertools.combinations(doms.index, 2):
            if i in to_remove or j in to_remove:
                continue
            dom_a = doms.loc[i, :]
            dom_b = doms.loc[j, :]
            if _is_overlapping(dom_a, dom_b):
                to_remove.add(j)
    pfam = pfam.drop(to_remove)
    return pfam


def load_DNA_binding_domains(add_additional_domains=True):
    """

    Note: the DBD names here are not exactly the same as the PFam domain or
    clan names.

    TODO: I'm not sure if adding the additional domains is a good idea,
    since some are e.g. non-DNA binding ZFs or domains of unknown function

    """
    dbd = pd.read_csv(
        DATA_DIR / "internal/a2_final_list_of_dbd_pfam_and_names_ZF_marked.txt",
        sep="\t",
    )
    clans = load_pfam_clans()
    dbd["clan"] = dbd["pfam"].map(clans)
    if not add_additional_domains:
        return dbd

    # hand adding missing domain see issue #61
    dbd = pd.concat(
        [
            dbd,
            pd.DataFrame(
                data=[
                    {
                        "dbd": "BTD",
                        "pfam": "PF09270",
                        "clan": clans.get("PF09270", np.nan),
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

    dbd_clans = {
        "CL0361",  # C2H2-ZF
        "CL0012",  # Histone (mostly DNA binding...)
        "CL0274",  # WRKY-GCM1
        "CL0114",  # HMG-box
        "CL0081",  # MBD-like
        "CL0073",  # P53-like
        "CL0407",  # TATA-Binding Protein like
        "CL0018",  # bZIP
    }
    if not all(c in dbd["clan"].unique() for c in dbd_clans):
        raise UserWarning()
    pfam_ac_to_name = (
        pd.read_csv(DATA_DIR / "external/Pfam-A.clans.tsv", sep="\t", header=None)
        .set_index(0)[4]
        .to_dict()
    )
    for pfam_id, clan_id in clans.items():
        if clan_id not in dbd_clans:
            continue
        if pfam_id not in dbd["pfam"].values:
            dbd = pd.concat(
                [
                    dbd,
                    pd.DataFrame(
                        data=[
                            {
                                "dbd": pfam_ac_to_name[pfam_id],
                                "pfam": pfam_id,
                                "clan": clan_id,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

    # TODO: figure out why that is there (PF02892)
    dbd = dbd.drop([56])  # duplicate entry in table

    return dbd


@functools.lru_cache()
def load_dbd_accessions():
    dbd = load_DNA_binding_domains()
    return set(dbd["pfam"].values).union(
        {"C2H2_ZF_array_" + str(i) for i in range(2, 30)}
    )


@cache_with_pickle
def load_annotated_6k_collection(
    path_6k_gtf=DATA_DIR / "internal/c_6k_unique_acc_aligns.gtf",
    path_6k_fa=DATA_DIR / "internal/j2_6k_unique_isoacc_and_nt_seqs.fa",
    path_gencode_aa_seq=DATA_DIR / "external/gencode.v30.pc_translations.fa",
    path_MANE_select=DATA_DIR / "external/MANE.GRCh38.v0.95.summary.txt",
    path_APPRIS=DATA_DIR
    / "external/APPRIS-annotations_human_GRCh38.p13_ensembl104.tsv",
    path_effector_domains=DATA_DIR
    / "external/Soto-et-al_MolCell_2022_Supplementary-tables.xlsx",
    path_disorder=DATA_DIR / "processed/TFiso1_disorder-and-ss_from-alphafold.tsv",
):
    """ """
    import pyranges

    # note that pyranges switches the indexing to python 0-indexed, half-open interval
    algn = pyranges.read_gtf(path_6k_gtf).df
    algn = algn.loc[algn["Start"] < algn["End"], :]  # filter out problematic rows
    nt_seq = {r.name: r for r in SeqIO.parse(path_6k_fa, format="fasta")}
    pfam = load_pfam_domains_6k()
    clones = load_valid_isoform_clones()
    algn = algn.loc[algn["transcript_id"].isin(clones["clone_acc"].unique()), :]
    nt_seq = {k: v for k, v in nt_seq.items() if k in clones["clone_acc"].unique()}
    tf_db = load_human_tf_db()
    hgnc_to_ensembl = tf_db.set_index("HGNC symbol")["Ensembl ID"].to_dict()
    ensembl_to_uniprot = pd.read_csv(
        DATA_DIR / "external/ensembl_to_uniprot_ids_human_tfs.tsv", sep="\t"
    )
    ensembl_to_uniprot = (
        ensembl_to_uniprot.groupby("ensembl_gene_id")["uniprot_ac"]
        .apply(lambda x: "/".join(x))
        .to_dict()
    )
    tf_family = load_tf_families()
    genes = {}
    for gene_name in algn["gene_id"].unique():
        if (
            gene_name == "PCGF6"
        ):  # has a 6nt insertion that doesn't map to reference genome
            continue
        orf_ids = algn.loc[algn["gene_id"] == gene_name, "transcript_id"].unique()
        missing = [orf for orf in orf_ids if orf not in nt_seq]
        if missing:
            raise ValueError(", ".join(missing) + " not in " + path_6k_fa)
        extra = [
            orf
            for orf in nt_seq.keys()
            if orf.split("|")[0] == gene_name and orf not in orf_ids
        ]
        if extra:
            warnings.warn(", ".join(extra) + " in fasta but not in gtf")
        isoforms = []
        for orf_id in orf_ids:
            exons = []
            columns = [
                "gene_id",
                "transcript_id",
                "Chromosome",
                "Strand",
                "Start",
                "End",
            ]
            for _i, row in algn.loc[
                algn["transcript_id"] == orf_id, columns
            ].iterrows():
                exons.append(isolib.Exon(*row.values))
            orf_name = orf_id.split("|")[0] + "-" + orf_id.split("|")[1].split("/")[0]
            isoforms.append(
                isolib.Cloned_Isoform(
                    clone_name=orf_name,
                    exons=exons,
                    clone_nt_seq=str(nt_seq[orf_id].seq),
                    clone_acc=orf_id,
                )
            )
        genes[gene_name] = isolib.Gene(
            gene_name,
            isoforms,
            ensembl_gene_id=hgnc_to_ensembl[gene_name],
            uniprot_ac=ensembl_to_uniprot.get(hgnc_to_ensembl[gene_name], None),
        )
        genes[gene_name].tf_family = tf_family[gene_name]
    pfam["gene_name"] = pfam["query name"].apply(lambda x: x.split("|")[0])
    for _i, row in pfam.iterrows():
        gene_name = row["gene_name"]
        iso_name = (
            row["query name"].split("|")[0]
            + "-"
            + row["query name"].split("|")[1].split("/")[0]
        )
        if gene_name not in genes or iso_name not in genes[gene_name]:
            continue
        genes[gene_name][iso_name].add_aa_seq_feature(
            category="Pfam_domain",
            name=row["target name"],
            accession=row["pfam_ac"],
            start=row["env_coord_from"] - 1,
            end=row["env_coord_to"],
            description=row["description of target"],
        )
    _make_c2h2_zf_arrays(genes)
    _add_dbd_flanks(genes)

    tfs_gencode = load_annotated_gencode_tfs(
        subset={g.ensembl_gene_id for g in genes.values()}
    )
    uncloned_orfs = defaultdict(list)
    for tf in genes.values():
        for gencode_isoform in tfs_gencode[tf.name].orfs:
            clone_match = False
            for cloned_isoform in tf.orfs:
                if gencode_isoform.exons == cloned_isoform.exons:
                    cloned_isoform.ensembl_transcript_ids = (
                        gencode_isoform.ensembl_transcript_ids
                    )
                    cloned_isoform.ensembl_protein_ids = (
                        gencode_isoform.ensembl_protein_ids
                    )
                    cloned_isoform.ensembl_transcript_names = (
                        gencode_isoform.ensembl_transcript_names
                    )
                    cloned_isoform.aa_seq_GENCODE = gencode_isoform.aa_seq
                    cloned_isoform.nt_seq_CDS_GENCODE = gencode_isoform.nt_seq
                    clone_match = True
            if not clone_match:
                uncloned_orfs[tf.name].append(gencode_isoform)
    for gene_name, isoforms in uncloned_orfs.items():
        if gene_name == "HSFY1":
            continue  # different strand isoforms???
        genes[gene_name].add_isoforms(isoforms)

    mane = pd.read_csv(path_MANE_select, sep="\t")
    mane_select = set(
        mane.loc[mane["MANE_status"] == "MANE Select", "Ensembl_nuc"]
        .str.slice(0, 15)
        .values
    )
    for tf in genes.values():
        if tf.ensembl_gene_id not in mane["Ensembl_Gene"].str.slice(0, 15).values:
            continue  # not all genes have a MANE select isoform
        for iso in tf.orfs:
            if hasattr(iso, "clone_acc") and iso.is_novel_isoform():
                iso.is_MANE_select_transcript = False
            else:
                iso.is_MANE_select_transcript = any(
                    t in mane_select for t in iso.ensembl_transcript_ids
                )

    appris = pd.read_csv(path_APPRIS, sep="\t")
    if appris["Transcript stable ID"].duplicated().any():
        raise UserWarning(
            "Unexpected duplicate ensembl transcript IDs in {}".format(path_APPRIS)
        )
    appris = appris.set_index("Transcript stable ID")["APPRIS annotation"].to_dict()

    def _consolidate_appris_annotations(annotations):
        return sorted(
            list(annotations), key=lambda x: int(x[-1]) - 99 * x.startswith("principle")
        )[0]

    for tf in genes.values():
        for iso in tf.orfs:
            if hasattr(iso, "clone_acc") and iso.is_novel_isoform():
                continue
            annotations = {
                appris[tid] for tid in iso.ensembl_transcript_ids if tid in appris
            }
            if len(annotations) > 0:
                iso.APPRIS_annotation = _consolidate_appris_annotations(annotations)

    reg_dom = pd.read_excel(path_effector_domains, sheet_name="Table S2")
    domain_type_full = {
        "AD": "Activation domain",
        "RD": "Repression domain",
        "Bif": "Bi-functional domain",
    }
    all_ensembl_gene_ids = {tf.ensembl_gene_id for tf in genes.values()}
    if not (
        reg_dom["TF name"].isin(genes.keys())
        == reg_dom["ENSEMBL gene ID"].isin(all_ensembl_gene_ids)
    ).all():
        print(
            reg_dom.loc[
                reg_dom["TF name"].isin(genes.keys())
                != reg_dom["ENSEMBL gene ID"].isin(all_ensembl_gene_ids),
                ["TF name", "ENSEMBL gene ID"],
            ]
        )
        raise UserWarning(
            "Problem with inconsistent gene names between effector domain file and cloned TFs"
        )
    for tf in genes.values():
        for _i, row in reg_dom.loc[reg_dom["TF name"] == tf.name, :].iterrows():
            desc = domain_type_full[row["Domain type"]]
            desc += "\nassay: " + row["Assay"]
            desc += "\nPMID: {}".format(row["Reference (PMID)"])
            if pd.notnull(row["Notes"]):
                desc += "Notes: " + row["Notes"]
            for iso in tf.orfs:
                if row["Sequence"] not in iso.aa_seq:
                    continue
                if len(re.findall("(?={})".format(row["Sequence"]), iso.aa_seq)) != 1:
                    raise UserWarning(
                        "Problem mapping effector domain: {} {}".format(row, iso)
                    )
                iso.add_aa_seq_feature(
                    category="effector_domain",
                    name=row["Domain type"],
                    accession=row["Effector domain ID"],
                    start=iso.aa_seq.find(row["Sequence"]),
                    end=iso.aa_seq.find(row["Sequence"]) + len(row["Sequence"]),
                    description=desc,
                )

    disorder = pd.read_csv(path_disorder, sep="\t")
    clones_with_disorder_data = set(disorder["clone_name"].unique())
    for tf in genes.values():
        for iso in tf.cloned_isoforms:
            if iso.name not in clones_with_disorder_data:
                print("missing disorder data for {}".format(iso.name))
                continue
            iso.disorder = disorder.loc[
                disorder["clone_name"] == iso.name, "is_disordered"
            ].values

    # sanity check
    for tf in genes.values():
        for iso in tf.cloned_isoforms:
            if hasattr(iso, "disorder"):
                if len(iso.disorder) != len(iso.aa_seq):
                    raise UserWarning(
                        "inconsistent amino acid sequence and disordered residues data for {}".format(
                            iso.name
                        )
                    )

    return genes


def _filter_gencode_gtf(out_file_path, genes_subset):
    print(
        "Filtering GENCODE .gtf file. Takes a few minutes but only needs to be done once."
    )
    import pyranges

    path_gencode_gtf = DATA_DIR / "external/gencode.v30.annotation.gtf"
    algn = pyranges.read_gtf(path_gencode_gtf, duplicate_attr=True)
    algn = algn[
        (algn.Feature == "CDS")
        & (algn.gene_type == "protein_coding")
        & (algn.transcript_type == "protein_coding")
    ]
    algn = algn[algn.tag.str.contains("basic")]
    algn = algn[algn.gene_id.str.replace(r"\..*", "", regex=True).isin(genes_subset)]
    algn.to_gtf(out_file_path)


def _filter_gencode_fasta(out_file_path, genes_subset):
    path_gencode_fa = DATA_DIR / "external/gencode.v30.pc_transcripts.fa"
    out_records = filter(
        lambda x: x.name.split("|")[1].split(".")[0] in genes_subset,
        SeqIO.parse(path_gencode_fa, format="fasta"),
    )
    SeqIO.write(out_records, out_file_path, "fasta")


def updated_gene_names():
    """We have two sources for HGNC symbols that were done at different times:
    GENCODE v30 and TF DB v1.01. We map the gencode to TF DB names.

    Note that this version of gencode can have the same gene name for different
    genes, e.g. ATF7 for ENSG00000170653 and ENSG00000267281
    """
    cached_file = Path(CACHE_DIR / "renamed_genes.tsv")
    if not cached_file.exists():
        import pyranges

        path_filtered_gencode_gtf = (
            CACHE_DIR / "filtered_CDS_PC_basic_TFs.gencode.v30.annotation.gtf"
        )
        algn = pyranges.read_gtf(path_filtered_gencode_gtf, duplicate_attr=True).df
        algn["gene_id"] = algn["gene_id"].str.replace(r"\..*", "", regex=True)
        algn["transcript_id"] = algn["transcript_id"].str.replace(
            r"\..*", "", regex=True
        )
        tf_db = pd.read_csv(DATA_DIR / "external/Human_TF_DB_v_1.01.csv")
        df = (
            pd.merge(algn, tf_db, how="inner", left_on="gene_id", right_on="Ensembl ID")
            .loc[:, ["Ensembl ID", "gene_name", "HGNC symbol"]]
            .rename(
                columns={
                    "gene_name": "HGNC symbol GENCODE v30",
                    "HGNC symbol": "HGNC symbol TF DB v1.01",
                }
            )
            .drop_duplicates()
        )
        df = df.loc[df["HGNC symbol GENCODE v30"] != df["HGNC symbol TF DB v1.01"], :]
        df.to_csv(cached_file, index=False, sep="\t")
    else:
        df = pd.read_csv(cached_file, sep="\t")
    return df.set_index("Ensembl ID")["HGNC symbol TF DB v1.01"].to_dict()


@cache_with_pickle
def load_annotated_gencode_tfs(
    subset=None,
    path_filtered_gencode_gtf=CACHE_DIR
    / "filtered_CDS_PC_basic_TFs.gencode.v30.annotation.gtf",
    path_filtered_gencode_fa=CACHE_DIR
    / "filtered_only_TFs.gencode.v30.pc_transcripts.fa",
    path_gencode_aa_seq=DATA_DIR / "external/gencode.v30.pc_translations.fa",
    path_MANE_select=DATA_DIR / "external/MANE.GRCh38.v0.95.summary.txt",
    path_APPRIS=DATA_DIR
    / "external/APPRIS-annotations_human_GRCh38.p13_ensembl104.tsv",
    path_effector_domains=DATA_DIR
    / "external/Soto-et-al_MolCell_2022_Supplementary-tables.xlsx",
):

    import pyranges  # this import is hidden as it's causing installation problems

    tf_ensembl_gene_ids = set(load_human_tf_db()["Ensembl ID"].values)
    tf_family = load_tf_families()
    # TODO: check for consistency of cached file with list of TF ensembl IDs
    if not path_filtered_gencode_gtf.exists():
        _filter_gencode_gtf(path_filtered_gencode_gtf, tf_ensembl_gene_ids)
    if not path_filtered_gencode_fa.exists():
        _filter_gencode_fasta(path_filtered_gencode_fa, tf_ensembl_gene_ids)
    # note that pyranges switches the indexing to python 0-indexed, half-open interval
    algn = pyranges.read_gtf(path_filtered_gencode_gtf, duplicate_attr=True).df
    algn["gene_id"] = algn["gene_id"].str.replace(r"\..*", "", regex=True)
    algn["transcript_id"] = algn["transcript_id"].str.replace(r"\..*", "", regex=True)
    # remove Y-chromosome copies of PAR region, leaving the X copies
    algn = algn.loc[~(algn["tag"].str.contains("PAR")), :]
    nt_seq = {
        r.name.split("|")[0].split(".")[0]: r
        for r in SeqIO.parse(path_filtered_gencode_fa, format="fasta")
    }
    if not algn["transcript_id"].isin(nt_seq).all():
        missing = algn.loc[~algn["transcript_id"].isin(nt_seq), "transcript_id"].values
        raise ValueError(", ".join(missing) + " not in " + path_filtered_gencode_fa)
    tf_db = load_human_tf_db()
    if tf_db["HGNC symbol"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    genes_to_rename = updated_gene_names()
    ensembl_to_uniprot = pd.read_csv(
        DATA_DIR / "external/ensembl_to_uniprot_ids_human_tfs.tsv", sep="\t"
    )
    ensembl_to_uniprot = (
        ensembl_to_uniprot.groupby("ensembl_gene_id")["uniprot_ac"]
        .apply(lambda x: "/".join(x))
        .to_dict()
    )
    tf_gene_names = set(algn["gene_name"].unique())
    valid_transcipts = set(algn["transcript_name"].unique())
    transcript_name_to_ensembl_gene_id = (
        algn.loc[:, ["transcript_name", "gene_id"]]
        .drop_duplicates()
        .set_index("transcript_name")["gene_id"]
    )
    aa_seqs = defaultdict(dict)
    duplicates = {}
    transcript_name_to_id = {}
    transcript_name_to_protein_id = {}
    for record in SeqIO.parse(path_gencode_aa_seq, "fasta"):
        protein_id, transcript_id = [
            x.split(".")[0] for x in record.name.split("|")[:2]
        ]
        transcript_name, gene_name = record.name.split("|")[5:7]
        transcript_name_to_id[transcript_name] = transcript_id
        transcript_name_to_protein_id[transcript_name] = protein_id
        if gene_name not in tf_gene_names or transcript_name not in valid_transcipts:
            continue
        aa_seqs[gene_name][transcript_name] = record.seq
    for gene_transcripts in aa_seqs.values():
        for transcript_name, seq in gene_transcripts.items():
            duplicates[transcript_name] = list(
                sorted([k for k, v in gene_transcripts.items() if v == seq])
            )
    unique_pc_transcripts = set([v[0] for v in duplicates.values()])

    genes = {}
    for gene_id in tqdm.tqdm(algn["gene_id"].unique()):
        if subset is not None and gene_id not in subset:
            continue
        if gene_id == "ENSG00000163602":
            continue  # RYBP has a sequencing error so don't have full CDS coordinates
        transcript_ids = algn.loc[algn["gene_id"] == gene_id, "transcript_id"].unique()
        isoforms = []
        for transcript_id in transcript_ids:
            transcript_name, gene_name = nt_seq[transcript_id].name.split("|")[4:6]
            if transcript_name not in unique_pc_transcripts:
                continue
            cds = _extract_region(nt_seq[transcript_id], "CDS", raise_error=True)
            utr5 = _extract_region(nt_seq[transcript_id], "UTR5", raise_error=False)
            utr3 = _extract_region(nt_seq[transcript_id], "UTR3", raise_error=False)
            exons = []
            columns = [
                "gene_id",
                "transcript_id",
                "Chromosome",
                "Strand",
                "Start",
                "End",
            ]
            for _i, row in algn.loc[
                algn["transcript_id"] == transcript_id, columns
            ].iterrows():
                exons.append(isolib.Exon(*row.values))
            isoforms.append(
                isolib.Isoform(
                    transcript_name,
                    exons,
                    CDS_nt_seq=cds,
                    aa_seq=str(aa_seqs[gene_name][transcript_name]),
                    UTR_5prime_nt_seq=utr5,
                    UTR_3prime_nt_seq=utr3,
                    ensembl_transcript_ids=[
                        transcript_name_to_id[x] for x in duplicates[transcript_name]
                    ],
                    ensembl_transcript_names=duplicates[transcript_name],
                    ensembl_protein_ids=[
                        transcript_name_to_protein_id[x]
                        for x in duplicates[transcript_name]
                    ],
                )
            )
        if gene_id in genes_to_rename:
            gene_name = genes_to_rename[gene_id]

        genes[gene_name] = isolib.Gene(
            gene_name,
            isoforms,
            ensembl_gene_id=gene_id,
            uniprot_ac=ensembl_to_uniprot.get(gene_id, None),
        )
        genes[gene_name].tf_family = tf_family[gene_name]

    pfam = load_pfam_domains_gencode()
    pfam["gene_name"] = pfam["query name"].apply(
        lambda x: "-".join(x.split("|")[0].split("-")[:-1])
    )
    pfam["ensembl_gene_id"] = (
        pfam["query name"]
        .apply(lambda x: x.split("|")[0])
        .map(transcript_name_to_ensembl_gene_id)
    )
    for _i, row in pfam.iterrows():
        gene_name = row["gene_name"]
        ensembl_gene_id = row["ensembl_gene_id"]
        if ensembl_gene_id in genes_to_rename:
            gene_name = genes_to_rename[ensembl_gene_id]
        transcript_name = row["query name"].split("|")[0]
        if gene_name not in genes or transcript_name not in genes[gene_name]:
            continue
        if transcript_name not in unique_pc_transcripts:
            continue
        genes[gene_name][transcript_name].add_aa_seq_feature(
            category="Pfam_domain",
            name=row["target name"],
            accession=row["pfam_ac"],
            start=row["env_coord_from"] - 1,
            end=row["env_coord_to"],
        )
    _make_c2h2_zf_arrays(genes)
    _add_dbd_flanks(genes)

    mane = pd.read_csv(path_MANE_select, sep="\t")
    mane_select = set(
        mane.loc[mane["MANE_status"] == "MANE Select", "Ensembl_nuc"]
        .str.slice(0, 15)
        .values
    )
    for tf in genes.values():
        if tf.ensembl_gene_id not in mane["Ensembl_Gene"].str.slice(0, 15).values:
            continue  # not all genes have a MANE select isoform
        for iso in tf.orfs:
            if hasattr(iso, "clone_acc") and iso.is_novel_isoform():
                continue
            iso.is_MANE_select_transcript = any(
                t in mane_select for t in iso.ensembl_transcript_ids
            )

    appris = pd.read_csv(path_APPRIS, sep="\t")
    if appris["Transcript stable ID"].duplicated().any():
        raise UserWarning(
            "Unexpected duplicate ensembl transcript IDs in {}".format(path_APPRIS)
        )
    appris = appris.set_index("Transcript stable ID")["APPRIS annotation"].to_dict()

    def _consolidate_appris_annotations(annotations):
        return sorted(
            list(annotations), key=lambda x: int(x[-1]) - 99 * x.startswith("principal")
        )[0]

    for tf in genes.values():
        for iso in tf.orfs:
            if hasattr(iso, "clone_acc") and iso.is_novel_isoform():
                continue
            annotations = {
                appris[tid] for tid in iso.ensembl_transcript_ids if tid in appris
            }
            if len(annotations) > 0:
                iso.APPRIS_annotation = _consolidate_appris_annotations(annotations)

    reg_dom = pd.read_excel(path_effector_domains, sheet_name="Table S2")
    reg_dom["TF name"] = (
        reg_dom["ENSEMBL gene ID"].map(genes_to_rename).fillna(reg_dom["TF name"])
    )
    domain_type_full = {
        "AD": "Activation domain",
        "RD": "Repression domain",
        "Bif": "Bi-functional domain",
    }
    all_ensembl_gene_ids = {tf.ensembl_gene_id for tf in genes.values()}
    if not (
        reg_dom["TF name"].isin(genes.keys())
        == reg_dom["ENSEMBL gene ID"].isin(all_ensembl_gene_ids)
    ).all():
        print(
            reg_dom.loc[
                reg_dom["TF name"].isin(genes.keys())
                != reg_dom["ENSEMBL gene ID"].isin(all_ensembl_gene_ids),
                ["TF name", "ENSEMBL gene ID"],
            ]
        )
        raise UserWarning(
            "Problem with inconsistent gene names between effector domain file and cloned TFs"
        )
    for tf in genes.values():
        for _i, row in reg_dom.loc[reg_dom["TF name"] == tf.name, :].iterrows():
            desc = domain_type_full[row["Domain type"]]
            desc += "\nassay: " + row["Assay"]
            desc += "\nPMID: {}".format(row["Reference (PMID)"])
            if pd.notnull(row["Notes"]):
                desc += "Notes: " + row["Notes"]
            for iso in tf.orfs:
                if row["Sequence"] not in iso.aa_seq:
                    continue
                if len(re.findall("(?={})".format(row["Sequence"]), iso.aa_seq)) != 1:
                    raise UserWarning(
                        "Problem mapping effector domain: {} {}".format(row, iso)
                    )
                iso.add_aa_seq_feature(
                    category="effector_domain",
                    name=row["Domain type"],
                    accession=row["Effector domain ID"],
                    start=iso.aa_seq.find(row["Sequence"]),
                    end=iso.aa_seq.find(row["Sequence"]) + len(row["Sequence"]),
                    description=desc,
                )

    return genes


def _extract_region(rec, region, raise_error=True):
    """
    Get subset of nucleotide sequence that correpsonds to different
    regions (e.g. CDS, 3'UTR) when bounds of regions are specified
    in desciption in fasta file, like this:
    >HES4-204|HES4|667|UTR5:1-8|CDS:9-578|UTR3:579-667|
    """
    bounds = [x for x in rec.name.split("|") if x.startswith(region + ":")]
    if len(bounds) == 0:
        if raise_error:
            raise UserWarning("Missing {} coordinates\n".format(region) + rec.name)
        else:
            return None
    if len(bounds) > 1:
        raise UserWarning(
            "More than one set of {} coordinates\n".format(region) + rec.name
        )
    start, stop = bounds[0][len(region) + 1 :].split("-")
    return str(rec.seq)[int(start) - 1 : int(stop)]


def _make_c2h2_zf_arrays(tfs, MAX_NUM_AA_C2H2_ZF_SEPERATION=10):
    clans = load_pfam_clans()
    C2H2_ZF_PFAM_CLAN_AC = "CL0361"
    c2h2_zf_pfam_ids = {k for k, v in clans.items() if v == C2H2_ZF_PFAM_CLAN_AC}
    for gene in tfs.values():
        for orf in gene.orfs:
            zfs = [d for d in orf.aa_seq_features if d.accession in c2h2_zf_pfam_ids]
            if len(zfs) < 2:
                continue
            zfs = sorted(zfs, key=lambda x: x.start)
            array = [zfs[0]]
            for zf in zfs[1:]:
                if zf.start <= (array[-1].end - 1) + MAX_NUM_AA_C2H2_ZF_SEPERATION:
                    array.append(zf)
                else:
                    if len(array) > 1:
                        orf.add_aa_seq_feature(
                            "ZF_array",
                            "C2H2_ZF_array_" + str(len(array)),
                            "C2H2_ZF_array_" + str(len(array)),
                            array[0].start,
                            array[-1].end,
                            description="{} C2H2 zinc fingers".format(len(array)),
                        )
                        for dom in array:
                            orf.remove_aa_seq_feature(dom.accession, dom.start, dom.end)
                    array = [zf]
            if len(array) > 1:
                orf.add_aa_seq_feature(
                    "ZF_array",
                    "C2H2_ZF_array_" + str(len(array)),
                    "C2H2_ZF_array_" + str(len(array)),
                    array[0].start,
                    array[-1].end,
                    description="{} C2H2 zinc fingers".format(len(array)),
                )
                for dom in array:
                    orf.remove_aa_seq_feature(dom.accession, dom.start, dom.end)


def _add_dbd_flanks(genes):
    for gene in genes.values():
        for isoform in gene.orfs:
            for dbd in isoform.dna_binding_domains:
                start_n = max(0, dbd.start - 15)
                end_n = dbd.start
                start_c = dbd.end
                end_c = min(len(isoform.aa_seq), dbd.end + 15)
                if dbd.start > 0:
                    isoform.add_aa_seq_feature(
                        "DBD_flank",
                        dbd.accession + "_flank_N",
                        "N_DBD_flank",
                        start_n,
                        end_n,
                    )
                if dbd.end < len(isoform.aa_seq):
                    isoform.add_aa_seq_feature(
                        "DBD_flank",
                        dbd.accession + "_flank_C",
                        "C_DBD_flank",
                        start_c,
                        end_c,
                    )
