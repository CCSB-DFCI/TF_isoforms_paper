from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from ucsc.api import Sequence

from .utils import DATA_DIR
from .protein_data import load_human_tf_db, load_cofactors, load_signaling_genes


def load_valid_isoform_clones():
    """

    TODO: this is including PCGF6 which has an insertion relative
    to the reference geneome, so that it's not consistent with the
    load_annotated_TFiso1 set. Check the PCGF6 sequences and probably
    remove from here.

    """
    df = pd.read_csv(
        DATA_DIR / "processed/valid-unique-isoform-clones_2021-07-20.tsv", sep="\t"
    )
    df["clone_name"] = df["clone_acc"].map(
        lambda x: x.split("|")[0] + "-" + x.split("|")[1].split("/")[0]
    )
    df = df.rename(columns={"gene": "gene_symbol"})
    df = df.loc[df["gene_symbol"] != "PCGF6", :]
    return df


# TODO: this is a mess because there are circular function dependencies here
def define_valid_isoform_clones():
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
    df["in_y1h"] = df["clone_acc"].isin(y1h["clone_acc"])
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
        iso_annot.set_index("clone_acc")["gc_match"] == 0
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
    require_at_least_two_partners=False,
    include_additional_paralog_pairs=True,
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
        include_additional_paralog_pairs (bool, optional): Whether to include
        pairs that did not come from the screens but were also tested because
        a paralog (or non-paralog control) gene interacted with that partner.
        Defaults to True.

    Returns:
        pandas.DataFrame

    """
    y2h = load_full_y2h_data_including_controls()
    if include_additional_paralog_pairs:
        categories = {
            "tf_isoform_ppis",
            "tf_paralog_ppis",
            "paralog_with_PDI",
            "non_paralog_control",
        }
    else:
        categories = {"tf_isoform_ppis"}
    ppi = y2h.loc[
        y2h["category"].isin(categories),
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


def _load_y2h_paralogs_additional_data():
    """Pairs tested in Y2H for the paralogs data, in addition to isoform pairs."""
    y2h = load_full_y2h_data_including_controls()
    pairs = load_paralog_pairs_tested_in_y2h()
    y2h_paralog = y2h.loc[
        y2h["category"].isin(
            ["tf_paralog_ppis", "paralog_with_PDI", "non_paralog_control"]
        ),
        :,
    ].copy()
    pair_map = defaultdict(set)
    for _i, row in pairs.iterrows():
        a, b = row["gene_symbol_a"], row["gene_symbol_b"]
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


def load_full_y2h_data_including_controls(
    add_partner_cateogories=False,
    remove_keratin_associated_proteins=True,
):
    """

    add_partner_cateogories (bool): include columns categorizing the interaction
        partner, e.g. co-factor, signaling protein, etc.
    remove_keratin_associated_proteins (bool): keratin associated proteins and
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


def load_paralog_pairs_tested_in_y2h(filter_for_valid_clones=True):
    """Pairs of TF gene paralogs and non-paralogs that were tested in Y2H pairwise tests.

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
        .loc[:, ["tf1", "tf2", "is_paralog_pair"]]
        .rename(
            columns={
                "tf1": "gene_symbol_a",
                "tf2": "gene_symbol_b",
                "AAseq_identity%": "aa_pct_seq_identity",
            }
        )
    )
    df = df.drop_duplicates()
    if filter_for_valid_clones:
        valid_clones = load_valid_isoform_clones()
        df = df.loc[
            df["gene_symbol_a"].isin(valid_clones["gene_symbol"])
            & df["gene_symbol_b"].isin(valid_clones["gene_symbol"]),
            :,
        ]
    if (df["gene_symbol_a"] == df["gene_symbol_b"]).any():
        raise ValueError("Same gene twice, should be two different genes")
    return df


def load_ppi_partner_categories():
    # table with gene_symbol_partner, category and cofactor_type
    tfdb = load_human_tf_db()
    tf_genes = set(tfdb["HGNC symbol"].unique())
    cof = load_cofactors()
    cofactor_genes = set(cof["gene_symbol"].unique())
    signaling_genes = load_signaling_genes()

    df = pd.DataFrame(
        data=load_full_y2h_data_including_controls()["db_gene_symbol"].unique(),
        columns=["gene_symbol_partner"],
    )
    df["category"] = "other"
    df.loc[df["gene_symbol_partner"].isin(tf_genes), "category"] = "TF"
    df.loc[df["gene_symbol_partner"].isin(cofactor_genes), "category"] = "cofactor"
    df.loc[df["gene_symbol_partner"].isin(signaling_genes), "category"] = "signaling"
    df = pd.merge(
        df,
        cof.loc[:, ["gene_symbol", "cofactor_type"]],
        left_on="gene_symbol_partner",
        right_on="gene_symbol",
        how="left",
    ).drop(columns=["gene_symbol"])

    return df


def _load_old_ppi_partner_categories():
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
    df = df.rename(columns={"partner": "gene_symbol_partner"})
    return df


def load_n2h_ppi_validation_data():
    df = pd.read_csv(DATA_DIR / "internal/N2H_results.tsv", sep="\t")
    df = df.dropna(subset=["score_pair", "score_empty-N1", "score_empty-N2"])
    if df.duplicated().any():
        raise UserWarning("unexpected duplicates")
    df["NLR"] = df["score_pair"] / df[["score_empty-N1", "score_empty-N2"]].max(axis=1)
    df["log2 NLR"] = df["NLR"].apply(np.log2)
    df["score_pair_log10"] = df["score_pair"].apply(np.log10)
    df["score_empty-N1_log10"] = df["score_empty-N1"].apply(np.log10)
    df["score_empty-N2_log10"] = df["score_empty-N2"].apply(np.log10)
    return df


def load_y1h_pdi_data(add_missing_data=False, include_pY1H_data=True):
    """Yeast one-hybrid results for the TF isoforms.

    Args:
        add_missing_data (bool, optional): include rows containing False for
            isoforms that were tested but where there were no hits for any
            isoform of that gene. Defaults to False.
        include_pY1H_data (bool, optional): Add data from additional experiment
            performed from the single cases in the paired yeast one-hybrid test.
            Defaults to True.

    Returns:
        pandas.DataFrame: yeast one-hybrid data, one row for each isoform,
        one column for each bait
    """
    df_long = pd.read_excel(
        DATA_DIR / "internal/TF isoforms eY1H calls 14JUL23.xlsx",
        sheet_name="List format",
    )
    if df_long["Bait"].isnull().any():
        raise UserWarning("unexpected missing values")
    if df_long["isoform ID"].isnull().any():
        raise UserWarning("unexpected missing values")
    df_long["Binary"] = df_long["Binary"].map(
        {"yes": True, "no": False, "inconclusive": np.nan}
    )

    df = (
        pd.concat(
            [df_long["isoform ID"], pd.get_dummies(df_long["Bait"])],
            axis=1,
        )
        .groupby(["isoform ID"])
        .sum()
        > 0
    ).reset_index()
    for _i, row in df_long.loc[df_long["Binary"] == False, :].iterrows():
        df.loc[(df["isoform ID"] == row["isoform ID"]), row["Bait"]] = False
    for _i, row in df_long.loc[df_long["Binary"].isnull(), :].iterrows():
        df.loc[(df["isoform ID"] == row["isoform ID"]), row["Bait"]] = np.nan
    clones = load_valid_isoform_clones()
    if clones["clone_name"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    clone_acc = clones.set_index("clone_name")["clone_acc"]
    gene_symbol = clones.set_index("clone_name")["gene_symbol"]
    df["gene_symbol"] = df["isoform ID"].map(gene_symbol)
    df["clone_acc"] = df["isoform ID"].map(clone_acc)
    df = df.drop(columns=["isoform ID"])
    df = df.loc[:, list(df.columns[-2:]) + list(df.columns[:-2])]

    # HACK adding data that is missing from excel file
    df.loc[df["gene_symbol"] == "TCF4", "HS1597"] = True

    df[df.columns[2:]] = df[df.columns[2:]].astype("boolean")
    zeros = pd.read_csv(DATA_DIR / "internal/a2_juan_isoforms_wo_pdi.tsv", sep="\t")
    zeros = zeros.rename(columns={"tf": "gene_symbol", "unique_acc": "clone_acc"})
    zeros = zeros.loc[~zeros["clone_acc"].isin(df["clone_acc"].values), :]
    zeros = pd.concat(
        [zeros, pd.DataFrame(data=False, index=zeros.index, columns=df.columns[2:])],
        axis=1,
    )
    zeros[zeros.columns[2:]] = zeros[zeros.columns[2:]].astype("boolean")
    df = pd.concat([df, zeros], axis=0, sort=False).reset_index(drop=True)

    if add_missing_data:
        tested = pd.read_excel(DATA_DIR / "internal/TF iso ey1h seq confirmed.xlsx")
        tested["Coordinate"] = tested["Coordinate"].apply(
            lambda x: x.split("-")[0].zfill(2)
            + x.split("-")[1][0]
            + x.split("-")[1][1:].zfill(2)
        )
        clones = load_valid_isoform_clones()
        if tested["Coordinate"].duplicated().any():
            raise UserWarning("unexpected duplicates")
        clones["Coordinate"] = clones["clone_acc"].str.slice(-5)
        if clones["Coordinate"].duplicated().any():
            raise UserWarning("unexpected duplicates")
        tested = pd.merge(tested, clones, how="left", on=["Coordinate"])
        tested = tested.loc[
            (tested["sequence verified?"] == "yes")
            & tested["clone_acc"].notnull()
            & ~tested["clone_acc"].isin(df["clone_acc"].values),
            ["gene_symbol", "clone_acc"],
        ].copy()
        tested.loc[:, df.columns[2:]] = False
        df = pd.concat([df, tested])

    df[df.columns[2:]] = df[df.columns[2:]].astype("boolean")
    if include_pY1H_data:
        pY1H = load_additional_PDI_data_from_unpaired_cases_in_paired_Y1H_experiment()
        df = pd.merge(df, pY1H, on=["gene_symbol", "clone_acc"], how="outer")
    df = df.sort_values(["gene_symbol", "clone_acc"])
    if df["clone_acc"].isnull().any():
        raise UserWarning("unexpected missing values")
    # HACK GATA2|3/4|12A02 & GATA2|4/4|11A12 have duplicate rows with no hits
    df = df.drop_duplicates()
    if df["clone_acc"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    return df


def load_Y1H_DNA_bait_sequences():
    cache_path = DATA_DIR / "processed/Y1H_DNA_baits.fa"
    if cache_path.exists():
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(cache_path, format="fasta")}
    df = pd.read_excel("../data/external/Fuxman-Bass-et-al_Cell_2015_Table-S8.xlsx")
    if df["Vista Enhancer ID"].duplicated().any():
        raise UserWarning("unexpected duplicate enhancer IDs")
    if df["Sequence"].duplicated().any():
        raise UserWarning("unexpected duplicate sequences")
    baits = df.set_index("Vista Enhancer ID")["Sequence"].to_dict()

    df = pd.read_excel(
        "../data/external/Fuxman-Bass-et-al_Cell_2015_Table-S10.xlsx",
        sheet_name="Noncoding mutations",
    )
    # reference genome version: GRCh37.p12

    def get_dna_sequence_from_coords(s):
        chrom, pos = s.split(":")
        chrom = "chr" + chrom[3:]  # to avoid problems from capital C
        start, stop = map(int, pos.replace(",", "").split("-"))
        try:
            ret = Sequence.get(genome="hg19", chrom=chrom, start=start - 1, end=stop)
        except:
            print("failed for: " + s)
            raise
        return ret.dna.upper()

    df["seq"] = df["Region amplified"].apply(get_dna_sequence_from_coords)
    baits = {**baits, **df.set_index("Mutation ID")["seq"].to_dict()}
    baits = {k.upper(): v for k, v in baits.items()}
    with open(cache_path, "w") as f:
        SeqIO.write(
            (SeqRecord(id=k, description="", seq=Seq(v)) for k, v in baits.items()),
            f,
            "fasta",
        )
    return baits


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
    df["gene_symbol"] = df["isoform"].apply(lambda x: x.split("-")[0])
    df["clone_acc"] = df["isoform"].map(clones)

    # kaia had to change this line - py version issue?
    # df = df.pivot(index=["gene_symbol", "clone_acc"], columns="bait gene", values="interaction")
    df = df.set_index(["gene_symbol", "clone_acc"]).pivot(columns="bait gene")[
        "interaction"
    ]

    df = df.astype("boolean")
    df = df.reset_index()
    df = df.loc[df["clone_acc"].notnull(), :]
    return df


def load_PDI_luciferase_validation_experiment():
    df = pd.read_excel(
        DATA_DIR / "internal/2019-08-28 Luciferase results.xlsx",
        sheet_name="Use for analyses",
    )
    df = df.rename(columns={"TF_Name": "gene_symbol", "unique_acc": "clone_acc"})
    df = df.iloc[:, 0:15]
    df["Y1H_positive"] = df["Interaction?"].map({"yes": True, "no": False})
    if df.isnull().any().any():
        raise UserWarning("Unexpected missing values")
    if df.duplicated(subset=["Bait", "clone_acc"]).any():
        raise UserWarning("unexpected duplicates")
    return df


def load_m1h_activation_data(add_missing_data=False):
    """ """
    df = pd.read_csv(DATA_DIR / "internal/a_m1h_final_table.tsv", sep="\t")
    df = df.rename(columns={"gene": "gene_symbol", "pos_acc": "clone_acc"})
    for column in df.columns:
        if column.startswith("M1H_rep"):
            df[column] = np.log2(df[column])
    if add_missing_data:
        isoforms = load_valid_isoform_clones()
        df = pd.merge(df, isoforms, on=["gene_symbol", "clone_acc"], how="outer")
    df = df.sort_values(["gene_symbol", "clone_acc"])
    return df
