from pathlib import Path

import pandas as pd

from ccsblib import paros_connection
from ccsblib import huri

from data_loading import load_full_y2h_data_including_controls


def isoform_clones():
    """The subset of TF isoform clones that have passed the stringent
       annotation process (at-length AA match to GENCODE, or approved by
       Gloria/GENCODE team).

    Returns:
        set(str): clone accession IDs

    """
    df = pd.read_sql(
        """SELECT symbol as gene,
                               unique_acc as clone_acc,
                               dup_idx
                          FROM tf_screen.iso6k_annotation
                         ORDER BY gene, clone_acc;""",
        paros_connection(),
    )
    return df


def tf_isoform_y2h_screen_results():
    """There were two screens performed:

    The cloned TF isoforms as AD-fusions against DB-fusions of:
    (1) ORFeome 9.1
    (2) Subset of TFs and co-factors

    Returns:
        pandas.DataFrame: for each pair, was it found in the first and second screens

    """
    qry = """SELECT ad_orf_id,
                    db_orf_id,
                    pool_name
               FROM swimseq_run.INGS_IST a,
                    swimseq_run.NGS_POOL b
              WHERE a.pool_id = b.pool_id
                AND a.pool_id in (787, 788)
                AND ist_score>=0.2;"""
    df = pd.read_sql(qry, paros_connection())
    df = (
        pd.get_dummies(df, columns=["pool_name"])
        .groupby(["ad_orf_id", "db_orf_id"])
        .sum()
        > 0
    ).reset_index()
    df = df.rename(
        columns={
            "pool_name_ds20180213_TFisoS04": "in_orfeome_screen",
            "pool_name_ds20180213_TFisoS05": "in_focussed_screen",
        }
    )
    if not (df["in_orfeome_screen"] | df["in_focussed_screen"]).all():
        raise UserWarning("Something went wrong...")
    return df


def _annotate_ppi_source_info(df):
    def load_orf_id_to_tf_gene():
        qry = """SELECT orf_id,
                        symbol AS gene_symbol
                FROM tf_screen.iso6k_sequences;"""
        df = pd.read_sql(qry, paros_connection())
        df = df.drop_duplicates()
        if df["orf_id"].duplicated().any():
            raise UserWarning("Unexpected duplicate ORF IDs")
        return df.set_index("orf_id")["gene_symbol"]

    screen = tf_isoform_y2h_screen_results()
    orf_id_to_tf_gene = load_orf_id_to_tf_gene()
    screen["ad_gene_symbol"] = screen["ad_orf_id"].map(orf_id_to_tf_gene)
    screen = screen.dropna(subset=["ad_gene_symbol"]).drop(columns=["ad_orf_id"])
    screen = screen.groupby(["ad_gene_symbol", "db_orf_id"]).any().reset_index()
    df = pd.merge(df, screen, how="left", on=["ad_gene_symbol", "db_orf_id"])
    df["in_orfeome_screen"] = df["in_orfeome_screen"].fillna(False)
    df["in_focussed_screen"] = df["in_focussed_screen"].fillna(False)
    hiu = huri.load_nw_hi_union(id_type="orf_id")
    id_map = huri.load_id_map("orf_id", "hgnc_symbol", via="ensembl_gene_id")
    hiu_pairs = set(
        (
            pd.merge(hiu, id_map, how="inner", left_on="orf_id_a", right_on="orf_id")[
                ["hgnc_symbol", "orf_id_b"]
            ].drop_duplicates()
        ).apply(lambda x: x["hgnc_symbol"] + "-" + str(x["orf_id_b"]), axis=1)
    )
    hiu_pairs = hiu_pairs.union(
        (
            pd.merge(hiu, id_map, how="inner", left_on="orf_id_b", right_on="orf_id")[
                ["hgnc_symbol", "orf_id_a"]
            ].drop_duplicates()
        ).apply(lambda x: x["hgnc_symbol"] + "-" + str(x["orf_id_a"]), axis=1)
    )
    df["in_hi_union"] = (df["ad_gene_symbol"] + "-" + df["db_orf_id"].astype(str)).isin(
        hiu_pairs
    )
    litbm = huri.load_nw_lit_bm_17(id_type="orf_id")
    litbm_pairs = set(
        (
            pd.merge(litbm, id_map, how="inner", left_on="orf_id_a", right_on="orf_id")[
                ["hgnc_symbol", "orf_id_b"]
            ].drop_duplicates()
        ).apply(lambda x: x["hgnc_symbol"] + "-" + str(x["orf_id_b"]), axis=1)
    )
    litbm_pairs = litbm_pairs.union(
        (
            pd.merge(litbm, id_map, how="inner", left_on="orf_id_b", right_on="orf_id")[
                ["hgnc_symbol", "orf_id_a"]
            ].drop_duplicates()
        ).apply(lambda x: x["hgnc_symbol"] + "-" + str(x["orf_id_a"]), axis=1)
    )
    df["in_lit_bm"] = (df["ad_gene_symbol"] + "-" + df["db_orf_id"].astype(str)).isin(
        litbm_pairs
    )
    return df


def isoform_and_paralog_y2h_data():
    """
    - NS: sequencing failed
    - NC: no call (e.g., mis-spotting)
    - AA: autoactivator

    """
    qry_a = """SELECT a.category,
                      a.ad_orf_id,
                      b.unique_acc AS ad_clone_acc,
                      a.ad_symbol AS ad_gene_symbol,
                      a.db_orf_id,
                      c.symbol AS db_gene_symbol,
                      a.final_score AS score,
                      d.standard_batch,
                      d.retest_pla,
                      d.retest_pos
                 FROM tf_screen.tf_isoform_final AS a
                 LEFT JOIN tf_screen.iso6k_sequences AS b
                   ON a.ad_orf_id = b.orf_id
                 LEFT JOIN horfeome_annotation_gencode27.orf_class_map_ensg AS c
                   ON a.db_orf_id = c.orf_id
                 LEFT JOIN tf_screen.retest AS d
                   ON a.retest_id = d.retest_id;"""
    df_a = pd.read_sql(qry_a, paros_connection())

    # remove duplicate ORF for gene DDX39B, where sequencing mostly failed
    df_a = df_a.loc[df_a["db_orf_id"] != 3677, :]

    df_a["category"] = df_a["category"].map(
        {
            "ppi": "tf_isoform_ppis",
            "ng_stem_cell_factor": "tf_isoform_ppis",
            "rrs": "rrs_isoforms",
            "litbm": "lit_bm_isoforms",
        }
    )
    qry_b = """SELECT a.simple_category AS category,
                      a.ad_orf_id,
                      b.unique_acc AS ad_clone_acc,
                      a.ad_symbol AS ad_gene_symbol,
                      a.db_orf_id,
                      c.symbol AS db_gene_symbol,
                      a.final_score AS score,
                      d.standard_batch,
                      d.retest_pla,
                      d.retest_pos
                 FROM tf_screen.paralog_final AS a
                 LEFT JOIN tf_screen.iso6k_sequences AS b
                   ON a.ad_orf_id = b.orf_id
                 LEFT JOIN horfeome_annotation_gencode27.orf_class_map_ensg AS c
                   ON a.db_orf_id = c.orf_id
                 LEFT JOIN tf_screen.retest AS d
                   ON a.retest_id = d.retest_id;"""
    df_b = pd.read_sql(qry_b, paros_connection())
    df_b["category"] = df_b["category"].map(
        {
            "paralog": "tf_paralog_ppis",
            "PDI_PPI": "paralog_with_PDI",
            "nonparalog": "non_paralog_control",
            "rrs": "rrs_paralogs",
            "litbm": "lit_bm_paralogs",
        }
    )
    df = pd.concat([df_a, df_b])
    # drop cases where single orf ID mapped to multiple gene symbols
    df = df.drop_duplicates(["category", "ad_orf_id", "db_orf_id"])
    # drop interaction partner ORFs whose sequence does not map to an ensembl gene
    df = df.dropna(subset=["db_gene_symbol"])

    non_protein_coding_genes = """RBMY2FP
                                  REXO1L6P
                                  SUMO1P1
                                  SLCO4A1-AS1
                                  ZNF213-AS1
                                  AC023509.3
                                  KRT8P41
                                  LINC01658""".split()
    df = df.loc[~df["db_gene_symbol"].isin(non_protein_coding_genes), :]

    # remove unexpected duplicates, some of which have inconsistent scores
    # see Issue #28
    df = df.loc[
        ~(
            df.duplicated(subset=["ad_clone_acc", "db_gene_symbol"], keep=False)
            & (df["category"] == "paralog_with_PDI")
        ),
        :,
    ]

    df = _annotate_ppi_source_info(df)
    return df


def n2h_ppi_validation_data():
    qry = """select a.test_orf_ida, a.test_orf_idb, a.test_pla, a.test_pos, b.score
    from tf_validation.validation_cp AS a
    left join tf_validation.mn2h_scoring AS b
    on a.standard_batch = b.standard_batch
    and a.test_pla_full = b.plate
    and a.test_pos = b.well
    where a.standard_batch = 'TFv02';"""
    df = pd.read_sql(qry, paros_connection())
    qry = """select orf_id1, orf_id2, source
    from tf_validation.validation_source
    where standard_batch = 'TFv02';"""
    source = pd.read_sql(qry, paros_connection())
    source["pair"] = (
        source[["orf_id1", "orf_id2"]].min(axis=1).astype(str)
        + "_"
        + source[["orf_id1", "orf_id2"]].max(axis=1).astype(str)
    )
    df["pair"] = (
        df[["test_orf_ida", "test_orf_idb"]].min(axis=1).astype(str)
        + "_"
        + df[["test_orf_ida", "test_orf_idb"]].max(axis=1).astype(str)
    )
    df = pd.merge(df, source.loc[:, ["pair", "source"]], how="left", on="pair")
    df.loc[
        (df["test_orf_ida"] == 0) & (df["test_orf_idb"] != 0), "source"
    ] = "empty-N1_control"
    df.loc[
        (df["test_orf_idb"] == 0) & (df["test_orf_ida"] != 0), "source"
    ] = "empty-N2_control"
    df.loc[df["test_pos"].isin(["B02", "C02"]), "source"] = "SKP1-SKP2_control"
    df.loc[df["test_pos"].isin(["D02", "E02"]), "source"] = "SKP1-BTRC_control"
    if df["source"].isnull().any():
        raise UserWarning("something wrong")
    df = df.loc[~df["source"].isin(["SKP1-SKP2_control", "SKP1-BTRC_control"]), :]

    # remove low liquid pairs
    missing = pd.read_excel(
        "../one_off/TFv02_empty_wells_after_cp.xlsx",
        sheet_name="Empty destination wells",
    )
    missing["Destination_plate"] = missing["Destination_plate"].str.slice(0, -3)
    missing["Destination_well"] = missing["Destination_well"].apply(
        lambda x: x[0] + str(x[1:]).zfill(2)
    )
    df.loc[
        (df["test_pla"].astype(str).str.zfill(3) + df["test_pos"]).isin(
            (
                missing["Destination_plate"].str.slice(-3) + missing["Destination_well"]
            ).values
        ),
        "score",
    ] = np.nan

    # remove low DNA count pairs
    qry = """
    select test_orf_id, vector , concentration_working_copy 
    from tf_validation.validation_nodes_final 
    where concentration_working_copy < 8 
    and node_pla_dna_working_copy  not like '%Low%'  ;
    """
    low_dna = pd.read_sql(qry, paros_connection())
    df.loc[
        df["test_orf_ida"].isin(
            low_dna.loc[low_dna["vector"] == "N1", "test_orf_id"].values
        ),
        "score",
    ] = np.nan
    df.loc[
        df["test_orf_idb"].isin(
            low_dna.loc[low_dna["vector"] == "N2", "test_orf_id"].values
        ),
        "score",
    ] = np.nan

    # match empty N1/N2 controls
    tmp_df = pd.merge(
        df.loc[~df["source"].isin(["empty-N1_control", "empty-N2_control"]), :],
        df.loc[
            df["source"] == "empty-N1_control",
            ["test_orf_idb", "test_pla", "test_pos", "score"],
        ],
        on=["test_orf_idb", "test_pla"],
        how="left",
        suffixes=("_pair", "_empty-N1"),
    )
    df = pd.merge(
        tmp_df,
        df.loc[
            df["source"] == "empty-N2_control",
            ["test_orf_ida", "test_pla", "test_pos", "score"],
        ],
        on=["test_orf_ida", "test_pla"],
        how="left",
    ).rename(columns={"test_pos": "test_pos_empty-N2", "score": "score_empty-N2"})

    # 23 -> plate 1
    # 24 -> plate 2
    # 25 -> plate 1, empty-N2
    # 26 -> plate 1, empty-N1
    # 27 -> plate 2, empty-N2
    # 28 -> plate 2, empty-N1
    qry = """select plate, well, score
    from tf_validation.mn2h_scoring AS b
    where standard_batch = 'TFv02'
    and plate between 'TFv02N2H_023' and 'TFv02N2H_028';"""
    prs = pd.read_sql(qry, paros_connection())
    prs = pd.concat(
        [
            (
                prs.loc[
                    prs["plate"].isin(["TFv02N2H_023", "TFv02N2H_025", "TFv02N2H_026"]),
                    :,
                ]
                .pivot_table(index="well", columns="plate", values="score")
                .reset_index()
                .rename(
                    columns={
                        "TFv02N2H_023": "score_pair",
                        "TFv02N2H_025": "score_empty-N2",
                        "TFv02N2H_026": "score_empty-N1",
                    }
                )
                .assign(plate=1)
            ),
            (
                prs.loc[
                    prs["plate"].isin(["TFv02N2H_024", "TFv02N2H_027", "TFv02N2H_028"]),
                    :,
                ]
                .pivot_table(index="well", columns="plate", values="score")
                .reset_index()
                .rename(
                    columns={
                        "TFv02N2H_024": "score_pair",
                        "TFv02N2H_027": "score_empty-N2",
                        "TFv02N2H_028": "score_empty-N1",
                    }
                )
                .assign(plate=2)
            ),
        ]
    )
    qry = """select plate, well, orf_ida as test_orf_ida, orf_idb as test_orf_idb, category as source
            from n2h.prsrrs
            where is_v2 = 1;"""
    prs = pd.merge(
        prs, pd.read_sql(qry, paros_connection()), how="inner", on=["plate", "well"]
    )

    prs["pair"] = (
        prs[["test_orf_ida", "test_orf_idb"]].min(axis=1).astype(str)
        + "_"
        + prs[["test_orf_ida", "test_orf_idb"]].max(axis=1).astype(str)
    )
    prs = prs.rename(columns={"well": "test_pos_pair"})
    prs["test_pos_empty-N1"] = prs["test_pos_pair"]
    prs["test_pos_empty-N2"] = prs["test_pos_pair"]
    prs["test_pla"] = prs["plate"].map({1: 23, 2: 24})
    prs = prs.drop(columns=["plate"])
    df = pd.concat([df, prs], ignore_index=True)

    df["source"] = df["source"].apply(
        lambda x: {"isoform_matched_negatives": "isoform_negatives"}.get(x, x)
    )
    df["source"] = df["source"].map(
        {
            "whole_tf_genes": "vignettes",
            "isoform_positives": "isoform positives",
            "TF_RRS": "RRS - TF space specific",
            "litbm": "Lit-BM - TF space specific",
            "isoform_negatives": "isoform negatives",
            "huri_rrs": "RRS - from HuRI",
            "lit_bm_2013_rand250": "Lit-BM-13",
            "PRS": "PRS - hPRS-v2",
            "RRS": "RRS - hRRS-v2",
        }
    )

    qry = """
    select orf_id, symbol, unique_acc
    from tf_screen.iso6k_sequences;
    """
    iso_orf_id_to_gene = (
        pd.read_sql(qry, paros_connection()).set_index("orf_id")["symbol"].to_dict()
    )
    orf_id_to_iso_acc = (
        pd.read_sql(qry, paros_connection()).set_index("orf_id")["unique_acc"].to_dict()
    )
    qry = qry = """
    select orf_id, entrez_gene_symbol
    from hi_ref.master_ref;
    """
    orf_id_to_partner_gene = (
        pd.read_sql(qry, paros_connection())
        .set_index("orf_id")["entrez_gene_symbol"]
        .to_dict()
    )
    df["clone_acc"] = df["test_orf_idb"].map(orf_id_to_iso_acc)
    df["gene_symbol_tf"] = df["test_orf_idb"].map(iso_orf_id_to_gene)
    df["gene_symbol_partner"] = df["test_orf_ida"].map(orf_id_to_partner_gene)

    y2h = load_full_y2h_data_including_controls()
    y2h = y2h.loc[y2h["category"] == "tf_isoform_ppis", :]
    y2h_positives = y2h.loc[
        y2h["Y2H_result"] == True, ["ad_orf_id", "db_orf_id"]
    ].values
    y2h_positives = set(map(tuple, y2h_positives))
    y2h_negatives = y2h.loc[
        y2h["Y2H_result"] == False, ["ad_orf_id", "db_orf_id"]
    ].values
    y2h_negatives = set(map(tuple, y2h_negatives))

    # check positives / negatives are in Y2H dataset
    def in_y2h_positves(row):
        pair = (row["test_orf_idb"], row["test_orf_ida"])
        return pair in y2h_positives

    def in_y2h_negatives(row):
        pair = (row["test_orf_idb"], row["test_orf_ida"])
        return pair in y2h_negatives

    df = df.loc[
        ~((df["source"] == "isoform positives") & ~df.apply(in_y2h_positves, axis=1)), :
    ]
    df = df.loc[
        ~((df["source"] == "isoform negatives") & ~df.apply(in_y2h_negatives, axis=1)),
        :,
    ]

    df = df.drop_duplicates()

    return df


func_to_file_name = [
    (isoform_clones, "isoform_clones"),
    (tf_isoform_y2h_screen_results, "tf_isoform_y2h_screen"),
    (isoform_and_paralog_y2h_data, "y2h_data"),
    (n2h_ppi_validation_data, "N2H_results"),
]
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / "data/internal"


def write_all_tables():
    for func, filename in func_to_file_name:
        func().to_csv(DATA_DIR / (filename + ".tsv"), sep="\t", index=False)


def check_consistency():
    all_good = True
    for func, filename in func_to_file_name:
        a = func()
        b = pd.read_csv(
            DATA_DIR / (filename + ".tsv"),
            sep="\t",
            na_values=[""],
            keep_default_na=False,
        )
        if not a.equals(b):
            all_good = False
            print(filename + ": inconsitent with database")
    if all_good:
        print("All files consistent with database")
