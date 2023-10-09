import pandas as pd
from goatools import obo_parser

from .utils import DATA_DIR, cache_with_pickle


def load_tf_families():
    """From the Lambert et al. review in Cell 2018

    Returns:
        pandas.Series: HGNC gene symbol to TF family

    """
    tf_fam = load_human_tf_db()
    tf_fam = tf_fam.loc[:, ["HGNC symbol", "DBD"]]
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


@cache_with_pickle
def load_cofactors():
    """ """
    cof = pd.read_csv(
        DATA_DIR / "external/AnimalTFDB4_Homo_sapiens_TF_cofactors.txt", sep="\t"
    )
    if cof["Symbol"].duplicated().any():
        raise UserWarning("unexpected duplicates")
    cof = cof.rename(columns={"Symbol": "gene_symbol"})

    tfdb = load_human_tf_db()
    tf_genes = set(tfdb["HGNC symbol"].unique())
    cof = cof.loc[~cof["gene_symbol"].isin(tf_genes), :]

    go = load_gene_ontology()
    goa = load_gene_ontology_annotations()
    GO_TERM_ID_ACTIVATOR = "GO:0003713"  # transcription coactivator  activity
    GO_TERM_ID_COREPRESSOR = "GO:0003714"  # transcription corepressor activity
    coactivator_go_terms = (
        go[GO_TERM_ID_ACTIVATOR].get_all_lower().union({GO_TERM_ID_ACTIVATOR})
    )
    corepressor_go_terms = (
        go[GO_TERM_ID_COREPRESSOR].get_all_lower().union({GO_TERM_ID_COREPRESSOR})
    )
    coactivator_genes = set(
        goa.loc[goa["GO ID"].isin(coactivator_go_terms), "DB Object Symbol"].unique()
    )
    coactivator_genes = coactivator_genes.difference(tf_genes)
    corepressor_genes = set(
        goa.loc[goa["GO ID"].isin(corepressor_go_terms), "DB Object Symbol"].unique()
    )
    corepressor_genes = corepressor_genes.difference(tf_genes)
    coactivator_corepressor_genes = coactivator_genes.intersection(corepressor_genes)
    coactivator_genes = coactivator_genes.difference(coactivator_corepressor_genes)
    corepressor_genes = corepressor_genes.difference(coactivator_corepressor_genes)
    cof["cofactor_type"] = "unknown"
    cof.loc[cof["gene_symbol"].isin(coactivator_genes), "cofactor_type"] = "coactivator"
    cof.loc[cof["gene_symbol"].isin(corepressor_genes), "cofactor_type"] = "corepressor"
    cof.loc[
        cof["gene_symbol"].isin(coactivator_corepressor_genes), "cofactor_type"
    ] = "both"

    return cof


def load_gene_ontology():
    go = obo_parser.GODag(
        DATA_DIR / "external/go-basic.obo", optional_attrs={"relationship"}
    )
    return go


@cache_with_pickle
def load_gene_ontology_annotations():
    # copied from: https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
    GAF_FORMAT_COLUMNS = [
        "DB",
        "DB Object ID",
        "DB Object Symbol",
        "Qualifier",
        "GO ID",
        "DB:Reference (IDB:Reference)",
        "Evidence Code",
        "With (or) From",
        "Aspect",
        "DB Object Name",
        "DB Object Synonym (ISynonym)",
        "DB Object Type",
        "Taxon(Itaxon)",
        "Date",
        "Assigned By",
        "Annotation Extension",
        "Gene Product Form ID",
    ]
    goa = pd.read_csv(
        DATA_DIR / "external/goa_human.gaf",
        sep="\t",
        comment="!",
        low_memory=False,
        names=GAF_FORMAT_COLUMNS,
    )
    goa = goa.loc[
        (~goa["Qualifier"].str.startswith("NOT|")) & goa["DB Object Symbol"].notnull(),
        :,
    ]
    return goa


@cache_with_pickle
def load_signaling_genes():
    tfdb = load_human_tf_db()
    tf_genes = set(tfdb["HGNC symbol"].unique())
    cof = load_cofactors()
    cofactor_genes = set(cof["gene_symbol"].values)
    go = load_gene_ontology()
    goa = load_gene_ontology_annotations()

    GO_TERM_ID_SIGNALING = "GO:0023052"  # signaling
    signaling_go_terms = (
        go[GO_TERM_ID_SIGNALING].get_all_lower().union(GO_TERM_ID_SIGNALING)
    )
    signaling_genes = set(
        goa.loc[goa["GO ID"].isin(signaling_go_terms), "DB Object Symbol"].unique()
    )
    signaling_genes = signaling_genes.difference(tf_genes.union(cofactor_genes))

    # BUG! nan in signaling genes
    return signaling_genes
