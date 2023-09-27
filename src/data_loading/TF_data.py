import pandas as pd

from .utils import DATA_DIR


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
