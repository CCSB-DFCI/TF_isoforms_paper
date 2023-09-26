import itertools
import functools
from pathlib import Path

import numpy as np
import pandas as pd
import tqdm
from Bio.PDB.DSSP import make_dssp_dict
from Bio.Data.IUPACData import protein_letters_3to1

from .utils import DATA_DIR, CACHE_DIR


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


"""
TODO: find where load_isoforms_of_paralogs_pairs went and add this to the above
function

from Bio import Align

from data_loading import load_valid_isoform_clones, load_paralog_pairs, load_isoforms_of_paralogs_pairs

def calculate_aa_sequence_id():
    isoforms = load_valid_isoform_clones()
    pairs = load_paralog_pairs()
    pairs = load_isoforms_of_paralogs_pairs(pairs, isoforms)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'

    def pairwise_seq_id(row):
        alignment = aligner.align(row['aa_seq_a'], row['aa_seq_b'])[0].__str__().split()[1]
        return alignment.count('|') / len(alignment) * 100

    pairs['pct_aa_seq_id'] = pairs.apply(pairwise_seq_id, axis=1)

    (pairs.loc[:, ['clone_acc_a', 'clone_acc_b', 'pct_aa_seq_id']]
        .to_csv('../data/processed/paralog_non_paralog_seq_id.tsv',
                sep='\t',
                index=False))

"""


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


def load_pfam_domains_TFiso1():
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


def load_pfam_domains_horfeome():
    """Filter and format the Pfam domain matches in the human ORFeome.

    Returns:
        pandas.DataFrame: ORFeome annotated with Pfam domains

    """
    filepath = Path("../data/internal/horfeome_hmmscan_pfam.domtab")
    evalue_cutoff = 1e-5
    pfam = read_hmmer3_domtab(filepath)
    pfam = pfam.loc[
        pfam["E-value"] <= evalue_cutoff,
        [
            "query name",
            "target accession",
            "target name",
            "description of target",
            "ali_coord_from",
            "ali_coord_to",
            "tlen",
        ],
    ].rename(
        columns={
            "query name": "orf_id",
            "target accession": "pfam_accession",
            "target name": "domain_name",
            "description of target": "domain_description",
            "ali_coord_from": "start",
            "ali_coord_to": "stop",
            "tlen": "domain_length",
        }
    )
    pfam["pfam_accession"] = pfam["pfam_accession"].str.replace(r"\..*", "", regex=True)
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


def convert_old_id_to_new_id(old_id):
    if len(old_id.split("xxx")) != 3:
        raise UserWarning("Unrecognized old isoform clone ID")
    return old_id.split("xxx")[1]


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


def load_disorder_predictions():
    cache_path = Path("../data/processed/TFiso1_disorder-and-ss_from-alphafold.tsv")
    if cache_path.exists():
        return pd.read_csv(cache_path, sep="\t")
    dssp_dir = Path("../data/processed/dssp_alphafold")
    dfs = []
    for dssp_file_path in dssp_dir.iterdir():
        dssp = make_dssp_dict(dssp_file_path)
        dfs.append(
            pd.DataFrame(
                data=[
                    (dssp_file_path.stem, k[1][1], v[0], v[1], v[2])
                    for k, v in dssp[0].items()
                ],
                columns=["clone_name", "position", "aa", "secondary_structure", "ASA"],
            )
        )
    df = pd.concat(dfs, axis=0, ignore_index=True)
    # NOTE: the Davey analysis uses GGXGG whereas I think this paper is GXG
    # Wilke: Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635
    max_asa = {
        "ALA": 129.0,
        "ARG": 274.0,
        "ASN": 195.0,
        "ASP": 193.0,
        "CYS": 167.0,
        "GLN": 225.0,
        "GLU": 223.0,
        "GLY": 104.0,
        "HIS": 224.0,
        "ILE": 197.0,
        "LEU": 201.0,
        "LYS": 236.0,
        "MET": 224.0,
        "PHE": 240.0,
        "PRO": 159.0,
        "SER": 155.0,
        "THR": 172.0,
        "TRP": 285.0,
        "TYR": 263.0,
        "VAL": 174.0,
    }
    max_asa = {protein_letters_3to1[k.capitalize()]: v for k, v in max_asa.items()}
    df["RSA"] = df["ASA"] / df["aa"].map(max_asa)
    df["RSA"] = df["RSA"].clip(upper=1.0)
    WINDOW_SIZE_RESIDUES = 20
    DISORDER_WINDOW_RSA_CUTOFF = 0.5
    rsa_window_col = f"RSA_window_{WINDOW_SIZE_RESIDUES}"
    df[rsa_window_col] = (
        df.groupby("clone_name")["RSA"]
        .rolling(
            window=WINDOW_SIZE_RESIDUES * 2 + 1,
            min_periods=WINDOW_SIZE_RESIDUES + 1,
            center=True,
        )
        .mean()
        .rename(rsa_window_col)
        .droplevel("clone_name")
    )
    df["is_disordered"] = df[rsa_window_col] >= DISORDER_WINDOW_RSA_CUTOFF

    DISORDER_HELIX_LENGTH_CUTOFF = 20
    to_change = []
    for clone_name, df_clone in df.groupby("clone_name"):
        helix_count = 0
        for _i, row in df_clone.iterrows():
            if row["secondary_structure"] == "H":
                helix_count += 1
            else:
                if helix_count >= DISORDER_HELIX_LENGTH_CUTOFF:
                    for i in range(
                        row["position"] - 1, row["position"] - helix_count, -1
                    ):
                        to_change.append((clone_name, i))
                helix_count = 0
        if helix_count >= DISORDER_HELIX_LENGTH_CUTOFF:
            for i in range(row["position"], row["position"] - helix_count, -1):
                to_change.append((clone_name, i))
    to_change = (df["clone_name"] + "_" + df["position"].astype(str)).isin(
        {a + "_" + str(b) for a, b in to_change}
    )
    print(
        f"{to_change.sum()} ({to_change.mean():.0%}) aa in helices of length 20 aa or more"
    )
    print(
        f"{df.loc[to_change, 'is_disordered'].mean():.0%} of residues in long helices misclassified as disordered"
    )
    df.loc[to_change, "is_disordered"] = False

    df.to_csv(cache_path, index=False, sep="\t")
    return df
