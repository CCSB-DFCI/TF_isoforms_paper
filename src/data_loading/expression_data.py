import numpy as np
import pandas as pd

from .utils import cache_with_pickle, DATA_DIR
from .clones_and_assays_data import load_valid_isoform_clones
from .load_annotated_TFs import (
    load_annotated_gencode_tfs,
    load_annotated_TFiso1_collection,
)


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
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.isoforms]):
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
    tfs = load_annotated_TFiso1_collection()
    df.index = df.index.map(lambda x: _convert_to_joint_clone_and_ensembl_id(x, tfs))
    df = df.loc[df.index != "fail", :]
    df = df.groupby("UID").sum()
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.isoforms]):
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
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.isoforms]):
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
    tfs = load_annotated_TFiso1_collection()
    df.index = df.index.map(lambda x: _convert_to_joint_clone_and_ensembl_id(x, tfs))
    df = df.loc[df.index != "fail", :]
    df = df.groupby("UID").sum()
    if df.shape[0] != len([orf for tf in tfs.values() for orf in tf.isoforms]):
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
