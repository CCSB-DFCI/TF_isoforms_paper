"""
We use GENCODE v30 for transcripts but we use the TF db HGNC symbols
which leads to 13 cases where the gene and transcript names aren't consistent
E.g. ZZZ3 has transcripts AC118549.1-202 and AC118549.1-204
"""

from pathlib import Path
import warnings
import re
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
import tqdm

from . import isolib
from .utils import cache_with_pickle, DATA_DIR, CACHE_DIR
from .TF_data import load_human_tf_db, load_tf_families
from .clones_and_assays_data import load_valid_isoform_clones
from .sequence_features import (
    load_pfam_domains_TFiso1,
    load_pfam_domains_gencode,
    load_pfam_clans,
)


DIMERIZING_TF_FAMILIES = {
    "bHLH",
    "bZIP",
    "Nuclear receptor",
    "E2F",
    "CENPB",
    "MADS box",
    "Grainyhead",
    "SAND",
    "Rel",
    "EBF1",
    "STAT",
    "IRF",
    "RFX",
    "HSF",
    "p53",
    "AP-2",
    "GCM",
    "BED ZF",
    "MADF",
    "ARID/BRIGHT",
    "Myb/SANT",
    "SMAD",
}


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
            cds = _extract_mRNA_sequence_region_from_fasta(
                nt_seq[transcript_id], "CDS", raise_error=True
            )
            utr5 = _extract_mRNA_sequence_region_from_fasta(
                nt_seq[transcript_id], "UTR5", raise_error=False
            )
            utr3 = _extract_mRNA_sequence_region_from_fasta(
                nt_seq[transcript_id], "UTR3", raise_error=False
            )
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
    _add_Pfam_domains_gencode(
        genes,
        transcript_name_to_ensembl_gene_id=transcript_name_to_ensembl_gene_id,
        genes_to_rename=genes_to_rename,
        unique_pc_transcripts=unique_pc_transcripts,
    )
    _make_c2h2_zf_arrays(genes)
    _add_dbd_flanks(genes)
    _add_MANE_and_APPRIS_annoations(
        genes=genes, path_MANE_select=path_MANE_select, path_APPRIS=path_APPRIS
    )
    _add_effector_domains(genes)

    return genes


@cache_with_pickle
def load_annotated_TFiso1_collection(
    path_TFiso1_gtf=DATA_DIR / "internal/c_6k_unique_acc_aligns.gtf",
    path_TFiso1_fa=DATA_DIR / "internal/j2_6k_unique_isoacc_and_nt_seqs.fa",
    path_gencode_aa_seq=DATA_DIR / "external/gencode.v30.pc_translations.fa",
    path_MANE_select=DATA_DIR / "external/MANE.GRCh38.v0.95.summary.txt",
    path_APPRIS=DATA_DIR
    / "external/APPRIS-annotations_human_GRCh38.p13_ensembl104.tsv",
    path_disorder=DATA_DIR / "processed/TFiso1_disorder-and-ss_from-alphafold.tsv",
):
    """ """
    import pyranges

    # note that pyranges switches the indexing to python 0-indexed, half-open interval
    algn = pyranges.read_gtf(path_TFiso1_gtf).df
    algn = algn.loc[algn["Start"] < algn["End"], :]  # filter out problematic rows
    nt_seq = {r.name: r for r in SeqIO.parse(path_TFiso1_fa, format="fasta")}
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
            raise ValueError(", ".join(missing) + " not in " + path_TFiso1_fa)
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
    _add_Pfam_domains_TFiso1(genes)
    _make_c2h2_zf_arrays(genes)
    _add_dbd_flanks(genes)

    tfs_gencode = load_annotated_gencode_tfs(
        subset={g.ensembl_gene_id for g in genes.values()}
    )
    uncloned_orfs = defaultdict(list)
    for tf in genes.values():
        for gencode_isoform in tfs_gencode[tf.name].isoforms:
            clone_match = False
            for cloned_isoform in tf.isoforms:
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

    _add_MANE_and_APPRIS_annoations(
        genes=genes, path_MANE_select=path_MANE_select, path_APPRIS=path_APPRIS
    )
    _add_effector_domains(genes=genes)
    _add_disordered_regions(genes=genes, path_disorder=path_disorder)

    return genes


def _add_Pfam_domains_TFiso1(genes):
    pfam = load_pfam_domains_TFiso1()
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


def _add_Pfam_domains_gencode(
    genes, transcript_name_to_ensembl_gene_id, genes_to_rename, unique_pc_transcripts
):
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


def _add_MANE_and_APPRIS_annoations(genes, path_MANE_select, path_APPRIS):
    mane = pd.read_csv(path_MANE_select, sep="\t")
    mane_select = set(
        mane.loc[mane["MANE_status"] == "MANE Select", "Ensembl_nuc"]
        .str.slice(0, 15)
        .values
    )
    for tf in genes.values():
        if tf.ensembl_gene_id not in mane["Ensembl_Gene"].str.slice(0, 15).values:
            continue  # not all genes have a MANE select isoform
        for iso in tf.isoforms:
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
        for iso in tf.isoforms:
            if hasattr(iso, "clone_acc") and iso.is_novel_isoform():
                continue
            annotations = {
                appris[tid] for tid in iso.ensembl_transcript_ids if tid in appris
            }
            if len(annotations) > 0:
                iso.APPRIS_annotation = _consolidate_appris_annotations(annotations)


def _add_effector_domains(genes):
    path_Soto_et_al = (
        DATA_DIR / "external/Soto-et-al_MolCell_2022_Supplementary-tables.xlsx"
    )
    path_Tycko_et_al = DATA_DIR / "external/Tycko-et-al_Cell_2020_Table-S4.xlsx"
    path_DelRosso_et_al = (
        DATA_DIR / "external/DelRosso-et-al_Nature_2023_Supplementary-Table-2.xlsx"
    )

    reg_dom = pd.read_excel(path_Soto_et_al, sheet_name="Table S2")
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
            for iso in tf.isoforms:
                if row["Sequence"] not in iso.aa_seq:
                    continue
                if len(re.findall("(?={})".format(row["Sequence"]), iso.aa_seq)) != 1:
                    raise UserWarning(
                        "Problem mapping effector domain: {} {}".format(row, iso)
                    )
                iso.add_aa_seq_feature(
                    category="effector_domain",
                    name=row["Domain type"],
                    accession="Soto_" + row["Effector domain ID"],
                    start=iso.aa_seq.find(row["Sequence"]),
                    end=iso.aa_seq.find(row["Sequence"]) + len(row["Sequence"]),
                    description=desc,
                )

    tycko = pd.concat(
        [
            pd.read_excel(path_Tycko_et_al, sheet_name=sheet_name)
            for sheet_name in ["NucRepr_data", "NucAct_data", "Tiling Repressors"]
        ]
    )
    tycko = tycko.loc[tycko["Hit"], :]
    tycko["Sequence"] = tycko["Sequence"].fillna(tycko["Extended Domain sequence"])
    tycko["Domain type"] = "RD"
    tycko.loc[tycko["Avg Act"].notnull(), "Domain type"] = "AD"
    # TODO: fix this, it's very ugly at the moment
    tycko["domain_accession"] = (
        "Tycko_"
        + tycko["Domain type"]
        + "_"
        + tycko["Gene entry name"].fillna("")
        + "_"
        + tycko["Domain ID"].fillna("")
        + "_"
        + ("tile-" + tycko["Tile number"].astype(str)).fillna("")
    )

    delrosso = pd.concat(
        [
            pd.read_excel(path_DelRosso_et_al, sheet_name=sheet_name)
            for sheet_name in ["Activation Domains", "Repression Domains"]
        ]
    )
    delrosso["Domain type"] = delrosso["Domain type"].map(
        {"AD": "AD", "pEF": "RD", "PGK": "RD", "PGKandpEF": "RD"}
    )
    # I looked at the DelRosso data and the HGNC names are consistent with
    # what we're using. The UniProt ACs have some disagreements
    # For 3% of these domains, the start and end positions don't match
    # the sequnece lenght.
    # There are 19 out of 475 that don't map to our cloned isoforms.
    # Some could be becuase of missing the reference sequence, but
    # I looked at a couple and they are due to small sequence differences
    # (e.g. FOXD4L1). Ignoring for now.
    #  - Luke

    for tf in genes.values():  # skipping the gene name matching with this one
        for _i, row in tycko.iterrows():
            desc = "Tycko et al. Cell 2020"
            for iso in tf.isoforms:
                if row["Sequence"] not in iso.aa_seq:
                    continue
                if len(re.findall("(?={})".format(row["Sequence"]), iso.aa_seq)) != 1:
                    raise UserWarning(
                        "Problem mapping effector domain: {} {}".format(row, iso)
                    )
                iso.add_aa_seq_feature(
                    category="effector_domain",
                    name=row["Domain type"],
                    accession=row["domain_accession"],
                    start=iso.aa_seq.find(row["Sequence"]),
                    end=iso.aa_seq.find(row["Sequence"]) + len(row["Sequence"]),
                    description=desc,
                )
    for tf in genes.values():
        for _i, row in delrosso.loc[delrosso["HGNC symbol"] == tf.name, :].iterrows():
            desc = "DelRosso et al. Nature 2022"
            for iso in tf.isoforms:
                if row["Sequence"] not in iso.aa_seq:
                    continue
                if len(re.findall("(?={})".format(row["Sequence"]), iso.aa_seq)) != 1:
                    raise UserWarning(
                        "Problem mapping effector domain: {} {}".format(row, iso)
                    )
                iso.add_aa_seq_feature(
                    category="effector_domain",
                    name=row["Domain type"],
                    accession="DelRosso_{}_{}_{}".format(
                        row["Domain type"], row["HGNC symbol"], row["Domain"]
                    ),
                    start=iso.aa_seq.find(row["Sequence"]),
                    end=iso.aa_seq.find(row["Sequence"]) + len(row["Sequence"]),
                    description=desc,
                )


def _add_disordered_regions(genes, path_disorder):
    """
    NOTE: this is just implemented for the cloned isoforms, not the full
    GENCODE set.
    """
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


def _extract_mRNA_sequence_region_from_fasta(rec, region, raise_error=True):
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
    start, stop = bounds[0][(len(region) + 1) :].split("-")
    return str(rec.seq)[(int(start) - 1) : int(stop)]


def _make_c2h2_zf_arrays(tfs, MAX_NUM_AA_C2H2_ZF_SEPERATION=10):
    clans = load_pfam_clans()
    C2H2_ZF_PFAM_CLAN_AC = "CL0361"
    c2h2_zf_pfam_ids = {k for k, v in clans.items() if v == C2H2_ZF_PFAM_CLAN_AC}
    for gene in tfs.values():
        for orf in gene.isoforms:
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
        for isoform in gene.isoforms:
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
