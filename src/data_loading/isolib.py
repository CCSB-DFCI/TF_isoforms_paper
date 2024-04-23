"""
authors: Gloria Sheynkman, Luke Lambourne
"""

import itertools
import random
import collections

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import pandas as pd
import seaborn as sns
from Bio.Seq import Seq
from Bio import Align


# moving this function as it is useful in other areas of the code
def _coords_transform_aa_seq_to_alignment(i, alignment):
    if i > len(alignment.replace("I", "")):
        raise ValueError("position is not in isoform AA sequence")
    aa_seq_indices = [
        "" if c == "I" else len(alignment[:j].replace("I", ""))
        for j, c in enumerate(alignment)
    ]
    return aa_seq_indices.index(i)


class GenomicFeature:
    """A contigous region on the genome

    Attributes:
        chrom (str): chromosome identifier
        strand (str): must be '+'/'-'
        start/end (int): chromosomal coordinates, follows python indexing
            convention of 0-indexed half-open interval, start must be a
            lower number than end, i.e. it is not the 'start'/'end' in
            the reading direction if it's on the negative strand.

    """

    def __init__(self, chrom, strand, start, end):
        self.chrom = chrom
        if strand not in {"+", "-"}:
            msg = strand + " is invalid value for strand (must be +/-)"
            raise ValueError(msg)
        self.strand = strand
        if start >= end:
            raise ValueError("start must be before end")
        self.start = start
        self.end = end

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        return (
            (self.chrom == other.chrom)
            and (self.start == other.start)
            and (self.end == other.end)
            and (self.strand == other.strand)
        )

    def __ne__(self, other):
        return not self.__eq__(other)


class Pathogenic_Coding_SNP(GenomicFeature):
    def __init__(
        self, chrom, strand, position, nt_change, aa_change, disease, mutation_id
    ):
        self.nt_change = nt_change
        self.aa_change = aa_change
        self.disease = disease
        self.mutation_id = mutation_id
        GenomicFeature.__init__(
            self,
            chrom,
            strand,
            position,
            position + 1,
        )


class ProteinSequenceFeature:
    """A contiguous stretch of amino acids within a protein.

    category (str): general class of feature, e.g. 'Pfam_domain', 'ELM_motif'
    accession (str): identifier for feature, e.g. 'PF00170'
    name (str): display name, e.g. 'bZIP'
    start/end (int): position within the protein sequence, follows python indexing
        convention of 0-indexed half-open interval, start must be before end

    """

    def __init__(self, category, accession, name, start, end, description=None):
        if start >= end:
            raise ValueError("start must be before end")
        self.start = start
        self.end = end
        self.category = category
        self.name = name
        self.accession = accession
        self.description = description

    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        s = "{}: {} {} {}-{}".format(
            self.category, self.accession, self.name, self.start + 1, self.end
        )
        return s


class Gene(GenomicFeature):
    """A collection of protein-coding isoforms belonging to the same gene

    Attributes:
        name (str): gene symbol
        isoforms (list(isolib.Isoform)): protein coding isoforms of gene

    """

    # NOTE: if you change arguments here, also need to change in add_isoforms method
    def __init__(self, name, isoforms, ensembl_gene_id=None, uniprot_ac=None):
        chroms = [iso.chrom for iso in isoforms]
        strands = [iso.strand for iso in isoforms]
        msg = "{} - {}".format(name, ", ".join([iso.name for iso in isoforms]))
        if len(isoforms) == 0:
            raise ValueError("Need at least one isoform to define a gene\n" + msg)
        if len(set(chroms)) > 1:
            raise ValueError("All isoforms must be on same chromosome\n" + msg)
        if len(set(strands)) > 1:
            raise ValueError("All isoforms must be same strand\n" + msg)
        chrom = chroms[0]
        strand = strands[0]
        self.ensembl_gene_id = ensembl_gene_id
        self.uniprot_ac = uniprot_ac

        # de-duplicate exons
        exon_bounds = {(exon.start, exon.end) for iso in isoforms for exon in iso.exons}
        exon_bounds = list(sorted(exon_bounds, key=lambda x: x[0]))
        merged_exon_bounds = []
        for i in range(len(exon_bounds) - 1):
            start_a, end_a = exon_bounds[i]
            start_b, end_b = exon_bounds[i + 1]
            if end_a < start_b:
                merged_exon_bounds.append(exon_bounds[i])
            else:  # overlapping exons
                exon_bounds[i + 1] = (start_a, max(end_a, end_b))
        merged_exon_bounds.append(exon_bounds[-1])
        if strand == "-":
            merged_exon_bounds = merged_exon_bounds[::-1]
        self.exons = []
        for start, end in merged_exon_bounds:
            self.exons.append(Exon(name, name, chrom, strand, start, end))
        # set exon number by gene of each isoform
        for iso in isoforms:
            for isoform_exon in iso.exons:
                for i, gene_exon in enumerate(self.exons):
                    if (
                        isoform_exon.start >= gene_exon.start
                        and isoform_exon.end <= gene_exon.end
                    ):
                        isoform_exon.exon_number_on_gene = i + 1
                        break

        self.number_of_isoforms = len(isoforms)
        self.name = name

        self._isoforms = list(
            sorted(isoforms, key=lambda x: int(x.name.split("-")[-1]))
        )

        self._iso_dict = {iso.name: iso for iso in self._isoforms}
        self._pairwise_changes = {}
        self.pathogenic_coding_SNPs = []
        GenomicFeature.__init__(
            self,
            chrom,
            strand,
            min([iso.start for iso in isoforms]),
            max([iso.end for iso in isoforms]),
        )

    def add_isoforms(self, new_isoforms):
        # NOTE: this is a little dangerous. If new arguments are added to
        # __init__ then they need to be added here. I couldn't think of a
        # better way to do it...
        self.__init__(
            self.name,
            self.isoforms + new_isoforms,
            ensembl_gene_id=self.ensembl_gene_id,
            uniprot_ac=self.uniprot_ac,
        )

    def add_pathogenic_coding_SNP(
        self, position, nt_change, aa_change, disease, mutation_id
    ):
        self.pathogenic_coding_SNPs.append(
            Pathogenic_Coding_SNP(
                self.chrom,
                self.strand,
                position,
                nt_change,
                aa_change,
                disease,
                mutation_id,
            )
        )

    @property
    def cloned_isoforms(self):
        return [iso for iso in self.isoforms if hasattr(iso, "clone_acc")]

    @property
    def _cloned_isoforms(self):
        return [iso for iso in self._isoforms if hasattr(iso, "clone_acc")]

    @property
    def GENCODE_isoforms(self):
        return [iso for iso in self.isoforms if iso.ensembl_transcript_ids is not None]

    @property
    def _GENCODE_isoforms(self):
        return [iso for iso in self._isoforms if iso.ensembl_transcript_ids is not None]

    @property
    def MANE_select_isoform(self):
        for iso in self._isoforms:
            if not hasattr(iso, "is_MANE_select_transcript"):
                return None
            if iso.is_MANE_select_transcript:
                return iso
        # for small number of cases where transcript is not in list of isoforms
        # e.g. was added after that version of GENCODE
        return None

    @property
    def has_MANE_select_isoform(self):
        if not hasattr(self._isoforms[0], "is_MANE_select_transcript"):
            return False
        return any(iso.is_MANE_select_transcript for iso in self._isoforms)

    @property
    def APPRIS_isoforms(self):
        appris_isos = {}
        for iso in self._isoforms:
            if hasattr(iso, "APPRIS_annotation"):
                if iso.APPRIS_annotation is not None:
                    appris_isos[iso.APPRIS_annotation] = iso
        sorted_keys = sorted(
            appris_isos.keys(),
            key=lambda x: int(x[-1]) - 99 * x.startswith("principal"),
        )
        appris_isos = {k: appris_isos[k] for k in sorted_keys}
        return appris_isos

    @property
    def reference_isoform(self):
        if self.has_MANE_select_isoform:
            return self.MANE_select_isoform
        elif len(self.APPRIS_isoforms) > 0:
            return list(self.APPRIS_isoforms.values())[0]
        elif len(self._GENCODE_isoforms) > 0:
            return list(
                sorted(
                    self._GENCODE_isoforms, key=lambda x: len(x.aa_seq), reverse=True
                )
            )[0]
        else:
            return self._isoforms[0]

    @property
    def alternative_isoforms(self):
        return [iso for iso in self.isoforms if iso != self.reference_isoform]

    @property
    def cloned_APPRIS_isoforms(self):
        return {
            annot: iso
            for annot, iso in self.APPRIS_isoforms.items()
            if hasattr(iso, "clone_acc")
        }

    @property
    def cloned_reference_isoform(self):
        if self.has_MANE_select_isoform:
            if hasattr(self.MANE_select_isoform, "clone_acc"):
                return self.MANE_select_isoform
        if len(self.cloned_APPRIS_isoforms) > 0:
            return list(self.cloned_APPRIS_isoforms.values())[0]
        if any(not iso.is_novel_isoform() for iso in self._cloned_isoforms):
            return list(
                sorted(
                    [
                        iso
                        for iso in self._cloned_isoforms
                        if not iso.is_novel_isoform()
                    ],
                    key=lambda x: len(x.aa_seq),
                    reverse=True,
                )
            )[0]
        return list(
            sorted(self._cloned_isoforms, key=lambda x: len(x.aa_seq), reverse=True)
        )[0]

    @property
    def cloned_MANE_select_isoform(self):
        if self.has_MANE_select_isoform:
            return hasattr(self.MANE_select_isoform, "clone_acc")
        return None

    @property
    def isoforms(self):
        return self._order_isoforms_for_display()

    def _order_isoforms_for_display(self):
        order = []
        if len(self._cloned_isoforms) > 0:
            order.append(self.cloned_reference_isoform)
        for iso in sorted(
            self._cloned_isoforms, key=lambda x: len(x.aa_seq), reverse=True
        ):
            if iso.name == self.cloned_reference_isoform.name:
                continue
            if iso.is_novel_isoform():
                continue
            order.append(iso)
        for iso in sorted(
            self._cloned_isoforms, key=lambda x: len(x.aa_seq), reverse=True
        ):
            if iso.name == self.cloned_reference_isoform.name:
                continue
            if not iso.is_novel_isoform():
                continue
            order.append(iso)
        if not hasattr(self.reference_isoform, "clone_acc"):
            order.append(self.reference_isoform)
        for iso in sorted(self._isoforms, key=lambda x: len(x.aa_seq), reverse=True):
            if hasattr(iso, "clone_acc"):
                continue
            if iso.name == self.reference_isoform.name:
                continue
            order.append(iso)
        if len(order) != len(self._isoforms):
            raise UserWarning("bug in code" + str(self))
        return order

    def alternative_start(self, isoform_a, isoform_b):
        def _start_pos(iso):
            # NOTE: the exons just gives you the CDS positions
            return iso.exons[0].start if self.strand == "+" else iso.exons[0].end

        a = _start_pos(self._iso_dict[isoform_a])
        b = _start_pos(self._iso_dict[isoform_b])
        return a != b

    def alternative_stop(self, isoform_a, isoform_b):
        def _stop_pos(iso):
            # NOTE: the exons just gives you the CDS positions
            return iso.exons[-1].end if self.strand == "+" else iso.exons[-1].start

        a = _stop_pos(self._iso_dict[isoform_a])
        b = _stop_pos(self._iso_dict[isoform_b])
        return a != b

    def alternative_internal_exon(self, isoform_a, isoform_b):
        a = [e.exon_number_on_gene for e in self._iso_dict[isoform_a].exons]
        b = [e.exon_number_on_gene for e in self._iso_dict[isoform_b].exons]
        start = max(min(a), min(b))
        stop = min(max(a), max(b))
        internal_a = {e for e in a if e > start and e < stop}
        internal_b = {e for e in b if e > start and e < stop}
        return internal_a != internal_b

    def alternative_3prime_acceptor(self, isoform_a, isoform_b):
        for exon_a in self._iso_dict[isoform_a].exons:
            for exon_b in self._iso_dict[isoform_b].exons:
                if exon_a.exon_number_on_gene == exon_b.exon_number_on_gene:
                    if self.strand == "+":
                        if exon_a.start != exon_b.start:
                            return True
                    else:
                        if exon_a.end != exon_b.end:
                            return True
        return False

    def alternative_5prime_donor(self, isoform_a, isoform_b):
        for exon_a in self._iso_dict[isoform_a].exons:
            for exon_b in self._iso_dict[isoform_b].exons:
                if exon_a.exon_number_on_gene == exon_b.exon_number_on_gene:
                    if self.strand == "+":
                        if exon_a.end != exon_b.end:
                            return True
                    else:
                        if exon_a.start != exon_b.start:
                            return True
        return False

    def exon_skipping(self, isoform_a, isoform_b):
        a = [e.exon_number_on_gene for e in self._iso_dict[isoform_a].exons]
        b = [e.exon_number_on_gene for e in self._iso_dict[isoform_b].exons]
        start = max(min(a), min(b))
        stop = min(max(a), max(b))
        internal_a = {e for e in a if e > start and e < stop}
        internal_b = {e for e in b if e > start and e < stop}
        return len(internal_b) < len(internal_a) and internal_b.issubset(internal_a)

    def mutually_exclusive_exons(self, isoform_a, isoform_b):
        """NOTE: this is the not-strict MXE definition where only a pair of transcripts are considered"""
        a = [e.exon_number_on_gene for e in self._iso_dict[isoform_a].exons]
        b = [e.exon_number_on_gene for e in self._iso_dict[isoform_b].exons]
        start = max(min(a), min(b))
        stop = min(max(a), max(b))
        internal_a = {e for e in a if e > start and e < stop}
        internal_b = {e for e in b if e > start and e < stop}
        return (
            len(internal_a.difference(internal_b)) > 0
            and len(internal_b.difference(internal_a)) > 0
        )

    def intron_retention(self, isoform_a, isoform_b):
        exons_a = self._iso_dict[isoform_a].exons
        exons_b = self._iso_dict[isoform_b].exons
        if self.strand == "-":
            exons_a = list(reversed(exons_a))
            exons_b = list(reversed(exons_b))
        for exon_a in exons_a:
            for i in range(len(exons_b) - 1):
                if exon_a.start < exons_b[i].end and exon_a.end > exons_b[i + 1].start:
                    return True
        for exon_b in exons_b:
            for i in range(len(exons_a) - 1):
                if exon_b.start < exons_a[i].end and exon_b.end > exons_a[i + 1].start:
                    return True
        return False

    def splicing_categories(self, isoform_a, isoform_b):
        return {
            "gene_symbol": self.name,
            "reference isoform": isoform_a,
            "alternative isoform": isoform_b,
            "alternative N-terminal": self.alternative_start(isoform_a, isoform_b),
            "alternative C-terminal": self.alternative_stop(isoform_a, isoform_b),
            "alternative internal exon": self.alternative_internal_exon(
                isoform_a, isoform_b
            ),
            "alternative 5' splice site": self.alternative_5prime_donor(
                isoform_a, isoform_b
            ),
            "alternative 3' splice site": self.alternative_3prime_acceptor(
                isoform_a, isoform_b
            ),
            "exon skipping": self.exon_skipping(isoform_a, isoform_b),
            "mutually exclusive exons": self.mutually_exclusive_exons(
                isoform_a, isoform_b
            ),
            "intron retention": self.intron_retention(isoform_a, isoform_b),
        }

    def genomic_alignment_of_aa_seqs(self, subset=None):
        """genomic co-ordinates of translated regions"""
        all_iso_names = [iso.name for iso in self.isoforms]
        if subset is None:
            subset = all_iso_names
        isoforms = [iso for iso in self.isoforms if iso.name in subset]
        if len(isoforms) != len(subset):
            msg = "Missing isoforms: "
            msg += "/".join([s for s in subset if s not in all_iso_names])
            raise ValueError(msg)
        gene_coords = sorted(
            list(set([res.coords[1] for iso in isoforms for res in iso.residues]))
        )
        tracks = {}
        for iso in sorted(isoforms, key=lambda x: x.name):
            if iso.name not in subset:
                continue
            iso_aa = {res.coords[1]: res.aa for res in iso.residues}
            tracks[iso.name] = "".join([iso_aa.get(i, "-") for i in gene_coords])
        if self.strand == "-":
            tracks = {k: v[::-1] for k, v in tracks.items()}
        return tracks

    def pairwise_changes_relative_to_reference(self, ref_iso_name, alt_iso_name):
        """
            Try and recreate Gloria's string of MMMMIIIIMMMRRRRMMMMFFF etc etc..

        Doesn't currently do substitutions

        """
        key = (ref_iso_name, alt_iso_name)
        if key in self._pairwise_changes:  # check cache
            return self._pairwise_changes[key]
        alignment = ""
        ref_iter = iter(self._iso_dict[ref_iso_name].residues)
        alt_iter = iter(self._iso_dict[alt_iso_name].residues)
        ref_res = next(ref_iter)
        alt_res = next(alt_iter)
        while True:
            if not any(i in ref_res.coords for i in alt_res.coords):
                if (self.strand == "+" and ref_res.coords[2] < alt_res.coords[2]) or (
                    self.strand == "-" and ref_res.coords[2] > alt_res.coords[2]
                ):
                    alignment += "D"
                    try:
                        ref_res = next(ref_iter)
                    except StopIteration:
                        alignment += "I"  # since alt_iter was incremented
                        for _remaining in alt_iter:
                            alignment += "I"
                        break
                else:
                    alignment += "I"
                    try:
                        alt_res = next(alt_iter)
                    except StopIteration:
                        alignment += "D"  # since ref_iter was incremented
                        for _remaining in ref_iter:
                            alignment += "D"
                        break
            else:
                if ref_res.coords == alt_res.coords:
                    alignment += "M"
                elif any(i == j for i, j in zip(ref_res.coords, alt_res.coords)):
                    # Different exon junctions
                    if ref_res.aa == alt_res.aa:
                        alignment += "M"
                    else:
                        if (
                            (
                                ref_res.coords[0] < alt_res.coords[0]
                                and self.strand == "+"
                            )
                            or (
                                ref_res.coords[2] > alt_res.coords[2]
                                and self.strand == "+"
                            )
                            or (
                                ref_res.coords[0] > alt_res.coords[0]
                                and self.strand == "-"
                            )
                            or (
                                ref_res.coords[2] < alt_res.coords[2]
                                and self.strand == "-"
                            )
                        ):
                            alignment += "I"
                            try:
                                alt_res = next(alt_iter)
                                continue
                            except StopIteration:
                                alignment += "D"  # since ref_iter was incremented
                                for _remaining in ref_iter:
                                    alignment += "D"
                                break
                        elif (
                            (
                                ref_res.coords[0] > alt_res.coords[0]
                                and self.strand == "+"
                            )
                            or (
                                ref_res.coords[2] < alt_res.coords[2]
                                and self.strand == "+"
                            )
                            or (
                                ref_res.coords[0] < alt_res.coords[0]
                                and self.strand == "-"
                            )
                            or (
                                ref_res.coords[2] > alt_res.coords[2]
                                and self.strand == "-"
                            )
                        ):
                            alignment += "D"
                            try:
                                ref_res = next(ref_iter)
                                continue
                            except StopIteration:
                                alignment += "I"  # since alt_iter was incremented
                                for _remaining in alt_iter:
                                    alignment += "I"
                                break
                        else:  # Here just the middle nt in codon matched
                            msg = "Unexpected alignement issue between: "
                            msg += ref_iso_name + " and " + alt_iso_name
                            raise UserWarning(msg)
                else:  # Frameshift
                    if any(
                        ref_res.coords[i] == alt_res.coords[(i + 1) % 3]
                        for i in range(3)
                    ):
                        alignment += "F"
                    else:
                        alignment += "f"

                try:
                    ref_res = next(ref_iter)
                except StopIteration:
                    for _remaining in alt_iter:
                        alignment += "I"
                    break
                try:
                    alt_res = next(alt_iter)
                except StopIteration:
                    alignment += "D"  # since we incremented above already
                    for _remaining in ref_iter:
                        alignment += "D"
                    break
        self._pairwise_changes[key] = alignment  # cache result
        return alignment

    def aa_seq_disruption(self, ref_iso_name, alt_iso_name, domain_start, domain_end):
        """Get pairwise alignment of isoform protein sequences. Return fraction of
        domain and insertion
        """
        algn = self.pairwise_changes_relative_to_reference(ref_iso_name, alt_iso_name)
        if len(algn.replace("I", "")) != len(self._iso_dict[ref_iso_name].aa_seq):
            msg = "Something is wrong\n"
            msg += ref_iso_name + ", " + alt_iso_name
            raise UserWarning(msg)

        def _coords_transform_aa_seq_to_alignment(i, alignment):
            if i > len(alignment.replace("I", "")):
                raise ValueError("position is not in isoform AA sequence")
            aa_seq_indices = [
                "" if c == "I" else len(alignment[:j].replace("I", ""))
                for j, c in enumerate(alignment)
            ]
            return aa_seq_indices.index(i)

        start = _coords_transform_aa_seq_to_alignment(domain_start, algn)
        end = _coords_transform_aa_seq_to_alignment(domain_end - 1, algn) + 1
        return {
            "deletion": algn[start:end].count("D"),
            "insertion": algn[start:end].count("I"),
            "frameshift": algn[start:end].count("F") + algn[start:end].count("f"),
        }

    def aa_feature_disruption(self, ref_iso_name):
        """

        TODO: name this better

        Args:
            ref_iso_name (str): [description]

        Returns:
            pandas.DataFrame: [description]
        """
        results = []
        ref_iso = self._iso_dict[ref_iso_name]
        row = {"gene_symbol": self.name, "ref_iso": ref_iso_name}
        for aa_feature in ref_iso.aa_seq_features:
            for alt_iso_name, alt_iso in self._iso_dict.items():
                if alt_iso_name == ref_iso_name:
                    continue
                r = self.aa_seq_disruption(
                    ref_iso_name, alt_iso_name, aa_feature.start, aa_feature.end
                )
                row.update(
                    {
                        "alt_iso": alt_iso_name,
                        "accession": aa_feature.accession,
                        "category": aa_feature.category,
                        "start_in_ref_iso": aa_feature.start,
                        "end_in_ref_iso": aa_feature.end,
                        "length": aa_feature.end - aa_feature.start,
                    }
                )
                row.update(r)
                results.append(row.copy())
        results = pd.DataFrame(results)
        return results

    def null_fraction_per_aa_feature(self, ref_iso_name):
        """Fraction of aa features that would be affected if they were evenly
            distributed along the AA sequence of the reference isoform.

        Note: this only depends on the two isoforms and the length of the
        domain.

        Args:
            ref_iso_name (str)

        """
        results = []
        ref_iso = self._iso_dict[ref_iso_name]
        row = {"gene_symbol": self.name, "ref_iso": ref_iso_name}
        aa_feature_lengths = {x.end - x.start for x in ref_iso.aa_seq_features}
        for length in aa_feature_lengths:
            for alt_iso_name, alt_iso in self._iso_dict.items():
                if alt_iso_name == ref_iso_name:
                    continue
                row.update(
                    {
                        "alt_iso": alt_iso_name,
                        "length": length,
                    }
                )
                rs = self._null_feature_disruption(ref_iso_name, alt_iso_name, length)

                def is_disrupted_at_all(feature_alignment_count):
                    return any(v > 0 for v in feature_alignment_count.values())

                def is_disrupted_by_some_fraction(
                    feature_alignment_count, feature_length, fraction
                ):
                    """
                    This is not perfect, because the insertions could be either
                    in the middle of the domain or just between the first two
                    amino acids and this treats both cases equivalently.
                    """
                    return (
                        sum(feature_alignment_count.values())
                        >= feature_length * fraction
                    )

                row["null_fraction_affected_at_all"] = sum(
                    [int(is_disrupted_at_all(r)) for r in rs]
                ) / len(rs)
                for frac in [1.0, 0.9, 0.5, 0.1]:
                    row[f"null_fraction_affected_{frac * 100:.0f}pct"] = sum(
                        [
                            int(is_disrupted_by_some_fraction(r, length, frac))
                            for r in rs
                        ]
                    ) / len(rs)
                results.append(row.copy())
        results = pd.DataFrame(results)
        return results

    def _null_feature_disruption(self, ref_iso_name, alt_iso_name, feature_length):
        algn = self.pairwise_changes_relative_to_reference(ref_iso_name, alt_iso_name)
        len_ref_iso_aa_seq = len(self._iso_dict[ref_iso_name].aa_seq)
        if len(algn.replace("I", "")) != len_ref_iso_aa_seq:
            msg = "Something is wrong\n"
            msg += ref_iso_name + ", " + alt_iso_name
            raise UserWarning(msg)
        coords_ref_iso_aa_seq = [
            "" if c == "I" else len(algn[:j].replace("I", ""))
            for j, c in enumerate(algn)
        ]
        two_mers = {"deletion": [], "insertion": [], "frameshift": []}
        for i in range(len_ref_iso_aa_seq - 1):
            start = coords_ref_iso_aa_seq.index(i)
            end = coords_ref_iso_aa_seq.index(i + 1) + 1
            two_mer = algn[start:end]
            two_mers["deletion"].append(two_mer.count("D"))
            two_mers["insertion"].append(two_mer.count("I"))
            two_mers["frameshift"].append(two_mer.count("F") + two_mer.count("f"))
        results = [{k: sum(v[:feature_length]) for k, v in two_mers.items()}]
        for i in range((len_ref_iso_aa_seq - feature_length) - 1):
            r = {
                k: (v - two_mers[k][i]) + two_mers[k][i + feature_length]
                for k, v in results[-1].items()
            }
            results.append(r)
        return results

    def random_aa_seq_feature_shuffle(self, ref_iso_name, n, subset=None):
        """Pick random, non-overlapping, aa sequence locations for the features.

        Args:
            ref_iso_name (str):
            n (int): number of randomizations
            subset (set(str), optional): [description]. Defaults to None.
        """
        if len(self._iso_dict) == 1:
            return pd.DataFrame([])

        # inspired by: https://stackoverflow.com/questions/18641272
        def _non_overlapping_random_feature_positons(len_aa_seq, feature_and_len):
            if sum([v[1] for v in feature_and_len]) > len_aa_seq:
                raise UserWarning("Impossible: " + str(ref_iso_name))
            indices = range(len_aa_seq - sum(v[1] - 1 for v in feature_and_len))
            result = []
            offset = 0
            n = len(feature_and_len)
            for (acc, l), i in zip(
                random.sample(feature_and_len, n), sorted(random.sample(indices, n))
            ):
                i += offset
                result.append((acc, (i, i + l)))
                offset += l - 1
            return result

        results = []
        ref_iso = self._iso_dict[ref_iso_name]
        row = {"gene_symbol": self.name, "ref_iso": ref_iso_name}
        for i in range(n):
            row["random_sample"] = i
            rnd_pos = _non_overlapping_random_feature_positons(
                len(ref_iso.aa_seq),
                [
                    (f.accession, len(f))
                    for f in ref_iso.aa_seq_features
                    if subset is None or f.accession in subset
                ],
            )
            for accession, (rnd_start, rnd_end) in rnd_pos:
                for alt_iso_name, alt_iso in self._iso_dict.items():
                    if alt_iso_name == ref_iso_name:
                        continue
                    row.update({"alt_iso": alt_iso_name})
                    row.update({"accession": accession})
                    r = self.aa_seq_disruption(
                        ref_iso_name, alt_iso_name, rnd_start, rnd_end
                    )
                    row.update(r)
                    row.update({"length": rnd_end - rnd_start})
                    results.append(row.copy())
        results = pd.DataFrame(results)
        return results

    def disordered_fraction_of_different_regions(self, ref_iso_name, alt_iso_name):
        algn = self.pairwise_changes_relative_to_reference(ref_iso_name, alt_iso_name)
        if not hasattr(self[ref_iso_name], "disorder") or not hasattr(
            self[alt_iso_name], "disorder"
        ):
            return np.nan
        ref_iter = iter(self[ref_iso_name].disorder)
        alt_iter = iter(self[alt_iso_name].disorder)
        merged_disorder = []
        for pos in algn:
            if pos == "I":
                merged_disorder.append(next(alt_iter))
            elif pos == "D":
                merged_disorder.append(next(ref_iter))
            else:
                merged_disorder.append(next(ref_iter))
                next(alt_iter)

        return np.mean(
            [
                is_disordered
                for pos, is_disordered in zip(algn, merged_disorder)
                if pos != "M"
            ]
        )

    def _get_exon_colors(self):
        """ """
        merged_exon_bounds = [(exon.start, exon.end) for exon in self.exons]
        # colors_frame_1 = sns.cubehelix_palette(len(merged_exon_bounds),
        #                                       start=2.3,
        #                                       rot=0.9,
        #                                       light=0.7,
        #                                       dark=0.3,
        #                                       hue=1)
        # TODO: will this error if too many exons?
        colors_frame_1 = sns.color_palette(
            palette="husl", n_colors=len(merged_exon_bounds)
        )

        colors_frame_2 = sns.cubehelix_palette(
            len(merged_exon_bounds), start=0.8, rot=0.9, light=0.7, dark=0.3, hue=1
        )
        colors_frame_3 = sns.cubehelix_palette(
            len(merged_exon_bounds), start=2.9, rot=0.9, light=0.7, dark=0.3, hue=1
        )

        def _pick_representative_cat(s):
            if s == "":
                return ""
            # reverse so returns f in case of e.g. s = 'Mf'
            return collections.Counter(s[::-1]).most_common(1)[0][0]

        exon_colors = {}
        for iso in self.isoforms:
            algn = self.pairwise_changes_relative_to_reference(
                self.isoforms[0].name, iso.name
            )
            algn = algn.replace("D", "")
            split_algn = [
                algn[
                    sum(len(e) for e in iso.exons[:i])
                    // 3 : sum(len(e) for e in iso.exons[: i + 1])
                    // 3
                ]
                for i in range(len(iso.exons))
            ]
            change_cats = [_pick_representative_cat(s) for s in split_algn]
            for exon, cat in zip(iso.exons, change_cats):
                merge_exon_idx = exon.exon_number_on_gene - 1
                if cat == "f":
                    exon_colors[(iso.name, exon.start, exon.end)] = colors_frame_3[
                        merge_exon_idx
                    ]
                elif cat == "F":
                    exon_colors[(iso.name, exon.start, exon.end)] = colors_frame_2[
                        merge_exon_idx
                    ]
                else:
                    exon_colors[(iso.name, exon.start, exon.end)] = colors_frame_1[
                        merge_exon_idx
                    ]
        return exon_colors

    def exon_diagram(
        self,
        ax=None,
        intron_nt_space=30,
        height=0.5,
        show_uncloned_isoforms=True,
        show_domains=False,
        show_matched_transcripts=True,
        show_mane_and_appris_annotations=True,
        show_pathogenic_variants=False,
        domain_font_size=6,
        subtle_splice_font_size=6,
        subtle_splice_threshold=20,
    ):
        if ax is None:
            ax = plt.gca()

        isoforms = self.isoforms if show_uncloned_isoforms else self.cloned_isoforms
        expanded_exon_bounds = [(exon.start, exon.end) for exon in self.exons]
        ref_exon_bounds = []
        for start, end in expanded_exon_bounds:
            for exon in [exon for iso in isoforms for exon in iso.exons]:
                if exon.start >= start and exon.end <= end:
                    ref_exon_bounds.append((exon.start, exon.end))
                    break
            else:
                ref_exon_bounds.append((start, end))
        diff_exon_ends = {}
        for iso in isoforms:
            for exon in iso.exons:
                merged_start, merged_end = ref_exon_bounds[exon.exon_number_on_gene - 1]
                if exon.start != merged_start:
                    diff_exon_ends[exon.start] = merged_start - exon.start
                if exon.end != merged_end:
                    diff_exon_ends[exon.end] = exon.end - merged_end
        if self.strand == "-":
            expanded_exon_bounds = expanded_exon_bounds[::-1]
        mapped_exon_bounds = [expanded_exon_bounds[0]]
        for i in range(1, len(expanded_exon_bounds)):
            a = mapped_exon_bounds[i - 1][1] + intron_nt_space
            b = a + (expanded_exon_bounds[i][1] - expanded_exon_bounds[i][0])
            mapped_exon_bounds.append((a, b))

        def _map_position(
            pos, bounds_in=expanded_exon_bounds, bounds_out=mapped_exon_bounds
        ):
            if len(bounds_in) != len(bounds_out):
                raise ValueError("Invalid boundary mapping")
            for (start_in, stop_in), (start_out, stop_out) in zip(
                bounds_in, bounds_out
            ):
                if pos >= start_in and pos < stop_in:
                    return start_out + (pos - start_in)
            msg = "can't map position outside of exon boundaries\n"
            msg += "position: {}\nboundaries: {}\n".format(pos, bounds_in)
            raise ValueError(msg)

        exon_colors = self._get_exon_colors()
        xmin = _map_position(expanded_exon_bounds[0][0])
        xmax = _map_position(expanded_exon_bounds[-1][1] - 1)
        for i, iso in enumerate(isoforms):
            for exon in iso.exons:
                x_start = _map_position(exon.start)
                x_stop = _map_position(exon.end - 1)
                box = patches.Rectangle(
                    [x_start, i],
                    x_stop - x_start,
                    height,
                    lw=1,
                    ec="k",
                    fc=exon_colors[(iso.name, exon.start, exon.end)],
                    joinstyle="round",
                )
                ax.add_patch(box)

                # draw number of NT for small exon boundary changes
                if exon.start in diff_exon_ends:
                    num_nt_diff = diff_exon_ends[exon.start]
                    if abs(num_nt_diff) <= subtle_splice_threshold:
                        ax.text(
                            _map_position(exon.start),
                            i + height + 0.03,
                            "{}{} nt".format(
                                "+" if num_nt_diff > 0 else "Δ", abs(num_nt_diff)
                            ),
                            ha="left",
                            va="top",
                            fontsize=subtle_splice_font_size,
                        )
                if exon.end in diff_exon_ends:
                    num_nt_diff = diff_exon_ends[exon.end]
                    if abs(num_nt_diff) <= subtle_splice_threshold:
                        ax.text(
                            _map_position(exon.end - 1),
                            i + height + 0.03,
                            "{}{} nt".format(
                                "+" if num_nt_diff > 0 else "Δ", abs(num_nt_diff)
                            ),
                            ha="left",
                            va="top",
                            fontsize=subtle_splice_font_size,
                        )

            if not show_domains:
                continue
            for dom in iso.aa_seq_features:
                if dom.name.endswith("_DBD_flank"):
                    continue
                dom_x_start = _map_position(iso.residues[dom.start].coords[0])
                dom_x_stop = _map_position(iso.residues[dom.end - 1].coords[2])
                dom_x_center = (dom_x_stop - dom_x_start) / 2 + dom_x_start
                dom_x_len = abs(dom_x_stop - dom_x_start)
                if i != 0:  # TODO: make this an option argument?
                    break
                ax_x_range = abs(xmax - xmin)
                ax_y_range = len(isoforms) + 1
                # the 0.3 below is just tuned by hand
                y_height = (ax_x_range * 0.3) / ax_y_range
                ax.annotate(
                    "",
                    xy=(dom_x_start, i),
                    xycoords="data",
                    xytext=(dom_x_stop, i),
                    textcoords="data",
                    arrowprops=dict(
                        arrowstyle="<->, head_length=0, head_width=0.01",
                        connectionstyle="bar, fraction={}".format(
                            y_height / (dom_x_len)
                        ),
                        color="darkorange",
                    ),
                )
                ax.text(
                    dom_x_center,
                    i - (1 - height),
                    dom.name,
                    rotation=90,
                    va="bottom",
                    ha="center",
                    fontsize=domain_font_size,
                )
        ax.set_yticks([y + height / 2 for y in range(len(isoforms))])
        ax.set_yticklabels([iso.name for iso in isoforms])
        ax.yaxis.set_tick_params(length=0)
        if show_matched_transcripts:
            ax.set_yticklabels(
                [
                    iso.clone_name if hasattr(iso, "clone_name") else ""
                    for iso in isoforms
                ]
            )
            for iso, y_pos in zip(isoforms, ax.get_yticks()):
                if hasattr(iso, "clone_acc") and iso.is_novel_isoform():
                    text = "     " + "Novel isoform"
                else:
                    text = "     " + "/".join(iso.ensembl_transcript_names)
                    if show_mane_and_appris_annotations:
                        if (
                            hasattr(iso, "is_MANE_select_transcript")
                            and iso.is_MANE_select_transcript
                        ):
                            text += " – MANE select"
                        if (
                            hasattr(iso, "APPRIS_annotation")
                            and iso.APPRIS_annotation is not None
                        ):
                            text += " – APPRIS " + iso.APPRIS_annotation
                ax.text(
                    x=xmin if self.strand == "-" else xmax,
                    y=y_pos,
                    s=text,
                    ha="left",
                    va="center",
                )
        ax.set_xticks([])
        x_pad = intron_nt_space * 3
        # intron dotted lines
        plt.hlines(
            y=[i + height / 2 for i in range(len(isoforms))],
            xmin=[_map_position(iso.start) for iso in isoforms],
            xmax=[_map_position(iso.end - 1) for iso in isoforms],
            color="black",
            ls="dotted",
            lw=1.5,
            zorder=0,
        )

        if show_pathogenic_variants:
            muts = list(
                sorted(
                    self.pathogenic_coding_SNPs, key=lambda x: int(x.aa_change[3:-3])
                )
            )
            for idx_mut, mut in enumerate(muts):
                # color by disease
                # add link
                # sort diseases?
                # change color or bold if it's early stop
                if mut.aa_change.endswith("Ter"):
                    ax.scatter(
                        x=[_map_position(mut.start)], y=[-0.2], color="red", marker="H"
                    )
                elif mut.aa_change.startswith("Ter"):
                    if self.strand == "-":
                        stop_codon_pos = _map_position(mut.start + 3) - 1
                    else:
                        stop_codon_pos = _map_position(mut.start - 3) + 1
                    ax.scatter(
                        x=[stop_codon_pos],
                        y=[-0.2],
                        color="green",
                        marker="H",
                    )
                else:
                    ax.vlines(
                        x=_map_position(mut.start),
                        ymin=-0.1,
                        ymax=-0.5,
                        linewidth=1,
                        color="black",
                    )
            diseases = {mut.disease for mut in muts}
            for idx_disease, disease in enumerate(diseases):
                y_pos = 1 + (idx_disease / 20)
                ax.text(
                    s=disease,
                    y=y_pos,
                    x=-1 / 50,
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                )
                for idx_mut, mut in enumerate(
                    filter(lambda x: x.disease == disease, muts)
                ):
                    url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}".format(
                        mut.mutation_id
                    )
                    ax.text(
                        s=mut.aa_change,
                        y=y_pos,
                        x=idx_mut / 6,
                        transform=ax.transAxes,
                        va="bottom",
                        color="red" if mut.aa_change.endswith("Ter") else "black",
                        url=url,
                        bbox=dict(
                            url=url,
                            color="w",
                            alpha=0.01,
                        ),
                    )

        if self.strand == "-":
            ax.set_xlim(xmax + x_pad, xmin - x_pad)
        else:
            ax.set_xlim(xmin - intron_nt_space * 3, xmax + intron_nt_space * 3)
        ax.set_ylim(len(isoforms), height - 1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    def protein_diagram(
        self,
        ax=None,
        only_cloned_isoforms=True,
        isoform_order=None,
        domain_label_rotation=0,
        draw_legend=True,
        # remove_overlapping_domains=False  # I broke this by changing the domains to be tuples including the isoform index
        draw_vertical_lines_at_domains=False,
    ):
        remove_overlapping_domains = False
        if isoform_order is not None:
            isoforms = [self[iso_id] for iso_id in isoform_order]
        elif only_cloned_isoforms:
            isoforms = [iso for iso in self.isoforms if hasattr(iso, "clone_acc")]
        else:
            isoforms = self.isoforms

        def _remove_overlapping_domains(domains):
            """Keep longest domains"""
            non_overlapping = []
            for dom in sorted(domains, key=lambda x: len(x), reverse=True):
                if any(
                    dom.start < d.end and dom.end > d.start for d in non_overlapping
                ):
                    continue
                non_overlapping.append(dom)
            return non_overlapping

        if ax is None:
            ax = plt.gca()
        gs = gridspec.GridSpecFromSubplotSpec(
            len(isoforms), 1, subplot_spec=ax.get_subplotspec(), hspace=2
        )

        max_seq_len = max(len(iso.aa_seq) for iso in isoforms)

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -99999
        aligner.extend_gap_score = -1
        aligner.target_end_gap_score = 0
        aligner.query_end_gap_score = 0

        def _get_offset(algn):
            t = algn.__str__().splitlines()[0]
            q = algn.__str__().splitlines()[2]
            target_end_gap = len(t.partition(t.strip("-"))[0])
            query_end_gap = len(q.partition(q.strip("-"))[0])
            return query_end_gap - target_end_gap

        aa_seqs = [x.aa_seq for x in isoforms]
        offsets = [
            0,
        ]
        for i in range(1, len(aa_seqs)):
            alignments = [aligner.align(aa_seqs[j], aa_seqs[i])[0] for j in range(i)]
            i_best_alignment = max(
                range(len(alignments)), key=lambda x: alignments[x].score
            )
            offsets.append(
                _get_offset(alignments[i_best_alignment]) + offsets[i_best_alignment]
            )
        x_max = max(len(iso.aa_seq) + x for iso, x in zip(isoforms, offsets))

        exon_colors = self._get_exon_colors()

        # TODO: draw domains on further down isoforms if they're not on the top
        # (i.e. same accession and overlap in genomic coords)
        domains_to_draw = []
        for iso_idx, iso in enumerate(isoforms):
            for d in iso.aa_seq_features:
                if d.accession.endswith("_flank_N") or d.accession.endswith("_flank_C"):
                    continue  # don't show DBD flank regions
                if iso_idx == 0:
                    domains_to_draw.append((iso_idx, d))
                elif d.category == "effector_domain" and d.accession not in {
                    x[1].accession for x in domains_to_draw
                }:
                    domains_to_draw.append(
                        (iso_idx, d)
                    )  # only show domains for top isoform unless domain not already displayed

        if remove_overlapping_domains:  # BUG: THIS NO LONGER WORKS
            print("WARNING remove_overlapping_domains doensn't work anymore")
            domains_to_draw = _remove_overlapping_domains(domains_to_draw)
        for i, isoform in enumerate(isoforms):
            ax = plt.subplot(gs.new_subplotspec((i, 0), rowspan=1, colspan=1))
            seq_len = len(isoform.aa_seq)
            height = 1
            prot = patches.Rectangle(
                (0.5 + offsets[i], 0.5 - height / 2),
                width=seq_len,
                height=height,
                clip_on=False,
                facecolor=None,
                edgecolor="grey",
                linewidth=1.5,
            )
            ax.add_patch(prot)

            exon_pos = 0
            for exon in isoform.exons:
                n_aa_exon = (exon.end - exon.start) / 3
                n_aa_exon = min(n_aa_exon, seq_len - exon_pos)
                exon = patches.Rectangle(
                    (0.5 + offsets[i] + exon_pos, 0.5 - height / 2),
                    width=n_aa_exon,
                    height=height,
                    clip_on=False,
                    facecolor=exon_colors[(isoform.name, exon.start, exon.end)],
                    edgecolor=None,
                    linewidth=0,
                )
                ax.add_patch(exon)
                exon_pos += n_aa_exon

            for iso_idx, domain in domains_to_draw:
                start, stop = domain.start, domain.end
                if domain.category == "Pfam_domain":
                    if domain in isoform.dna_binding_domains:
                        dom_color = "blue"
                    else:
                        dom_color = "purple"
                elif domain.category == "effector_domain":
                    if domain.name == "AD":
                        dom_color = "green"
                    elif domain.name == "RD":
                        dom_color = "red"
                    else:
                        dom_color = "orange"
                else:
                    dom_color = "purple"
                if draw_vertical_lines_at_domains:
                    ax.axvline(
                        x=(start + 0.5) + offsets[iso_idx],
                        ymin=-1,
                        ymax=2,
                        clip_on=False,
                        linestyle="--",
                        linewidth=2,
                        color=dom_color,
                    )
                    ax.axvline(
                        x=(stop + 0.5) + offsets[iso_idx],
                        ymin=-1,
                        ymax=2,
                        clip_on=False,
                        linestyle="--",
                        linewidth=1.5,
                        color=dom_color,
                    )
                if i == iso_idx:
                    if domain.category == "effector_domain":
                        y_pos = (0.5 - height / 2) - height * 0.3
                        va = "top"
                    else:
                        y_pos = (0.5 + height / 2) + height * 0.3
                        va = "bottom"
                    ax.text(
                        x=(start + (stop - start) / 2) + offsets[i],
                        y=y_pos,
                        s=domain.name,
                        url="http://pfam.xfam.org/family/" + domain.accession,
                        ha="center",
                        va=va,
                        color=dom_color,
                        fontweight="bold",
                        fontsize=7,
                        rotation=domain_label_rotation,
                        bbox=dict(
                            url="https://pfam.xfam.org/family/" + domain.accession,
                            color="w",
                            alpha=0.01,
                        ),
                    )
                    x_ax_len = x_max - min(offsets)
                    if domain.category == "effector_domain":
                        y_pos = (0.5 - height / 2) - height * 0.2
                    else:
                        y_pos = (0.5 + height / 2) + height * 0.2
                    ax.axhline(
                        xmin=(start + (offsets[i] - min(offsets))) / x_ax_len,
                        xmax=(stop + (offsets[i] - min(offsets))) / x_ax_len,
                        y=y_pos,
                        clip_on=False,
                        linestyle="-",
                        linewidth=2.5,
                        color=dom_color,
                        alpha=0.8,
                    )

            ax.text(
                x=-0.5 + min(offsets),
                y=0.5,
                s=isoform.name,
                horizontalalignment="right",
                verticalalignment="center",
            )
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["bottom"].set_position(("data", 0.5 - height / 2))
            # domain_pos = [(d.start, d.end) for d in domains_to_draw]
            # xticks = [1] + [i for ij in domain_pos for i in ij] + [seq_len]
            xticks = [1, seq_len]
            ax.set_xticks([x + offsets[i] for x in xticks])
            ax.set_xticklabels([str(x) for x in xticks])
            ax.xaxis.set_tick_params(rotation=90)
            ax.set_yticks([])

            def _extract_pmid(s):
                if "PMID: " not in s:
                    return []
                return [
                    l[len("PMID: ") :] for l in s.splitlines() if l.startswith("PMID: ")
                ][0].split(", ")

            if draw_legend and i == 0:
                unique_domains = {
                    (d.name, d.description): d for _i, d in domains_to_draw
                }.values()
                i = 0
                for d in unique_domains:
                    if d.category == "Pfam_domain":
                        url = "https://pfam.xfam.org/family/" + d.accession
                    elif d.category == "effector_domain":
                        pmids = _extract_pmid(d.description)
                        if len(pmids) == 1:
                            url = "https://pubmed.ncbi.nlm.nih.gov/" + pmids[0]
                        elif len(pmids) > 1:
                            url = "https://pubmed.ncbi.nlm.nih.gov/?term=" + "+".join(
                                pmids
                            )
                    else:
                        url = ""
                    ax.text(
                        x=x_max + 0.5 + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.05,
                        y=ax.get_ylim()[1] - i * 0.5,
                        s=(
                            d.name + " – " + d.description
                            if d.description is not None
                            else d.name
                        ),
                        va="top",
                        ha="left",
                        fontsize=7,
                        url=url,
                        bbox=dict(
                            url=url, color="w", alpha=0.01
                        ),  # if you set alpha to 0 then you're not able to click the link
                    )
                    increment = (
                        len(d.description.splitlines()) + 1
                        if d.description is not None
                        else 2
                    )
                    i += increment
            ax.set_ylim(0, 1)
            ax.set_xlim(0.5 + min(offsets), x_max + 0.5)

    def iso_pairs(self):
        return list(itertools.combinations(self.isoforms, 2))

    def __getitem__(self, iso_id):
        if iso_id in self._iso_dict:
            return self._iso_dict[iso_id]
        for iso in self.isoforms:
            if hasattr(iso, "clone_acc"):
                if iso_id == iso.clone_acc:
                    return iso
            if iso.ensembl_transcript_names is not None:
                if iso_id in iso.ensembl_transcript_names:
                    return iso
                if iso_id in iso.ensembl_transcript_ids:
                    return iso
                if iso_id in iso.ensembl_protein_ids:
                    return iso
        raise KeyError(iso_id)

    def __contains__(self, iso_id):
        try:
            self[iso_id]
            return True
        except KeyError:
            return False

    def __repr__(self):
        s = "Gene: {}\n".format(self.name)
        s += "Isoforms: " + str([iso.name for iso in self.isoforms])
        return s


class Isoform(GenomicFeature):
    """Protein coding isoform of a gene

    Attributes:
        name (str): isoform name
        exons (list(isolib.Exon)): exons in isoform
        nt_seq (str): nucleotide sequence of CDS, excluding stop codon

    """

    def __init__(
        self,
        name,
        exons,
        CDS_nt_seq,
        aa_seq=None,
        UTR_5prime_nt_seq=None,
        UTR_3prime_nt_seq=None,
        ensembl_protein_ids=None,
        ensembl_transcript_ids=None,
        ensembl_transcript_names=None,
    ):
        self.name = name
        self.ensembl_protein_ids = ensembl_protein_ids
        self.ensembl_transcript_names = ensembl_transcript_names
        self.ensembl_transcript_ids = ensembl_transcript_ids

        if isinstance(CDS_nt_seq, str):
            if len(CDS_nt_seq) % 3 != 0:
                raise ValueError("CDS sequence length must be multiple of 3")
            stop_codons = {"TAG", "TAA", "TGA"}
            codons = [CDS_nt_seq[i : i + 3] for i in range(0, len(CDS_nt_seq), 3)]
            if any(codon in stop_codons for codon in codons[:-1]):
                raise ValueError("CDS sequnce must not contain early stop codon")
            self.CDS_nt_seq = CDS_nt_seq
            if CDS_nt_seq[-3:] in stop_codons:
                self.nt_seq = CDS_nt_seq[:-3]
            else:
                self.nt_seq = CDS_nt_seq
            if aa_seq is not None:
                self.aa_seq = aa_seq
            else:
                self.aa_seq = str(Seq(self.nt_seq).translate(to_stop=True))
            self.codons = [
                self.nt_seq[i : i + 3] for i in range(0, len(self.aa_seq) * 3, 3)
            ]
        else:
            raise NotImplementedError()
        if len(self.nt_seq) != len(self.aa_seq) * 3:
            raise UserWarning("Inconsistent DNA and AA sequences – " + self.name)

        if len(exons) == 0:
            raise ValueError(
                self.name + " - Need at least one exon to define an isoform"
            )
        chroms = [exon.chrom for exon in exons]
        strands = [exon.strand for exon in exons]
        if len(set(chroms)) > 1:
            raise ValueError(self.name + " - All exons must be on same chromosome")
        if len(set(strands)) > 1:
            raise ValueError(self.name + " - All exons must be same strand")
        len_all_exons = sum(len(e) for e in exons)
        if len_all_exons != len(self.nt_seq):
            msg = "exons not same length as CDS sequence minus stop\n"
            msg += "exons = {} nt, CDS = {} nt\n".format(
                len_all_exons, len(self.nt_seq)
            )
            msg += self.name
            raise UserWarning(msg)
        is_neg_strand = strands[0] == "-"
        self.exons = sorted(exons, key=lambda x: x.start, reverse=is_neg_strand)

        residues = []
        genomic_coords = [i for exon in self.exons for i in exon.genomic_coords()]
        codon_genomic_coords = [
            (a, b, c)
            for a, b, c in zip(
                genomic_coords[0::3], genomic_coords[1::3], genomic_coords[2::3]
            )
        ]
        for aa, codon, codon_coords in zip(
            self.aa_seq, self.codons, codon_genomic_coords
        ):
            residues.append(Residue(aa, codon, chroms[0], strands[0], codon_coords))
        self.residues = residues
        self.UTR_5prime_nt_seq = UTR_5prime_nt_seq
        self.UTR_3prime_nt_seq = UTR_3prime_nt_seq
        self.aa_seq_features = []
        GenomicFeature.__init__(
            self,
            chroms[0],
            strands[0],
            min([exon.start for exon in self.exons]),
            max([exon.end for exon in self.exons]),
        )

    def add_aa_seq_feature(
        self, category, accession, name, start, end, description=None
    ):
        if end > len(self.aa_seq):
            raise ValueError("Feature bounds outside protein")
        self.aa_seq_features.append(
            ProteinSequenceFeature(category, accession, name, start, end, description)
        )

    def remove_aa_seq_feature(self, accession, start, end):
        idx = [
            i
            for i, f in enumerate(self.aa_seq_features)
            if f.accession == accession and f.start == start and f.end == end
        ]
        if len(idx) == 0:
            raise ValueError("Feature not present")
        if len(idx) > 1:
            raise UserWarning("Unexpected duplicate domains – " + self.name)
        del self.aa_seq_features[idx[0]]

    @property
    def dna_binding_domains(self):
        # hiding import to avoid circular imports
        from data_loading import load_dbd_accessions

        dbd_acc = load_dbd_accessions()
        return [dom for dom in self.aa_seq_features if dom.accession in dbd_acc]

    def __repr__(self):
        if self.ensembl_protein_ids is not None:
            ids = " / ".join(
                [
                    "|".join(self.ensembl_transcript_names),
                    "|".join(self.ensembl_protein_ids),
                    "|".join(self.ensembl_transcript_ids),
                ]
            )
        else:
            ids = self.name
        s = "Isoform: {}\nlength: {} aa".format(ids, len(self.aa_seq))
        return s


class Cloned_Isoform(Isoform):
    def __init__(
        self,
        clone_name,
        exons,
        clone_nt_seq,
        clone_acc=None,
        aa_seq=None,
        UTR_5prime_nt_seq=None,
        UTR_3prime_nt_seq=None,
        ensembl_transcript_ids=None,
        ensembl_transcript_names=None,
        ensembl_protein_ids=None,
    ):
        self.clone_name = clone_name
        self.clone_acc = clone_acc
        self.clone_nt_seq = clone_nt_seq

        def _trim_clone_nt_seq_to_CDS_and_strip_stop_codon(clone_seq):
            """
            The clones can have untranslated regions at the 3' end, for example
            if there's a inclusion of an exon containing a stop codon, in between
            the primers. There is also sometimes a stop codon at the end and sometimes
            not. Also not all sequences have length a multiple of 3. This function
            normalizes the nucleotide sequence to just be the coding sequence without
            the stop codon.
            """
            stop_codons = {"TAG", "TAA", "TGA"}
            codons = [clone_seq[i : i + 3] for i in range(0, len(clone_seq), 3)]
            if any(codon in stop_codons for codon in codons):
                return "".join(
                    codons[
                        : min(
                            codons.index(stop) for stop in stop_codons if stop in codons
                        )
                    ]
                )
            elif len(clone_seq) % 3 == 0:
                return clone_seq
            else:
                return clone_seq[: -(len(clone_seq) % 3)]

        cds_nt_seq = _trim_clone_nt_seq_to_CDS_and_strip_stop_codon(clone_nt_seq)
        is_neg_strand = exons[0].strand == "-"
        ordered_exons = sorted(exons, key=lambda x: x.start, reverse=is_neg_strand)
        if sum(len(e) for e in exons) != len(cds_nt_seq):
            # removing non-coding sequence from exons
            trimmed_exons = []
            cummulative_n_nt = 0
            for exon in ordered_exons:
                if cummulative_n_nt < len(cds_nt_seq):
                    trimmed_exons.append(exon)
                cummulative_n_nt += len(exon)
            if is_neg_strand:
                trimmed_exons[-1].start = trimmed_exons[-1].end - (
                    len(cds_nt_seq) - sum(len(e) for e in trimmed_exons[:-1])
                )
            else:
                trimmed_exons[-1].end = trimmed_exons[-1].start + (
                    len(cds_nt_seq) - sum(len(e) for e in trimmed_exons[:-1])
                )
            if sum(len(e) for e in trimmed_exons) != len(cds_nt_seq):
                raise UserWarning(
                    self.name, " - problem with reducing exons to coding-only"
                )
        else:
            trimmed_exons = exons

        Isoform.__init__(
            self,
            name=clone_name,
            exons=trimmed_exons,
            CDS_nt_seq=cds_nt_seq,
            aa_seq=aa_seq,
            UTR_5prime_nt_seq=UTR_5prime_nt_seq,
            UTR_3prime_nt_seq=UTR_3prime_nt_seq,
            ensembl_protein_ids=ensembl_protein_ids,
            ensembl_transcript_ids=ensembl_transcript_ids,
            ensembl_transcript_names=ensembl_transcript_ids,
        )

    def is_novel_isoform(self):
        return self.ensembl_transcript_ids is None

    def aa_seq_differences_from_GENCODE(self):
        if self.is_novel_isoform():
            raise UserWarning("Novel isoform doesn't match to GENCODE")
        if self.aa_seq == self.aa_seq_GENCODE:
            return None
        msg = []
        for i, (aa_clone, aa_gc) in enumerate(zip(self.aa_seq, self.aa_seq_GENCODE)):
            if aa_clone != aa_gc:
                msg.append(aa_gc + str(i + 1) + aa_clone)
        return ", ".join(msg)

    def nt_seq_differences_from_GENCODE_CDS(self):
        if self.is_novel_isoform():
            raise UserWarning("Novel isoform doesn't match to GENCODE")
        if self.clone_nt_seq == self.nt_seq_CDS_GENCODE:
            return None
        msg = []
        for i, (nt_clone, nt_gc) in enumerate(
            zip(self.clone_nt_seq, self.nt_seq_CDS_GENCODE)
        ):
            if nt_clone != nt_gc:
                msg.append(str(i + 1) + nt_gc + ">" + nt_clone)
        return ", ".join(msg)

    def __repr__(self):
        s = "Clone acc: {}\n".format(self.clone_acc)
        s += Isoform.__repr__(self)
        return s


class Exon(GenomicFeature):
    """Exon within a protein coding isoform transcript

    Attributes:
        gene_id / transcript_id (str): associated gene and transcript
        chrom (str): chromosome identifier
        strand (str): must be '+'/'-'
        start/end (int): chromosomal coordinates, follows python indexing
            convention of 0-indexed half-open interval, start must be before end

    """

    def __init__(self, gene_id, transcript_id, chrom, strand, start, end):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        GenomicFeature.__init__(self, chrom, strand, start, end)

    @property
    def is_sym(self):
        """Determine if exon is symmetrical in terms of translation frame."""
        return True if self.len % 3 == 0 else False

    def genomic_coords(self):
        """TODO: maybe needs a better name"""
        if self.strand == "+":
            return range(self.start, self.end)
        elif self.strand == "-":
            return reversed(range(self.start, self.end))
        else:
            raise UserWarning("Strand not set correctly (+/-)")

    def __repr__(self):
        s = "Exon of {} {} nt\n".format(self.transcript_id, len(self))
        return s


class Residue:
    """An amino acid, with its position on the genome

    Attributes:
        aa (str): single-letter amino acid code
        codon (str): nucleotide triplet
        chrom (str): chromosome identifier
        strand (str): must be '+'/'-'
        positions ((int, int, int)): genomic coordinates of each of three nucleotides

    """

    def __init__(self, aa, codon, chrom, strand, coords):
        # TODO: check valid AA and codon
        self.aa = aa
        self.codon = codon
        self.chrom = chrom
        if strand not in {"+", "-"}:
            msg = strand + " is invalid value for strand (must be +/-)"
            raise ValueError(msg)
        self.strand = strand
        if len(coords) != 3 or len(set(coords)) != 3:
            raise ValueError("Invalid genomic coordinates for residue")
        self.coords = coords
