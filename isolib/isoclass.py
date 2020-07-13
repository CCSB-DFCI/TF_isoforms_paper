# title           :isoclass.py
# description     :Classes representing isoform-related objects.
# author          :Gloria Sheynkman
# date            :May 1st, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================

import itertools
import random

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from Bio.Seq import Seq


class GenomicFeature:
    """A contigous region on the genome

    Attributes:
        chrom (str): chromosome identifier
        strand (str): must be '+'/'-'
        start/end (int): chromosomal coordinates, follows python indexing
            convention of 0-indexed half-open interval, start must be before end

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


class ProteinSequenceFeature:
    """A contiguous stretch of amino acids within a protein.

        category (str): general class of feature, e.g. 'Pfam_domain', 'ELM_motif'
        accession (str): identifier for feature, e.g. 'PF00170'
        name (str): display name, e.g. 'bZIP'
        start/end (int): position within the protein sequence, follows python indexing
            convention of 0-indexed half-open interval, start must be before end

    """

    def __init__(self, category, accession, name, start, end):
        if start >= end:
            raise ValueError("start must be before end")
        self.start = start
        self.end = end
        self.category = category
        self.name = name
        self.accession = accession

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
        orfs (list(isolib.ORF)): protein coding isoforms of gene

    """

    def __init__(self, name, orfs):
        chroms = [orf.chrom for orf in orfs]
        strands = [orf.strand for orf in orfs]
        msg = '{} - {}'.format(name, ', '.join([orf.name for orf in orfs]))
        if len(orfs) == 0:
            raise ValueError("Need at least one ORF to define a gene\n" + msg)
        if len(set(chroms)) > 1:
            raise ValueError("All ORFs must be on same chromosome\n" + msg)
        if len(set(strands)) > 1:
            raise ValueError("All ORFs must be same strand\n" + msg)
        # de-duplicate exons
        gene_exons = {}
        for orf in orfs:
            for exon in orf.exons:
                gene_exons[(exon.start, exon.end)] = exon
        self.exons = sorted(
            gene_exons.values(), key=lambda x: x.start, reverse=(strands[0] == "-")
        )
        self.number_of_isoforms = len(orfs)
        self.name = name
        self.orfs = list(
            sorted(orfs, key=lambda x: int(x.name.split("-")[-1]))
        )
        self._orf_dict = {orf.name: orf for orf in self.orfs}
        self._pairwise_changes = {}
        GenomicFeature.__init__(
            self,
            chroms[0],
            strands[0],
            min([orf.start for orf in orfs]),
            max([orf.end for orf in orfs]),
        )

    def genomic_alignment_of_aa_seqs(self, subset=None):
        """genomic co-ordinates of translated regions"""
        all_orf_names = [orf.name for orf in self.orfs]
        if subset is None:
            subset = all_orf_names
        orfs = [orf for orf in self.orfs if orf.name in subset]
        if len(orfs) != len(subset):
            msg = "Missing ORFs: "
            msg += "/".join([s for s in subset if s not in all_orf_names])
            raise ValueError(msg)
        gene_coords = sorted(
            list(set([res.coords[1] for orf in orfs for res in orf.residues]))
        )
        tracks = {}
        for orf in sorted(orfs, key=lambda x: x.name):
            if orf.name not in subset:
                continue
            orf_aa = {res.coords[1]: res.aa for res in orf.residues}
            tracks[orf.name] = "".join([orf_aa.get(i, "-") for i in gene_coords])
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
        alignment = ''
        ref_iter = iter(self._orf_dict[ref_iso_name].residues)
        alt_iter = iter(self._orf_dict[alt_iso_name].residues)
        ref_res = next(ref_iter)
        alt_res = next(alt_iter)
        while True:
            if not any(i in ref_res.coords for i in alt_res.coords):
                if ((self.strand == '+' and ref_res.coords[2] < alt_res.coords[2])
                   or (self.strand == '-' and ref_res.coords[2] > alt_res.coords[2])):
                    alignment += 'D'
                    try:
                        ref_res = next(ref_iter)
                    except StopIteration:
                        alignment += 'I'  # since alt_iter was incremented
                        for _remaining in alt_iter:
                            alignment += 'I'
                        break
                else:
                    alignment += 'I'
                    try:
                        alt_res = next(alt_iter)
                    except StopIteration:
                        alignment += 'D'  # since ref_iter was incremented
                        for _remaining in ref_iter:
                            alignment += 'D'
                        break
            else:
                if ref_res.coords == alt_res.coords:
                    alignment += 'M'
                elif any(i == j for i, j in zip(ref_res.coords, alt_res.coords)):
                    # Different exon junctions
                    if ref_res.aa == alt_res.aa:
                        alignment += 'M'
                    else:
                        if ((ref_res.coords[0] < alt_res.coords[0] and self.strand == '+')
                           or (ref_res.coords[2] > alt_res.coords[2] and self.strand == '+')
                           or (ref_res.coords[0] > alt_res.coords[0] and self.strand == '-')
                           or (ref_res.coords[2] < alt_res.coords[2] and self.strand == '-')):
                            alignment += 'I'
                            try:
                                alt_res = next(alt_iter)
                                continue
                            except StopIteration:
                                alignment += 'D'  # since ref_iter was incremented
                                for _remaining in ref_iter:
                                    alignment += 'D'
                                break
                        elif ((ref_res.coords[0] > alt_res.coords[0] and self.strand == '+')
                              or (ref_res.coords[2] < alt_res.coords[2] and self.strand == '+')
                              or (ref_res.coords[0] < alt_res.coords[0] and self.strand == '-')
                              or (ref_res.coords[2] > alt_res.coords[2] and self.strand == '-')):
                            alignment += 'D'
                            try:
                                ref_res = next(ref_iter)
                                continue
                            except StopIteration:
                                alignment += 'I'  # since alt_iter was incremented
                                for _remaining in alt_iter:
                                    alignment += 'I'
                                break
                        else:  # Here just the middle nt in codon matched
                            msg = 'Unexpected alignement issue between: '
                            msg += ref_iso_name + ' and ' + alt_iso_name
                            raise UserWarning(msg)
                else:  # Frameshift
                    if any(ref_res.coords[i] == alt_res.coords[(i + 1) % 3] for i in range(3)):
                        alignment += 'F'
                    else:
                        alignment += 'f'

                try:
                    ref_res = next(ref_iter)
                except StopIteration:
                    for _remaining in alt_iter:
                        alignment += 'I'
                    break
                try:
                    alt_res = next(alt_iter)
                except StopIteration:
                    alignment += 'D'  # since we incremented above already
                    for _remaining in ref_iter:
                        alignment += 'D'
                    break
        self._pairwise_changes[key] = alignment   # cache result
        return alignment

    def aa_seq_disruption(self, ref_iso_name, alt_iso_name, domain_start, domain_end):
        """Get pairwise alignment of orf protein sequences. Return fraction of
           domain and insertion
        """
        algn = self.pairwise_changes_relative_to_reference(ref_iso_name,
                                                           alt_iso_name)
        if len(algn.replace("I", "")) != len(self._orf_dict[ref_iso_name].aa_seq):
            msg = 'Something is wrong\n'
            msg += ref_iso_name + ', ' + alt_iso_name
            raise UserWarning(msg)

        def _coords_transform_aa_seq_to_alignment(i, alignment):
            if i > len(alignment.replace("I", "")):
                raise ValueError("position is not in ORF AA sequence")
            aa_seq_indices = [
                "" if c == "I" else len(alignment[:j].replace("I", ""))
                for j, c in enumerate(alignment)
            ]
            return aa_seq_indices.index(i)

        start = _coords_transform_aa_seq_to_alignment(domain_start, algn)
        end = _coords_transform_aa_seq_to_alignment(domain_end - 1, algn) + 1
        return {'deletion': algn[start:end].count('D'),
                'insertion': algn[start:end].count('I'),
                'frameshift': algn[start:end].count('F') + algn[start:end].count('f')}

    def aa_feature_disruption(self, ref_iso_name):
        """

        Args:
            ref_iso_name (str): [description]

        Returns:
            [type]: [description]
        """
        results = []
        ref_iso = self._orf_dict[ref_iso_name]
        row = {'gene': self.name, 'ref_iso': ref_iso_name}
        for aa_feature in ref_iso.aa_seq_features:
            for alt_iso_name, alt_iso in self._orf_dict.items():
                if alt_iso_name == ref_iso_name:
                    continue
                row.update({'alt_iso': alt_iso_name})
                row.update({'accession': aa_feature.accession})
                r = self.aa_seq_disruption(
                    ref_iso_name, alt_iso_name, aa_feature.start, aa_feature.end
                )
                row.update(r)
                row.update({'length': aa_feature.end - aa_feature.start})
                results.append(row.copy())
        results = pd.DataFrame(results)
        return results

    def null_fraction_per_aa_feature(self, ref_iso_name):
        """Fraction of aa features that would be affected if they were evenly
            distributed along the AA sequence of the reference isoform.

        Args:
            ref_iso_name (str)

        """
        results = []
        cache = {}
        ref_iso = self._orf_dict[ref_iso_name]
        row = {'gene': self.name, 'ref_iso': ref_iso_name}
        for aa_feature in ref_iso.aa_seq_features:
            for alt_iso_name, alt_iso in self._orf_dict.items():
                if alt_iso_name == ref_iso_name:
                    continue
                row.update({'alt_iso': alt_iso_name,
                            'accession': aa_feature.accession})
                aa_feature_length = aa_feature.end - aa_feature.start
                if (alt_iso_name, aa_feature_length) in cache:
                    row['null_fraction_affected'] = cache[(alt_iso_name, aa_feature_length)]
                    results.append(row.copy())
                    continue
                rs = self._null_feature_disruption(ref_iso_name, alt_iso_name, aa_feature_length)

                def is_disrupted(feature_alignment_count):
                    return any(v > 0 for v in feature_alignment_count.values())

                fraction_affected = sum([int(is_disrupted(r)) for r in rs]) / len(rs)
                row['null_fraction_affected'] = fraction_affected
                results.append(row.copy())
                cache[(alt_iso_name, aa_feature_length)] = fraction_affected
        results = pd.DataFrame(results)
        return results

    def _null_feature_disruption(self, ref_iso_name, alt_iso_name, feature_length):
        algn = self.pairwise_changes_relative_to_reference(ref_iso_name,
                                                           alt_iso_name)
        len_ref_iso_aa_seq = len(self._orf_dict[ref_iso_name].aa_seq)
        if len(algn.replace("I", "")) != len_ref_iso_aa_seq:
            msg = 'Something is wrong\n'
            msg += ref_iso_name + ', ' + alt_iso_name
            raise UserWarning(msg)
        coords_ref_iso_aa_seq = [
            "" if c == "I" else len(algn[:j].replace("I", ""))
            for j, c in enumerate(algn)
        ]
        two_mers = {'deletion': [], 'insertion': [], 'frameshift': []}
        for i in range(len_ref_iso_aa_seq - 1):
            start = coords_ref_iso_aa_seq.index(i)
            end = coords_ref_iso_aa_seq.index(i + 1) + 1
            two_mer = algn[start:end]
            two_mers['deletion'].append(two_mer.count('D'))
            two_mers['insertion'].append(two_mer.count('I'))
            two_mers['frameshift'].append(two_mer.count('F') + two_mer.count('f'))
        results = [{k: sum(v[:feature_length]) for k, v in two_mers.items()}]
        for i in range((len_ref_iso_aa_seq - feature_length) - 1):
            r = {k: (v - two_mers[k][i]) + two_mers[k][i + feature_length] for k, v in results[-1].items()}
            results.append(r)
        return results

    def random_aa_seq_feature_shuffle(self, ref_iso_name, n, subset=None):
        """Pick random, non-overlapping, aa sequence locations for the features.

        Args:
            ref_iso_name (str):
            n (int): number of randomizations
            subset (set(str), optional): [description]. Defaults to None.
        """
        if len(self._orf_dict) == 1:
            return pd.DataFrame([])

        # inspired by: https://stackoverflow.com/questions/18641272
        def _non_overlapping_random_feature_positons(len_aa_seq, feature_and_len):
            if sum([v[1] for v in feature_and_len]) > len_aa_seq:
                raise UserWarning('Impossible: ' + str(ref_iso_name))
            indices = range(len_aa_seq - sum(v[1] - 1 for v in feature_and_len))
            result = []
            offset = 0
            n = len(feature_and_len)
            for (acc, l), i in zip(random.sample(feature_and_len, n),
                                   sorted(random.sample(indices, n))):
                i += offset
                result.append((acc, (i, i + l)))
                offset += l - 1
            return result

        results = []
        ref_iso = self._orf_dict[ref_iso_name]
        row = {'gene': self.name, 'ref_iso': ref_iso_name}
        for i in range(n):
            row['random_sample'] = i
            rnd_pos = _non_overlapping_random_feature_positons(len(ref_iso.aa_seq),
                                                               [(f.accession, len(f)) for f in ref_iso.aa_seq_features if subset is None or f.accession in subset])
            for accession, (rnd_start, rnd_end) in rnd_pos:
                for alt_iso_name, alt_iso in self._orf_dict.items():
                    if alt_iso_name == ref_iso_name:
                        continue
                    row.update({'alt_iso': alt_iso_name})
                    row.update({'accession': accession})
                    r = self.aa_seq_disruption(
                        ref_iso_name, alt_iso_name, rnd_start, rnd_end
                    )
                    row.update(r)
                    row.update({'length': rnd_end - rnd_start})
                    results.append(row.copy())
        results = pd.DataFrame(results)
        return results

    def exon_diagram(self,
                     intron_nt_space=30,
                     height=0.5,
                     draw_domains=True,
                     ax=None,
                     domain_font_size=6):
        if ax is None:
            ax = plt.gca()
        exon_bounds = [(exon.start, exon.end) for exon in self.exons]
        if self.strand == "-":
            exon_bounds = exon_bounds[::-1]
        merged_exon_bounds = []
        for i in range(len(exon_bounds) - 1):
            # they're sorted above so start_a <= start_b
            start_a, end_a = exon_bounds[i]
            start_b, end_b = exon_bounds[i + 1]
            if end_a < start_b:
                merged_exon_bounds.append(exon_bounds[i])
            else:
                exon_bounds[i + 1] = (start_a, max(end_a, end_b))
        merged_exon_bounds.append(exon_bounds[-1])
        mapped_exon_bounds = [merged_exon_bounds[0]]
        for i in range(1, len(merged_exon_bounds)):
            a = mapped_exon_bounds[i - 1][1] + intron_nt_space
            b = a + (merged_exon_bounds[i][1] - merged_exon_bounds[i][0])
            mapped_exon_bounds.append((a, b))

        def _map_position(
            pos, bounds_in=merged_exon_bounds, bounds_out=mapped_exon_bounds
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

        xmin = _map_position(merged_exon_bounds[0][0])
        xmax = _map_position(merged_exon_bounds[-1][1] - 1)
        for i, orf in enumerate(self.orfs):
            for exon in orf.exons:
                x_start = _map_position(exon.start)
                x_stop = _map_position(exon.end - 1)
                box = mpatches.Rectangle(
                    [x_start, i],
                    x_stop - x_start,
                    height,
                    lw=1,
                    ec="k",
                    fc="mediumpurple",
                    joinstyle="round",
                )
                ax.add_patch(box)
            if not draw_domains:
                continue
            for dom in orf.aa_seq_features:
                dom_x_start = _map_position(orf.residues[dom.start].coords[0])
                dom_x_stop = _map_position(orf.residues[dom.end - 1].coords[2])
                dom_x_center = (dom_x_stop - dom_x_start) / 2 + dom_x_start
                dom_x_len = abs(dom_x_stop - dom_x_start)
                if i != 0:  # TODO: make this an option argument?
                    break
                ax_x_range = abs(xmax - xmin)
                ax_y_range = len(self.orfs) + 1
                # the 0.3 below is just tuned by hand
                y_height = (ax_x_range * 0.3) / ax_y_range
                ax.annotate('',
                            xy=(dom_x_start, i),
                            xycoords='data',
                            xytext=(dom_x_stop, i),
                            textcoords='data',
                            arrowprops=dict(arrowstyle="<->, head_length=0, head_width=0.01",
                                            connectionstyle="bar, fraction={}".format(y_height / (dom_x_len)),
                                            color='darkorange',))
                ax.text(dom_x_center,
                        i - (1 - height),
                        dom.name,
                        rotation=90,
                        va='bottom',
                        ha='center',
                        fontsize=domain_font_size
                        )
        ax.set_yticks([y + height / 2 for y in range(len(self.orfs))])
        ax.set_yticklabels([orf.name for orf in self.orfs])
        ax.yaxis.set_tick_params(length=0)
        ax.set_xticks([])
        x_pad = intron_nt_space * 3
        # intron dotted lines
        plt.hlines(
            y=[i + height / 2 for i in range(len(self.orfs))],
            xmin=[_map_position(orf.start) for orf in self.orfs],
            xmax=[_map_position(orf.end - 1) for orf in self.orfs],
            color="black",
            ls="dotted",
            lw=1.5,
            zorder=0,
        )
        if self.strand == "-":
            ax.set_xlim(xmax + x_pad, xmin - x_pad)
        else:
            ax.set_xlim(xmin - intron_nt_space * 3, xmax + intron_nt_space * 3)
        ax.set_ylim(len(self.orfs), height - 1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    def orf_pairs(self):
        return list(itertools.combinations(self.orfs, 2))

    def __getitem__(self, orf_id):
        if orf_id in self._orf_dict:
            return self._orf_dict[orf_id]
        for orf in self.orfs:
            if orf_id == orf.clone_acc:
                return orf
        raise KeyError()

    def __contains__(self, orf_id):
        return orf_id in self._orf_dict

    def __repr__(self):
        s = "Gene: {}\n".format(self.name)
        s += "Isoforms: " + str([orf.name for orf in self.orfs])
        return s


class ORF(GenomicFeature):
    """Protein coding isoform of a gene

    Attributes:
        name (str): isoform name
        exons (list(isolib.Exon)): exons in isoform
        nt_seq (str): nucleotide sequence

    """

    def __init__(
        self,
        name,
        exons,
        nt_seq,
        aa_seq=None,
        ensembl_transcript_id=None,
        ensembl_protein_id=None,
        clone_acc=None,
    ):
        self.name = name
        self.ensembl_transcript_id = ensembl_transcript_id
        self.ensembl_protein_id = ensembl_protein_id
        self.clone_acc = clone_acc
        if len(exons) == 0:
            raise ValueError(self.name + " - Need at least one exon to define a gene")
        chroms = [exon.chrom for exon in exons]
        strands = [exon.strand for exon in exons]
        if len(set(chroms)) > 1:
            raise ValueError(self.name + " - All exons must be on same chromosome")
        if len(set(strands)) > 1:
            raise ValueError(self.name + " - All exons must be same strand")
        is_neg_strand = strands[0] == "-"
        self.exons = sorted(exons, key=lambda x: x.start, reverse=is_neg_strand)
        if aa_seq is not None:
            self.aa_seq = aa_seq
        if isinstance(nt_seq, str):
            self.nt_seq = nt_seq
            if aa_seq is None:
                self.aa_seq = str(Seq(self.nt_seq).translate(to_stop=True))
            self.codons = [
                self.nt_seq[i:i + 3] for i in range(0, len(self.aa_seq) * 3, 3)
            ]
        else:
            raise NotImplementedError()
        residues = []
        genomic_coords = [i for exon in self.exons for i in exon.genomic_coords()]
        codon_genomic_coords = [(a, b, c) for a, b, c in zip(genomic_coords[0::3],
                                                             genomic_coords[1::3],
                                                             genomic_coords[2::3])
                                ]
        len_all_exons = sum(len(e) for e in self.exons)
        if ((len_all_exons != len(self.aa_seq) * 3)
           and (len_all_exons != len(self.nt_seq))):
            msg = """Genome alignment issues for {}\n
                     {} exons with cumulative length: {}\n
                     length aa seq: {}\n
                     length nt seq: {}\n
                  """
            msg = msg.format(self.name,
                             len(self.exons),
                             len_all_exons,
                             len(self.aa_seq),
                             len(self.nt_seq))
            raise UserWarning(msg)
        for aa, codon, codon_coords in zip(
            self.aa_seq, self.codons, codon_genomic_coords
        ):
            residues.append(Residue(aa, codon, chroms[0], strands[0], codon_coords))
        self.residues = residues
        self.aa_seq_features = []
        GenomicFeature.__init__(
            self,
            chroms[0],
            strands[0],
            min([exon.start for exon in exons]),
            max([exon.end for exon in exons]),
        )

    def add_aa_seq_feature(self, category, accession, name, start, end):
        if end > len(self.aa_seq):
            raise ValueError("Feature bounds outside protein")
        self.aa_seq_features.append(
            ProteinSequenceFeature(category, accession, name, start, end)
        )

    def remove_aa_seq_feature(self, accession, start, end):
        idx = [i for i, f in enumerate(self.aa_seq_features)
               if f.accession == accession
               and f.start == start
               and f.end == end]
        if len(idx) == 0:
            raise ValueError('Feature not present')
        if len(idx) > 1:
            raise UserWarning('Unexpected duplicate domains')
        del self.aa_seq_features[idx[0]]


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
