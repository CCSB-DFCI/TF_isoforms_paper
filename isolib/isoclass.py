# title           :isoclass.py
# description     :Classes representing isoform-related objects.
# author          :Gloria Sheynkman
# date            :May 1st, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================

import itertools

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
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
        s = '{}: {} {} {}-{}'.format(self.category,
                                 self.accession,
                                 self.name,
                                 self.start + 1,
                                 self.end)
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
        if len(orfs) == 0:
            raise ValueError("Need at least one ORF to define a gene")
        if len(set(chroms)) > 1:
            raise ValueError("All ORFs must be on same chromosome")
        if len(set(strands)) > 1:
            raise ValueError("All ORFs must be same strand")
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
        self.orfs = list(sorted(orfs,
                                key=lambda x: int(x.name.split('|')[1].split('/')[0])))
        self._orf_dict = {orf.name: orf for orf in self.orfs}
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
            msg = 'Missing ORFs: '
            msg += '/'.join([s for s in subset if s not in all_orf_names])
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

    def aa_seq_disruption(self, ref_iso_name, alt_iso_name, domain_start, domain_end):
        """Get pairwise alignment of orf protein sequences. Return fraction of
           domain and insertion
        """
        algn = self.genomic_alignment_of_aa_seqs(subset=[ref_iso_name, alt_iso_name])

        def _coords_transform_alignment_to_aa_seq(i, alignment):
            if alignment[i] == '-':
                raise ValueError('position is not in ORF AA sequence')
            return len(alignment[:i+1].replace('-', ''))

        def _coords_transform_aa_seq_to_aligment(i, alignment):
            if i > len(alignment.replace('-', '')):
                raise ValueError('position is not in ORF AA sequence')
            aa_seq_indices = ['' if c == '-' else len(alignment[:j].replace('-', '')) for j, c in enumerate(alignment)]
            return aa_seq_indices.index(i)

        start = _coords_transform_aa_seq_to_aligment(domain_start,
                                                     algn[ref_iso_name])
        end = _coords_transform_aa_seq_to_aligment(domain_end,
                                                   algn[ref_iso_name])
        n_aa_insert = algn[ref_iso_name][start:end].count('-')
        n_aa_del = algn[alt_iso_name][start:end].count('-')
        # todo: process frameshift, insertion / deletion / substitution
        return (n_aa_insert, n_aa_del)

    def aa_feature_disruption(self, ref_iso_name):
        results = {}
        ref_iso = self._orf_dict[ref_iso_name]
        for aa_feature in ref_iso.aa_seq_features:
            for alt_iso_name, alt_iso in self._orf_dict.items():
                if alt_iso_name == ref_iso_name:
                    continue
                r = self.aa_seq_disruption(ref_iso_name,
                                           alt_iso_name,
                                           aa_feature.start,
                                           aa_feature.end)
                results[(ref_iso_name, alt_iso_name, aa_feature.accession)] = r
        return results

    def exon_diagram(self, intron_nt_space=30, height=0.5, ax=None):
        if ax is None:
            ax = plt.gca()
        exon_bounds = [(exon.start, exon.end) for exon in self.exons]
        if self.strand == '-':
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

        def _map_position(pos,
                          bounds_in=merged_exon_bounds,
                          bounds_out=mapped_exon_bounds):
            if len(bounds_in) != len(bounds_out):
                raise ValueError('Invalid boundary mapping')
            for (start_in, stop_in), (start_out, stop_out) in zip(
                bounds_in, bounds_out
            ):
                if pos >= start_in and pos < stop_in:
                    return start_out + (pos - start_in)
            msg = "can't map position outside of exon boundaries\n"
            msg += 'position: {}\nboundaries: {}\n'.format(pos, bounds_in)
            raise ValueError(msg)

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
        ax.set_yticks([y + height / 2 for y in range(len(self.orfs))])
        ax.set_yticklabels([orf.name for orf in self.orfs])
        ax.yaxis.set_tick_params(length=0)
        ax.set_xticks([])
        xmin = _map_position(merged_exon_bounds[0][0])
        xmax = _map_position(merged_exon_bounds[-1][1] - 1)
        x_pad = intron_nt_space * 3
        # intron dotted lines
        plt.hlines(y=[i + height / 2 for i in range(len(self.orfs))],
                   xmin=[_map_position(orf.start) for orf in self.orfs],
                   xmax=[_map_position(orf.end - 1) for orf in self.orfs],
                   color='black',
                   ls="dotted",
                   lw=1.5,
                   zorder=0)
        if self.strand == '-':
            ax.set_xlim(xmax + x_pad,
                        xmin - x_pad)
        else:
            ax.set_xlim(xmin - intron_nt_space * 3,
                        xmax + intron_nt_space * 3)
        ax.set_ylim(len(self.orfs), height - 1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    def orf_pairs(self):
        return list(itertools.combinations(self.orfs, 2))

    def __getitem__(self, orf_id):
        return self._orf_dict[orf_id]

    def __contains__(self, orf_id):
        return orf_id in self._orf_dict

    def __repr__(self):
        s = 'Gene: {}\n'.format(self.name)
        s += 'Isoforms: ' + str([orf.name for orf in self.orfs])
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
        ensembl_transcipt_id=None,
        ensembl_protein_id=None,
        orf_id=None,
    ):
        self.name = name
        self.ensembl_transcipt_id = ensembl_transcipt_id
        self.ensembl_protein_id = ensembl_protein_id
        self.orf_id = orf_id
        if len(exons) == 0:
            raise ValueError("Need at least one exon to define a gene")
        chroms = [exon.chrom for exon in exons]
        strands = [exon.strand for exon in exons]
        if len(set(chroms)) > 1:
            raise ValueError("All exons must be on same chromosome")
        if len(set(strands)) > 1:
            raise ValueError("All exons must be same strand")
        is_neg_strand = strands[0] == "-"
        self.exons = sorted(exons, key=lambda x: x.start, reverse=is_neg_strand)
        if isinstance(nt_seq, str):
            self.nt_seq = nt_seq
            self.aa_seq = str(Seq(self.nt_seq).translate(to_stop=True))
            self.codons = [
                self.nt_seq[i: i + 3] for i in range(0, len(self.aa_seq) * 3, 3)
            ]
        else:
            raise NotImplementedError()
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
            raise ValueError('Feature bounds outside protein')
        self.aa_seq_features.append(ProteinSequenceFeature(category,
                                                           accession,
                                                           name,
                                                           start,
                                                           end))


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
