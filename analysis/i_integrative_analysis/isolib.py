"""
authors: Gloria Sheynkman, Luke Lambourne
"""
import itertools
import random
import collections

from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import pandas as pd
import seaborn as sns
from Bio.Seq import Seq
from Bio import Align


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
        orfs (list(isolib.Isoform)): protein coding isoforms of gene

    """

    def __init__(self, name, orfs):
        chroms = [orf.chrom for orf in orfs]
        strands = [orf.strand for orf in orfs]
        msg = '{} - {}'.format(name, ', '.join([orf.name for orf in orfs]))
        if len(orfs) == 0:
            raise ValueError("Need at least one isoform to define a gene\n" + msg)
        if len(set(chroms)) > 1:
            raise ValueError("All isoforms must be on same chromosome\n" + msg)
        if len(set(strands)) > 1:
            raise ValueError("All isoforms must be same strand\n" + msg)
        chrom = chroms[0]
        strand = strands[0]

        # de-duplicate exons
        exon_bounds = {(exon.start, exon.end) for iso in orfs for exon in iso.exons}
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
        for iso in orfs:
            for isoform_exon in iso.exons:
                for i, gene_exon in enumerate(self.exons):
                    if isoform_exon.start >= gene_exon.start and isoform_exon.end <= gene_exon.end:
                        isoform_exon.exon_number_on_gene = i + 1
                        break

        self.number_of_isoforms = len(orfs)
        self.name = name

        self.orfs = list(
            sorted(orfs, key=lambda x: int(x.name.split("-")[-1]))
        )

        self._orf_dict = {orf.name: orf for orf in self.orfs}
        self._pairwise_changes = {}
        GenomicFeature.__init__(
            self,
            chrom,
            strand,
            min([orf.start for orf in orfs]),
            max([orf.end for orf in orfs]),
        )

    def alternative_start(self, isoform_a, isoform_b):

        def _start_pos(iso):
            return iso.exons[0].start if self.strand == '+' else iso.exons[0].end

        a = _start_pos(self._orf_dict[isoform_a])
        b = _start_pos(self._orf_dict[isoform_b])
        return a != b

    def alternative_stop(self, isoform_a, isoform_b):

        def _stop_pos(iso):
            # TODO: is this correct? as there might be untranslated region in exon?
            return iso.exons[-1].end if self.strand == '+' else iso.exons[-1].start

        a = _stop_pos(self._orf_dict[isoform_a])
        b = _stop_pos(self._orf_dict[isoform_b])
        return a != b

    def alternative_internal_exon(self, isoform_a, isoform_b):
        a = [e.exon_number_on_gene for e in self._orf_dict[isoform_a].exons]
        b = [e.exon_number_on_gene for e in self._orf_dict[isoform_b].exons]
        start = min(min(a), min(b))
        stop = min(max(a), max(b))
        internal_a = {e for e in a if e > start and e < stop}
        internal_b = {e for e in b if e > start and e < stop}
        return internal_a != internal_b

    def alternative_3prime(self, isoform_a, isoform_b):
        for exon_a in self._orf_dict[isoform_a].exons:
            for exon_b in self._orf_dict[isoform_b].exons:
                if exon_a.exon_number_on_gene == exon_b.exon_number_on_gene:
                    if self.strand == '+':
                        if exon_a.start != exon_b.start:
                            return True
                    else:
                        if exon_a.end != exon_b.end:
                            return True
        return False

    def alternative_5prime(self, isoform_a, isoform_b):
        for exon_a in self._orf_dict[isoform_a].exons:
            for exon_b in self._orf_dict[isoform_b].exons:
                if exon_a.exon_number_on_gene == exon_b.exon_number_on_gene:
                    if self.strand == '+':
                        if exon_a.end != exon_b.end:
                            return True
                    else:
                        if exon_a.start != exon_b.start:
                            return True
        return False

    def exon_skipping(self, isoform_a, isoform_b):
        raise NotImplementedError()

    def exon_switching(self, isoform_a, isoform_b):
        raise NotImplementedError()

    def genomic_alignment_of_aa_seqs(self, subset=None):
        """genomic co-ordinates of translated regions"""
        all_orf_names = [orf.name for orf in self.orfs]
        if subset is None:
            subset = all_orf_names
        orfs = [orf for orf in self.orfs if orf.name in subset]
        if len(orfs) != len(subset):
            msg = "Missing isoforms: "
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
                raise ValueError("position is not in isoform AA sequence")
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

    def _get_exon_colors(self):
        """
        """
        merged_exon_bounds = [(exon.start, exon.end) for exon in self.exons]
        #colors_frame_1 = sns.cubehelix_palette(len(merged_exon_bounds),
        #                                       start=2.3,
        #                                       rot=0.9,
        #                                       light=0.7,
        #                                       dark=0.3,
        #                                       hue=1)
        # TODO: will this error if too many exons?
        colors_frame_1 = sns.color_palette(palette='bright',
                                           n_colors=len(merged_exon_bounds))

        colors_frame_2 = sns.cubehelix_palette(len(merged_exon_bounds),
                                               start=0.8,
                                               rot=0.9,
                                               light=0.7,
                                               dark=0.3,
                                               hue=1)
        colors_frame_3 = sns.cubehelix_palette(len(merged_exon_bounds),
                                               start=2.9,
                                               rot=0.9,
                                               light=0.7,
                                               dark=0.3,
                                               hue=1)

        def _pick_representative_cat(s):
            if s == '':
                return ''
            # reverse so returns f in case of e.g. s = 'Mf'
            return collections.Counter(s[::-1]).most_common(1)[0][0]

        exon_colors = {}
        for orf in self.orfs:
            algn = self.pairwise_changes_relative_to_reference(self.orfs[0].name, orf.name)
            algn = algn.replace('D', '')
            split_algn = [algn[sum(len(e) for e in orf.exons[:i]) // 3:
                               sum(len(e) for e in orf.exons[:i + 1]) // 3] for i in range(len(orf.exons))]
            change_cats = [_pick_representative_cat(s) for s in split_algn]
            for exon, cat in zip(orf.exons, change_cats):
                merge_exon_idx = exon.exon_number_on_gene - 1
                if cat == 'f':
                    exon_colors[(orf.name, exon.start, exon.end)] = colors_frame_3[merge_exon_idx]
                elif cat == 'F':
                    exon_colors[(orf.name, exon.start, exon.end)] = colors_frame_2[merge_exon_idx]
                else:
                    exon_colors[(orf.name, exon.start, exon.end)] = colors_frame_1[merge_exon_idx]
        return exon_colors

    def exon_diagram(self,
                     intron_nt_space=30,
                     height=0.5,
                     draw_domains=False,
                     ax=None,
                     domain_font_size=6,
                     subtle_splice_font_size=6,
                     subtle_splice_threshold=20):
        if ax is None:
            ax = plt.gca()

        merged_exon_bounds = [(exon.start, exon.end) for exon in self.exons]
        diff_exon_ends = {}
        for orf in self.orfs:
            for exon in orf.exons:
                merged_start, merged_end = merged_exon_bounds[exon.exon_number_on_gene - 1]
                if exon.start != merged_start:
                    diff_exon_ends[exon.start] = exon.start - merged_start
                if exon.end != merged_end:
                    diff_exon_ends[exon.end] = merged_end - exon.end
        if self.strand == "-":
            merged_exon_bounds = merged_exon_bounds[::-1]
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

        exon_colors = self._get_exon_colors()
        xmin = _map_position(merged_exon_bounds[0][0])
        xmax = _map_position(merged_exon_bounds[-1][1] - 1)
        for i, orf in enumerate(self.orfs):
            for exon in orf.exons:
                x_start = _map_position(exon.start)
                x_stop = _map_position(exon.end - 1)
                box = patches.Rectangle(
                    [x_start, i],
                    x_stop - x_start,
                    height,
                    lw=1,
                    ec="k",
                    fc=exon_colors[(orf.name, exon.start, exon.end)],
                    joinstyle="round",
                )
                ax.add_patch(box)
    
                # draw number of NT for small exon boundary changes
                if exon.start in diff_exon_ends:
                    num_nt_diff = diff_exon_ends[exon.start]
                    if num_nt_diff <= subtle_splice_threshold:
                        ax.text(_map_position(exon.start),
                                i + height + 0.03,
                                '{} nt'.format(num_nt_diff),
                                ha='left',
                                va='top',
                                fontsize=subtle_splice_font_size)
                if exon.end in diff_exon_ends:
                    num_nt_diff = diff_exon_ends[exon.end]
                    if num_nt_diff <= subtle_splice_threshold:
                        ax.text(_map_position(exon.end - 1),
                                i + height + 0.03,
                                '{} nt'.format(num_nt_diff),
                                ha='left',
                                va='top',
                                fontsize=subtle_splice_font_size)

            if not draw_domains:
                continue
            for dom in orf.aa_seq_features:
                if dom.name.endswith('_DBD_flank'):
                    continue
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

    def protein_diagram(self, ax=None, isoform_order=None, protein_color='pink', domain_label_rotation=0):

        def _remove_overlapping_domains(domains):
            """Keep longest domains"""
            non_overlapping = []
            for dom in sorted(domains, key=lambda x: len(x), reverse=True):
                if any(dom.start < d.end and dom.end > d.start for d in non_overlapping):
                    continue
                non_overlapping.append(dom)
            return non_overlapping

        if ax is None:
            ax = plt.gca()
        gs = gridspec.GridSpecFromSubplotSpec(len(self.orfs), 1,
                                              subplot_spec=ax.get_subplotspec(),
                                              hspace=2)
        max_seq_len = max(len(iso.aa_seq) for iso in self.orfs)

        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -99999
        aligner.extend_gap_score = -1
        aligner.target_end_gap_score = 0
        aligner.query_end_gap_score = 0

        def _get_offset(algn):
            t = algn.__str__().splitlines()[0]
            q = algn.__str__().splitlines()[2]
            target_end_gap = len(t.partition(t.strip('-'))[0])
            query_end_gap = len(q.partition(q.strip('-'))[0])
            return query_end_gap - target_end_gap

        aa_seqs = [x.aa_seq for x in self.orfs]
        offsets = [0,]
        for i in range(1, len(aa_seqs)):
            alignments = [aligner.align(aa_seqs[j], aa_seqs[i])[0] for j in range(i)]
            i_best_alignment = max(range(len(alignments)), key=lambda x: alignments[x].score)
            offsets.append(_get_offset(alignments[i_best_alignment]) + offsets[i_best_alignment])
        offsets
        x_max = max(len(iso.aa_seq) + x for iso, x in zip(self.orfs, offsets))

        exon_colors = self._get_exon_colors()

        domains = [d for d in self.orfs[0].aa_seq_features
                    if not (d.accession.endswith('_flank_N') or
                            d.accession.endswith('_flank_C'))]
        domains = _remove_overlapping_domains(domains)
        if isoform_order is not None:
            isoforms = [self[iso_id] for iso_id in isoform_order]
        else:
            isoforms = self.orfs
        for i, isoform in enumerate(isoforms):
            ax = plt.subplot(gs.new_subplotspec((i, 0), rowspan=1, colspan=1))
            seq_len = len(isoform.aa_seq)
            ax.set_ylim(0, 1)
            ax.set_xlim(0.5 - min(offsets), x_max + 0.5)
            height = 1
            prot = patches.Rectangle((0.5 + offsets[i],
                                      0.5 - height / 2),
                                     width=seq_len,
                                     height=height,
                                     clip_on=False,
                                     #facecolor=protein_color,
                                     facecolor=None,
                                     edgecolor='grey',
                                     linewidth=1.5)
            ax.add_patch(prot)

            exon_pos = 0
            for exon in isoform.exons:
                n_aa_exon = (exon.end - exon.start) / 3
                n_aa_exon = min(n_aa_exon,
                                seq_len - exon_pos)
                exon = patches.Rectangle((0.5 + offsets[i] + exon_pos,
                                          0.5 - height / 2),
                                         width=n_aa_exon,
                                         height=height,
                                         clip_on=False,
                                         facecolor=exon_colors[(isoform.name,
                                                                exon.start,
                                                                exon.end)],
                                         edgecolor=None,
                                         linewidth=0)
                ax.add_patch(exon)
                exon_pos += n_aa_exon

            for domain in domains:
                start, stop = domain.start, domain.end
                ax.axvline(x=(start - 0.5) + offsets[0],
                           ymin=-1,
                           ymax=2,
                           clip_on=False,
                           linestyle='--',
                           linewidth=2,
                           color='green')
                ax.axvline(x=(stop + 0.5) + offsets[0],
                           ymin=-1,
                           ymax=2,
                           clip_on=False,
                           linestyle='--',
                           linewidth=1.5,
                           color='green')
                if i == 0:
                    ax.text(x=(start + (stop - start) / 2) + offsets[i],
                            y=1.2,
                            s=domain.name,
                            url='http://pfam.xfam.org/family/' + domain.accession,
                            ha='center',
                            va='bottom',
                            color='green',
                            fontweight='bold',
                            fontsize=7,
                            rotation=domain_label_rotation,)

            ax.text(x=-0.5,
                    y=0.5,
                    s=isoform.name,
                    horizontalalignment='right',
                    verticalalignment='center')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['bottom'].set_position(('data', 0.5 - height / 2))
            #domain_pos = [(d.start, d.end) for d in domains]
            #xticks = [1] + [i for ij in domain_pos for i in ij] + [seq_len]
            xticks = [1, seq_len]
            ax.set_xticks([x + offsets[i] for x in xticks])
            ax.set_xticklabels([str(x) for x in xticks])
            ax.xaxis.set_tick_params(rotation=90)
            ax.set_yticks([])

    def orf_pairs(self):
        return list(itertools.combinations(self.orfs, 2))

    def __getitem__(self, orf_id):
        if orf_id in self._orf_dict:
            return self._orf_dict[orf_id]
        for orf in self.orfs:
            if orf_id == orf.clone_acc:
                return orf
        raise KeyError(orf_id)

    def __contains__(self, orf_id):
        return orf_id in self._orf_dict

    def __repr__(self):
        s = "Gene: {}\n".format(self.name)
        s += "Isoforms: " + str([orf.name for orf in self.orfs])
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
        ensembl_transcript_names=None
    ):
        self.name = name
        self.ensembl_protein_ids = ensembl_protein_ids
        self.ensembl_transcript_names = ensembl_transcript_names
        self.ensembl_transcript_ids = ensembl_transcript_ids


        if isinstance(CDS_nt_seq, str):
            if len(CDS_nt_seq) % 3 != 0:
                raise ValueError('CDS sequence length must be multiple of 3')
            stop_codons = {'TAG', 'TAA', 'TGA'}
            codons = [CDS_nt_seq[i:i + 3] for i in range(0, len(CDS_nt_seq), 3)]
            if any(codon in stop_codons for codon in codons[:-1]):
                raise ValueError('CDS sequnce must not contain early stop codon')
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
                self.nt_seq[i:i + 3] for i in range(0, len(self.aa_seq) * 3, 3)
            ]
        else:
            raise NotImplementedError()
        if len(self.nt_seq) != len(self.aa_seq) * 3:
            raise UserWarning('Inconsistent DNA and AA sequences – ' + self.name)

        if len(exons) == 0:
            raise ValueError(self.name + " - Need at least one exon to define an isoform")
        chroms = [exon.chrom for exon in exons]
        strands = [exon.strand for exon in exons]
        if len(set(chroms)) > 1:
            raise ValueError(self.name + " - All exons must be on same chromosome")
        if len(set(strands)) > 1:
            raise ValueError(self.name + " - All exons must be same strand")
        len_all_exons = sum(len(e) for e in exons)
        if len_all_exons != len(self.nt_seq):
            msg = 'exons not same length as CDS sequence minus stop\n'
            msg += 'exons = {} nt, CDS = {} nt\n'.format(len_all_exons, len(self.nt_seq))
            msg += self.name
            raise UserWarning(msg)
        is_neg_strand = strands[0] == "-"
        self.exons = sorted(exons, key=lambda x: x.start, reverse=is_neg_strand)

        residues = []
        genomic_coords = [i for exon in self.exons for i in exon.genomic_coords()]
        codon_genomic_coords = [(a, b, c) for a, b, c in zip(genomic_coords[0::3],
                                                             genomic_coords[1::3],
                                                             genomic_coords[2::3])
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

    @property
    def dna_binding_domains(self):
        # hiding import to avoid circular imports
        from data_loading import load_dbd_accessions
        dbd_acc = load_dbd_accessions()
        return [dom for dom in self.aa_seq_features if dom.accession in dbd_acc]

    def __repr__(self):
        if self.ensembl_protein_ids is not None:
            ids = ' / '.join(['|'.join(self.ensembl_transcript_names),
                            '|'.join(self.ensembl_protein_ids),
                            '|'.join(self.ensembl_transcript_ids)])
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
            stop_codons = {'TAG', 'TAA', 'TGA'}
            codons = [clone_seq[i:i + 3] for i in range(0, len(clone_seq), 3)]
            if any(codon in stop_codons for codon in codons):
                return ''.join(codons[:min(codons.index(stop) for stop in stop_codons if stop in codons)])
            elif len(clone_seq) % 3 == 0:
                return clone_seq
            else:
                return clone_seq[:-(len(clone_seq) % 3)]

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
                trimmed_exons[-1].start = trimmed_exons[-1].end - (len(cds_nt_seq) - sum(len(e) for e in trimmed_exons[:-1]))
            else:
                trimmed_exons[-1].end = trimmed_exons[-1].start + (len(cds_nt_seq) - sum(len(e) for e in trimmed_exons[:-1]))
            if sum(len(e) for e in trimmed_exons) != len(cds_nt_seq):
                raise UserWarning(self.name, ' - problem with reducing exons to coding-only')
        else:
            trimmed_exons = exons

        Isoform.__init__(self,
                         name=clone_name,
                         exons=trimmed_exons,
                         CDS_nt_seq=cds_nt_seq,
                         aa_seq=aa_seq,
                         UTR_5prime_nt_seq=UTR_5prime_nt_seq,
                         UTR_3prime_nt_seq=UTR_3prime_nt_seq,
                         ensembl_protein_ids=ensembl_protein_ids,
                         ensembl_transcript_ids=ensembl_transcript_ids,
                         ensembl_transcript_names=ensembl_transcript_ids)


    def is_novel_isoform(self):
        return self.ensembl_transcript_ids is None

    def __repr__(self):
        s = 'Clone acc: {}\n'.format(self.clone_acc)
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
