import csv
from collections import Counter
from itertools import islice
from mist.apps import utils
from mist.lib.constants import WTA, TARGETED


class LibraryMetrics:
    """
    LibraryMetrics keeps track of the metrics pertaining to a single library, used within the
    AnotateReads node. It also contains the helper functions, including (1) combine_metrics, used to sum up the
    metrics dictionary from multiple LibraryMetrics instances; and (2) parse_read_file, which calculates the
    library metrics directly from R1 and R2 files.
    """

    def __init__(self, trueno, assay, umi_option, filter_stats_fps=None, quality_metrics_stats_fps=None):
        self.trueno = trueno
        self.assay = assay
        self.umi_option = umi_option
        self.total_reads = 0

        if filter_stats_fps is not None:
            self.filter_stats_fps = filter_stats_fps
        else:
            self.filter_stats_fps = []
        if quality_metrics_stats_fps is not None:
            self.quality_metrics_stats_fps = quality_metrics_stats_fps
        else:
            self.quality_metrics_stats_fps = []

        # targeted
        self.start_pos = 0
        self.align_len = 0
        self.num_mapped = 0
        self.tag_start_pos = 0
        self.tag_align_len = 0

        # wta
        self.not_uniq = 0
        self.ambiguous = 0
        self.no_feature = 0
        self.not_aligned = 0
        self.mapped_to_genes = 0
        self.valid = 0
        self.total_phix = 0
        self.isCell = 0
        self.isCell_not_uniq = 0
        self.isCell_ambiguous = 0
        self.isCell_no_feature = 0
        self.isCell_not_aligned = 0
        self.isCell_mapped_to_genes = 0

        # trueno
        self.total_tag_reads = 0
        self.tag_is_cell = 0
        self.tag_valid = 0
        self.tag_counts_each = {}

    def parse_row(self, read1, read2):
        """determine the validity of the read pair, keeping track of the overall metrics simultaneously

        Args:
            read1, read2: a single, corresponding line from the R1 and R2 file

        Returns:
            whether the read pair is valid

        """

        try:
            cell_label, cell_label_mismatch, umi, poly_t_check_passing = read1
        except ValueError:
            raise ValueError("Unable to parse read 1: {}".format(read1))
        else:
            poly_t_check_passing = {"T": True, "F": False}[poly_t_check_passing]
        if self.assay == TARGETED:
            try:
                gene, cigar, start_pos_passing, length_match_passing, is_phiX, abumi = read2
            except ValueError:
                raise ValueError("Unable to parse read 2: {}".format(read2))
            else:
                start_pos_passing = {"1": True, "0": False, "x": "skip/unaligned"}[start_pos_passing]
                length_match_passing = {"1": True, "0": False, "*": "skip/unaligned"}[length_match_passing]
        elif self.assay == WTA:
            try:
                gene, score, status, read2_three_prime, transcript_name, pos4, fragment_length, abumi = read2
            except ValueError:
                raise ValueError("Unable to parse read 2: {}".format(read2))

        self.total_reads += 1

        is_tag = gene.endswith('stAbO')
        # Determine R2 status without considering R1
        if self.assay == TARGETED:
            if start_pos_passing is True:
                self.start_pos += 1
                if is_tag is True:
                    self.tag_start_pos += 1

            if length_match_passing is True:
                self.align_len += 1
                if is_tag is True:
                    self.tag_align_len += 1
                if start_pos_passing is True:
                    self.num_mapped += 1

        if self.assay == 'WTA':
            if gene != '*':  # counts uniquely aligned to genes and multialigned to unique genes
                self.mapped_to_genes += 1
            elif status.startswith('alignment'):
                self.not_uniq += 1
            elif status.startswith('ambiguous'):
                # only includes ambiguously mapped reads with multiple gene features, matching htseq definition
                self.ambiguous += 1
            elif status == 'x':
                self.not_aligned += 1
            else:
                # match htseq definition
                self.no_feature += 1

        # skip reads mapped to phix
        if gene == 'phiX174':
            self.total_phix += 1
            self.mapped_to_genes -= 1
            return False

        if is_tag:
            self.total_tag_reads += 1
            try:
                self.tag_counts_each[gene] += 1
            except KeyError:
                self.tag_counts_each[gene] = 1

        # only look at cells with proper cell label
        if self.assay == WTA and 'x' not in cell_label:
            self.isCell += 1
            if gene != '*':  # counts uniquely aligned to genes and multialigned to unique genes
                self.isCell_mapped_to_genes += 1
            elif status.startswith('alignment'):
                self.isCell_not_uniq += 1
            elif status.startswith('ambiguous'):
                # only includes ambiguously mapped reads with multiple gene features, matching htseq definition
                self.isCell_ambiguous += 1
            elif status == 'x':
                self.isCell_not_aligned += 1
            else:
                # match htseq definition
                self.isCell_no_feature += 1
            if is_tag is True:
                self.tag_is_cell += 1
        elif self.assay != WTA and 'x' not in cell_label:
            self.isCell += 1
            if is_tag is True:
                self.tag_is_cell += 1
        else:
            return False

        # only valid if MI is correct length and PolyT exists

        if self.assay == WTA:
            # only include gene uniquely mapped to genes (htseq definition)
            # or multialigned to unique genes
            # in the future row[7] not in ['x', 'intergenic', 'no_feature'] or row[7].startswith('ambiguous')
            # can be further annotated
            if 'x' in cell_label \
                    or len(umi) < 8 \
                    or (status != '1' and 'MAcorrected' not in status):
                return False
        else:
            # meet start position and minimum alignment length
            if start_pos_passing is not True or length_match_passing is not True:
                return False

        if not poly_t_check_passing:
            return False

        if abumi:
            if self.umi_option == 1:
                umi = abumi
            elif self.umi_option == 2:
                umi = umi + abumi

        if 'N' in umi:
            return False

        self.valid += 1
        if is_tag:
            self.tag_valid += 1

        return True

    def calculate_relative_metrics(self):
        """determine relative metrics

        Returns:
            dictionary that contains the read stats and, if applicable, the filtering stats and tag stats

        """

        statistics = {}

        percentage_cells = 100.0 * self.isCell / self.total_reads
        percentage_valid = 100.0 * self.valid / self.total_reads
        percentage_reads_phix = 100.0 * self.total_phix / self.total_reads

        read_stats = [self.total_reads,
                      self.valid, percentage_valid,
                      self.isCell, percentage_cells,
                      self.total_phix, percentage_reads_phix]

        if self.assay == 'WTA':
            percentage_not_uniq = 100.0 * self.not_uniq / self.total_reads
            percentage_ambiguous = 100.0 * self.ambiguous / self.total_reads
            percentage_no_feature = 100.0 * self.no_feature / self.total_reads
            percentage_not_aligned = 100.0 * self.not_aligned / self.total_reads
            percentage_mapped_to_genes = 100.0 * self.mapped_to_genes / self.total_reads
            percentage_isCell_not_uniq = 100.0 * self.isCell_not_uniq / self.total_reads
            percentage_isCell_ambiguous = 100.0 * self.isCell_ambiguous / self.total_reads
            percentage_isCell_no_feature = 100.0 * self.isCell_no_feature / self.total_reads
            percentage_isCell_not_aligned = 100.0 * self.isCell_not_aligned / self.total_reads
            percentage_isCell_mapped_to_genes = 100.0 * self.isCell_mapped_to_genes / self.total_reads

            read_stats.extend([self.mapped_to_genes, percentage_mapped_to_genes,
                               self.not_aligned, percentage_not_aligned,
                               self.not_uniq, percentage_not_uniq,
                               self.no_feature, percentage_no_feature,
                               self.ambiguous, percentage_ambiguous,
                               self.isCell_mapped_to_genes, percentage_isCell_mapped_to_genes,
                               self.isCell_not_aligned, percentage_isCell_not_aligned,
                               self.isCell_not_uniq, percentage_isCell_not_uniq,
                               self.isCell_no_feature, percentage_isCell_no_feature,
                               self.isCell_ambiguous, percentage_isCell_ambiguous])

        else:
            percentage_start_pos = 100.0 * self.start_pos / self.total_reads
            percentage_align_len = 100.0 * self.align_len / self.total_reads
            percentage_mapped = 100.0 * self.num_mapped / self.total_reads

            read_stats.extend([self.num_mapped, percentage_mapped,
                               self.start_pos, percentage_start_pos,
                               self.align_len, percentage_align_len])

        statistics['read_stats'] = read_stats

        if self.trueno:
            try:
                percentage_cells_tag = 100.0 * self.tag_is_cell / self.total_tag_reads
                percentage_valid_tag = 100.0 * self.tag_valid / self.total_tag_reads
            except ZeroDivisionError:
                percentage_cells_tag = percentage_valid_tag = 0

            tag_stats = [self.total_tag_reads,
                         self.tag_valid, percentage_valid_tag,
                         self.tag_is_cell, percentage_cells_tag]

            if self.assay != 'WTA':
                try:
                    percentage_start_pos_tag = 100.0 * self.tag_start_pos / self.total_tag_reads
                    percentage_align_len_tag = 100.0 * self.tag_align_len / self.total_tag_reads
                except ZeroDivisionError:
                    percentage_start_pos_tag = percentage_align_len_tag = 0

                tag_stats.extend([self.tag_start_pos, percentage_start_pos_tag,
                                  self.tag_align_len, percentage_align_len_tag])

            statistics['tag_stats'] = tag_stats
            statistics['tag_counts_each'] = self.tag_counts_each
            statistics['total_reads'] = self.total_reads
            statistics['total_tag_reads'] = self.total_tag_reads

        if self.filter_stats_fps:
            stats = []
            for file in self.filter_stats_fps:
                with utils.quick_gzip_open(file) as f:
                    stats.append([float(x) for x in f.readline().split(',')])

            # add up all counts
            stats = [sum(metric) for metric in zip(*stats)]
            total_num_reads = stats[0]
            # compute percentages
            pct = []
            for count in islice(stats, 1, None):
                pct.append(100 - (count * 100.0 / total_num_reads))
            filter_stats = [total_num_reads, stats[1], pct[0], stats[2], pct[1], stats[3], pct[2], stats[4], pct[3]]

            statistics['filter_stats'] = filter_stats

        if self.quality_metrics_stats_fps:
            picard_quality_metric_dicts = []
            combined_picard_quality_metrics = Counter()
            total_read_length = 0
            for quality_metrics_fp in self.quality_metrics_stats_fps:
                with utils.quick_gzip_open(quality_metrics_fp) as quality_metrics_file:
                    quality_metrics_reader = csv.DictReader(quality_metrics_file)
                    quality_metrics = next(quality_metrics_reader)
                    quality_metrics = {key: float(val) for key, val in quality_metrics.items()}
                    total_read_length += float(quality_metrics['TOTAL_READS']) * float(
                        quality_metrics['READ_LENGTH'])
                    picard_quality_metric_dicts.append(quality_metrics)
                    combined_picard_quality_metrics.update(quality_metrics)

            ratio_q30_bases = float(combined_picard_quality_metrics['PF_Q30_BASES']) / float(combined_picard_quality_metrics['PF_BASES'])
            pct_q30_bases = ratio_q30_bases * 100.0
            combined_picard_quality_metrics['PCT_Q30_BASES'] = pct_q30_bases

            combined_picard_quality_metrics['READ_LENGTH'] = int(
                total_read_length / combined_picard_quality_metrics['TOTAL_READS'])

            statistics['picard_stats'] = combined_picard_quality_metrics

        return statistics


    @classmethod
    def parse_read_file(cls, r1s, r2s, trueno, assay, filter_stats_fp=None, quality_metrics_fp=None):

        """

        Args:
            R1s and R2s: paths, alphabetically sorted, whose metrics will be ascertained

        Returns: the metrics for the R1s and R2s, for testing purposes

        """

        lib_metrics = cls(trueno, assay, filter_stats_fp, quality_metrics_fp)

        read1f = utils.csv_input(r1s)
        read2f = utils.csv_input(r2s)

        for r1, r2 in zip(read1f, read2f):
            lib_metrics.parse_row(r1, r2)

        stats = lib_metrics.calculate_relative_metrics()

        return stats

    @classmethod
    def combine_metrics(cls, libraries, metric_name, inverse_ratios=False):
        """

        Args:
            libraries: library files, including R1, R2 and filtering stats
            metric_name: the metric to be combined, which should be formatted thusly:
                [total reads, metric_1, pct_metric_1, metric_2, pct_metric_2...]
            inverse_ratios: if true, take the ratio of the total to the total; if false,
                take (1 - the ratio of the subset to the subset)

        Returns: Combined metric (in list form) from the specified libraries

        """

        metrics = [library[metric_name] for library in libraries if metric_name in library]
        cols = zip(*metrics)

        assert all(len(metrics[0]) == len(other_metric) for other_metric in metrics)

        combined_stats = []
        for i, col in enumerate(cols):
            if i > 0 and i % 2 == 0:  # Every other column, starting with the third, is a percentage
                total_reads = combined_stats[0]
                previous_metric = combined_stats[-1]
                try:
                    if inverse_ratios is False:
                        ratio = (float(previous_metric) / float(total_reads)) * 100
                    else:
                        ratio = (1.0 - float(previous_metric) / float(total_reads)) * 100
                except ZeroDivisionError:
                    raise ZeroDivisionError(
                        "Total reads is 0, please double check! "
                        "One possible error is sample-tag-version "
                        "is used for runs without using sample tags."
                    )
                combined_stats.append(round(ratio, ndigits=2))
            else:
                sum_of_metrics = sum(col)
                combined_stats.append(sum_of_metrics)
        return combined_stats

    @classmethod
    def combine_picard_metrics(cls, libraries):
        """

        Args:
            qual_metric_dicts: library files, including R1, R2 and filtering stats

        Returns: Combined picard metrics

        """

        qual_metric_dicts = [library['picard_stats'] for library in libraries if 'picard_stats' in library]

        combined_picard_quality_metrics = Counter()
        total_read_length = 0
        for picard_quality_metrics_dict in qual_metric_dicts:
            combined_picard_quality_metrics.update(picard_quality_metrics_dict)
            total_read_length += float(picard_quality_metrics_dict['TOTAL_READS']) * float(picard_quality_metrics_dict['READ_LENGTH'])

        # Recalculate Pct_Q30_Bases
        ratio_q30_bases = float(combined_picard_quality_metrics['PF_Q30_BASES']) / float(combined_picard_quality_metrics['PF_BASES'])
        pct_q30_bases = ratio_q30_bases * 100.0
        combined_picard_quality_metrics['PCT_Q30_BASES'] = pct_q30_bases

        # Recalculate READ_LENGTH
        combined_picard_quality_metrics['READ_LENGTH'] = int(total_read_length/combined_picard_quality_metrics['TOTAL_READS'])

        return combined_picard_quality_metrics

