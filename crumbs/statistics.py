# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of seq_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with seq_crumbs. If not, see <http://www.gnu.org/licenses/>.

# pylint: disable=C0111

from __future__ import division
from collections import Counter
import operator
import re

from crumbs.settings import get_setting
from crumbs.iterutils import rolling_window
from crumbs.utils import approx_equal
from crumbs.seq import get_str_seq, get_length, get_int_qualities


LABELS = {'title': 'histogram', 'xlabel': 'values',
          'ylabel': 'count', 'minimum': 'minimum',
          'maximum': 'maximum', 'average': 'average',
          'variance': 'variance', 'sum': 'sum',
          'items': 'items', 'quartiles': 'quartiles'}


class IntCounter(Counter):
    '''This is a subclass  of the Counter on python collections.

    It adds some statistical functionalities to the class'''

    def __init__(self, *args, **kwargs):
        'It inialices a '
        super(IntCounter, self).__init__(*args, **kwargs)
        self.labels = LABELS.copy()

    @property
    def min(self):
        'Get the minimun value'
        return min(self.keys())

    @property
    def max(self):
        'Get the maximun value'
        return max(self.keys())

    @property
    def count(self):
        'It returns the count of the values stored in the array'
        return sum(self.values())

    @property
    def median(self):
        'It calculates the median of the values appended'
        quotient, remainder = divmod(self.count, 2)
        if remainder == 0:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            return (val1 + val2) / 2
        else:
            return self._get_value_for_index(quotient)

    @property
    def sum(self):
        'It gets the sum of the values'
        sum_ = 0
        for index, value in self.items():
            sum_ += (index * value)
        return int(sum_)

    @property
    def average(self):
        'It calculates the average'
        count = self.count
        sum_ = self.sum
        return sum_ / count

    @property
    def variance(self):
        'It gets the variance of the values'
        mean = self.average
        sum_ = 0
        for index, counts in self.items():
            sum_ += ((index - mean) ** 2) * counts
        return sum_ / self.count

    @property
    def quartiles(self):
        'It returns the quartiles'
        num_items = self.count
        if num_items < 4:
            msg = 'At least 4 values are required to calculate the quartiles'
            raise RuntimeError(msg)
        # quartile 1
        quotient, remainder = divmod(num_items + 1, 4)
        if not remainder:
            quartile1 = self._get_value_for_index(quotient - 1)
        else:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            quartile1 = (val1 + val2) / 2
        # quartile 3
        quotient, remainder = divmod((num_items + 1) * 3, 4)
        if not remainder:
            quartile3 = self._get_value_for_index(quotient - 1)
        else:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            quartile3 = (val1 + val2) / 2
        return quartile1, self.median, quartile3

    @property
    def irq(self):
        'It gets the interquartile range'
        # pylint: disable=W0612
        quart1 = self.quartiles[0]
        quart3 = self.quartiles[2]
        return quart3 - quart1

    @property
    def outlier_limits(self):
        'It returns the intercuartile'
        # pylint: disable=W0612
        quart1 = self.quartiles[0]
        quart3 = self.quartiles[2]

        iqr = self.irq
        limit_distance = round(iqr * 1.5)

        start = int(quart1 - limit_distance)
        end = int(quart3 + limit_distance)
        return (start, end)

    def _get_value_for_index(self, position):
        '''It takes a position and it returns the value for the given index'''
        cum_count = 0
        for index in sorted(self.keys()):
            count = self[index]
            cum_count += count
            if position <= cum_count - 1:
                return index
        else:
            if position >= cum_count:
                raise IndexError('You asked for an index beyond the scope')
            return index

    def _calculate_dist_range(self, min_, max_, outlier_threshold):
        'it calculates the range for the histogram'
        if ((min_ is not None or max_ is not None) and
            outlier_threshold is None):
            msg = 'You can not pass max, min and outlier_threslhosld to '
            msg += 'calculate distribution range'
            raise ValueError(msg)

        if min_ is None:
            min_ = self.min
        if max_ is None:
            max_ = self.max

        if outlier_threshold:
            left_limit = self.count * outlier_threshold / 100
            rigth_limit = self.count - left_limit
            left_value = self._get_value_for_index(left_limit)
            rigth_value = self._get_value_for_index(rigth_limit)

            if min_ < left_value:
                min_ = left_value
            if max_ > rigth_value:
                max_ = rigth_value
        return min_, max_

    def calculate_bin_edges(self, min_, max_, n_bins=None):
        'It calculates the bin_edges'
        min_bins = get_setting('MIN_BINS')
        max_bins = get_setting('MAX_BINS')
        if n_bins is None:
            num_values = max_ - min_
            if num_values == 0:
                n_bins = 1
            elif num_values < min_bins:
                n_bins = num_values
            else:
                n_bins = int(self.count / get_setting('MEAN_VALUES_IN_BIN'))
                if n_bins < min_bins:
                    n_bins = min_bins
                if n_bins > max_bins:
                    n_bins = max_bins
                if n_bins > num_values:
                    n_bins = num_values

        # now we can calculate the bin edges
        distrib_span = max_ - min_ if max_ != min_ else 1

        if distrib_span % n_bins:
            distrib_span = distrib_span + n_bins - (distrib_span % n_bins)
        bin_span = distrib_span // n_bins
        bin_edges = [min_ + bin_ * bin_span for bin_ in range(n_bins + 1)]
        return bin_edges

    def calculate_distribution(self, bins=None, min_=None, max_=None,
                               outlier_threshold=None):
        'It returns an histogram with the given range and bin'
        distrib = []
        min_, max_ = self._calculate_dist_range(min_, max_, outlier_threshold)
        if min_ is None or max_ is None:
            return None
        bin_edges = self.calculate_bin_edges(min_, max_, bins)
        for bin_index, left_edge in enumerate(bin_edges):
            try:
                rigth_edge = bin_edges[bin_index + 1]
            except IndexError:
                break
            sum_values = 0

            for index2 in sorted(self.keys()):
                value = self[index2]
                if index2 > rigth_edge:
                    break

                elif (left_edge <= index2  and index2 < rigth_edge or
                     left_edge <= index2 and index2 == max_):
                    sum_values += value

            distrib.append(sum_values)
        return {'counts': distrib, 'bin_limits': bin_edges}

    def update_labels(self, labels):
        'It prepares the labels for output files'
        self.labels.update(labels)

    def count_relative_to_value(self, value, comparison):
        'It counts the ints greater, equal, etc, relative to the given value.'
        return sum([counts for val, counts in self.items()
                                                    if comparison(val, value)])

    def __add__(self, other):
        'Add counts from two counters.'
        counter_python = super(IntCounter, self).__add__(other)
        return self.__class__(counter_python)

    def __str__(self):
        'It writes some basic stats of the values'
        if self.count != 0:
            labels = self.labels
            # now we write some basic stats
            format_num = lambda x: '{:,d}'.format(x) if isinstance(x, int) \
                                                                else '%.2f' % x
            text = '{}: {}\n'.format(labels['minimum'], format_num(self.min))
            text += '{}: {}\n'.format(labels['maximum'], format_num(self.max))
            text += '{}: {}\n'.format(labels['average'],
                                      format_num(self.average))

            if labels['variance'] is not None:
                text += '{}: {}\n'.format(labels['variance'],
                                          format_num(self.variance))
            if labels['sum'] is not None:
                text += '{}: {}\n'.format(labels['sum'],
                                          format_num(self.sum))
            if labels['items'] is not None:
                text += '{}: {}\n'.format(labels['items'], self.count)
            text += '\n'

            distrib = self.calculate_distribution()
            text += draw_histogram_ascii(distrib['bin_limits'],
                                         distrib['counts'])
            return text
        return ''


def draw_histogram_ascii(bin_limits, counts):
    'It draws an ASCII histogram'

    fill_char = '*'

    assert len(bin_limits) == len(counts) + 1
    # pylint: disable=W0108
    number_to_str = lambda n: '{:d}'.format(n) if isinstance(n, int) else \
                                                            '{:.2f}'.format(n)

    # we gather all bin limits and we calculate the longest number
    bin_start = None
    bin_end = bin_limits[0]
    max_ndigits = len(number_to_str(bin_end))
    max_count_ndigits = 0
    bins = []
    for bin_limit, cnt in zip(bin_limits[1:], counts):
        bin_start, bin_end = bin_end, bin_limit
        n_digits = len(number_to_str(bin_end))
        if max_ndigits < n_digits:
            max_ndigits = n_digits
        n_digits = len(number_to_str(cnt))
        if max_count_ndigits < n_digits:
            max_count_ndigits = n_digits
        bins.append((bin_start, bin_end))

    limit_fmt_int = '{:>' + str(max_ndigits) + 'd}'
    limit_fmt_float = '{:>' + str(max_ndigits) + '.5f}'
    limit_to_padded_str = lambda n: limit_fmt_int.format(n) \
                           if isinstance(n, int) else limit_fmt_float.format(n)

    count_fmt = '{:>' + str(max_count_ndigits) + 'd}'
    count_to_padded_str = lambda n: count_fmt.format(n)

    result = []
    for bin_, cnt in zip(bins, counts):
        line = ''
        line += '['
        line += limit_to_padded_str(bin_[0])
        line += ' , '
        line += limit_to_padded_str(bin_[1])
        line += '[ ('
        line += count_to_padded_str(cnt)
        line += '): '
        result.append(line)

    # pylint: disable=W0141
    max_count = max(counts)
    max_header_len = max(map(len, result))
    max_hist_width = get_setting('MAX_WIDTH_ASCII_PLOT') - max_header_len
    counts_ratio = max_hist_width / max_count

    result2 = []
    for line, cnt in zip(result, counts):
        line += fill_char * int(cnt * counts_ratio)
        line += '\n'
        result2.append(line)

    return ''.join(result2)


class IntBoxplot(object):
    'It represents a Box and whisker plot'
    def __init__(self):
        'The init'
        self.counts = {}

    def append(self, category, value):
        'It appends a value to the distribution corresponding to a category'
        counts = self.counts
        try:
            cat_counts = counts[category]
        except KeyError:
            counts[category] = IntCounter()
            cat_counts = counts[category]
        cat_counts[value] += 1

    @property
    def aggregated_array(self):
        'It returns the IntSummarizedArray of all appended values.'
        aggregated_dist = None
        for cat_count in self.counts.viewvalues():
            if aggregated_dist is None:
                aggregated_dist = cat_count
            else:
                aggregated_dist += cat_count
        return aggregated_dist

    @property
    def ascii_plot(self):
        'It returns an string with an ASCII representation'
        distributions = self.counts
        categories = sorted(distributions.viewkeys())
        distrib_descriptions = {}
        for category in categories:
            distrib = distributions[category]
            min_ = distrib.min
            max_ = distrib.max
            try:
                quart1, median, quart3 = distrib.quartiles
                distrib_descriptions[category] = {'min': min_, 'max': max_,
                                            'quart1': quart1, 'median': median,
                                            'quart3': quart3}
            except RuntimeError:
                distrib_descriptions[category] = None

        min_value = min((d['min'] for d in distrib_descriptions.values()
                                                             if d is not None))
        max_value = max((d['max'] for d in distrib_descriptions.values()
                                                             if d is not None))

        len_str_float = lambda x: len('{:.1f}'.format(x))

        # how wide is the widest number of each type (min, max, median, etc)
        widths = {}
        for distrib in distrib_descriptions.values():
            if  distrib is None:
                continue
            for key in distrib.keys():
                value = len_str_float(distrib[key])
                if key not in widths or widths[key] < value:
                    widths[key] = value

        widths['labels'] = max([len(str(d)) for d in distrib_descriptions])

        category_format = '{:>' + str(widths['labels']) + 's}'
        header_format = category_format + ':{:>' + str(widths['min']) + '.1f}'
        header_format += ',{:>' + str(widths['quart1']) + '.1f}'
        header_format += ',{:>' + str(widths['median']) + '.1f}'
        header_format += ',{:>' + str(widths['quart3']) + '.1f}'
        header_format += ',{:>' + str(widths['max']) + '.1f}'

        header_len = None
        headers, distribs = [], []
        for category in categories:
            distrib = distrib_descriptions[category]
            if distrib is not None:
                header = header_format.format(str(category), distrib['min'],
                                          distrib['quart1'], distrib['median'],
                                          distrib['quart3'], distrib['max'])
                header += ' '
                if header_len is None or header_len < len(header):
                    header_len = len(header)
            else:
                header = None
            headers.append(header)
            distribs.append(distrib)

        plot_width = (get_setting('MAX_WIDTH_ASCII_PLOT') - widths['labels'] -
                      len(str(max_value)) - header_len)
        val_per_pixel = (max_value - min_value) / plot_width
        if val_per_pixel == 0:
            val_per_pixel = max_value / plot_width

        to_axis_scale = lambda x: int((x - min_value) / val_per_pixel)

        result = ''

        for distrib, header in zip(distribs, headers):
            if distrib is not None:
                line = header
                min_ = to_axis_scale(distrib['min'])
                line += ' ' * (min_ - 1)
                line += '<'  # min
                quart1 = to_axis_scale(distrib['quart1'])
                line += '-' * ((quart1 - min_) - 1)
                line += '['  # quartil 1
                median = to_axis_scale(distrib['median'])
                line += '=' * ((median - quart1) - 1)
                line += '|'  # median 1
                quart3 = to_axis_scale(distrib['quart3'])
                line += '=' * ((quart3 - median) - 1)
                line += ']'  # quart 2
                max_ = to_axis_scale(distrib['max'])
                line += '-' * ((max_ - quart3) - 1)
                line += '>'  # max
            else:
                line = category_format.format(str(category), 0, 0, 0, 0, 0)
            line += '\n'
            result += line

        axis1 = ' ' * header_len + '+' + '-' * (plot_width - 2) + '+' + '\n'
        axis2 = ' ' * header_len + str(min_value)
        axis2 += ' ' * (plot_width - len(str(min_value)) - len(str(max_value))
                       + 1)
        axis2 += str(max_value) + '\n'

        result += axis1
        result += axis2
        return result

    def __nonzero__(self):
        'It returns True if the object holds any counts'
        return bool(self.counts)


class NuclFreqsPlot(object):
    'It represents the base frequencies along the read lengts'
    def __init__(self,
                 count_up_to_base=get_setting('DEF_PLOT_FREQS_UP_TO_BASE')):
        'The init'
        self.counts = {}
        self.count_up_to_base = count_up_to_base

    def append(self, base_index, nucleotide):
        'It appends a value to the distribution corresponding to a category'
        count_up_to_base = self.count_up_to_base
        if count_up_to_base is not None and base_index > count_up_to_base:
            return
        counts = self.counts
        try:
            cat_counts = counts[base_index]
        except KeyError:
            counts[base_index] = Counter()
            cat_counts = counts[base_index]
        nucleotide = nucleotide.upper()
        if nucleotide not in ('A', 'C', 'T', 'G'):
            nucleotide = 'N'
        cat_counts[nucleotide] += 1

    @property
    def ascii_plot(self):
        'It plots columns with the nucleotide frequencies'
        nucls = ('A', 'C', 'G', 'T', 'N')
        plot_nucls = ('a', 'C', 'g', 'T', 'n')
        locs = self.counts.keys()
        if not locs:
            return ''
        loc_max = max(locs)
        loc_min = min(locs)
        loc_width = len(str(loc_max))
        loc_fmt = '{:>' + str(loc_width) + 'd}'
        counts = self.counts

        def _header_for_nucl(loc):
            'It returns the header for the given position'
            header = loc_fmt.format(loc) + ' ('
            count = counts.get(loc, {})
            freqs = [count.get(n, 0) for n in nucls]
            tot_bases = sum(freqs)
            freqs = [f / tot_bases for f in freqs]
            freq_strs = ['{}: {:.2f}'.format(n, f)
                                                 for n, f in zip(nucls, freqs)]
            header += ', '.join(freq_strs) + ') | '
            return header, freqs

        header_len = len(_header_for_nucl(0)[0])
        plot_width = get_setting('MAX_WIDTH_ASCII_PLOT') - header_len
        val_per_pixel = 1 / plot_width
        plot = ''
        for loc in range(loc_min, loc_max + 1):
            header, freqs = _header_for_nucl(loc)
            assert approx_equal(sum(freqs), 1)
            line = header

            remainder_freqs = [float(re.sub('\d\.', '0.', str(f)))
                                                                for f in freqs]
            round_freqs = [int(round(f / val_per_pixel)) for f in freqs]

            pixels_remaining = plot_width - sum(round_freqs)

            if pixels_remaining > 0:
                add_to_freq = remainder_freqs.index(max(remainder_freqs))
                round_freqs[add_to_freq] += (plot_width - sum(round_freqs))
            elif pixels_remaining < 0:
                add_to_freq = remainder_freqs.index(min(remainder_freqs))
                round_freqs[add_to_freq] -= (plot_width - sum(round_freqs))
            assert approx_equal(sum(round_freqs), plot_width)
            line += ''.join([n * f for f, n in zip(round_freqs, plot_nucls)])
            line += '\n'
            plot += line
        return plot


class KmerCounter(object):
    'It counts kmers in the given sequences'
    def __init__(self, kmer_size=get_setting('DEFAULT_KMER_SIZE')):
        'initiator'
        self._kmer_size = kmer_size
        self._counter = Counter()

    def count_seq(self, serie):
        'It adds the kmers of the given iterable/serie'
        for kmer in rolling_window(serie, self._kmer_size):
            self._counter[kmer] += 1

    @property
    def values(self):
        'It returns the values of the counter'
        return iter(self._counter.viewvalues())

    def most_common(self, num_items):
        'return most common kmers with their counts'
        return self._counter.most_common(num_items)


def _calculate_rawscore(string):
    'It returns a non-normalized dustscore'
    triplet_counts = Counter()
    for triplet in rolling_window(string, 3):
        # It should do something with non ATCG, but we sacrifice purity for
        # speed. Maybe we should reconsider this
        triplet_counts[triplet.upper()] += 1

    return sum(tc * (tc - 1) * 0.5 for tc in triplet_counts.viewvalues())


def calculate_dust_score(seq):
    '''It returns the dust score.

    From: "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA
    Sequences"
    doi:10.1089/cmb.2006.13.1028

    and re-implemented from PRINSEQ
    '''
    seq = get_str_seq(seq)
    length = len(seq)
    if length == 3:
        return 0
    if length <= 5:
        return None

    windowsize = get_setting('DUST_WINDOWSIZE')
    windowstep = get_setting('DUST_WINDOWSTEP')

    dustscores = []
    if length > windowsize:
        windows = 0
        for seq_in_win in rolling_window(seq, windowsize, windowstep):
            score = _calculate_rawscore(seq_in_win)
            dustscores.append(score / (windowsize - 2))
            windows += 1
        remaining_seq = seq[windows * windowstep:]
    else:
        remaining_seq = seq

    if remaining_seq > 5:
        length = len(remaining_seq)
        score = _calculate_rawscore(remaining_seq)
        dustscore = score / (length - 3) * (windowsize - 2) / (length - 2)
        dustscores.append(dustscore)

    # max score should be 100 not 31
    dustscore = sum(dustscores) / len(dustscores) * 100 / 31
    return dustscore


def calculate_nx(int_counter, percentage):
    '''It calcalutes N50, N90 etc.

    Modified from wikipedia:
    Given a set of sequences of varying lengths, the N50 length is defined as
    the length N for which half of all bases in the sequences are in a sequence
    of length L >= N
    '''
    if not int_counter:
        return None
    total_len = int_counter.sum
    len_desired = int(total_len * percentage / 100)
    accum_len = 0
    for length, n_seqs in sorted(int_counter.viewitems(),
                                 key=operator.itemgetter(0), reverse=True):
        accum_len += length * n_seqs
        if accum_len >= len_desired:
            return length


def calculate_sequence_stats(seqs, kmer_size=None, do_dust_stats=False,
                             nxs=None):
    'It calculates some stats for the given seqs.'
    # get data
    lengths = IntCounter()
    quals_per_pos = IntBoxplot()
    nucl_freq = NuclFreqsPlot()
    kmer_counter = KmerCounter(kmer_size) if kmer_size else None
    dustscores = IntCounter()
    for seq in seqs:
        lengths[get_length(seq)] += 1
        try:
            quals = get_int_qualities(seq)
        except AttributeError:
            quals = []
        for index, qual in enumerate(quals):
            quals_per_pos.append(index + 1, qual)
        str_seq = get_str_seq(seq)
        for index, nucl in enumerate(str_seq):
            nucl_freq.append(index, nucl)
        if kmer_counter is not None:
            kmer_counter.count_seq(str_seq)
        if do_dust_stats:
            dustscore = calculate_dust_score(seq)
            if dustscore is not None:
                dustscores[int(dustscore)] += 1

    lengths.update_labels({'sum': 'tot. residues', 'items': 'num. seqs.'})

    # length distribution
    lengths_srt = 'Length stats and distribution.\n'
    lengths_srt += '------------------------------\n'
    nxs = sorted(nxs) if nxs else []
    for nx in sorted(nxs):
        lengths_srt += 'N{:d}: {:d}\n'.format(nx, calculate_nx(lengths, nx))
    lengths_srt += str(lengths)
    lengths_srt += '\n'

    # agregate quals
    if quals_per_pos:
        quals = quals_per_pos.aggregated_array
        quals.update_labels({'sum': None, 'items': 'tot. base pairs'})

        q30 = quals.count_relative_to_value(30, operator.ge) / quals.count
        q30 *= 100

        q20 = quals.count_relative_to_value(20, operator.ge) / quals.count
        q20 *= 100

        # qual distribution
        qual_str = 'Quality stats and distribution.\n'
        qual_str += '-------------------------------\n'
        qual_str += 'Q20: {:.2f}\n'.format(q20)
        qual_str += 'Q30: {:.2f}\n'.format(q30)
        qual_str += str(quals)
        qual_str += '\n'

        # qual per position boxplot
        qual_boxplot = 'Boxplot for quality per position.\n'
        qual_boxplot += '---------------------------------\n'
        qual_boxplot += quals_per_pos.ascii_plot
        qual_boxplot += '\n'
    else:
        qual_str = ''
        qual_boxplot = ''

    # nucl freqs
    freq_str = 'Nucleotide frequency per position.\n'
    freq_str += '----------------------------------\n'
    freq_str += nucl_freq.ascii_plot
    freq_str += '\n'

    # kmer_distriubution
    kmer_str = ''
    if kmer_counter is not None:
        kmers = IntCounter(kmer_counter.values)
        if kmers:
            kmers.update_labels({'sum': None, 'items': 'num. kmers'})
            kmer_str = 'Kmer distribution\n'
            kmer_str += '-----------------\n'
            kmer_str += str(kmers)
            kmer_str += '\n'
            kmer_str += 'Most common kmers:\n'
            for kmer, number in kmer_counter.most_common(20):
                kmer_str += '\t{}: {}\n'.format(kmer, number)

    dust_str = ''
    if dustscores:
        dustscores.update_labels({'sum': None, 'items': 'num. seqs.'})
        dust_str = 'Dustscores stats and distribution.\n'
        dust_str += '----------------------------------\n'
        dust7 = (dustscores.count_relative_to_value(7, operator.gt) /
                 dustscores.count)
        dust_str += '% above 7 (low complexity): {:.2f}\n'.format(dust7)
        dust_str += str(dustscores)
        dust_str += '\n'

    return {'length': lengths_srt,
            'quality': qual_str,
            'nucl_freq': freq_str,
            'qual_boxplot': qual_boxplot,
            'kmer': kmer_str,
            'dustscore': dust_str}


def count_seqs(seqs):
    'It counts the number of sequences and the total length.'
    num_seqs = 0
    total_len = 0
    for seq in seqs:
        total_len += get_length(seq)
        num_seqs += 1

    return {'num_seqs': num_seqs,
            'total_length': total_len}


class BestItemsKeeper(object):
    '''It keeps a sorted list with the best items added.

    The interface is similar to a list, but with add and update
    '''
    def __init__(self, num_items, initializer=None, key=None, reverse=False):
        '''The initiator

        num_items: The number of items to keep (int)
        initializer: An iterator to initiate the list of best items
        key: Analogous to the key in sorted
        '''
        self._num_items = num_items
        if key is None:
            key = lambda x: x
        self._key = key
        self._best_items = []
        self._best0_key = None
        self._reverse = reverse
        if initializer is not None:
            self.update(initializer)

    def add(self, item):
        best = self._best_items
        key = self._key

        if len(best) >= self._num_items:
            remove_last_item = self._best0_key < key(item)
            if self._reverse:
                remove_last_item = not(remove_last_item)
            if remove_last_item:
                best.pop()
            else:
                return
        self._insort(item)

    def update(self, items):
        for item in items:
            self.add(item)

    def _insort(self, item):
        'Insert item inside the sorted list.'

        items = self._best_items
        key = self._key
        reverse = self._reverse
        # range in which to look for the insertion place
        low, high = 0, len(items)

        item_val = key(item)

        while low < high:
            mid = (low + high) // 2
            insert = item_val > key(items[mid])
            if reverse:
                insert = not(insert)
            if insert:
                high = mid
            else:
                low = mid + 1
        items.insert(low, item)
        self._best0_key = key(items[0])

    def __getitem__(self, index):
        return self._best_items[index]

    def __eq__(self, other):
        return self._best_items == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return str(self._best_items)
