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


from crumbs.settings import MAX_WIDTH_ASCII_PLOT


def draw_histogram(bin_limits, counts):
    'It draws an ASCII histogram'

    fill_char = '*'

    assert len(bin_limits) == len(counts) + 1

    # pylint: disable=W0108
    number_to_str = lambda n: '{:d}'.format(n)

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

    limit_fmt = '{:>' + str(max_ndigits) + 'd}'
    limit_to_padded_str = lambda n: limit_fmt.format(n)

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
    max_hist_width = MAX_WIDTH_ASCII_PLOT - max_header_len
    counts_ratio = max_hist_width / max_count

    result2 = []
    for line, cnt in zip(result, counts):
        line += fill_char * int(cnt * counts_ratio)
        line += '\n'
        result2.append(line)

    return ''.join(result2)
