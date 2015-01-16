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

import re
from itertools import izip_longest, chain

from toolz import first

from crumbs.exceptions import (PairDirectionError, InterleaveError,
                               ItemsNotSortedError)
from crumbs.seq.seqio import write_seqs
from crumbs.seq.seq import get_title, get_name
from crumbs.utils.tags import FWD, REV
from crumbs.utils.file_utils import flush_fhand
from crumbs.iterutils import sorted_items, group_in_packets_fill_last
from crumbs.collectionz import KeyedSet


def _parse_pair_direction_and_name(seq):
    'It parses the description field to get the name and the pair direction'
    return _parse_pair_direction_and_name_from_title(get_title(seq))


def _parse_pair_direction_and_name_from_title(title):
    'It guesses the direction from the title line'
    reg_exps = ['(.+)[/|\\\\](\d+)', '(.+)\s(\d+):.+', '(.+)\.(\w)\s?']
    for reg_exp in reg_exps:
        match = re.match(reg_exp, title)
        if match:
            name, direction = match.groups()
            if direction in ('1', 'f'):
                direction = FWD
            elif direction in ('2', 'r'):
                direction = REV
            else:
                raise PairDirectionError('unknown direction descriptor')
            return name, direction
    raise PairDirectionError('Unable to detect the direction of the seq')


def _get_paired_and_orphan(reads, ordered, max_reads_memory, temp_dir):
    if ordered:
        sorted_reads = reads
    else:
        def _key(seq):
            return get_title(seq)
        sorted_reads = sorted_items(reads, _key, max_reads_memory, temp_dir)
    return group_pairs_by_name(sorted_reads)


def match_pairs(reads, out_fhand, orphan_out_fhand, out_format, ordered=True,
                check_order_buffer_size=0, max_reads_memory=None,
                temp_dir=None):
    '''It matches the seq pairs in an iterator and splits the orphan seqs.'''
    counts = 0
    check_order_buffer = KeyedSet()
    for pair in _get_paired_and_orphan(reads, ordered, max_reads_memory,
                                       temp_dir):
        if len(pair) == 1:
            write_seqs(pair, orphan_out_fhand, out_format)
            try:
                name = _parse_pair_direction_and_name(pair[0])[0]
            except PairDirectionError:
                name = get_name(pair[0])
            if ordered and counts < check_order_buffer_size:
                counts += 1
                if not check_order_buffer.check_add(name):
                    msg = 'Reads are not ordered by pairs.Use unordered option'
                    raise ItemsNotSortedError(msg)
            elif ordered and counts >= check_order_buffer_size:
                if name in check_order_buffer:
                    msg = 'Reads are not ordered by pairs.Use unordered option'
                    raise ItemsNotSortedError(msg)
        elif len(pair) == 2:
            write_seqs(pair, out_fhand, out_format)
    flush_fhand(orphan_out_fhand)
    flush_fhand(out_fhand)


def _check_name_and_direction_match(*seqs):
    'It fails if the names do not match or if the directions are equal'
    n_seqs = len(seqs)
    names = set()
    directions = set()
    for seq in seqs:
        name, direction = _parse_pair_direction_and_name(seq)
        names.add(name)
        directions.add(direction)

    if len(names) > 1:
        msg = 'The read names from a pair do not match: %s'
        msg %= ','.join(names)
        raise InterleaveError(msg)
    if len(directions) != n_seqs:
        msg = 'A pair has repeated directions: ' + first(names)
        raise InterleaveError(msg)


def interleave_pairs(seqs1, seqs2, skip_checks=False):
    '''A generator that interleaves the paired reads found in two iterators.

    It will fail if forward and reverse reads do not match in both sequence
    iterators.
    '''
    for seq1, seq2 in izip_longest(seqs1, seqs2, fillvalue=None):
        if not skip_checks:
            if seq1 is None or seq2 is None:
                msg = 'The files had a different number of sequences'
                raise InterleaveError(msg)
            _check_name_and_direction_match(seq1, seq2)
        if seq1 is not None:
            yield seq1
        if seq2 is not None:
            yield seq2


def deinterleave_pairs(seqs, out_fhand1, out_fhand2, out_format):
    '''It splits a sequence iterator with alternating paired reads in two.

    It will fail if forward and reverse reads are not alternating.
    '''
    for pair in group_pairs(seqs, n_seqs_in_pair=2):
        write_seqs((pair[0],), out_fhand1, out_format)
        write_seqs((pair[1],), out_fhand2, out_format)
    out_fhand1.flush()
    out_fhand2.flush()


def group_pairs_by_name(seqs, all_pairs_same_n_seqs=False):
    paired_seqs = []
    prev_name = None
    n_seqs_per_pair = None
    for seq in iter(seqs):
        try:
            name = _parse_pair_direction_and_name(seq)[0]
        except PairDirectionError:
            name = None
        if name is None or (paired_seqs and name != prev_name):
            if all_pairs_same_n_seqs:
                if n_seqs_per_pair is None:
                    n_seqs_per_pair = len(paired_seqs)
                elif n_seqs_per_pair != len(paired_seqs):
                    msg = 'Pair had different number of reads: '
                    msg += prev_name
                    raise InterleaveError(msg)
            if paired_seqs:
                yield paired_seqs
                paired_seqs = []
        paired_seqs.append(seq)
        prev_name = name
    if paired_seqs:
        if all_pairs_same_n_seqs and n_seqs_per_pair is not None:
            if n_seqs_per_pair != len(paired_seqs):
                msg = 'Pair had different number of reads: '
                msg += prev_name
                raise InterleaveError(msg)
        yield paired_seqs


def _get_first_pair_by_name(seqs):
    paired_seqs = []
    prev_name = None
    for seq in iter(seqs):
        name = _parse_pair_direction_and_name(seq)[0]
        if prev_name is None:
            prev_name = name
        if name != prev_name:
            return paired_seqs, seq
        else:
            paired_seqs.append(seq)
    return None, None


def group_pairs(seqs, n_seqs_in_pair=None, check_all_same_n_seqs=True,
                check_name_matches=True):

    seqs = iter(seqs)
    if n_seqs_in_pair is None:
        first_pair, next_read = _get_first_pair_by_name(seqs)
        if first_pair is None:
            n_seqs_in_pair = None
        else:
            yield first_pair
            n_seqs_in_pair = len(first_pair)
            seqs = chain([next_read], seqs)

    if n_seqs_in_pair == 1:
        # No need to check anything, a pair cannot have less than one read
        # or more than one name
        check_all_same_n_seqs = False
        check_name_matches = False

    if n_seqs_in_pair:
        pairs = group_in_packets_fill_last(seqs, packet_size=n_seqs_in_pair)
        for pair in pairs:
            pair = filter(lambda seq: seq is not None, pair)
            if check_all_same_n_seqs and n_seqs_in_pair != len(pair):
                msg = 'The last pair has fewer reads'
                raise InterleaveError(msg)
            if check_name_matches:
                _check_name_and_direction_match(*pair)
            yield pair
