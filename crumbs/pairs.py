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
from itertools import izip_longest

from crumbs.utils.optional_modules import _index
from crumbs.exceptions import (MaxNumReadsInMem, PairDirectionError,
                               InterleaveError, MalformedFile,
                               ItemsNotSortedError)
from crumbs.seqio import write_seqs
from crumbs.settings import get_setting
from crumbs.third_party.index import FastqRandomAccess, index
from crumbs.seq import get_title, SeqWrapper, get_name
from crumbs.utils.file_formats import get_format, remove_multiline
from crumbs.utils.tags import FWD, REV, SEQRECORD
from crumbs.utils.file_utils import flush_fhand
from crumbs.iterutils import sorted_items
from crumbs.collectionz import OrderedSet


def _parse_pair_direction_and_name(seq):
    'It parses the description field to get the name and the pair direction'
    return _parse_pair_direction_and_name_from_title(get_title(seq))


def _parse_pair_direction_and_name_from_title(title):
    'It gueses the direction from the title line'
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


def _index_seq_file(fpath, file_format=None):
    '''It indexes a seq file using Biopython index.

    It uses the title line line as the key and not just the id.
    '''
    if file_format is None:
        file_format = get_format(open(fpath))

    file_format = remove_multiline(file_format)

    # pylint: disable W0212
    # we monkey patch to be able to index using the whole tile line and not
    # only the id. We need it because in a pair end file sequences with the
    # same id could be found
    accessor = _index._FormatToRandomAccess
    old_accessor = accessor.copy()
    accessor['fastq'] = FastqRandomAccess
    accessor['astq-sanger'] = FastqRandomAccess
    accessor['fastq-solexa'] = FastqRandomAccess
    accessor['fastq-illumina'] = FastqRandomAccess

    file_index = index(fpath, format=file_format)

    _index._FormatToRandomAccess = old_accessor

    return file_index


def _match_pairs_from_sorted_reads(sorted_reads):
    prev_reads = []
    prev_reads_name = None
    prev_reads_directions = []
    #{'read': None, 'name': None, 'direction': None}
    for read in sorted_reads:
        title = get_title(read)
        try:
            name, direction = _parse_pair_direction_and_name_from_title(title)
        except PairDirectionError:
            name, direction = None, None

        if direction is not None and prev_reads_name is None:
            prev_reads = [read]
            prev_reads_name = name
            prev_reads_directions = [direction]
            continue
        elif direction is None:
            if prev_reads:
                yield prev_reads
                prev_reads = []
                prev_reads_name = None
                prev_reads_directions = []
            yield [read]
        elif name == prev_reads_name:
            if direction == REV:
                prev_reads.append(read)
                prev_reads_directions.append(direction)
            else:
                prev_reads.insert(0, read)
                prev_reads_directions.insert(0, direction)
        else:
            if prev_reads:
                yield prev_reads
            prev_reads = [read]
            prev_reads_name = name
            prev_reads_directions = [direction]
    if prev_reads:
        yield prev_reads


def _get_paired_and_orphan(reads, ordered, max_reads_memory, temp_dir):
    if ordered:
        sorted_reads = reads
    else:
        def _key(seq):
            return get_title(seq)
        sorted_reads = sorted_items(reads, _key, max_reads_memory, temp_dir)
    return _match_pairs_from_sorted_reads(sorted_reads)


def match_pairs(reads, out_fhand, orphan_out_fhand, out_format, ordered=True,
                check_order_buffer_size=0, max_reads_memory=None,
                temp_dir=None):
    '''It matches the seq pairs in an iterator and splits the orphan seqs.
    It assumes that sequences are already sorted'''
    counts = 0
    check_order_buffer = OrderedSet()
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
    orphan_out_fhand.flush()
    out_fhand.flush()


def _get_paired_and_orphan_old(index_):
    'It guesses the paired and the orphan seqs'
    fwd_reads = {}
    rev_reads = {}

    orphan_titles = []
    for title in index_.iterkeys():
        try:
            name, direction = _parse_pair_direction_and_name_from_title(title)
        except PairDirectionError:
            orphan_titles.append(title)
            continue

        if direction in (FWD, REV):
            reads_info = fwd_reads if direction == FWD else rev_reads
            reads_info[name] = title
        else:
            raise RuntimeError('Unknown direction for read.')

    fwd_names = set(fwd_reads.viewkeys())
    rev_names = set(rev_reads.viewkeys())

    paired = fwd_names.intersection(rev_names)

    paired_titles = []
    for paired_name in paired:
        paired_titles.append(fwd_reads[paired_name])
        paired_titles.append(rev_reads[paired_name])

    fwd_orphans = fwd_names.difference(rev_names)
    orphan_titles.extend(fwd_reads[orphan] for orphan in fwd_orphans)
    del fwd_orphans

    rev_orphans = rev_names.difference(fwd_names)
    orphan_titles.extend(rev_reads[orphan] for orphan in rev_orphans)

    return paired_titles, orphan_titles


def match_pairs_unordered_old(seq_fpath, out_fhand, orphan_out_fhand, out_format):
    'It matches the seq pairs in an iterator and splits the orphan seqs'
    indx = _index_seq_file(seq_fpath)
    paired, orphans = _get_paired_and_orphan(indx)

    # write paired
    write_seqs((SeqWrapper(SEQRECORD, indx[title], None) for title in paired),
                out_fhand, out_format)

    # orphans
    write_seqs((SeqWrapper(SEQRECORD, indx[title], None) for title in orphans),
               orphan_out_fhand, out_format)


def match_pairs_old(seqs, out_fhand, orphan_out_fhand, out_format,
                memory_limit=get_setting('DEFAULT_SEQS_IN_MEM_LIMIT')):
    'It matches the seq pairs in an iterator and splits the orphan seqs'
    buf_fwd = {'index': {}, 'items': []}
    buf_rev = {'index': {}, 'items': []}
    buf1, buf2 = buf_rev, buf_fwd  # for the all orphan case
    for seq in seqs:
        try:
            seq_name, direction = _parse_pair_direction_and_name(seq)
        except PairDirectionError:
            write_seqs([seq], orphan_out_fhand, out_format)
            continue

        # buf1 -> buffer for the reads with the same orientation as the
        # current one
        # buf2 -> buffer for the reads with the reverse orientation as the
        # current one

        if direction == FWD:
            buf1 = buf_fwd
            buf2 = buf_rev
        else:
            buf1 = buf_rev
            buf2 = buf_fwd

        try:
            matching_seq_index = buf2['index'][seq_name]
        except KeyError:
            matching_seq_index = None

        if matching_seq_index is None:
            # add to buff
            buf1['items'].append(seq)
            buf1['index'][seq_name] = len(buf1['items']) - 1
            # check mem limit
            sum_items = len(buf2['items'] + buf1['items'])
            if memory_limit is not None and sum_items >= memory_limit:
                error_msg = 'There are too many consecutive non matching seqs'
                error_msg += ' in your input. We have reached the memory limit'
                error_msg += '. Are you sure that the reads are sorted and '
                error_msg += 'interleaved?. You could try with the unordered'
                error_msg += ' algorith'
                raise MaxNumReadsInMem(error_msg)
        else:
            # write seqs from buffer1
            orphan_seqs = buf2['items'][:matching_seq_index]
            matching_seq = buf2['items'][matching_seq_index]
            write_seqs(orphan_seqs, orphan_out_fhand, out_format)
            write_seqs([matching_seq, seq], out_fhand, out_format)
            # fix buffer 1
            if matching_seq_index != len(buf2['items']) - 1:
                msg = 'The given files are not sorted (ordered) and '
                msg = 'interleaved. You could try with the unordered algorithm'
                raise MalformedFile(msg)
            buf2 = {'index': {}, 'items': []}
            # writes seqs from buffer 2 and fix buffer2
            write_seqs(buf1['items'], orphan_out_fhand, out_format)
            buf1 = {'index': {}, 'items': []}

        if direction == FWD:
            buf_fwd = buf1
            buf_rev = buf2
        else:
            buf_rev = buf1
            buf_fwd = buf2

    else:
        orphan_seqs = buf1['items'] + buf2['items']
        write_seqs(orphan_seqs, orphan_out_fhand, out_format)

    orphan_out_fhand.flush()
    flush_fhand(out_fhand)


def _check_name_and_direction_match(seq1, seq2):
    'It fails if the names do not match or if the direction are equal'
    name1, direction1 = _parse_pair_direction_and_name(seq1)
    name2, direction2 = _parse_pair_direction_and_name(seq2)
    if name1 != name2:
        msg = 'The reads from the two files do not match: {}, {}'
        msg = msg.format(name1, name2)
        raise InterleaveError(msg)
    if direction1 == direction2:
        msg = 'Two paired reads have the same direction: {}, {}'
        msg = msg.format(name1 + ' ' + direction1,
                         name2 + ' ' + direction2)
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
    while True:
        try:
            seq1 = seqs.next()
        except StopIteration:
            seq1 = None
        try:
            seq2 = seqs.next()
        except StopIteration:
            seq2 = None
        if seq1 is None:
            break  # we have consumed the input iterator completely
        if seq2 is None:
            msg = 'The file had an odd number of sequences'
            raise InterleaveError(msg)
        _check_name_and_direction_match(seq1, seq2)
        write_seqs([seq1], out_fhand1, out_format)
        write_seqs([seq2], out_fhand2, out_format)
    out_fhand1.flush()
    out_fhand2.flush()


def group_seqs_in_pairs(seqs):
    '''It generates lists of paired reads.

    The paired reads should be interleaved.
    '''
    paired_seqs = []
    prev_name = None
    for seq in iter(seqs):
        name = _parse_pair_direction_and_name(seq)[0]
        if prev_name and name != prev_name:
            yield paired_seqs
            paired_seqs = []
        paired_seqs.append(seq)
        prev_name = name
    if paired_seqs:
        yield paired_seqs
