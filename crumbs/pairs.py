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

from Bio.SeqIO import _index

from crumbs.exceptions import (MaxNumReadsInMem, PairDirectionError,
                               InterleaveError)
from crumbs.utils.tags import FWD, REV
from crumbs.seqio import write_seqrecords, _remove_one_line
from crumbs.settings import get_setting
from crumbs.third_party.index import FastqRandomAccess, index
from crumbs.utils.seq_utils import guess_format, get_title


def _parse_pair_direction_and_name(seq):
    'It parses the description field to get the name and the pair direction'
    return _parse_pair_direction_and_name_from_title(get_title(seq))


def _parse_pair_direction_and_name_from_title(title):
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
        file_format = guess_format(open(fpath))

    file_format = _remove_one_line(file_format)

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


def _get_paired_and_orphan(index_):
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


def match_pairs_unordered(seq_fpath, out_fhand, orphan_out_fhand, out_format):
    'It matches the seq pairs in an iterator and splits the orphan seqs'
    index_ = _index_seq_file(seq_fpath)
    paired, orphans = _get_paired_and_orphan(index_)

    # write paired
    write_seqrecords((index_[title] for title in paired), out_fhand,
                     out_format)

    # orphans
    write_seqrecords((index_[title] for title in orphans), orphan_out_fhand,
                     out_format)


def match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                memory_limit=get_setting('DEFAULT_SEQS_IN_MEM_LIMIT')):
    'It matches the seq pairs in an iterator and splits the orphan seqs'
    buf_fwd = {'index': {}, 'items': []}
    buf_rev = {'index': {}, 'items': []}
    buf1, buf2 = buf_fwd, buf_rev   # for the all orphan case
    for seq in seqs:
        try:
            seq_name, direction = _parse_pair_direction_and_name(seq)
        except PairDirectionError:
            write_seqrecords([seq], orphan_out_fhand, out_format)
            continue

        if direction == FWD:
            buf1 = buf_rev
            buf2 = buf_fwd
        else:
            buf1 = buf_fwd
            buf2 = buf_rev

        try:
            matching_seq_index = buf1['index'][seq_name]
        except KeyError:
            matching_seq_index = None

        if matching_seq_index is None:
            # add to buff
            buf2['items'].append(seq)
            buf2['index'][seq_name] = len(buf2['items']) - 1
            # check mem limit
            sum_items = len(buf1['items'] + buf2['items'])
            if memory_limit is not None and sum_items >= memory_limit:
                error_msg = 'There are too many consecutive non matching seqs'
                error_msg += ' in your input. We have reached the memory limit'
                raise MaxNumReadsInMem(error_msg)
        else:
            # write seqs from buffer1
            orphan_seqs = buf1['items'][:matching_seq_index]
            matching_seq = buf1['items'][matching_seq_index]
            write_seqrecords(orphan_seqs, orphan_out_fhand, out_format)
            write_seqrecords([matching_seq, seq], out_fhand, out_format)
            # fix buffers 1
            buf1['items'] = buf1['items'][matching_seq_index + 1:]
            buf1['index'] = {s: i for i, s in enumerate(buf1['items'])}

            # writes seqs from buffer 2 and fix buffer2
            write_seqrecords(buf2['items'], orphan_out_fhand, out_format)
            buf2['items'] = []
            buf2['index'] = {}
    else:
        orphan_seqs = buf1['items'] + buf2['items']
        write_seqrecords(orphan_seqs, orphan_out_fhand, out_format)

    orphan_out_fhand.flush()
    out_fhand.flush()


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
        write_seqrecords([seq1], out_fhand1, out_format)
        write_seqrecords([seq2], out_fhand2, out_format)
    out_fhand1.flush()
    out_fhand2.flush()
