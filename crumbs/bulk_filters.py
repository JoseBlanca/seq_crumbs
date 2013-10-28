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

from subprocess import Popen, PIPE
import sys
from itertools import imap
from tempfile import NamedTemporaryFile

from crumbs.seq import get_str_seq, get_title, get_str_qualities
from crumbs.pairs import group_seqs_in_pairs
from crumbs.seqio import read_seqs
from crumbs.utils.tags import SEQITEM, GUESS_FORMAT
from crumbs.utils.file_utils import flush_fhand
from crumbs.iterutils import group_in_packets_fill_last
from crumbs.seq import SeqItem, SeqWrapper
from crumbs.exceptions import UnknownFormatError


def _tabbed_pairs_equal(pair1, pair2):
    pair1 = pair1.split('\t')
    pair2 = pair2.split('\t')
    if len(pair1) != len(pair2):
        return False

    seqs_in_pair1 = pair1[1::3]
    seqs_in_pair2 = pair2[1::3]
    for read1, read2 in zip(seqs_in_pair1, seqs_in_pair2):
        if read1 != read2:
            return False
    return True


def _seqitem_pairs_equal(pair1, pair2):
    if len(pair1) != len(pair2):
        return False
    elif len(pair1) == 2:
        for read1, read2 in zip(pair1, pair2):
            if not get_str_seq(read1) == get_str_seq(read2):
                return False
        return True
    else:
        return get_str_seq(pair1) == get_str_seq(pair2)

def _seq_to_tabbed_str(read):
    title = '@' + get_title(read)
    title = title.replace('\t', ' ')
    sequence = get_str_seq(read)
    quality = get_str_qualities(read)
    return '\t'.join([title, sequence, quality])


def _tabbed_pair_to_seqs_seqitem(pair_line, file_format):
    pair_items = pair_line.rstrip().split('\t')

    if 'fastq' not in file_format:
        raise NotImplemented('Not implemented yet for fasta')
    seqs = []
    for read_items in group_in_packets_fill_last(pair_items, 3):
        title = read_items[0]
        lines = [read_items[0]+'\n', read_items[1]+'\n', '+\n', read_items[2]+'\n']
        seqitem = SeqItem(title.split()[0][1:], lines)
        seq = SeqWrapper(SEQITEM, seqitem, file_format)
        seqs.append(seq)
    return seqs


def _pair_to_tabbed_str(pair):
    reads = [_seq_to_tabbed_str(read) for read in pair]
    reads.append('\n')
    return '\t'.join(reads)


def _read_pairs(in_fhands):
    seqs = read_seqs(in_fhands, prefered_seq_classes=[SEQITEM])
    seqs = list(seqs)
    pairs = group_seqs_in_pairs(seqs)
    return pairs


def _convert_fastq_to_tabbed_pairs(in_fhands, out_fhand):
    'It converts fastq files to one line per pair format'
    pairs = _read_pairs(in_fhands)
    pairs_in_one_line = imap(_pair_to_tabbed_str, pairs)
    for pair in pairs_in_one_line:
        try:
            out_fhand.write(pair)
        except IOError, error:
            # The pipe could be already closed
            if not 'Broken pipe' in str(error):
                raise
    flush_fhand(out_fhand)


CONVERT_FASTQ_TO_LINES_SCRIPT = '''
from sys import argv, stdout
from crumbs.bulk_filters import _convert_fastq_to_tabbed_pairs
in_fpaths = argv[1:]
in_fhands = [open(fpath) for fpath in in_fpaths]
_convert_fastq_to_tabbed_pairs(in_fhands, out_fhand=stdout)
'''


def filter_duplicates(in_fpaths, in_format=GUESS_FORMAT):
    if not in_fpaths:
        raise ValueError('At least one input fpath is required')
    for fpath in in_fpaths:
        fhand = open(fpath)
        if fhand.next() == '':
            raise UnknownFormatError('The file is empty')
        fhand.close()
    file_format = in_format
    to_lines_script = NamedTemporaryFile(prefix='fastq_to_lines.',
                                         suffix='.py')
    to_lines_script.write(CONVERT_FASTQ_TO_LINES_SCRIPT)
    to_lines_script.flush()
    to_lines_cmd = [sys.executable, to_lines_script.name]
    to_lines_cmd.extend(in_fpaths)

    to_lines = Popen(to_lines_cmd, stdout=PIPE)
    sort = Popen(['sort', '-k', '2,5', '-T', '-u'],
                 stdin=to_lines.stdout, stdout=PIPE)

    prev_pair = None
    for pair_line in sort.stdout:
        if prev_pair == None:
            duplicated = False
        else:
            duplicated = _tabbed_pairs_equal(pair_line, prev_pair)
        if not duplicated:
            yield _tabbed_pair_to_seqs_seqitem(pair_line,
                                               file_format=file_format)
        prev_pair = pair_line

    to_lines_script.close()
    #print to_lines.returncode
    #print sort.returncode
