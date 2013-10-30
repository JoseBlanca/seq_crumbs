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

from crumbs.seq import (get_str_seq, get_title, get_str_qualities,
                        get_file_format)
from crumbs.pairs import group_seqs_in_pairs
from crumbs.seqio import read_seqs
from crumbs.utils.tags import SEQITEM
from crumbs.iterutils import group_in_packets_fill_last
from crumbs.seq import SeqItem, SeqWrapper


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
    else:
        for read1, read2 in zip(pair1, pair2):
            if not get_str_seq(read1) == get_str_seq(read2):
                return False
        return True


def _seq_to_tabbed_str(read):
    if 'fastq' in get_file_format(read):
        title = '@' + get_title(read)
        title = title.replace('\t', ' ')
        sequence = get_str_seq(read)
        quality = get_str_qualities(read)
        return '\t'.join([title, sequence, quality])
    elif 'fasta' in get_file_format(read):
        title = '>' + get_title(read)
        title = title.replace('\t', ' ')
        sequence = get_str_seq(read)
        return '\t'.join([title, sequence])
    else:
        raise ValueError('Format not supported')


def _tabbed_pair_to_seqs_seqitem(pair_line, file_format):
    pair_items = pair_line.rstrip().split('\t')
    seqs = []
    if 'fasta' in file_format:
        for read_items in group_in_packets_fill_last(pair_items, 2):
            title = '>' + read_items[0][1:]
            lines = [title + '\n', read_items[1] + '\n']
            seqitem = SeqItem(title.split()[0][1:], lines)
            seq = SeqWrapper(SEQITEM, seqitem, file_format)
            seqs.append(seq)
    elif 'fastq' in file_format:
        for read_items in group_in_packets_fill_last(pair_items, 3):
            title = read_items[0]
            lines = [read_items[0] + '\n', read_items[1] + '\n',
                     '+\n', read_items[2] + '\n']
            seqitem = SeqItem(title.split()[0][1:], lines)
            seq = SeqWrapper(SEQITEM, seqitem, file_format)
            seqs.append(seq)
    else:
        raise ValueError('Format not supported')
    return seqs


def _pair_to_tabbed_str(pair):
    reads = [_seq_to_tabbed_str(read) for read in pair]
    reads.append('\n')
    return '\t'.join(reads)


def _read_pairs(in_fhands, paired_reads):
    seqs = read_seqs(in_fhands, prefered_seq_classes=[SEQITEM])
    if paired_reads:
        pairs = group_seqs_in_pairs(seqs)
    else:
        pairs = ([seq] for seq in seqs)
    return pairs


def _convert_fastq_to_tabbed_pairs(in_fhands, paired_reads):
    'It converts fastq files to one line per pair format'
    pairs = _read_pairs(in_fhands, paired_reads)
    pairs_in_one_line = imap(_pair_to_tabbed_str, pairs)
    for pair in pairs_in_one_line:
        yield pair


CONVERT_FASTQ_TO_LINES_SCRIPT = '''
from sys import argv, stdout
from crumbs.bulk_filters import _convert_fastq_to_tabbed_pairs
in_fpaths = argv[1:]
in_fhands = [open(fpath) for fpath in in_fpaths]
_convert_fastq_to_tabbed_pairs(in_fhands, out_fhand=stdout)
'''
UNIQUE_SEQITEM_SCRIPT = '''
import sys
from crumbs.bulk_filters import _unique_and_to_pairs
from crumbs.seqio import write_seqs
in_fhand = sys.stdin
out_format = sys.argv[1]
filtered_pairs = _unique_and_to_pairs(in_fhand, out_format)
for pair in filtered_pairs:
    write_seqs(pair, sys.stdout)
sys.stdout.flush()
'''


def _unique_and_to_pairs(in_fhand, file_format):
    prev_pair = None
    for pair_line in in_fhand:
        if prev_pair == None:
            duplicated = False
        else:
            duplicated = _tabbed_pairs_equal(pair_line, prev_pair)
        if not duplicated:
            yield _tabbed_pair_to_seqs_seqitem(pair_line,
                                               file_format=file_format)
        prev_pair = pair_line


def filter_duplicates(in_fhands, out_fhand, paired_reads, out_format='fastq'):
    '''It filters exact duplicated sequences even with different
    qualities or names. The output is given in fastq format by default,
    but allows also fasta format'''
    if not in_fhands:
        raise ValueError('At least one input fhand is required')
    if paired_reads:
        keys = '3,7'
    else:
        keys = '3'
    sort = Popen(['sort', '-k', keys, '-T', '-u'],
                 stdin=PIPE, stdout=PIPE)

    uniq_and_to_fastq_script = NamedTemporaryFile()
    uniq_and_to_fastq_script.write(UNIQUE_SEQITEM_SCRIPT)
    uniq_and_to_fastq_script.flush()
    cmd = [sys.executable, uniq_and_to_fastq_script.name, str(out_format)]
    uniq_and_to_fastq = Popen(cmd, stdin=sort.stdout, stdout=out_fhand)

    for line_pair in _convert_fastq_to_tabbed_pairs(in_fhands, paired_reads):
        sort.stdin.write(line_pair)
    sort.stdin.close()
    sort.wait()
    uniq_and_to_fastq.wait()
    uniq_and_to_fastq_script.close()
    assert not sort.returncode
    assert not uniq_and_to_fastq.returncode
