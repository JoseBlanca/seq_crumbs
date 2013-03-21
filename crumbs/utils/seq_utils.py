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

import re
import itertools
from multiprocessing import Pool

from Bio.Seq import Seq

from crumbs.utils.tags import UPPERCASE, LOWERCASE, SWAPCASE
from crumbs.seq import get_description, get_name


# pylint: disable=R0903
# pylint: disable=C0111


def replace_seq_same_length(seqrecord, seq_str):
    'It replaces the str with another of equal length keeping the annots.'
    annots = seqrecord.letter_annotations
    seqrecord.letter_annotations = {}
    alphabet = seqrecord.seq.alphabet
    seqrecord.seq = Seq(seq_str, alphabet)
    seqrecord.letter_annotations = annots
    return seqrecord


def uppercase_length(string):
    'It returns the number of uppercase characters found in the string'
    return len(re.findall("[A-Z]", string))


def get_uppercase_segments(string):
    '''It detects the unmasked regions of a sequence

    It returns a list of (start, end) tuples'''
    start = 0
    for is_upper, group in itertools.groupby(string, lambda x: x.isupper()):
        group = list(group)
        end = start + len(group) - 1
        if is_upper:
            yield start, end
        start = end + 1


class ChangeCase(object):
    'It changes the sequence case.'

    def __init__(self, action):
        'The initiator'
        if action not in (UPPERCASE, LOWERCASE, SWAPCASE):
            msg = 'Action should be: uppercase, lowercase or invertcase'
            raise ValueError(msg)
        self.action = action

    def __call__(self, seqrecords):
        'It changes the case of the seqrecords.'
        action = self.action
        processed_seqs = []
        for seqrecord in seqrecords:
            str_seq = str(seqrecord.seq)
            if action == UPPERCASE:
                str_seq = str_seq.upper()
            elif action == LOWERCASE:
                str_seq = str_seq.lower()
            elif action == SWAPCASE:
                str_seq = str_seq.swapcase()
            else:
                raise NotImplementedError()
            seqrecord = replace_seq_same_length(seqrecord, str_seq)
            processed_seqs.append(seqrecord)
        return processed_seqs


def append_to_description(seqrecord, text):
    'it appends the text to the seqrecord description'
    desc = get_description(seqrecord)
    if desc in (None, get_name(seqrecord), '<unknown description>'):
        desc = ''
    desc += text
    seqrecord.object.description = desc


class _FunctionRunner(object):
    'a class to join all the mapper functions in a single function'
    def __init__(self, map_functions):
        'Class initiator'
        self.map_functions = map_functions

    def __call__(self, seq_packet):
        'It runs all the map_functions for each seq_packet '
        processed_packet = seq_packet
        for map_function in self.map_functions:
            processed_packet = map_function(processed_packet)
        return processed_packet


def process_seq_packets(seq_packets, map_functions, processes=1,
                        keep_order=True):
    'It processes the SeqRecord packets'
    if processes > 1:
        workers = Pool(processes=processes)
        mapper = workers.imap if keep_order else workers.imap_unordered
    else:
        workers = None
        mapper = itertools.imap
    run_functions = _FunctionRunner(map_functions)

    seq_packets = mapper(run_functions, seq_packets)

    return seq_packets, workers
