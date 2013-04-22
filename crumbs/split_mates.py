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

# similar software SffToCA

from crumbs.blast import BlasterForFewSubjects
from crumbs.seqio import write_seqs
from crumbs.settings import get_setting
from crumbs.utils.tags import NUCL, SEQITEM
from crumbs.seq import (assing_kind_to_seqs, get_name, slice_seq, get_length,
                        copy_seq, SeqItem)

# pylint: disable=R0903


class MatePairSplitter(object):
    'It splits the input sequences with the provided linkers.'

    def __init__(self, linkers=None):
        'The initiator'
        if linkers is None:
            linkers = get_setting('LINKERS')
            linkers = [SeqItem(str(i), '>%d\n%s\n' % (i, l)) for i, l in enumerate(linkers)]
            linkers = assing_kind_to_seqs(SEQITEM, linkers, 'fasta')
        self.linkers = list(linkers)

    def __call__(self, seqs):
        'It splits a list of sequences with the provided linkers'
        seq_fhand = write_seqs(seqs, file_format='fasta')
        seq_fhand.flush()

        min_identity = 87.0
        min_len = 13
        filters = [{'kind': 'min_length', 'min_num_residues': min_len,
                    'length_in_query': False, 'filter_match_parts': True},
                   {'kind': 'score_threshold', 'score_key': 'identity',
                   'min_score': min_identity}]

        matcher = BlasterForFewSubjects(seq_fhand.name, self.linkers,
                                        program='blastn', filters=filters,
                                        params={'task': 'blastn-short'},
                                        elongate_for_global=True,
                                        seqs_type=NUCL)
        new_seqs = []
        for seq in seqs:
            segments = matcher.get_matched_segments_for_read(get_name(seq))
            if segments is not None:
                split_seqs = self._split_by_mate_linker(seq, segments)
            else:
                split_seqs = [seq]
            for seq in split_seqs:
                new_seqs.append(seq)
        return new_seqs

    def _split_by_mate_linker(self, seq, (segments, is_partial)):
        'It splits the seqs using segments'

        if not segments:
            return [copy_seq(seq)]

        elongated_match = is_partial
        if len(segments) == 1:
            segment_start = segments[0][0]
            segment_end = segments[0][1]
            if segment_start == 0:
                new_seq = slice_seq(seq, segment_end + 1, None)
                return [new_seq]
            elif segment_end == get_length(seq):
                new_seq = slice_seq(seq, None, segment_start)
                return [new_seq]
            else:
                new_seq1 = slice_seq(seq, None, segment_start)
                new_seq2 = slice_seq(seq, segment_end + 1, None)
                if elongated_match:
                    name = get_name(seq) + '_pl'
                    new_seq1 = copy_seq(new_seq1, name=name + '.part1')
                    new_seq2 = copy_seq(new_seq2, name=name + '.part2')
                else:
                    name = get_name(seq)
                    new_seq1 = copy_seq(new_seq1, name=name + r'\1')
                    new_seq2 = copy_seq(new_seq2, name=name + r'\2')
                return [new_seq1, new_seq2]
        else:
            seqs = []
            counter = 1
            seq_start = 0
            for segment_start, segment_end in segments:
                if segment_start == 0:
                    continue
                new_seq = slice_seq(seq, seq_start, segment_start)
                seq_name = get_name(seq) + '_mlc.part{0:d}'.format(counter)
                new_seq = copy_seq(new_seq, name=seq_name)
                seqs.append(new_seq)
                counter += 1
                seq_start = segment_end + 1
            else:
                if segment_end != get_length(seq) + 1:
                    new_seq = slice_seq(seq, segment_end + 1, None)
                    name = get_name(seq) + '_mlc.part{0:d}'.format(counter)
                    new_seq = copy_seq(new_seq, name=name)
                    seqs.append(new_seq)
            return seqs
