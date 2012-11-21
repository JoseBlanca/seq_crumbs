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
from crumbs.seqio import write_seqrecords

from crumbs.settings import LINKERS
from crumbs.utils.tags import NUCL

# pylint: disable=R0903


class MatePairSplitter(object):
    'It splits the input sequences with the provided linkers.'

    def __init__(self, linkers=None):
        'The initiator'
        self.linkers = LINKERS if linkers is None else linkers

    def __call__(self, seqs):
        'It splits a list of sequences with the provided linkers'
        seq_fhand = write_seqrecords(seqs, file_format='fasta')
        seq_fhand.flush()

        min_identity = 87.0
        min_len = 17
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
        for seqrec in seqs:
            segments = matcher.get_matched_segments_for_read(seqrec.id)
            if segments is not None:
                split_seqs = self._split_by_mate_linker(seqrec, segments)
            else:
                split_seqs = [seqrec]
            for seq in split_seqs:
                new_seqs.append(seq)
        return new_seqs

    def _split_by_mate_linker(self, seqrec, (segments, is_partial)):
        'It splits the seqs using segments'

        if not segments:
            return [seqrec]

        elongated_match = is_partial
        if len(segments) == 1:
            segment_start = segments[0][0]
            segment_end = segments[0][1]
            if segment_start == 0:
                new_seqrec = seqrec[segment_end + 1:]
                return [new_seqrec]
            elif segment_end == len(seqrec) - 1:
                new_seqrec = seqrec[:segment_start]
                return [new_seqrec]
            else:
                new_seqrec1 = seqrec[:segment_start]
                new_seqrec2 = seqrec[segment_end + 1:]
                id_ = seqrec.id
                if elongated_match:
                    id_ = seqrec.id + '_pl'
                    new_seqrec1.id = id_ + '.part1'
                    new_seqrec2.id = id_ + '.part2'
                else:
                    new_seqrec1.id = id_ + r'\1'
                    new_seqrec2.id = id_ + r'\2'
                return [new_seqrec1, new_seqrec2]
        else:
            seqrecords = []
            counter = 1
            seq_start = 0
            for segment_start, segment_end in segments:
                if segment_start == 0:
                    continue
                seqrecord = seqrec[seq_start:segment_start]
                seqrecord.id = seqrec.id + '_mlc.part{0:d}'.format(counter)
                seqrecords.append(seqrecord)
                counter += 1
                seq_start = segment_end + 1
            else:
                if segment_end != len(seqrec) + 1:
                    seqrecord = seqrec[segment_end + 1:]
                    seqrecord.id = seqrec.id + '_mlc.part{0:d}'.format(counter)
                    seqrecords.append(seqrecord)
            return seqrecords
