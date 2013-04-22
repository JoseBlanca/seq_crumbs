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

from __future__ import division
from array import array

from crumbs.utils.optional_modules import SffIterator

# pylint: disable=R0913


def _min_left_clipped_seqs(sff_fhand, trim, min_left_clip):
    'It generates sequences (as tuples) given a path to a SFF file.'
    for record in SffIterator(sff_fhand, trim=False):
        annots = record.annotations
        clip_qual = annots['clip_qual_left']
        clip_adapt = annots['clip_adapter_left']
        clip = max(min_left_clip, clip_qual, clip_adapt)
        seq = record.seq
        if trim:
            record.annotations = {}
            record = record[clip:]
        else:
            annots['clip_qual_left'] = clip
            annots['clip_adapter_left'] = clip
            seq = seq[:clip].lower() + seq[clip:].upper()
            quals = record.letter_annotations['phred_quality']
            record.letter_annotations = {}
            record.seq = seq
            dict.__setitem__(record._per_letter_annotations,
                             "phred_quality", quals)
        yield record


class SffExtractor(object):
    'This class extracts the reads from an SFF file'
    def __init__(self, sff_fhands, trim=False, min_left_clip=0,
                 nucls_to_check=50, max_nucl_freq_threshold=0.5):
        'It inits the class'
        self.fhands = sff_fhands
        self.trim = trim
        self.min_left_clip = min_left_clip

        # checking
        self.nucls_to_check = nucls_to_check
        self.max_nucl_freq_threshold = max_nucl_freq_threshold
        self.nucl_counts = {}

    @property
    def seqs(self):
        'It yields all sequences'
        for fhand in self.fhands:
            self._prepare_nucl_counts(fhand.name)
            if not self.min_left_clip:
                seqs = SffIterator(fhand, trim=self.trim)
            else:
                seqs = _min_left_clipped_seqs(fhand, self.trim,
                                              self.min_left_clip)
            for record in seqs:
                self._update_nucl_counts(str(record.seq), fhand.name)
                yield record

    def _prepare_nucl_counts(self, fpath):
        'It prepares the structure to store the nucleotide counts'
        counts = {'A': array('L', [0] * self.nucls_to_check),
                  'T': array('L', [0] * self.nucls_to_check),
                  'C': array('L', [0] * self.nucls_to_check),
                  'G': array('L', [0] * self.nucls_to_check)}
        self.nucl_counts[fpath] = counts

    def _update_nucl_counts(self, seq, fpath):
        'Given a seq (as a string) it updates the nucleotide counts'
        seq = seq[:self.nucls_to_check]
        counts = self.nucl_counts
        for index, nucl in enumerate(seq):
            try:
                counts[fpath][nucl][index] += 1
            except KeyError:
                pass  # we do not count the lowercase letters

    @property
    def clip_advice(self):
        'It checks how many positions have a high max nucl freq.'
        advices = {}
        for fhand in self.fhands:
            fpath = fhand.name
            counts = self.nucl_counts[fpath]
            treshold = self.max_nucl_freq_threshold
            pos_above_threshold = 0
            seq_above_threshold = ''
            index = 0
            for index in range(self.nucls_to_check):
                num_nucls = [counts['A'][index], counts['T'][index],
                             counts['C'][index], counts['G'][index]]
                tot_nucls = sum(num_nucls)
                if not tot_nucls:
                    continue
                freq_nucls = [i / tot_nucls for i in num_nucls]
                above_threshold = [i >= treshold for i in freq_nucls]
                if any(above_threshold):
                    pos_above_threshold += 1
                    seq_above_threshold += _get_nucl_with_max_freq('ATCG',
                                                                   freq_nucls)
                else:
                    break
            if pos_above_threshold:
                if self.trim:
                    # number of nucleotides to remove next time, the ones
                    # that we have detected plus the ones already removed
                    advice = index + self.min_left_clip, seq_above_threshold
                else:
                    advice = index, seq_above_threshold
            else:
                advice = None
            advices[fpath] = advice
        return advices


def _get_nucl_with_max_freq(nucls, freq_nucls):
    'It returns the nucleotide with the maximum frequency'
    max_ = None
    for index, freq in enumerate(freq_nucls):
        if max_ is None or max_ < freq:
            max_ = freq
            nucl = nucls[index]
    return nucl
