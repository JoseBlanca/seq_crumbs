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

from itertools import compress
import os.path

from crumbs.annotation import (PolyaAnnotator, EstscanOrfAnnotator,
                               BlastAnnotator)
from crumbs.utils.tags import REVERSE


class TranscriptOrientator(object):
    '''This class orientates the transcripts

    It can take into account: poly-A, ORFs and blast matches'''

    def __init__(self, polya_params=None, estscan_params=None,
                 blast_params=None):
        'Initiatilzes the class'

        self._annotators = self._create_pipeline(polya_params,
                                                   estscan_params,
                                                   blast_params)

    def _create_pipeline(self, polya_params, estscan_params, blast_params):
        'It creates the annotation pipeline'
        # pylint: disable=W0142
        annotators = []
        if polya_params:
            annotator = PolyaAnnotator(**polya_params)
            annotators.append({'annotator': annotator,
                               'feat_selector': self._polya_selector})
        if estscan_params:
            annotator = EstscanOrfAnnotator(**estscan_params)
            annotators.append({'annotator': annotator,
                               'feat_selector': self._orf_selector})

        if blast_params:
            for blast_param in blast_params:
                annotator = BlastAnnotator(**blast_param)
                annotators.append({'annotator': annotator,
                                   'feat_selector': self._match_part_selector})
        return annotators

    @staticmethod
    def _select_features_by_type(features, kind):
        'it selects features by type'
        return [feat for feat in features if feat.type == kind]

    def _polya_selector(self, features, blastdb=None):
        'It selects the polya'
        # pylint: disable=W0613
        feats = self._select_features_by_type(features, 'polyA_sequence')
        if feats:
            return feats[0]

    def _orf_selector(self, features, blastdb=None):
        'it returns the longest feature'
        # pylint: disable=W0613
        features = self._select_features_by_type(features, 'ORF')
        if not features:
            return None

        lengths = [len(feat) for feat in features]
        return features[lengths.index(max(lengths))]

    def _match_part_selector(self, features, blastdb):
        'it return the match_part with the best e-value'
        features = self._select_features_by_type(features, 'match_part')
        features = [f for f in features if f.qualifiers['blastdb'] == blastdb]

        if not features:
            return None
        scores = [feat.qualifiers['score'] for feat in features]
        return features[scores.index(min(scores))]

    @staticmethod
    def _guess_orientations(seqrecords, feat_selector, blastdb):
        '''It returns the orientation of the annotated transcripts.'''
        orientations = []
        for seqrecord in seqrecords:
            feature = feat_selector(seqrecord.features, blastdb=blastdb)
            orientation = None if feature is None else feature.strand
            orientations.append(orientation)
        return orientations

    def __call__(self, seqrecords):
        'It makes the job'
        orientations = None

        for annotator in self._annotators:
            if orientations:
                to_annalyze = [not o for o in orientations]
                seqrecords_to_annalyze = list(compress(seqrecords,
                                                       to_annalyze))
            else:
                orientations = [None] * len(seqrecords)
                seqrecords_to_annalyze = seqrecords

            feat_selector = annotator['feat_selector']
            annotator = annotator['annotator']
            try:
                blastdb = os.path.basename(annotator.blastdb)
            except AttributeError:
                blastdb = None
            annot_seqrecords = annotator(seqrecords_to_annalyze)
            annot_strands = self._guess_orientations(annot_seqrecords,
                                                     feat_selector,
                                                     blastdb=blastdb)

            analyzed_seqs_index = 0
            for index, orientation in enumerate(orientations):
                if orientation is None:
                    orientations[index] = annot_strands[analyzed_seqs_index]
                    analyzed_seqs_index += 1

        # Now we reverse the seqs that we have guess that are reversed
        reorientated_seqrecords = []
        for orientation, seqrecord in zip(orientations, seqrecords):
            if orientation == REVERSE:
                seqrecord = seqrecord.reverse_complement()
            reorientated_seqrecords.append(seqrecord)
        return reorientated_seqrecords
