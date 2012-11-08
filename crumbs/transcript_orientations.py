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
from random import choice
import os.path

from crumbs.annotation import (PolyaAnnotator, EstscanOrfAnnotator,
                               BlastAnnotator)
from crumbs.utils.tags import REVERSE, FORWARD


class TranscriptOrientator(object):
    '''This class orientates the transcript seqs looking into the polia, orf
    and blast matches against other databases. In this order'''

    def __init__(self, polya_params=None, estscan_params=None,
                 blast_params=None):
        'Initiatilzes the class'

        self._annotators = self._create_pipeline(polya_params,
                                                   estscan_params,
                                                   blast_params)
        # analises

    def _create_pipeline(self, polya_params, estscan_params, blast_params):
        'It creates the annotation pipeline'
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
        for feature in features:
            if feature.type == kind:
                yield feature

    def _polya_selector(self, features, db=None):
        'It selects the polya'
        feats = list(self._select_features_by_type(features, 'polyA_sequence'))
        if feats:
            return feats[0]

    def _orf_selector(self, features, db=None):
        'it returns the longest feature'
        features = list(self._select_features_by_type(features, 'ORF'))
        if not features:
            return None

        lengths = [len(feat) for feat in features]
        return features[lengths.index(max(lengths))]

    def _match_part_selector(self, features, db):
        'it return the match_part with the best e-value'
        features = list(self._select_features_by_type(features, 'match_part'))
        features = [f for f in features if f.qualifiers['blastdb'] == db]

        if not features:
            return None
        scores = [feat.qualifiers['score'] for feat in features]
        return features[scores.index(min(scores))]

    def _guess_orientations(self, seqrecords, feat_selector, db):
        '''It returns the orientation of the transcripts.'''
        orientations = []
        for seqrecord in seqrecords:
            feature = feat_selector(seqrecord.features, db=db)
            if feature is not None:
                orientation = feature.strand
            else:
                orientation = None
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
                db = os.path.basename(annotator.blastdb)
            except AttributeError:
                db = None
            annot_seqrecords = annotator(seqrecords_to_annalyze)
            annot_strands = self._guess_orientations(annot_seqrecords,
                                                     feat_selector, db=db)

            counter = 0
            for index, orientation in enumerate(orientations):
                if orientation is None:
                    orientations[index] = annot_strands[counter]
                    counter += 1
        reorientated_seqrecords = []
        for orientation, seqrecord in zip(orientations, seqrecords):
            if orientation == None:
                orientation = choice([REVERSE, FORWARD])
            if orientation == REVERSE:
                seqrecord = seqrecord.reverse_complement()
            reorientated_seqrecords.append(seqrecord)
        return reorientated_seqrecords

