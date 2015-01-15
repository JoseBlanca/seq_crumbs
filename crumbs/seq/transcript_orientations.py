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
from crumbs.utils.seq_utils import append_to_description
from crumbs.utils.tags import SEQRECORD
from crumbs.seq import SeqWrapper


class TranscriptOrientator(object):
    '''This class orientates the transcripts

    It can take into account: poly-A, ORFs and blast matches'''

    def __init__(self, polya_params=None, estscan_params=None,
                 blast_params=None):
        self._polya_params = polya_params
        self._estscan_params = estscan_params
        self._blast_params = blast_params

        self._annotators = self._create_pipeline()

    def _create_pipeline(self):
        'It creates the annotation pipeline'
        # pylint: disable=W0142
        annotators = []
        if self._polya_params:
            annotators.append({'name': 'polyA'})
        if self._estscan_params:
            annotators.append({'name': 'estscan_orf'})
        if self._blast_params:
            for blast_param in self._blast_params:
                annotators.append({'name': 'blast',
                                   'blastdb': blast_param['blastdb']})
        return annotators

    @staticmethod
    def _select_features_by_type(features, kind):
        'it selects features by type'
        return [feat for feat in features if feat.type == kind]

    def _polya_selector(self, features):
        'It selects the polya'
        # pylint: disable=W0613
        feats = self._select_features_by_type(features, 'polyA_sequence')
        if feats:
            return feats[0]

    def _orf_selector(self, features):
        'it returns the longest feature'
        # pylint: disable=W0613
        features = self._select_features_by_type(features, 'ORF')
        if not features:
            return None

        lengths = [len(feat) for feat in features]
        return features[lengths.index(max(lengths))]

    def _match_part_selector(self, features, blastdb):
        'it return the match_part with the best e-value'
        blastdb = os.path.basename(blastdb)
        features = self._select_features_by_type(features, 'match_part')
        features = [f for f in features if f.qualifiers['blastdb'] == blastdb]

        if not features:
            return None
        scores = [feat.qualifiers['score'] for feat in features]
        return features[scores.index(min(scores))]

    def _guess_orientations(self, seqs, annotator_name, blastdb):
        '''It returns the orientation of the annotated transcripts.'''
        orientations = []
        for seq in seqs:
            if annotator_name == 'polyA':
                feature = self._polya_selector(seq.object.features)
            elif annotator_name == 'estscan_orf':
                feature = self._orf_selector(seq.object.features)
            elif  annotator_name == 'blast':
                feature = self._match_part_selector(seq.object.features,
                                                    blastdb=blastdb)
            else:
                raise NotImplementedError('This annotator type not supported')

            orientation = None if feature is None else feature.strand
            orientations.append(orientation)
        return orientations

    def _get_annotator(self, annotator_name, blastdb):
        'It prepares and returns the annotator'
        if annotator_name == 'polyA':
            annotator = PolyaAnnotator(**self._polya_params)
        elif annotator_name == 'estscan_orf':
            annotator = EstscanOrfAnnotator(**self._estscan_params)
        elif annotator_name == 'blast':
            blast_param = None
            for blast_param_ in self._blast_params:
                if blastdb == blast_param_['blastdb']:
                    blast_param = blast_param_
                    break
            annotator = BlastAnnotator(**blast_param)
        else:
            raise NotImplementedError('This annotator type not supported')
        return annotator

    def __call__(self, seqs):
        'It orientates seqs, that should have a SeqRecord in it'
        orientations = None
        orientation_log = [None] * len(seqs)

        for annotator in self._annotators:
            if orientations:
                to_annalyze = [not o for o in orientations]
                seqs_to_analyze = list(compress(seqs, to_annalyze))
            else:
                orientations = [None] * len(seqs)
                seqs_to_analyze = seqs

            annotator_name = annotator['name']
            blastdb = annotator.get('blastdb', None)
            annotator = self._get_annotator(annotator_name, blastdb)

            annot_seqrecords = annotator(seqs_to_analyze)
            annot_strands = self._guess_orientations(annot_seqrecords,
                                                     annotator_name,
                                                     blastdb=blastdb)

            if blastdb:
                annotator_name += ' ' + os.path.basename(blastdb)

            analyzed_seqs_index = 0
            for index, orientation in enumerate(orientations):
                if orientation is None:
                    orientations[index] = annot_strands[analyzed_seqs_index]
                    if annot_strands[analyzed_seqs_index] == -1:  # reverse
                        orientation_log[index] = annotator_name
                    analyzed_seqs_index += 1
        # Now we reverse the seqs that we have guess that are reversed
        reorientated_seqrecords = []
        for orientation, seq, reason in zip(orientations, seqs,
                                            orientation_log):
            if orientation == -1:
                rev_seqrecord = seq.object.reverse_complement(id=True,
                                                              description=True,
                                                              annotations=True,
                                                              features=True,
                                                              dbxrefs=True,
                                                              name=True)
                seq = SeqWrapper(SEQRECORD, rev_seqrecord, None)
                # we mark the reason why it has been reversed
                text = '(reversed because of: {})'.format(reason)
                append_to_description(seq, text)

            reorientated_seqrecords.append(seq)
        return reorientated_seqrecords
