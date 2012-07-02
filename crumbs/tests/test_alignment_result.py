'''It tests the representation of the results from programs like blast,
ssaha2, etc. that align one sequence against a database.'''

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


# pylint: disable=R0201
# pylint: disable=R0904

import unittest
import os
from StringIO import StringIO
from tempfile import NamedTemporaryFile
import math


from crumbs.alignment_result import (BlastParser, TabularBlastParser,
                                     alignment_results_scores, ExonerateParser,
                                     filter_alignments, covered_segments,
                                     TextBlastParser)


from crumbs.tests.utils import TEST_DATA_DIR


def floats_are_equal(num1, num2, ratio=0.01):
    'Given two numbers it returns True if they are similar'
    if num1 == 0.0:
        if num2 == 0.0:
            return True
        else:
            return False
    log1 = math.log(float(num1))
    log2 = math.log(float(num2))
    return abs(log1 - log2) < ratio


def _check_sequence(sequence, expected):
    'It matches a sequence against an expected result'
    if 'name' in expected:
        assert sequence['name'] == expected['name']
    if 'description' in expected:
        assert sequence['description'] == expected['description']
    if 'length' in expected:
        assert sequence['length'] == expected['length']


def _check_match_part(match_part, expected):
    'It matches a match_part against an expected result'
    if 'query_start' in expected:
        assert match_part['query_start'] == expected['query_start']
    if 'query_end' in expected:
        assert match_part['query_end'] == expected['query_end']
    if 'query_strand' in expected:
        assert match_part['query_strand'] == expected['query_strand']
    if 'subject_start' in expected:
        assert match_part['subject_start'] == expected['subject_start']
    if 'subject_end' in expected:
        assert match_part['subject_end'] == expected['subject_end']
    if 'subject_strand' in expected:
        assert match_part['subject_strand'] == expected['subject_strand']
    for key in expected['scores']:
        assert floats_are_equal(match_part['scores'][key],
                                 expected['scores'][key])


def _check_blast(blast, expected):
    'It matches a blast results against the expected result'
    if 'query' in expected:
        _check_sequence(blast['query'], expected['query'])
    if 'matches' in expected:
        for match_index, expt_match in enumerate(expected['matches']):
            real_match = blast['matches'][match_index]
            if 'subject' in expt_match:
                _check_sequence(real_match['subject'],
                                expt_match['subject'])
            if 'match_parts' in expt_match:
                for match_part_index, expt_match_part in \
                                        enumerate(expt_match['match_parts']):
                    bl_match_part = real_match['match_parts'][match_part_index]
                    _check_match_part(bl_match_part, expt_match_part)
            if 'scores' in expt_match:
                for key in expt_match['scores']:
                    assert floats_are_equal(real_match['scores'][key],
                                             expt_match['scores'][key])
            if 'start' in expt_match:
                assert real_match['start'] == expt_match['start']
            if 'end' in expt_match:
                assert real_match['end'] == expt_match['end']
            if 'subject_start' in expt_match:
                assert real_match['subject_start'] == expt_match['subject_start']
            if 'subject_end' in expt_match:
                assert real_match['subject_end'] == expt_match['subject_end']


class BlastParserTest(unittest.TestCase):
    'It test the blast parser'

    def test_blast_parser(self):
        'It test the blast parser'
        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast.xml'))
        parser = BlastParser(fhand=blast_file)

        expected_results = [
            {'query':{'name':'cCL1Contig2',
                      'description':"<unknown description>",
                      'length':1924},
             'matches':[
                 {'subject':{'name':'chr18',
                             'description':'No definition line found',
                             'length':19691255},
                  'scores':{'expect':4.60533e-35},
                  'match_parts':[{'query_start':276, 'query_end':484,
                                  'query_strand':-1,
                                  'subject_start':477142,
                                  'subject_end':477350,
                                  'subject_strand':1,
                                  'scores':{'expect':    4.60533e-35,
                                            'similarity':84.2,
                                            'identity':  84.2}
                                 }],
                 }
             ]
            },
            {'query':{'name':'cCL1Contig3',
                      'description':"<unknown description>",
                      'length':629},
            },
            {}, {}
        ]
        n_blasts = 0
        for index, blast in enumerate(parser):
            _check_blast(blast, expected_results[index])
            n_blasts += 1
        assert n_blasts == 4

        #with the subject id given in the xml blast
        expected_results = [
            {'query':{'name':'cCL1Contig2',
                      'description':'<unknown description>',
                      'length':1924}}, {}, {}, {}]
        parser = BlastParser(fhand=blast_file)
        for index, blast in enumerate(parser):
            _check_blast(blast, expected_results[index])

        # Check using def as acceion in all the blasts
        #It changes depending on the blast output format. depends on version
        blast_file = open(os.path.join(TEST_DATA_DIR, 'melon_tair.xml'))
        parser = BlastParser(fhand=blast_file)
        assert parser.next()['matches'][0]['subject']['name'] == 'tair1'

    def test_blast_tab_parser(self):
        'It test the blast tabular parser'
        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast.tab'))
        parser = TabularBlastParser(fhand=blast_file)

        expected_results = [
            {'query':{'name':'primer'},
             'matches':[
                 {'subject':{'name':'seq_with_primer2'},
                  'scores':{'expect': 6e-10},
                  'match_parts':[{'query_start':6, 'query_end':27,
                                  'subject_start':15,
                                  'subject_end':36,
                                  'scores':{'expect':    6e-10,
                                            'identity':  100.0}
                                 }],
                 },
                {'subject':{'name':'seq_with_primer'},
                  'scores':{'expect': 6e-10},
                  'match_parts':[{'query_start':6, 'query_end':27,
                                  'subject_start':15,
                                  'subject_end':36,
                                  'scores':{'expect':    6e-10,
                                            'identity':  100.0}
                                 },
                                 {'query_start':8, 'query_end':27,
                                  'subject_start':182,
                                  'subject_end':201,
                                  'scores':{'expect':    1e-8,
                                            'identity':  100.0}
                                 }
                  ],
                 }
            ],
            },
            {'query':{'name':'primer2'},
            },
            {}, {}
        ]
        n_blasts = 0
        for index, blast in enumerate(parser):
            _check_blast(blast, expected_results[index])
            n_blasts += 1
        assert n_blasts == 2

    def test_blast_no_result(self):
        'It test that the xml output can be and empty string'
        blast_file = NamedTemporaryFile()
        blasts = BlastParser(fhand=blast_file)

        filters = [{'kind': 'best_scores',
                    'score_key': 'expect',
                    'max_score': 1e-4,
                    'score_tolerance': 10
                   }]
        filt_b = filter_alignments(blasts, config=filters,)
        try:
            filt_b.next()
            self.fail()
        except StopIteration:
            pass

    def test_blast_text_parser(self):
        'It test the blast text parser'
        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast.blast'))
        parser = TextBlastParser(fhand=blast_file)
        alignments = list(parser)
        expected_results = [
            {'query':{'name':'seq_with_primer', 'length': 301},
             'matches':[
                 {'subject':{'name':'seq_with_primer', 'length':301},
                  'scores':{'expect': 1e-159},
                  'match_parts':[{'query_start':0, 'query_end':300,
                                  'subject_start':0,
                                  'subject_end':300,
                                  'scores':{'expect':   1e-159,
                                            'identity': 100.0}
                                 }],
                 },
                 {'subject':{'name':'seq_with_primer2', 'length':262},
                  'scores':{'expect': 5e-126},
                  'match_parts':[{'query_start':0, 'query_end':281,
                                  'subject_start':0,
                                  'subject_end':261,
                                  'scores':{'expect':   5e-126,
                                            'identity': 92.0,
                                            'score':    432}
                                 },
                                 {'query_start':282, 'query_end':300,
                                  'subject_start':18,
                                  'subject_end':36,
                                  'scores':{'expect':   6e-05,
                                            'identity': 94.0,
                                            'score':    30.1}
                                 }],
                 },
             ]
            },
            {'query':{'name':'seq_with_primer2', 'length': 262},
             'matches':[
                 {'subject':{'name':'seq_with_primer', 'length':301},
                  'scores':{'expect': 5e-126},
                  'match_parts':[{'query_start':0, 'query_end':261,
                                  'subject_start':0,
                                  'subject_end':281,
                                  'scores':{'expect':   5e-126,
                                            'identity': 92.0}
                                 },
                                 {'query_start':18, 'query_end':36,
                                  'subject_start':282,
                                  'subject_end':300,
                                  'scores':{'expect':   6e-5,
                                            'identity': 94.0}
                                 },
                                ],
                 },
                 {'subject':{'name':'seq_with_primer2', 'length':262},
                  'scores':{'expect': 2e-138},
                  'match_parts':[{'query_start':0, 'query_end':261,
                                  'subject_start':0,
                                  'subject_end':261,
                                  'scores':{'expect':   2e-138,
                                            'identity': 100.0}
                                 },
                                ],
                 }
             ],
            }
        ]

        n_blasts = 0
        for index, blast in enumerate(alignments):
            _check_blast(blast, expected_results[index])
            n_blasts += 1
        assert n_blasts == 2

        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast2.blast'))
        parser = TextBlastParser(fhand=blast_file)
        alignments = list(parser)

        expected_results = [
            {'query':{'name':'arabi', 'length': 456},
             'matches':[
                 {'subject':{'name':'AT1G55265.1',
                             'description': 'a gene',
                             'length':693},
                  'scores':{'expect': 0.0},
                  'match_parts':[{'query_start':0, 'query_end':455,
                                  'subject_start':237,
                                  'subject_end':692,
                                  'scores':{'expect':   0.0,
                                            'identity': 100.0}
                                 }],
                 },
              ]
             }
        ]

        n_blasts = 0
        for index, blast in enumerate(alignments):
            _check_blast(blast, expected_results[index])
            n_blasts += 1
        assert n_blasts == 1


def _summarize_matches(parser):
    '''Given a alignment result parser it returns a dict with the matches for
    each query'''
    summary = {}
    for result in parser:
        query_name = result['query']['name']
        matches = result['matches']
        summary[query_name] = matches
    return summary


def _check_match_summary(match_summary, expected):
    '''Given a match summary it checks that the correct number of hits
    remain after a match filtering'''
    for query_name in expected:
        if expected[query_name] == 0:
            assert query_name not in match_summary
        else:
            assert len(match_summary[query_name]) == expected[query_name]


EXONERATE_OUTPUT = '''
" (peio@oiz: ~/work_in/tmp/exonerate)Command line: [exonerate --model affine:local query.fasta adaptators.fasta --showalignment False --showvulgar False --ryo cigar_like:%S %ql %tl\n]
Hostname: [oiz]
cigar_like:prueb2 7 31 + adaptor1 5 29 + 120 55 34 95
cigar_like:prueba 0 29 + adaptor2 5 34 + 135 54 35 95
cigar_like:prueba 27 54 + adaptor2 5 32 + 136 54 35 95
cigar_like:prueba 0 28 + adaptor3 5 32 + 126 54 35 95
cigar_like:prueba 27 54 + adaptor3 5 33 + 131 54 35 95
-- completed exonerate analysis
 '''


class ExonerateParserTest(unittest.TestCase):
    '''It tests the exonertae parser output - Cigar output '''
    @staticmethod
    def test_exonerate_parser():
        'It tests exonerate parser'
        parser = ExonerateParser(fhand=StringIO(EXONERATE_OUTPUT))

        expected_results = [
            {'query':{'name':'prueba', 'length':54},
             'matches':[
                 {'subject':{'name':'adaptor2', 'length':35},
                  'scores':{'score':136},
                  'start':0,
                  'end':53,
                  'match_parts':[{'query_start':0, 'query_end':28,
                                  'query_strand':1,
                                  'subject_start':5,
                                  'subject_end':33,
                                  'subject_strand':1,
                                  'scores':{'score':135, 'similarity':95}},
                                  {'query_start':27, 'query_end':53,
                                  'query_strand':1,
                                  'subject_start':5,
                                  'subject_end':31,
                                  'subject_strand':1,
                                  'scores':{'score':136, 'similarity':95}},
                                  ]}],
                 },
            {}, {}]

        n_exonerates = 0
        for index, exonerate in enumerate(parser):
            _check_blast(exonerate, expected_results[index])
            n_exonerates += 1
        assert n_exonerates == 2


class AlignmentFilters(unittest.TestCase):
    'It tests the filtering and mapping of alignment structures'
    def test_best_score_mapper(self):
        'We keep the matches with the best scores'
        filter1 = {'kind': 'best_scores',
                   'score_key': 'expect',
                   'max_score': 1e-3,
                   'score_tolerance': 10
                   }

        align1 = {'matches': [{'scores':{'expect':1e-4},
                               'start':0,
                               'end':100,
                               'subject_start':0,
                               'subject_end':100,
                               'match_parts':[{'scores':{'expect':1e-4},
                                              'query_start':0, 'query_end':10,
                                              'subject_start':0,
                                              'subject_end':10,
                                             },
                                             {'scores':{'expect':5e-4},
                                              'query_start':30, 'query_end':40,
                                              'subject_start':30,
                                              'subject_end':40,
                                             },
                                             {'scores':{'expect':1e-3},
                                              'query_start':50, 'query_end':60,
                                              'subject_start':50,
                                              'subject_end':60,
                                             },
                                             {'scores':{'expect':1e-2},
                                             'query_start':80, 'query_end':100,
                                              'subject_start':80,
                                              'subject_end':100,
                                             }
                                            ],
                               },
                               {'scores':{'expect':1e-3},
                               'match_parts':[{'scores':{'expect':1e-3}}],
                               },
                               {'scores':{'expect':1e-2},
                                'match_parts':[{'scores':{'expect':1e-2}}],
                               }
                             ]
                 }
        align2 = {'matches': [{'scores':{'expect':1e-2},
                               'match_parts':[{'scores':{'expect':1e-2}}],
                              }]}
        alignments = [align1, align2]
        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter1]))
        expected_align1 = {'matches': [{'scores':{'expect':1e-4},
                               'start':0,
                               'end':40,
                               'subject_start':0,
                               'subject_end':40,
                               'match_parts':[{'scores':{'expect':1e-4},
                                              'query_start':0, 'query_end':10,
                                              'subject_start':0,
                                              'subject_end':10,
                                             },
                                             {'scores':{'expect':5e-4},
                                              'query_start':30, 'query_end':40,
                                              'subject_start':30,
                                              'subject_end':40,
                                             },
                                            ],
                               },
                             ]
                 }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1

    def test_min_score_mapper(self):
        'We keep the matches with the scores above the threshold'
        filter1 = {'kind': 'score_threshold',
                   'score_key': 'score',
                   'min_score': 100,
                   }

        align1 = {'matches': [{'scores':{'score':400},
                               'start':0,
                               'end':100,
                               'subject_start':0,
                               'subject_end':100,
                               'match_parts':[{'scores':{'score':400},
                                              'query_start':0, 'query_end':10,
                                              'subject_start':0,
                                              'subject_end':10,
                                             },
                                             {'scores':{'score':300},
                                              'query_start':30, 'query_end':40,
                                              'subject_start':30,
                                              'subject_end':40,
                                             },
                                             {'scores':{'score':50},
                                              'query_start':50, 'query_end':60,
                                              'subject_start':50,
                                              'subject_end':60,
                                             },
                                             {'scores':{'score':40},
                                             'query_start':80, 'query_end':100,
                                              'subject_start':80,
                                              'subject_end':100,
                                             }
                                            ],
                               },
                               {'scores':{'score':20},
                               'match_parts':[{'scores':{'score':20}}],
                               },
                               {'scores':{'score':90},
                                'match_parts':[{'scores':{'score':90}}],
                               }
                             ]
                 }
        align2 = {'matches': [{'scores':{'score':20},
                               'match_parts':[{'scores':{'score':20}}],
                              }]}
        alignments = [align1, align2]
        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter1]))
        expected_align1 = {'matches': [{'scores':{'score':400},
                               'start':0,
                               'end':40,
                               'subject_start':0,
                               'subject_end':40,
                               'match_parts':[{'scores':{'score':400},
                                              'query_start':0, 'query_end':10,
                                              'subject_start':0,
                                              'subject_end':10,
                                             },
                                             {'scores':{'score':300},
                                              'query_start':30, 'query_end':40,
                                              'subject_start':30,
                                              'subject_end':40,
                                              },
                                            ],
                               },
                             ]
                 }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1

    def test_max_score_mapper(self):
        filter1 = {'kind': 'best_scores',
                   'score_key': 'expect',
                   'max_score': 1e-3
                   }

        align1 = {'matches': [{'scores':{'expect':1e-4},
                               'start':0,
                               'end':100,
                               'subject_start':0,
                               'subject_end':100,
                               'match_parts':[{'scores':{'expect':1e-4},
                                              'query_start':0, 'query_end':10,
                                              'subject_start':0,
                                              'subject_end':10,
                                             },
                                             {'scores':{'expect':5e-4},
                                              'query_start':30, 'query_end':40,
                                              'subject_start':30,
                                              'subject_end':40,
                                             },
                                             {'scores':{'expect':1e-3},
                                              'query_start':50, 'query_end':60,
                                              'subject_start':50,
                                              'subject_end':60,
                                             },
                                             {'scores':{'expect':1e-2},
                                             'query_start':80, 'query_end':100,
                                              'subject_start':80,
                                              'subject_end':100,
                                             }
                                            ],
                               },
                               {'scores':{'expect':1e-3},
                               'match_parts':[{'scores':{'expect':1e-3}}],
                               },
                               {'scores':{'expect':1e-2},
                                'match_parts':[{'scores':{'expect':1e-2}}],
                               }
                             ]
                 }
        align2 = {'matches': [{'scores':{'expect':1e-2},
                               'match_parts':[{'scores':{'expect':1e-2}}],
                              }]}
        alignments = [align1, align2]
        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter1]))
        expected_align1 = {'matches': [{'scores':{'expect':1e-4},
                               'start':0,
                               'end':60,
                               'subject_start':0,
                               'subject_end':60,
                               'match_parts':[{'scores':{'expect':1e-4},
                                              'query_start':0, 'query_end':10,
                                              'subject_start':0,
                                              'subject_end':10,
                                             },
                                             {'scores':{'expect':5e-4},
                                              'query_start':30, 'query_end':40,
                                              'subject_start':30,
                                              'subject_end':40,
                                             },
                                             {'scores':{'expect':1e-3},
                                              'query_start':50, 'query_end':60,
                                              'subject_start':50,
                                              'subject_end':60,
                                             },
                                            ],
                               },
                               {'scores':{'expect':1e-3},
                               'match_parts':[{'scores':{'expect':1e-3}}],
                               },
                             ]
                 }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1

    def test_min_length_mapper(self):
        'We can filter the matches according to their length'
        filter1 = {'kind': 'min_length',
                   'min_num_residues': 100,
                   'length_in_query': True,
                   }

        align1 = {'matches': [{'match_parts':[{'query_start':0,
                                               'query_end':100,
                                               'subject_start':0,
                                               'subject_end':100},
                                              {'query_start':0,
                                               'query_end':50,
                                               'subject_start':0,
                                               'subject_end':50}, ]},
                              {'match_parts':[{'query_start':0, 'query_end':50,
                                               'subject_start':0,
                                               'subject_end':100, }]},
                             ]
                 }
        alignments = [align1]

        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter1]))
        expected_align1 = {'matches': [{'start':0, 'end':100,
                                        'subject_start':0,
                                        'subject_end':100, },
                                      ]
                          }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1
        assert len(filtered_alignments[0]['matches'][0]['match_parts']) == 2

        # Now filtering every match part
        filter1 = {'kind': 'min_length',
                   'min_num_residues': 100,
                   'length_in_query': True,
                   'filter_match_parts': True
                   }
        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter1]))
        expected_align1 = {'matches': [{'start':0, 'end':100,
                                        'subject_start':0,
                                        'subject_end':100, },
                                      ]
                          }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1
        assert len(filtered_alignments[0]['matches'][0]['match_parts']) == 1

        filter_ = {'kind': 'min_length',
                    'min_num_residues': 100,
                    'length_in_query': False,
                    }
        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter_]))
        expected_align1 = {'matches': [{'start':0, 'end':100,
                                        'subject_start':0,
                                        'subject_end':100, },
                                        {'start':0, 'end':50,
                                        'subject_start':0,
                                        'subject_end':100, },
                                      ]
                          }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1

        filter_ = {'kind': 'min_length',
                  'min_percentage': 90,
                  'length_in_query': True,
                 }
        align1 = {'query': {'length': 100},
                  'matches': [{'match_parts':[{'query_start':0, 'query_end':90,
                                               'subject_start':0,
                                               'subject_end':100, }]},
                              {'match_parts':[{'query_start':0, 'query_end':50,
                                               'subject_start':0,
                                               'subject_end':100, }]},
                             ]
                 }
        alignments = [align1]

        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter_]))
        expected_align1 = {'matches': [{'start':0, 'end':90,
                                        'subject_start':0,
                                        'subject_end':100, },
                                      ]
                          }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1

        filter_ = {'kind': 'min_length',
                   'min_percentage': 90,
                   'length_in_query': False,
                  }
        align1 = {'matches': [{'subject': {'length': 100},
                               'match_parts':[{'query_start':0,
                                               'query_end':100,
                                               'subject_start':0,
                                               'subject_end':90, }]},
                              {'subject': {'length': 100},
                               'match_parts':[{'query_start':0,
                                               'query_end':100,
                                               'subject_start':0,
                                               'subject_end':89, }]},
                             ]
                 }
        alignments = [align1]

        filtered_alignments = list(filter_alignments(alignments,
                                                     config=[filter_]))
        expected_align1 = {'matches': [{'start':0, 'end':100,
                                        'subject_start':0, 'subject_end':90, },
                                      ]
                          }
        _check_blast(filtered_alignments[0], expected_align1)
        assert len(filtered_alignments) == 1

    def test_no_filter(self):
        'It test the blast parser'
        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast.xml'))
        parser = BlastParser(fhand=blast_file)
        match_summary = _summarize_matches(parser)
        #lcl|2_0 cCL1Contig2
        #lcl|3_0 cCL1Contig3
        #lcl|4_0 cCL1Contig4
        #lcl|5_0 cCL1Contig5
        expected = {'cCL1Contig2': 3, 'cCL1Contig3': 1,
                    'cCL1Contig4': 5, 'cCL1Contig5': 8}
        _check_match_summary(match_summary, expected)

    def test_best_scores_filter(self):
        'We can keep the hits with the bests expects'
        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast.xml'))
        filters = [{'kind': 'best_scores',
                    'score_key': 'expect',
                    'max_score': 1e-4,
                    'score_tolerance': 10
                   }]
        expected = {'cCL1Contig2': 2, 'cCL1Contig3': 1,
                     'cCL1Contig4': 1, 'cCL1Contig5': 2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = filter_alignments(blasts, config=filters)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

    def test_min_scores_filter(self):
        'We can keep the hits scores above the given one'
        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast.xml'))

        #with evalue
        filters = [{'kind': 'score_threshold',
                    'score_key': 'expect',
                    'max_score': 1e-34,
                   }]
        expected = {'cCL1Contig2': 2, 'cCL1Contig3': 0,
                     'cCL1Contig4': 2, 'cCL1Contig5': 2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = filter_alignments(blasts, config=filters)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with similartiry
        filters = [{'kind': 'score_threshold',
                    'score_key': 'similarity',
                    'min_score': 92,
                   }]
        expected = {'cCL1Contig2': 0, 'cCL1Contig3': 0,
                     'cCL1Contig4': 1, 'cCL1Contig5': 2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = filter_alignments(blasts, config=filters)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

    def test_min_length_filter(self):
        'We can keep the hits length above the given one'
        blast_file = open(os.path.join(TEST_DATA_DIR, 'blast.xml'))

        #with the min length given in base pairs
        filters = [{'kind': 'min_length',
                    'min_num_residues': 500,
                    'length_in_query': True
                   }]
        expected = {'cCL1Contig2': 3, 'cCL1Contig3': 0,
                     'cCL1Contig4': 1, 'cCL1Contig5': 1}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = filter_alignments(blasts, config=filters)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with the min length given in query
        filters = [{'kind': 'min_length',
                    'min_percentage': 70,
                    'length_in_query': True
                   }]
        expected = {'cCL1Contig2': 0, 'cCL1Contig3': 0,
                     'cCL1Contig4': 2, 'cCL1Contig5': 0}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = filter_alignments(blasts, config=filters)
        match_summary = _summarize_matches(filtered_blasts)
        #print match_summary
        _check_match_summary(match_summary, expected)

        #with the min length given in subject %
        filters = [{'kind': 'min_length',
                    'min_percentage': 0.002,
                    'length_in_query': False
                   }]
        expected = {'cCL1Contig2': 3, 'cCL1Contig3': 0,
                     'cCL1Contig4': 1, 'cCL1Contig5': 2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = filter_alignments(blasts, config=filters)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)


class AlignmentSearchSimilDistribTest(unittest.TestCase):
    'It test that we can calculate the distribution of similarity'

    @staticmethod
    def test_alignment_results_scores():
        'It tests that we can get all the scores from a result'
        query = {'length': 25, 'name': 'query'}
        subject1 = {'length': 32, 'name': 'hola'}
        subject2 = {'length': 32, 'name': 'query'}
        match_part1 = {'scores': {'similarity': 90.0},
                       'query_start': 10,
                       'query_end': 20,
                       'query_strand': 1,
                       'subject_start': 0,
                       'subject_end': 10,
                       'subject_strand': 1
                      }
        match_part2 = {'scores': {'similarity': 80.0},
                       'query_start': 10,
                       'query_end': 20,
                       'query_strand': 1,
                       'subject_start': 0,
                       'subject_end': 10,
                       'subject_strand': 1
                      }
        result1 = {'query': query,
                   'matches':
                        [{'subject': subject1,
                          'start': 10, 'end': 20,
                          'scores': {'expect': 0.01},
                          'match_parts': [match_part1]
                         },
                         {'subject': subject2,
                          'start': 10, 'end': 20,
                          'scores': {'expect': 0.01},
                          'match_parts': [match_part2]
                         }]}
        results = [result1]
        #filtering the query-subj with the same names
        scores = alignment_results_scores(results, ['similarity'])
        assert scores == [90.0]
        #without the filtering
        scores = alignment_results_scores(results, ['similarity'],
                                          filter_same_query_subject=False)
        assert scores == [90.0, 80.0]

        #we can ask for several scores at once
        scores = alignment_results_scores(results,
                                          ['similarity', 'similarity'],
                                          filter_same_query_subject=False)
        assert scores == [[90.0, 80.0], [90.0, 80.0]]


#class WaterTests(unittest.TestCase):
#    'Water related tests'
#    @staticmethod
#    def test_build_water_relations():
#        '''it test the function that makes the relations between two sequences
#         using a markx10 format file'''
#       seq = 'ATGGCTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCA'
#       seq += 'AGCTAGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTTTTATGTA'
#       seq += 'CTGTTTTNACTCGCANGACCAACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAG'
#       seq += 'GGCNTGAAGGTGTGCCCACCACTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGA'
#       seq += 'TATGAGTAACGAGCAATTGGGAAAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCT'
#       seq += 'GCATTGAATTCGACATTCACAGTGGATTCGTTTACCGTGAGACCCACAGGTCACCAGGATAC'
#       seq += 'TTCGATGGACGCTACTGGACCATGTGGAAGCTGCCCATGTTTGGCTGCACCGAT'
#
#       seq2 = 'ATGGCTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCTGCTCAA'
#       seq2 += 'GCTAGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTNACTCGCA'
#       seq2 += 'NGACCAACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAGGGCNTGAAGGTGTG'
#       seq2 += 'CCCACCACTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGATATGAGTAACGAG'
#       seq2 += 'CAATTGGGAAAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCTGCATTGAATTCG'
#       seq2 += 'ACATTCACAGTGGATTCGTTTACCGTGAGACCCACAGGTCACCAGGATACTTCGATGGACG'
#       seq2 += 'CTACTGGACCATGTGGAAGCTGCCCATGTTTGGCTGCACCGAT'
#
#        subject_seq = {'length': len(seq), 'name': 'subject'}
#        query_seq = {'length': len(seq2), 'name': 'query'}
#
#        subject_fhand = temp_fasta_file(subject_seq)
#        parameters = {'subject': subject_fhand.name}
#        aligner = create_runner(tool='water', parameters=parameters)
#        result_fhand = aligner(query_seq)['water']
#        relations = build_relations_from_aligment(result_fhand,
#                                                query_name=query_seq.name,
#                                                subject_name=subject_seq.name)
#        assert relations == {'query': [(0, 50), (51, 112), (113, 409)],
#                             'subject': [(0, 50), (52, 113), (129, 425)]}


class MergeMatchesTests(unittest.TestCase):
    'Tests overlaping match_parts merging'

    def test_match_parts_merging(self):
        '''Merging match parts'''
        # ---   ---
        #  ---
        match_part1 = {'query_start': 80, 'query_end': 100,
                       'subject_start': 180, 'subject_end': 200}
        match_part2 = {'query_start': 90, 'query_end': 110,
                       'subject_start': 190, 'subject_end': 210}
        match_part3 = {'query_start': 190, 'query_end': 200,
                       'subject_start': 290, 'subject_end': 300}
        mparts = [match_part1, match_part2, match_part3]
        covered_segs = covered_segments(mparts)
        assert covered_segs == [(80, 110), (190, 200)]

        match_part2 = {'query_start': 90, 'query_end': 110,
                       'subject_start': 190, 'subject_end': 210}
        match_part3 = {'query_start': 190, 'query_end': 200,
                       'subject_start': 290, 'subject_end': 300}
        covered_segs = covered_segments([match_part2, match_part3])
        assert covered_segs == [(90, 110), (190, 200)]

        # q 0---10
        # s 0---10
        #   q 5---15
        #   s 15--25

        match_part1 = {'query_start': 0, 'query_end': 10,
                       'subject_start': 0, 'subject_end': 10}
        match_part2 = {'query_start': 5, 'query_end': 15,
                       'subject_start': 15, 'subject_end': 25}
        mparts = [match_part1, match_part2]
        covered_segs = covered_segments(mparts)
        assert covered_segs == [(0, 15)]

        covered_segs = covered_segments(mparts, in_query=False)
        assert covered_segs == [(0, 10), (15, 25)]

        match_part1 = {'query_start': 1, 'query_end': 10,
                       'subject_start': 1, 'subject_end': 10}
        match_part2 = {'query_start': 11, 'query_end': 20,
                       'subject_start': 11, 'subject_end': 20}
        match_part3 = {'query_start': 30, 'query_end': 40,
                       'subject_start': 30, 'subject_end': 40}
        mparts = [match_part1, match_part2, match_part3]
        covered_segs = covered_segments(mparts)
        assert covered_segs == [(1, 20), (30, 40)]

        match_part1 = {'query_start': 1, 'query_end': 10,
                       'subject_start': 1, 'subject_end': 10}
        match_part2 = {'query_start': 10, 'query_end': 20,
                       'subject_start': 10, 'subject_end': 20}
        mparts = [match_part1, match_part2]
        covered_segs = covered_segments(mparts)
        assert covered_segs == [(1, 20)]

        match_part1 = {'query_start': 80, 'query_end': 85,
                       'subject_start': 180, 'subject_end': 185}
        match_part2 = {'query_start': 90, 'query_end': 110,
                       'subject_start': 190, 'subject_end': 210}
        match_part3 = {'query_start': 190, 'query_end': 200,
                       'subject_start': 290, 'subject_end': 300}
        mparts = [match_part1, match_part2, match_part3]
        covered_segs = covered_segments(mparts, merge_segments_closer=10)
        assert covered_segs == [(80, 110), (190, 200)]
        covered_segs = covered_segments(mparts)
        assert covered_segs == [(80, 85), (90, 110), (190, 200)]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'BlastParserTest.test_blast_tab_parser']
    unittest.main()
