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

from crumbs.utils.segments_utils import (get_longest_segment, get_all_segments,
                                         get_complementary_segments,
                                         get_longest_complementary_segment)


class SegmentsTest(unittest.TestCase):
    'It tests the segments functions'

    @staticmethod
    def test_get_longest_section():
        'It gets the longest section from a list of sections'

        segments = [(0, 3), (10, 34)]
        assert (10, 34) == get_longest_segment(segments)

        segments = [(0, 3), (10, 13)]
        segment = get_longest_segment(segments)
        assert segment == (0, 3) or segment == (10, 13)

    @staticmethod
    def test_get_all_segments():
        'Give a list of discontinious segments we get all segments'
        segments = get_all_segments([(0, 10), (15, 20)], 30)
        assert segments == [((0, 10), True), ((11, 14), False),
                            ((15, 20), True), ((21, 29), False)]
        segments = get_all_segments([(15, 20)], 30)
        assert segments == [((0, 14), False),
                            ((15, 20), True), ((21, 29), False)]
        segments = get_all_segments([(15, 29)], 30)
        assert segments == [((0, 14), False), ((15, 29), True)]

    @staticmethod
    def test_non_matched():
        'Given a list of segments we get the complementary matches'
        segments = get_complementary_segments([(0, 10), (15, 20)], 30)
        assert segments == [(11, 14), (21, 29)]

    @staticmethod
    def test_get_longest_complementary():
        'It test that we can get the longest complementary segment'
        segments = [(0, 200)]
        result = get_longest_complementary_segment(segments, seq_len=200)
        assert result is None


if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
