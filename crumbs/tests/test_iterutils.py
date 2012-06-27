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

import unittest

from crumbs.iterutils import sample, length

# pylint: disable=R0201
# pylint: disable=R0904


class IterutilsTest(unittest.TestCase):
    'It tests the iterator tools.'

    def check_sampled_items(self, items, sampled_items, num_sampled_items):
        'It checks that the sample is valid'
        assert len(sampled_items) == num_sampled_items

        # no repeated ones
        sampled_items = set(sampled_items)
        assert len(sampled_items) == num_sampled_items

        assert sampled_items.issubset(items)

    def test_sample(self):
        'We can sample an iterator'
        length = 100
        items = range(length)

        num_sampled_items = 10
        for num_sampled_items in (10, 90):
            sampled_items = list(sample(items, length, num_sampled_items))
            self.check_sampled_items(items, sampled_items, num_sampled_items)

    def test_length(self):
        'We can count an iterator'
        items = xrange(10)
        assert length(items) == 10


if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
