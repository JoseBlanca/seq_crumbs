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
from crumbs.collectionz import OrderedSet


class TestCollections(unittest.TestCase):

    def test_ordered_set(self):
        in_list = [1, 2, 3, 4, 5, 8, 10]
        not_in_list = [6, 9, 11, 13]
        ordered_set = OrderedSet(in_list)
        for item in in_list:
            assert item in ordered_set
        assert ordered_set.check_add(7)
        assert ordered_set.list == [1, 2, 3, 4, 5, 7, 8, 10]
        assert not ordered_set.check_add(2)
        assert ordered_set.list == [1, 2, 3, 4, 5, 7, 8, 10]
        assert ordered_set.check_add(0)
        assert ordered_set.list == [0, 1, 2, 3, 4, 5, 7, 8, 10]
        for item in not_in_list:
            assert item not in ordered_set


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestCollections']
    unittest.main()
