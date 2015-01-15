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
# pylint: disable=W0402
# pylint: disable=C0111

import unittest

from crumbs.utils.sam import (bit_tags_to_int_flag, int_flag_to_bit_tags,
                              IS_PAIRED, IS_IN_PROPER_PAIR)


class FlagTests(unittest.TestCase):
    def test_create_flag(self):
        assert bit_tags_to_int_flag([IS_PAIRED]) == 1
        assert bit_tags_to_int_flag([IS_IN_PROPER_PAIR]) == 2
        assert bit_tags_to_int_flag([IS_PAIRED, IS_IN_PROPER_PAIR]) == 3

    def test_flag_to_binary(self):
        assert not int_flag_to_bit_tags(0)
        assert IS_PAIRED in int_flag_to_bit_tags(1)
        assert IS_IN_PROPER_PAIR in int_flag_to_bit_tags(2)
        assert IS_PAIRED in int_flag_to_bit_tags(1 | 2)
        assert IS_IN_PROPER_PAIR in int_flag_to_bit_tags(1 | 2)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'ComplexityFilterTest']
    unittest.main()