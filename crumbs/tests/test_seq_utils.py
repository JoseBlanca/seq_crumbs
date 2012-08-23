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

from crumbs.seq_utils import uppercase_length


class UppercaseLengthTest(unittest.TestCase):
    'It tests the uppercase character count'
    def test_uppercase_length(self):
        'It counts the number of uppercase letters in a string'
        assert uppercase_length('aCTaGGt') == 4
        assert uppercase_length('acagt') == 0

if __name__ == "__main__":
#    import sys;sys.argv = ['', 'TestPool']
    unittest.main()
