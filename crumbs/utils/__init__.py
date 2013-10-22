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

from sys import version_info

if version_info[0] < 3 or (version_info[0] == 3 and version_info[1] < 3):
    # until python 3.3 the standard file module has no support for
    # wrapping file object and required to open a new file
    # bz2file is a backport of the python 3.3 std library module
    try:
        from bz2file import BZ2File
    except ImportError:
        pass
else:
    from bz2 import BZ2File


def approx_equal(float_a, float_b, tolerance=0.01):
    'an aproach to compare two floating numbers'
    tol = tolerance
    return (abs(float_a - float_b) / max(abs(float_a), abs(float_b))) < tol
