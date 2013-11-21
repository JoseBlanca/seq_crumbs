# Copyright 2013 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of seq_crumbs.
# vcf_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# vcf_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with vcf_crumbs. If not, see <http://www.gnu.org/licenses/>.

import unittest
from tempfile import NamedTemporaryFile

from vcf_crumbs.utils import compress_with_bgzip, uncompress_gzip


class CompressTest(unittest.TestCase):
    def test_bgzip_compression(self):
        orig = 'hola\ncaracola\n'
        orig_fhand = NamedTemporaryFile()
        orig_fhand.write(orig)
        orig_fhand.flush()

        compressed_fhand = NamedTemporaryFile(suffix='.gz')
        compress_with_bgzip(orig_fhand, compressed_fhand)

        compressed_fhand.seek(0)
        compressed = compressed_fhand.read()
        orig_fhand.seek(0)
        assert orig_fhand.read() == orig

        uncompressed_fhand = NamedTemporaryFile()
        uncompress_gzip(compressed_fhand, uncompressed_fhand)
        compressed_fhand.seek(0)
        assert compressed_fhand.read() == compressed

        uncompressed_fhand.seek(0)
        assert uncompressed_fhand.read() == orig

if __name__ == "__main__":
#     import sys;sys.argv = ['', 'FilterTest.test_close_to_filter']
    unittest.main()
