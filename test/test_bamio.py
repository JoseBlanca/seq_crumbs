#!/usr/bin/env python

# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of bam_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with bam_crumbs. If not, see <http://www.gnu.org/licenses/>.

import unittest
from tempfile import NamedTemporaryFile
from bam_crumbs.bamio import sam_to_bam, bam_to_sam

SAM = '''@SQ	SN:SGN-U576692	LN:1714
@SQ	SN:SGN-U572743	LN:833
SGN-E221403	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATYAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0	RG:Z:g2
SGN-E221664	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0	RG:Z:g1
'''


class BamIOTest(unittest.TestCase):
    'It test bam io functions'

    def test_bamio_converters(self):
        'test converters'
        sam_fhand = NamedTemporaryFile(suffix='.sam')
        sam_fhand.write(SAM)
        sam_fhand.flush()
        bam_fhand = NamedTemporaryFile(suffix='.bam')

        sam_to_bam(sam_fhand, bam_fhand)

        sam_fhand2 = NamedTemporaryFile(suffix='.sam')
        bam_to_sam(bam_fhand, sam_fhand2)

        result = open(sam_fhand2.name).read()
        assert result == SAM


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
