import unittest
from tempfile import NamedTemporaryFile

import pysam

from bam_crumbs.coord_transforms import ReadRefCoord


SAM = '''@HD\tVN:1.3\tSO:coordinate
@SQ\tSN:ref\tLN:45
r001\t0\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTAAGATAAAGGATACTG\t*
r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t*
r003\t0\tref\t9\t30\t5H6M\t*\t0\t0\tAGCTAA\t*\tNM:i:1
r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCTTCAGC\t*
r005\t16\tref\t29\t30\t6H5M\t*\t0\t0\tTAGGC\t*\tNM:i:0
r008\t16\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGCCAT\t*
r009\t16\tref\t7\t30\t8M2I5M1D1M\t=\t37\t39\tAATCTATTttCTACCA\t*
r010\t0\tref\t9\t30\t3H6M1P1I4M\t*\t0\t0\tAGATAAGGATA\t*
'''


class BasePositionTest(unittest.TestCase):
    def test_base_position(self):
        #                    1111  111111222222222233333333334444444
        # Coor     01234567890123  456789012345678901234567890123456
        # ref      AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT**
        #                01234567890123 456
        # +r001          TAAGATAAAGGATA*CTG
        #               012345678 90123
        # +r002         aaaAGATAA*GGATA
        #                  012345
        # +r003       gcctaAGCTAA
        #                           012345              67890
        # +r004                     ATAGCT..............TCAGC
        #                                        43210
        # -r005                            ttagctTAGGC
        #                                                876543210
        # -r008                                          CAGCGCCAT
        #                543210987654321 0
        # -r009          AATCTATTTTCTACC*A
        #               012345678 90123
        # +r010         aaaAGATAA*GGATA

        fhand = NamedTemporaryFile(suffix='.sam')
        fhand.write(SAM)
        fhand.flush()
        sam = pysam.AlignmentFile(fhand.name)
        alig_reads = list(sam)
        coords = [ReadRefCoord(read, sam) for read in alig_reads]
        blocks = coords[0].blocks
        assert blocks == [(6, 13, 0, 7), (None, None, 8, 9), (14, 17, 10, 13),
                          (18, 18, None, None), (19, 21, 14, 16)]

        blocks = coords[4].blocks
        assert blocks == [(28, 32, 4, 0)]

        blocks = coords[5].blocks
        assert blocks == [(36, 44, 8, 0)]

        blocks = coords[6].blocks
        assert blocks == [(6, 13, 15, 8), (None, None, 7, 6), (14, 18, 5, 1),
                          (19, 19, None, None), (20, 20, 0, 0)]

    def test_guess_base_position(self):
        fhand = NamedTemporaryFile(suffix='.sam')
        fhand.write(SAM)
        fhand.flush()
        sam = pysam.AlignmentFile(fhand.name)
        coords = [ReadRefCoord(read, sam) for read in sam]
        # out of the reads
        assert coords[0].get_read_pos(('ref', 3)) is None
        assert coords[0].get_read_pos_counting_from_end(('ref', 3)) is None

        # inside a deletion
        assert coords[0].get_read_pos(('ref', 18)) is None
        assert coords[0].get_read_pos_counting_from_end(('ref', 18)) is None

        # before insertion
        assert coords[0].get_read_pos(('ref', 12)) == 6
        assert coords[0].get_read_pos_counting_from_end(('ref', 12)) == -11
        # after insertion
        assert coords[0].get_read_pos(('ref', 16)) == 12
        assert coords[0].get_read_pos_counting_from_end(('ref', 16)) == -5
        # reverse
        assert coords[6].get_read_pos(('ref', 9)) == 12
        assert coords[6].get_read_pos_counting_from_end(('ref', 9)) == -4
        # reverse
        assert coords[6].get_read_pos(('ref', 16)) == 3
        assert coords[6].get_read_pos_counting_from_end(('ref', 16)) == -13
        # soft_clip
        assert coords[1].get_read_pos(('ref', 6)) is None
        assert coords[1].get_read_pos_counting_from_end(('ref', 6)) is None

        assert coords[1].get_read_pos(('ref', 16)) == 12
        assert coords[1].get_read_pos_counting_from_end(('ref', 16)) == -2

        # TODO hard_clip
        assert coords[7].get_read_pos(('ref', 16)) == 9
        assert coords[7].get_read_pos_counting_from_end(('ref', 16)) == -2

        sam = pysam.AlignmentFile(fhand.name)
        coords = [ReadRefCoord(read, sam, True) for read in sam]
        assert coords[7].get_read_pos(('ref', 16)) == 12
        assert coords[7].get_read_pos_counting_from_end(('ref', 16)) == -2

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'CalmdTest']
    unittest.main()
