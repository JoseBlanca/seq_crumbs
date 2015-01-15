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
import os.path
from tempfile import NamedTemporaryFile
from subprocess import check_output, CalledProcessError

from crumbs.sff_extract import SffExtractor
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.bin_utils import BIN_DIR

# pylint: disable=R0201
# pylint: disable=R0904


class SffExtractTest(unittest.TestCase):
    'It tests the SFF extraction functionality.'

    def test_seqs_in_sff(self):
        'It extracts an sff file'

        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')

        # No trim
        seqs = SffExtractor([open(sff_fpath, 'rb')]).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('tcagGGTCTACATGTTGGTTAACCCGTACTGAT')

        # trimmed
        seqs = SffExtractor([open(sff_fpath, 'rb')], trim=True).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('GGTCTACATGTTGGTTAACCCGTACTGATTTGA')

        # min left_clip and no trim
        seqs = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=5).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('tcaggGTCTACATGTTGGTTAACCCGTACTGAT')

        # min left_clip and trim
        seqs = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=5,
                            trim=True).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('GTCTACATGTTGGTTAACCCGTACTGAT')

        # empty file
        empty_fhand = NamedTemporaryFile()
        try:
            seqs = SffExtractor([open(empty_fhand.name, 'rb')]).seqs
            list(seqs)
            self.fail('ValueError expected.')
        except ValueError:
            pass

        # Wrong file type
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>a_seq\n' + 'ACTG' * 30 + '\n')
        fasta_fhand.flush()
        try:
            seqs = SffExtractor([open(fasta_fhand.name, 'rb')]).seqs
            list(seqs)
        except ValueError:
            pass

    def test_check_nucl_counts(self):
        'It checks that the nucleotide freqs are all below the given threshold'
        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')

        extractor = SffExtractor([open(sff_fpath, 'rb')])
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert extractor.clip_advice[sff_fpath] == (5, 'A')

        extractor = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=4,
                            trim=False)
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert extractor.clip_advice[sff_fpath] == (5, 'A')

        extractor = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=4,
                            trim=True)
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert extractor.clip_advice[sff_fpath] == (5, 'A')

        extractor = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=5,
                            trim=True)
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert not extractor.clip_advice[sff_fpath]


class SffExtractBinTest(unittest.TestCase):
    'It tests the sff_extract binary'

    def test_extract_sff(self):
        'It tests the sff_extract binary'
        sff_bin = os.path.join(BIN_DIR, 'sff_extract')
        assert 'usage' in check_output([sff_bin, '-h'])
        assert 'usage' in check_output([sff_bin])

        # version string
        stderr = NamedTemporaryFile()
        check_output([sff_bin, '--version'], stderr=stderr)
        assert 'version' in open(stderr.name).read()

        # clipping warning
        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')
        cmd = [sff_bin, sff_fpath]
        stderr = NamedTemporaryFile()
        try:
            check_output(cmd, stderr=stderr)
            self.fail()
        except CalledProcessError:
            assert 'Countermeasures' in open(stderr.name).read()

        # warning avoided
        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')
        cmd = [sff_bin, '--max_percentage', '70.0', sff_fpath]
        stderr = NamedTemporaryFile()
        assert 'Countermeasures' not in check_output(cmd, stderr=stderr)

        # min left clip
        cmd = [sff_bin, '--min_left_clip', '5', sff_fpath]
        assert 'Countermeasures' not in check_output(cmd)

        # clip
        cmd = [sff_bin, '--min_left_clip', '5', '--clip', sff_fpath]
        stdout = check_output(cmd)
        assert '@E3MFGYR02JWQ7T\nGTCTACATGTTG' in stdout

        # file does not exist
        cmd = [sff_bin, 'no_file']
        stderr = NamedTemporaryFile()
        try:
            check_output(cmd, stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            assert 'A file was not found: no_file' in open(stderr.name).read()

        # create traceinfofile
        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')
        trace_fhand = NamedTemporaryFile()
        cmd = [sff_bin, sff_fpath, '-x', trace_fhand.name]
        stderr = NamedTemporaryFile()
        try:
            # it fails because of the clipping warning
            check_output(cmd, stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            pass
        xml = open(trace_fhand.name).read()
        xml_part = "<trace_name>E3MFGYR02JWQ7T</trace_name>\n\t\t"
        xml_part += "<clip_quality_left>5</clip_quality_left>"
        assert xml_part in xml


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'SffExtractBinTest.test_extract_sff']
    unittest.main()
