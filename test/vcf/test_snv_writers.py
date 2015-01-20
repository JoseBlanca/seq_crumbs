
import unittest
from os.path import join as pjoin
from io import StringIO
from io import BytesIO
from os import remove
from tempfile import NamedTemporaryFile
from subprocess import check_output

from vcf import Reader

from crumbs.vcf.snv import VCFReader
from crumbs.utils.file_utils import compress_with_bgzip, index_vcf_with_tabix
from crumbs.vcf.writers import (IlluminaWriter, _replace_snvs_with_iupac,
                                RQTLWriter, DEF_PHYS_TO_GENET_DIST,
                                write_parent_checker)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111

VCF_PATH = pjoin(TEST_DATA_DIR, 'sample.vcf.gz')
REF_PATH = pjoin(TEST_DATA_DIR, 'sample_ref.fasta')
VCF_INDEL_PATH = pjoin(TEST_DATA_DIR, 'sample_indel.vcf.gz')


class BareBoneSNP(object):
    def __init__(self, alleles, start, end, is_sv=False):
        self.alleles = alleles
        self.POS = start + 1
        self.start = start
        self.end = end
        allele_lens = [len(allele) for allele in alleles]
        self.is_sv = is_sv
        if is_sv:
            # sv
            self.is_snp = False
            self.is_deletion = False
            self.is_indel = False
        elif len(set(allele_lens)) == 1:
            # snp
            self.is_snp = True
            self.is_deletion = False
            self.is_indel = False
        elif allele_lens[0] == 1:
            # insertion
            self.is_snp = False
            self.is_deletion = False
            self.is_indel = True
        else:
            # deletion
            self.is_snp = False
            self.is_deletion = True
            self.is_indel = True
        self.is_sv = is_sv


class IlluminaWriterTest(unittest.TestCase):
    def test_illumina_writer(self):
        # 01234
        # 1234567890
        # CCTGATTT-A
        # TAACGA
        #   -  C -A
        vcf = '''##fileformat=VCFv4.1
#CHROM POS ID REF ALT QUAL FILTER INFO
ref 1 . C T 10 PASS .
ref 2 . CT CA,C 10 PASS .
ref 3 . T A 10 PASS .
ref 4 . G C 10 PASS .
ref 5 . A G 10 PASS .
ref 6 . T A,C 10 PASS .
ref 7 . TT T 10 PASS .
ref 8 . T TA 10 PASS .
ref 10 . A C 10 PASS .
'''

        vcf = vcf.replace(' ', '\t')
        vcf_fhand = NamedTemporaryFile(suffix='.vcf')
        vcf_fhand.write(vcf)
        vcf_fhand.flush()
        vcf_compressed = NamedTemporaryFile(suffix='.vcf.gz')
        compress_with_bgzip(vcf_fhand, vcf_compressed)
        index_vcf_with_tabix(vcf_compressed.name)

        ref_fhand = NamedTemporaryFile(suffix='.fasta')
        ref_fhand.write('>ref\nACTGATTTA\n')
        ref_fhand.flush()

        out_fhand1 = StringIO()
        writer = IlluminaWriter(ref_fhand.name, out_fhand1, min_length=0,
                                vcf_fpath=vcf_compressed.name)
        for snp in Reader(filename=vcf_compressed.name):
            writer.write(snp)

        # With no SNPs converted to IUPAC around
        out_fhand2 = StringIO()
        writer = IlluminaWriter(ref_fhand.name, out_fhand2, min_length=0)
        for snp in Reader(filename=vcf_compressed.name):
            writer.write(snp)

        remove(vcf_compressed.name + '.tbi')
        expected = u'CHROM\tPOS\tID\tseq\n'
        expected += u'ref\t1\t.\t[C/T]*WSRHT-^A\n'
        expected += u'ref\t2\t.\tY[CT/CA/C]SRHT-^A\n'
        expected += u'ref\t3\t.\tYC[T/A]SRHT-^A\n'
        expected += u'ref\t4\t.\tY*W[G/C]RHT-^A\n'
        expected += u'ref\t5\t.\tY*WS[A/G]HT-^A\n'
        expected += u'ref\t6\t.\tY*WSR[T/A/C]T-^A\n'
        expected += u'ref\t7\t.\tY*WSRH[TT/T]A\n'
        expected += u'ref\t8\t.\tY*WSRHT[T/TA]A\n'
        expected += u'ref\t10\t.\tY*WSRHT-^A[A/C]\n'
        assert expected == out_fhand1.getvalue()

        expected = u'CHROM\tPOS\tID\tseq\nref\t1\t.\t[C/T]'
        assert expected in out_fhand1.getvalue()

    def test_iupac_formater(self):
        seq = 'ACTGA'
        snp = BareBoneSNP(['T', 'A'], 3, 4)
        iupac_seq = _replace_snvs_with_iupac(seq, [snp], seq_offset=1)
        assert iupac_seq == 'ACWGA'

        snp = BareBoneSNP(['T', 'A', 'G'], 3, 4)
        iupac_seq = _replace_snvs_with_iupac(seq, [snp], seq_offset=1)
        assert iupac_seq == 'ACDGA'

        snp = BareBoneSNP(['T', 'A', 'G', 'C'], 3, 4)
        iupac_seq = _replace_snvs_with_iupac(seq, [snp], seq_offset=1)
        assert iupac_seq == 'ACNGA'

        snp = BareBoneSNP(['CT', 'C'], 2, 4)
        iupac_seq = _replace_snvs_with_iupac(seq, [snp], seq_offset=1)
        assert iupac_seq == 'AC-GA'

        snp = BareBoneSNP(['C', 'CC'], 2, 3)
        iupac_seq = _replace_snvs_with_iupac(seq, [snp], seq_offset=1)
        assert iupac_seq == 'AC^TGA'

    def test_errors(self):
        # 01234
        # 1234567890
        # CCTGATTT-A
        # TAACGA
        #   -  C -A
        vcf = '''##fileformat=VCFv4.1
#CHROM POS ID REF ALT QUAL FILTER INFO
ref 1 . C T 10 PASS .
ref 2 . CT CA,C 10 PASS .
ref 3 . T A 10 PASS .
ref 4 . G C 10 PASS .
ref 5 . A G 10 PASS .
ref 6 . T A,C 10 PASS .
ref 7 . TT T 10 PASS .
ref 8 . T TA 10 PASS .
ref 10 . A C 10 PASS .
'''

        vcf = vcf.replace(' ', '\t')
        vcf_fhand = NamedTemporaryFile(suffix='.vcf')
        vcf_fhand.write(vcf)
        vcf_fhand.flush()
        vcf_compressed = NamedTemporaryFile(suffix='.vcf.gz')
        compress_with_bgzip(vcf_fhand, vcf_compressed)
        index_vcf_with_tabix(vcf_compressed.name)

        ref_fhand = NamedTemporaryFile(suffix='.fasta')
        ref_fhand.write('>ref\nACTGATTTA\n')
        ref_fhand.flush()

        out_fhand = StringIO()
        writer = IlluminaWriter(ref_fhand.name, out_fhand,
                                vcf_fpath=vcf_compressed.name)
        snps = Reader(filename=vcf_compressed.name)
        snp = snps.next()
        try:
            writer.write(snp)
            self.fail('NotEnoughAdjacentSequenceError expected')
        except IlluminaWriter.NotEnoughAdjacentSequenceError:
            pass

    def test_run_binary(self):
        binary = pjoin(BIN_DIR, 'write_snps_for_illumina')
        assert 'usage' in check_output([binary, '-h'])

        reference = pjoin(TEST_DATA_DIR, 'sample_ref.fasta')
        vcf = pjoin(TEST_DATA_DIR, 'sample.vcf.gz')

        cmd = [binary, '-r', reference, '-m', '0', vcf]
        stdout = check_output(cmd)
        assert 'GAAAT[A/C]AA' in stdout


class RQTLWriterTest(unittest.TestCase):
    VCF_HEADER = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=mysnpprogram
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Read Depth">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Read Depth">
'''
    vcf = '''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
20\t11\t.\tG\tA\t29\tPASS\tNS=3\tGT\t./.\t./.\t./.\t./.\t./.\t./.
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/.\t0/.\t0/.\t0/1'''

    expected = '''phenot,,,1,2,3,4,5,6
1,20,1.100000e-05,NA,NA,NA,NA,NA,NA
2,20,1.400000e-05,AA,BB,B-,A-,A-,AB
'''

    def test_rqtl_writer(self):
        vcf = StringIO(unicode(self.VCF_HEADER + self.vcf))
        snps = list(VCFReader(vcf).parse_snvs())

        fhand = StringIO()
        writer = RQTLWriter(fhand, phys_to_genet_dist=DEF_PHYS_TO_GENET_DIST)
        for snp in snps:
            writer.write(snp)
        assert fhand.getvalue() == self.expected

    def test_bin(self):
        vcf = self.VCF_HEADER + self.vcf
        vcf_fhand = NamedTemporaryFile(suffix='.vcf')
        vcf_fhand.write(vcf)
        vcf_fhand.flush()

        binary = pjoin(BIN_DIR, 'write_snps_for_rqtl')
        assert 'usage' in check_output([binary, '-h'])

        cmd = [binary, vcf_fhand.name]
        stdout = check_output(cmd)
        assert stdout == self.expected

        binary = pjoin(BIN_DIR, 'write_snps_for_rqtl')
        assert 'usage' in check_output([binary, '-h'])

        cmd = [binary, '-d', '100000', vcf_fhand.name]
        stdout = check_output(cmd)
        expected = '''phenot,,,1,2,3,4,5,6
1,20,1.100000e-04,NA,NA,NA,NA,NA,NA
2,20,1.400000e-04,AA,BB,B-,A-,A-,AB
'''
        assert stdout == expected


class ParentCheckerWriterTest(unittest.TestCase):

    VCF_HEADER = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=mysnpprogram
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Read Depth">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Read Depth">
'''
    vcf = '''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
20\t11\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/.\t1/.\t0/0\t0/0\t1/1\t1/1
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t./.\t./.\t1/1\t0/1\t0/1\t0/1
20\t15\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t./.\t1/1\t1/1
20\t16\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t1/1\t1/1\t0/0\t0/0
20\t17\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t2/2\t2/2\t1/1\t1/1
20\t18\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t1/1
'''
    expected = '''ID\t\20_11\t\20_16\20_17
S3\tA\tA\tA
S4\tA\tA\tA
S5\tB\tB\tB
S6\tB\tB\tB
'''

    def test_parentChecker(self):
        vcf = StringIO(unicode(self.VCF_HEADER + self.vcf))

        results_fhand = BytesIO()
        phys_fhand = BytesIO()

        write_parent_checker(vcf, parents_a=['S1'], parents_b=['S2'],
                             genos_fhand=results_fhand,
                             phys_map_fhand=phys_fhand)
        print results_fhand.getvalue()
        print phys_fhand.getvalue()


if __name__ == "__main__":
    import sys;sys.argv = ['', 'ParentCheckerWriterTest.test_parentChecker']
    unittest.main()
