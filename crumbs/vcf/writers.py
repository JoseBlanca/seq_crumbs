from __future__ import division

from array import array

from Bio import SeqIO
from vcf import Reader

from crumbs.vcf.ab_coding import ABCoder, DEF_AB_CODER_THRESHOLD
from crumbs.vcf.snv import get_or_create_id


# Missing docstring
# pylint: disable=C0111

IUPAC_CODES = {('A', 'C'): 'M',
               ('A', 'G'): 'R',
               ('A', 'T'): 'W',
               ('C', 'G'): 'S',
               ('C', 'T'): 'Y',
               ('G', 'T'): 'K',
               ('A', 'C', 'G'): 'V',
               ('A', 'C', 'T'): 'H',
               ('A', 'G', 'T'): 'D',
               ('C', 'G', 'T'): 'B',
               ('A', 'C', 'G', 'T'): 'N'}
IUPAC_CODES = {tuple(sorted(alleles)): code
               for alleles, code in IUPAC_CODES.items()}


def _get_unicode_alleles(snv):
    return [unicode(allele) for allele in snv.alleles]


def _replace_snvs_with_iupac(str_seq, snvs, seq_offset):
    # TODO add parameter that accepts a function to be run fore
    # every SNP, to select which ones get put as IUPAC, e.g. maf filter
    seq = list(str_seq)
    for snv in snvs:
        start = snv.start - seq_offset
        if start < 0:
            continue
        end = snv.end - seq_offset
        if end > len(seq):
            continue
        alleles = tuple(sorted(_get_unicode_alleles(snv)))
        if snv.is_snp:
            iupac_snv = IUPAC_CODES[alleles]
            seq[start] = iupac_snv
        elif snv.is_indel:
            if snv.is_deletion:
                for loc in range(start + 1, end):
                    seq[loc] = '-'
            elif len(alleles) > 2:
                # complex event
                seq[start] = '*'
            else:
                seq[start] = seq[start] + '^'
        else:
            seq[start] = seq[start] + '*'
    return u''.join(seq)


def _build_snv_section(snv):
    alleles = '/'.join(_get_unicode_alleles(snv))
    return '[' + alleles + ']'


class IlluminaWriter(object):
    '''It writes the SNPs in Illumina format

    ref_fpath should be in fasta format and it has to have a name attribute.
    min_maf controls the SNPs reported in the adjacent segments as IUPAC codes.
    '''

    # TODO add extra error classes
    # TODO include the error classes inside this class to easy access
    class NotEnoughAdjacentSequenceError(Exception):
        pass

    def __init__(self,  ref_fpath, out_fhand, length=60, vcf_fpath=None,
                 min_length=None):
        ''''It inits.

        The vcf will be used to replace in the reference sequence the SNPs
        around the SNP of interest with IUPAC codes
        '''
        self._sep = u'\t'
        self._len = length
        if min_length is None:
            min_length = length
        if min_length > length:
            msg = 'Minimum length must be smaller than required length'
            raise ValueError(msg)
        self._min_len = min_length

        self._ref_seqs = SeqIO.index(ref_fpath, format='fasta')

        if vcf_fpath:
            self._snvs = Reader(filename=vcf_fpath)
        else:
            self._snvs = None
        self._out_fhand = out_fhand
        out_fhand.write(u'CHROM\tPOS\tID\tseq\n')
        self._prev_chrom = None

    def write(self, snv):
        chrom_name = snv.CHROM

        prev_chrom = self._prev_chrom
        if prev_chrom is None or prev_chrom.name != chrom_name:
            chrom = self._ref_seqs[chrom_name]
            self._prev_chrom = chrom
        else:
            chrom = prev_chrom

        length = self._len
        min_len = self._min_len

        snv_start = snv.start   # 0 based
        snv_end = snv.end       # 1 based
        desired_start = snv_start - length  # desired segment start
        end = snv_end + length      # desired segment end
        chrom_seq = chrom.seq
        first_segment = unicode(chrom_seq[desired_start:snv_start])

        if len(first_segment) < min_len:
            msg = "Not enough sequence in 3'. ID: %s, POS: %d, CHROM: %s"
            msg %= (snv.ID, snv.POS, snv.CHROM)
            raise self.NotEnoughAdjacentSequenceError(msg)

        if self._snvs:
            real_start = snv_start - len(first_segment)
            close_snvs = self._snvs.fetch(chrom.name, start=real_start,
                                          end=snv_start)
            first_segment = _replace_snvs_with_iupac(first_segment, close_snvs,
                                                     seq_offset=real_start)

        snv_segment = _build_snv_section(snv)
        second_segment = unicode(chrom_seq[snv_end:end])
        if len(second_segment) < min_len:
            msg = "Not enough sequence in 5'. ID: %s, POS: %d, CHROM: %s"
            msg %= (snv.ID, snv.POS, snv.CHROM)
            raise self.NotEnoughAdjacentSequenceError(msg)

        if self._snvs:
            real_end = snv_end + len(second_segment)
            close_snvs = self._snvs.fetch(chrom.name, start=snv_end,
                                          end=real_end)
            second_segment = _replace_snvs_with_iupac(second_segment,
                                                      close_snvs,
                                                      seq_offset=snv_end)

        out_fhand = self._out_fhand
        sep = self._sep
        out_fhand.write(unicode(snv.CHROM))
        out_fhand.write(sep)
        out_fhand.write(unicode(snv.POS))
        out_fhand.write(sep)
        snp_id = snv.ID
        if snp_id is None:
            snp_id = u'.'
        out_fhand.write(snp_id)
        out_fhand.write(sep)
        out_fhand.write(first_segment)
        out_fhand.write(snv_segment)
        out_fhand.write(second_segment)
        out_fhand.write(u'\n')

    def flush(self):
        self._out_fhand.flush()

    def close(self):
        self._out_fhand.close()


DEF_PHYS_TO_GENET_DIST = 1000000


class RQTLWriter(object):
    geno_coding = {('A', 'A'): 'AA',
                   ('A', '-'): 'A-',
                   ('-', 'A'): 'A-',
                   ('B', 'B'): 'BB',
                   ('B', '-'): 'B-',
                   ('-', 'B'): 'B-',
                   ('A', 'B'): 'AB',
                   ('B', 'A'): 'AB'}

    def __init__(self, out_fhand, phys_to_genet_dist=None):
        self._fhand = out_fhand
        self._snv_count = 0
        self.sep = ','
        self._first_snp = True
        self._phys_to_genet_dist = phys_to_genet_dist

    def _write_header_line(self, snv):
        line = ['phenot', '', '']
        line.extend(unicode(idx) for idx in range(1, len(snv.calls) + 1))
        self._fhand.write(self.sep.join(line))
        self._fhand.write(unicode('\n'))

    def _calls_to_ab(self, snv):
        if self._first_snp:
            self._write_header_line(snv)
            self._first_snp = False

        genotype = []
        ab_coding = {'.': '-'}

        for call in snv.record.samples:
            if not call.called:
                ab_geno = 'NA'
            else:
                geno = call.gt_alleles
                ab_geno = []
                for allele in geno:
                    if allele not in ab_coding:
                        # We do not know the phase of the parents so
                        # AB won't be related with the parent phases
                        ab_allele = chr(len(ab_coding) + 64)
                        ab_coding[allele] = ab_allele
                    else:
                        ab_allele = ab_coding[allele]
                    ab_geno.append(ab_allele)
                # we have to destroy the phase of the individual
                try:
                    ab_geno = self.geno_coding[tuple(ab_geno)]
                except KeyError:
                    msg = 'Until now, there is only support for bialleic SNPs'
                    raise NotImplementedError(msg)
            genotype.append(ab_geno)
        return genotype

    def write(self, snv):
        line = []
        self._snv_count += 1
        line.append(unicode(self._snv_count))
        line.append(unicode(snv.chrom))

        pos = snv.pos + 1
        if self._phys_to_genet_dist is not None:
            pos /= self._phys_to_genet_dist
        else:
            msg = 'The genetic map availability has not been implemented yet'
            raise NotImplementedError(msg)
        pos = '%.6e' % (pos,)
        line.append(pos)

        ab_genotypes = self._calls_to_ab(snv)
        line.extend(ab_genotypes)

        self._fhand.write(self.sep.join(line))
        self._fhand.write(unicode('\n'))

    def flush(self):
        self._fhand.flush()

    def close(self):
        self._fhand.close()


def _code_to_one_letter(geno):
    if geno is None or None in geno:
        return '-'
    if len(set(geno)) == 1:
        return list(set(geno))[0]
    else:
        return 'H'


def write_parent_checker(vcf_fhand, parents_a, parents_b, genos_fhand,
                         phys_map_fhand=None,
                         coder_threshold=DEF_AB_CODER_THRESHOLD):
    sep = '\t'
    coder = ABCoder(vcf_fhand, parents_a, parents_b, threshold=coder_threshold)
    samples = coder.offspring
    snp_ids = []
    snp_genos = []
    coding = 'ascii'
    for snp, genos in coder.recode_genotypes(samples):
        snp_id = get_or_create_id(snp).encode(coding)
        snp_ids.append(snp_id)
        geno_array = array('c')
        for geno in genos.values():
            geno_array.append(_code_to_one_letter(geno))
        snp_genos.append(geno_array)
        if phys_map_fhand is not None:
            phys_map_fhand.write(snp_id)
            phys_map_fhand.write(sep)
            phys_map_fhand.write(str(snp.POS))
            phys_map_fhand.write(sep)
            phys_map_fhand.write(snp.CHROM.encode(coding))
            phys_map_fhand.write('\n')

    genos_fhand.write('ID')
    genos_fhand.write(sep)
    genos_fhand.write(sep.join(snp_ids))
    genos_fhand.write('\n')

    for sample_idx, sample in enumerate(samples):
        genos_fhand.write(sample.encode(coding))
        genos_fhand.write(sep)
        to_write = sep.join(snp_genos[snp_idx][sample_idx] for snp_idx in range(len(snp_ids)))
        genos_fhand.write(to_write)
        genos_fhand.write('\n')

    return coder
