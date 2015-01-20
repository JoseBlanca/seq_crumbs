from __future__ import division

from collections import namedtuple, Counter, OrderedDict
from itertools import imap
from array import array

from matplotlib.figure import Figure

from vcf import Reader as pyvcfReader

from crumbs.iterutils import RandomAccessIterator
from crumbs.vcf.ld import _count_biallelic_haplotypes
from crumbs.vcf.filters import _print_figure

DEF_AB_CODING_WIN = 31
MORE_THAN_2_CODINGS = 'More than 2 codings'
NOT_ENOUGH_SUPPORT = 'Not enough support'
ENOUGH_SUPPORT = 'Enough support'
NO_INFO = 'No information'
DEF_AB_CODER_THRESHOLD = 0.9

AlleleCoding = namedtuple('AlleleCoding', ['A', 'B'])


class GetCoding(object):
    def __init__(self, parents_a, parents_b):
        self.parents_a = parents_a
        self.parents_b = parents_b

    @staticmethod
    def _get_unique_allele(genotypes, samples):
        gts = [genotypes.get(sample, None) for sample in samples]
        alleles = set(allele for gt in gts if gt for allele in gt if allele != '.')
        if not alleles or len(alleles) > 1:
            return None
        elif len(alleles) == 1:
            return list(alleles)[0]

    def __call__(self, snp):
        gts = {call.sample: call.gt_alleles for call in snp.samples if call.called}
        allele_a = self._get_unique_allele(gts, self.parents_a)
        allele_b = self._get_unique_allele(gts, self.parents_b)
        if allele_a is None and allele_b is None:
            coding = None
        elif (allele_a is not None and allele_b is not None):
            if allele_a != allele_b:
                coding = AlleleCoding(allele_a, allele_b)
            else:
                coding = None
        else:
            coding = None
        return coding


class ABCoder(object):
    def __init__(self, vcf_fhand, parents_a, parents_b,
                 threshold=DEF_AB_CODER_THRESHOLD,
                 offspring=None, window=DEF_AB_CODING_WIN):
        self._reader = pyvcfReader(vcf_fhand)
        self.parents_a = parents_a
        self.parents_b = parents_b
        self._offspring = offspring
        self.window = window
        self.threshold = threshold
        self.log = Counter()
        self.indexes = array('f')

    @property
    def offspring(self):
        if self._offspring is not None:
            return self._offspring
        off = set(self._reader.samples)
        offspring = off.difference(self.parents_a).difference(self.parents_b)
        offspring = list(sorted(offspring))
        self._offspring = offspring
        return offspring

    def _deduce_coding(self, snp_and_coding, snp1_calls, snp2_idxs):
        votes = Counter()
        offspring = self.offspring
        for snp2_idx in snp2_idxs:
            snp2, coding2 = snp_and_coding[snp2_idx]
            if coding2 is None:
                continue
            snp2_calls = [snp2.genotype(sample) for sample in offspring]
            haplos = _count_biallelic_haplotypes(snp1_calls, snp2_calls,
                                                 return_alleles=True)
            if haplos is None:
                continue
            else:
                haplo_cnt, alleles = haplos
            alleles_in_major_haplo = {alleles.b: alleles.a,
                                      alleles.B: alleles.A}
            if haplo_cnt is None:
                continue

            if (coding2.A not in alleles_in_major_haplo or
               coding2.B not in alleles_in_major_haplo):
                # The offspring alleles in snp2 do not match the alleles
                # in the parents
                continue
            # This is allele A in snp 1
            allele1A = alleles_in_major_haplo[coding2.A]
            # This is allele B in snp 1
            allele1B = alleles_in_major_haplo[coding2.B]
            voted_coding1 = AlleleCoding(allele1A, allele1B)

            recomb_rate = (haplo_cnt.aB + haplo_cnt.Ab) / sum(haplo_cnt)
            weight = 2 * (0.5 - recomb_rate) if recomb_rate < 0.5 else 0
            votes[voted_coding1] += weight
        if not votes or sum(votes.values()) == 0:
            deduced_coding1 = None
            self.log[NO_INFO] += 1
        elif len(votes) > 2:
            deduced_coding1 = None
            self.log[MORE_THAN_2_CODINGS] += 1
        else:
            deduced_coding1 = votes.most_common(1)[0][0]
            index = votes[deduced_coding1] / sum(votes.values())
            self.indexes.append(index)
            if index < self.threshold:
                deduced_coding1 = None
                self.log[NOT_ENOUGH_SUPPORT] += 1
            else:
                self.log[ENOUGH_SUPPORT] += 1
        if deduced_coding1 is None:
            return None
        return {deduced_coding1.A: 'A', deduced_coding1.B: 'B'}

    def _map_to_ab(self, call, coding):
        if not call.called:
            return None
        try:
            mapped = [coding[allele] for allele in call.gt_alleles]
        except KeyError:
            msg = 'An allele is not found in the deduced AB coding: '
            msg += str(allele)
            raise RuntimeError(msg)
        return mapped

    def recode_genotypes(self, samples=None):
        get_coding = GetCoding(self.parents_a, self.parents_b)

        def mapper(snp):
            return snp, get_coding(snp)

        win = self.window
        snp_and_coding = RandomAccessIterator(imap(mapper, self._reader),
                                              rnd_access_win=win)
        offspring = self.offspring
        for idx, (snp1, coding1) in enumerate(snp_and_coding):
            snp1_calls = [snp1.genotype(sample) for sample in offspring]
            start = idx - win
            if start < 0:
                start = 0

            snp2_idxs = []
            for snp2_idx in range(start, start + win + 1):
                try:
                    snp2_chrom = snp_and_coding[snp2_idx][0].CHROM
                except IndexError:
                    continue
                if snp2_chrom == snp1.CHROM:
                    snp2_idxs.append(snp2_idx)

            coding1 = self._deduce_coding(snp_and_coding, snp1_calls,
                                          snp2_idxs)
            if coding1 is None:
                # We haven't manage to deduce the AB coding for this snp
                continue
            coding1['.'] = '.'
            if samples is None:
                calls = snp1.samples
            else:
                calls = [snp1.genotype(sample) for sample in samples]
            recoded = OrderedDict((call.sample, self._map_to_ab(call, coding1)) for call in calls)
            yield snp1, recoded

    def plot_index_hist(self, fhand):
        fig = Figure()
        axes = fig.add_subplot(111)
        axes.hist(self.indexes, fill=True, log=True, bins=20, rwidth=1)
        axes.axvline(x=self.threshold)
        axes.set_xlabel('support index for parental coding')
        axes.set_ylabel('num. SNPs')
        _print_figure(axes, fig, fhand, plot_legend=False)

    def write_log(self, fhand):
        tot = sum(self.log.values())
        fhand.write('%d SNPs divided in:\n' % tot)
        enough = self.log[ENOUGH_SUPPORT]
        fhand.write('Enough support: %d (%.1f%%)\n' % (enough,
                                                       enough / tot * 100))
        not_enough = self.log[NOT_ENOUGH_SUPPORT]
        fhand.write('Not enough support: %d (%.1f%%)\n' % (not_enough,
                                                           not_enough / tot * 100))
        no_info = self.log[NO_INFO]
        fhand.write('No information: %d (%.1f%%)\n' % (no_info,
                                                       no_info / tot * 100))
