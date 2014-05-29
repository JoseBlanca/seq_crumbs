from __future__ import division

from vcf.model import make_calldata_tuple, _Call

from vcf_crumbs.annotation import choose_samples
from vcf_crumbs.vcf_stats import get_call_data, GQ


def remove_low_quality_gt(record, gq_threshold, vcf_variant):
    ''''Filters by genotype quality and sets to uncalled those that do not meet
    the requierements'''
    for index_, call in enumerate(record):
        call_data = get_call_data(call, vcf_variant)
        if call_data[GQ] is None or call_data[GQ] < gq_threshold:
            calldata = make_calldata_tuple(record.FORMAT.split(':'))
            sampdat = []
            for format_def in record.FORMAT.split(':'):
                if format_def == 'GT':
                    value = None
                else:
                    value = getattr(call.data, format_def)
                sampdat.append(value)
            call = _Call(record, call.sample, calldata(*sampdat))
            call_data = get_call_data(call, vcf_variant)
            record.samples[index_] = call
    return record


class GenotypesInSamplesFilter(object):
    'Filter by genotypes in specific samples'

    def __init__(self, genotypes, samples=None, n_min_called=None):
        self.genotypes = genotypes
        self.genotypes.append(None)
        self.samples = samples
        self.n_min_called = n_min_called

    def __call__(self, snv):
        chosen_samples = choose_samples(snv, self.samples)
        n_called = 0
        for call in chosen_samples:
            if call.gt_type not in self.genotypes:
                return False
            if call.gt_type is not None:
                n_called += 1
        if self.n_min_called is None or n_called >= self.n_min_called:
            return True
        else:
            return False


class AlleleNumberFilter(object):
    'Filter by number of different alleles'

    def __init__(self, n_alleles, samples=None):
        self.n_alleles = n_alleles
        self.samples = samples

    def __call__(self, snv):
        chosen_samples = choose_samples(snv, self.samples)
        alleles = set()
        for call in chosen_samples:
            if call.gt_bases is None:
                continue
            if call.phased:
                bases = call.gt_bases.split('|')
            else:
                bases = call.gt_bases.split('/')
            for allele in bases:
                alleles.add(allele)
        return len(snv.alleles) == self.n_alleles


def get_n_missing_genotypes(record):
    n_missing = 0
    for call in record:
        if not call.called:
            n_missing += 1
    return n_missing


class MissingGenotypesFilter(object):
    'Filter by maximim number of missing genotypes'

    def __init__(self, max_missing_genotypes, samples=None):
        self.max_missing_genotypes = max_missing_genotypes
        self.samples = samples

    def __call__(self, snv):
        chosen_samples = choose_samples(snv, self.samples)
        n_missing = get_n_missing_genotypes(chosen_samples)
        return n_missing <= self.max_missing_genotypes
