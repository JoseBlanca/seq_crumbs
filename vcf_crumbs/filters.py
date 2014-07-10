from __future__ import division

from crumbs.iterutils import group_in_packets


# TODO. This is not a filter, it's a mapper because sets calls to
# not called
def remove_low_quality_gt(record, gq_threshold, vcf_variant):
    ''''Filters by genotype quality and sets to uncalled those that do not meet
    the requierements'''
    for index_, call in enumerate(record):
        call_data = get_call_data(call, vcf_variant)
        gq = call_data.get(GQ, None)
        if gq is None or call_data[GQ] < gq_threshold:
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


PASSED = 'passed'
FILTERED_OUT = 'filtered_out'


def group_in_filter_packets(items, items_per_packet):
    for packet in group_in_packets(items, items_per_packet):
        yield {PASSED: packet, FILTERED_OUT: []}


class _BaseFilter(object):
    def __init__(self, reverse=False):
        self.reverse = reverse

    def _setup_checks(self, filterpacket):
        pass

    def _do_check(self, seq):
        raise NotImplementedError()

    def __call__(self, filterpacket):
        self._setup_checks(filterpacket)
        reverse = self.reverse
        items_passed = []
        filtered_out = filterpacket[FILTERED_OUT][:]
        for item in filterpacket[PASSED]:
            passed = self._do_check(item)
            if reverse:
                passed = not passed
            if passed:
                items_passed.append(item)
            else:
                filtered_out.append(item)

        return {PASSED: items_passed, FILTERED_OUT: filtered_out}


class CallRateFilter(_BaseFilter):
    'Filter by the min. number of genotypes called'

    def __init__(self, min_calls=None, min_call_rate=None, reverse=False):
        super(CallRateFilter, self).__init__(reverse=reverse)
        if min_calls is not None and min_call_rate is not None:
            msg = 'Both min_calls and min_call rate cannot be given'
            msg += 'at the same time'
            raise ValueError(msg)
        elif min_calls is None and min_call_rate is None:
            msg = 'min_calls or call_rate should be given'
            raise ValueError(msg)
        self.min_calls = min_calls
        self.min_call_rate = min_call_rate

    def _do_check(self, snv):
        if self.min_calls:
            if snv.num_called >= self.min_calls:
                return True
            else:
                return False
        else:
            if snv.call_rate >= self.min_call_rate:
                return True
            else:
                return False


class BiallelicFilter(_BaseFilter):
    'Filter the biallelic SNPs'

    def __init__(self, reverse=False):
        super(BiallelicFilter, self).__init__(reverse=reverse)

    def _do_check(self, snv):
        if len(snv.alleles) == 2:
            return True
        else:
            return False


class IsSNPFilter(_BaseFilter):
    def _do_check(self, snv):
        return snv.is_snp


class GenotypeQualFilter(_BaseFilter):
    def __init__(self, min_qual, reverse=False):
        super(GenotypeQualFilter, self).__init__(reverse=reverse)
        self.min_qual = min_qual

    def _do_check(self, snv):
        qual = snv.qual
        if qual is None:
            return False
        else:
            return qual >= self.min_qual
