from __future__ import division

from crumbs.iterutils import group_in_packets

# Missing docstring
# pylint: disable=C0111

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


class ObsHetFilter(_BaseFilter):
    def __init__(self, min_het=None, max_het=None, remove_nd=True,
                 reverse=False):
        super(ObsHetFilter, self).__init__(reverse=reverse)
        if min_het is None and max_het is None:
            msg = 'At least one value should be given for the min or max het'
            raise ValueError(msg)
        self.min_het = min_het
        self.max_het = max_het
        self.remove_nd = remove_nd

    def _do_check(self, snv):
        min_het = self.min_het
        max_het = self.max_het
        het = snv.obs_het
        if het is None and self.remove_nd:
            return False
        if min_het is not None and het < min_het:
            return False
        if max_het is not None and het > max_het:
            return False
        return True
