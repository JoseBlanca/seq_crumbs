from __future__ import division

from collections import OrderedDict

from crumbs.iterutils import group_in_packets

from vcf_crumbs.snv import VCFReader, VCFWriter

# Missing docstring
# pylint: disable=C0111
# Too few pulic methods
# pylint: disable=R0903

PASSED = 'passed'
FILTERED_OUT = 'filtered_out'

SNPS_PER_FILTER_PACKET = 50


def group_in_filter_packets(items, items_per_packet):
    for packet in group_in_packets(items, items_per_packet):
        yield {PASSED: packet, FILTERED_OUT: []}


def _write_log(log_fhand, tot_snps, passed_snps):
    log_fhand.write('SNVs processed: ' + str(tot_snps) + '\n')
    good_snps = sum(passed_snps.values())
    msg = 'SNVs passsed: ' + str(good_snps) + '\n'
    log_fhand.write(msg)
    msg = 'SNVs filtered out: ' + str(tot_snps - good_snps) + '\n'
    log_fhand.write(msg)
    log_fhand.write('Number of SNVs that passed each filter\n')
    for filter_, count in passed_snps.items():
        msg = filter_ + ': ' + str(count) + '\n'
        log_fhand.write(msg)
    log_fhand.flush()


def filter_snvs(in_fhand, out_fhand, filters, filtered_fhand=None,
                template_fhand=None, log_fhand=None, reader_kwargs=None):
    if reader_kwargs is None:
        reader_kwargs = {}
    reader = VCFReader(in_fhand, **reader_kwargs)

    template_reader = reader if template_fhand is None else VCFReader(template_fhand)
    writer = VCFWriter(out_fhand, template_reader=template_reader)
    if filtered_fhand:
        filtered_writer = VCFWriter(filtered_fhand,
                                    template_reader=template_reader)
    else:
        filtered_writer = None

    packets = group_in_filter_packets(reader.parse_snvs(),
                                      SNPS_PER_FILTER_PACKET)
    tot_snps = 0
    passed_snps = OrderedDict()
    for packet in packets:
        tot_snps += len(packet[PASSED]) + len(packet[FILTERED_OUT])
        for filter_ in filters:
            packet = filter_(packet)
            filter_name = filter_.__class__.__name__
            if filter_name not in passed_snps:
                passed_snps[filter_name] = 0
            passed_snps[filter_name] += len(packet[PASSED])

        for snv in packet[PASSED]:
            writer.write_snv(snv)
        if filtered_writer:
            for snv in packet[FILTERED_OUT]:
                filtered_writer.write_snv(snv)

    if log_fhand:
        _write_log(log_fhand, tot_snps, passed_snps)

    writer.flush()


class _BaseFilter(object):
    def __init__(self, samples_to_consider=None, reverse=False):
        self.reverse = reverse
        self.samples_to_consider = samples_to_consider

    def _setup_checks(self, filterpacket):
        pass

    def _do_check(self, seq):
        raise NotImplementedError()

    def __call__(self, filterpacket):
        self._setup_checks(filterpacket)
        reverse = self.reverse
        items_passed = []
        filtered_out = filterpacket[FILTERED_OUT][:]
        samples_to_consider = self.samples_to_consider
        for snv in filterpacket[PASSED]:
            if samples_to_consider is not None:
                snv_to_check = snv.filter_calls_by_sample(samples=samples_to_consider,
                                                          reverse=False)
            else:
                snv_to_check = snv

            passed = self._do_check(snv_to_check)
            if reverse:
                passed = not passed
            if passed:
                items_passed.append(snv)
            else:
                filtered_out.append(snv)

        return {PASSED: items_passed, FILTERED_OUT: filtered_out}


class CallRateFilter(_BaseFilter):
    'Filter by the min. number of genotypes called'

    def __init__(self, min_calls=None, min_call_rate=None, reverse=False,
                 samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(CallRateFilter, self).__init__(**parent_kwargs)
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

    def __init__(self, reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(BiallelicFilter, self).__init__(**parent_kwargs)

    def _do_check(self, snv):
        if len(snv.alleles) == 2:
            return True
        else:
            return False


class IsSNPFilter(_BaseFilter):
    def _do_check(self, snv):
        return snv.is_snp


class SnvQualFilter(_BaseFilter):
    def __init__(self, min_qual, reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(SnvQualFilter, self).__init__(**parent_kwargs)
        self.min_qual = min_qual

    def _do_check(self, snv):
        qual = snv.qual
        if qual is None:
            return False
        else:
            return qual >= self.min_qual


class ObsHetFilter(_BaseFilter):
    def __init__(self, min_het=None, max_het=None, remove_nd=True,
                 reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(ObsHetFilter, self).__init__(**parent_kwargs)
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


class MafFilter(_BaseFilter):
    def __init__(self, min_maf=None, max_maf=None, remove_nd=True,
                 reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(MafFilter, self).__init__(**parent_kwargs)
        if min_maf is None and max_maf is None:
            msg = 'At least one value should be given for the min or max het'
            raise ValueError(msg)
        self.min_maf = min_maf
        self.max_maf = max_maf
        self.remove_nd = remove_nd

    def _do_check(self, snv):
        min_maf = self.min_maf
        max_maf = self.max_maf
        maf = snv.maf
        if maf is None and self.remove_nd:
            return False
        if min_maf is not None and maf < min_maf:
            return False
        if max_maf is not None and maf > max_maf:
            return False
        return True
