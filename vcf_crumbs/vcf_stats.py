from __future__ import division


from collections import Counter
from array import array


from  vcf import Reader

from crumbs.seq import get_name, get_length
from crumbs.seqio import read_seqs
from crumbs.statistics import IntCounter


VARSCAN = 'VarScan'
GATK = 'gatk'
FREEBAYES = 'freebayes'

# taken from pyvcf
HOM_REF = 0
HET = 1
HOM_ALT = 2

DP = 'deep'
ACS = 'alt_counts'
RC = 'ref_count'
GQ = 'genotype_quality'
GT = 'genotype'


def get_snpcaller_name(reader):
    metadata = reader.metadata
    if 'source' in metadata:
        if 'VarScan2' in metadata['source']:
            return VARSCAN
        elif 'freebayes' in metadata['source'][0].lower():
            return FREEBAYES
    if 'UnifiedGenotyper' in metadata:
        return GATK
    raise NotImplementedError('Can not get snp caller of the vcf file')


def get_call_data(call, vcf_variant):
    data = call.data
    gt = data.GT
    try:
        gq = data.GQ
    except AttributeError:
        gq = None
    dp = data.DP
    if vcf_variant == GATK:
        alleles = [int(a) for a in gt.split('/')]
        rd = 0
        ad = []
        for allele_count, allele in zip(data.AD, alleles):
            if allele == 0:
                rd = allele_count
            else:
                ad.append(allele_count)
    elif vcf_variant == VARSCAN:
        rd = data.RD
        ad = [data.AD]
    elif vcf_variant == FREEBAYES:
        rd = data.RO
        ad = data.AO
        if isinstance(ad, int):
            ad = [ad]
    else:
        raise NotImplementedError('Not using one of the supported snp callers')
    calldata = {GT: gt, GQ: gq, DP: dp, RC: rd, ACS: ad}
    return calldata


def calculate_maf(snp, vcf_variant):
    total_ad = 0
    total_rd = 0
    mafs = {}
    for call in snp.samples:
        if call.called:
            calldata = get_call_data(call, vcf_variant)
            rd = calldata[RC]
            ad = sum(calldata[ACS])
            if rd + ad == 0:
                #freebayes writes some call data although it has no read counts
                # for this sample. We have to pass those
                continue
            mafs[call.sample] = max([rd, ad]) / sum([rd, ad])
            total_ad += ad
            total_rd += rd
    values = [total_ad, total_rd]
    total = sum(values)
    if not total:
        return None
    maf = max(values) / total
    mafs['all'] = maf
    return mafs


class VcfStats(object):
    def __init__(self, vcf_path, gq_threshold=None, selected_samples=None):
        self.reader = Reader(filename=vcf_path)
        self.vcf_variant = get_snpcaller_name(self.reader)
        self._gq_threslhold = 0 if gq_threshold is None else gq_threshold
        self._samples = set()
        self._selected_samples = selected_samples
        self._snps_per_chromo = Counter()
        # mafs
        self._mafs = {}
        self._variable_gt_per_snp = IntCounter()
        #het by sample
        self._het_by_sample = {}

        int_code = 'I'
        float_code = 'f'
        self._genotype_qualities = IntCounter()

        self._call_data = {HOM_REF: {'x': array(int_code), 'y':
                                     array(int_code),
                                     'value': array(float_code)},
                          HET: {'x': array(int_code), 'y': array(int_code),
                                'value': array(float_code)},
                          HOM_ALT: {'x': array(int_code), 'y': array(int_code),
                                    'value': array(float_code)}
                         }
        self._missing_calls_prc = []
        self._calculate()

    def _calculate(self):
        vcf_variant = self.vcf_variant
        gq_threshold = self._gq_threslhold
        chrom_counts = self._snps_per_chromo
        mafs = self._mafs
        call_datas = self._call_data
        samples = self._samples
        het_by_sample = self._het_by_sample
        gqs = self._genotype_qualities
        variable_gt_per_snp = self._variable_gt_per_snp
        missing_calls_prc = self._missing_calls_prc
        selected_samples = self._selected_samples

        for snp in self.reader:
            chrom_counts[snp.CHROM] += 1

            #maf
            maf_per_sample = calculate_maf(snp, vcf_variant=vcf_variant)
            if maf_per_sample:
                for sample, maf in maf_per_sample.items():
                    if sample not in mafs:
                        mafs[sample] = IntCounter()
                    if maf:
                        maf = int(maf * 100)
                        mafs[sample][maf] += 1

            #sample_data
            num_diff_to_ref_gts = 0
            num_called = 0
            total_calls = 0
            for call in snp.samples:
                sample_name = call.sample
                if (selected_samples is not None and
                    sample_name not in selected_samples):
                    continue
                total_calls += 1
                samples.add(sample_name)
                if sample_name not in het_by_sample:
                    het_by_sample[sample_name] = IntCounter({'num_gt': 0,
                                                             'num_het': 0})
                if call.called:
                    calldata = get_call_data(call, vcf_variant)
                    gq = calldata[GQ]
                    call_datas[call.gt_type]['x'].append(calldata[RC])
                    call_datas[call.gt_type]['y'].append(sum(calldata[ACS]))
                    call_datas[call.gt_type]['value'].append(gq)
                    het_by_sample[sample_name]['num_gt'] += 1
                    gqs[int(gq)] += 1
                    if (call.is_het and gq >= gq_threshold):
                        het_by_sample[sample_name]['num_het'] += 1

                    if (gq >= gq_threshold and call.gt_type == 2):
                        num_diff_to_ref_gts += 1
                    if gq >= gq_threshold:
                        num_called += 1
            missing_calls_prc.append(num_called * 100 / total_calls)

            if snp.num_called:
                if num_called != 0:
                    perc_variable = int((num_diff_to_ref_gts / num_called) * 100)
                    variable_gt_per_snp[perc_variable] += 1

    @property
    def snps_per_chromosome(self):
        return self._snps_per_chromo

    @property
    def mafs(self):
        return self._mafs

    @property
    def variable_gt_per_snp(self):
        return self._variable_gt_per_snp

    @property
    def call_data(self):
        return self._call_data

    @property
    def genotype_qualities(self):
        return self._genotype_qualities

    @property
    def het_by_sample(self):
        return self._het_by_sample

    @property
    def samples(self):
        return list(self._samples)

    @property
    def missing_calls_prc(self):
        return self._missing_calls_prc


def get_data_from_vcf(vcf_path, gq_threshold=0):
    reader = Reader(filename=vcf_path)
    vcf_variant = get_snpcaller_name(reader)
    typecode = 'I'
    value_typecode = 'f'
    data = {'samples': set(),
            'snps_per_chromo': Counter(),
            'maf_per_snp': [],
            'variable_gt_per_snp': array('f'),
            'het_by_sample': {},
            'genotype_qualities': array('f'),
            'missing_calls_prc': [],
            'call_data': {HOM_REF: {'x': array(typecode), 'y': array(typecode),
                                    'value': array(value_typecode)},
                          HET: {'x': array(typecode), 'y': array(typecode),
                                'value': array(value_typecode)},
                          HOM_ALT: {'x': array(typecode), 'y': array(typecode),
                                    'value': array(value_typecode)}
                         }
           }

    chrom_counts = data['snps_per_chromo']
    mafs = data['maf_per_snp']
    call_datas = data['call_data']
    samples = data['samples']
    het_by_sample = data['het_by_sample']
    gqs = data['genotype_qualities']
    missing_calls_prc = data['missing_calls_prc']
    for snp in reader:
        chrom_counts[snp.CHROM] += 1
        maf = calculate_maf(snp, vcf_variant=vcf_variant)
        if maf is not None:
            mafs.append(maf)
        #sample_data
        num_diff_to_ref_gts = 0
        num_called = 0
        total_calls = 0
        for call in snp.samples:
            total_calls += 1
            sample_name = call.sample
            if sample_name not in het_by_sample:
                het_by_sample[sample_name] = IntCounter({'num_gt': 0,
                                                         'num_het': 0})
            samples.add(sample_name)
            if call.called:
                calldata = get_call_data(call, vcf_variant)
                gq = calldata[GQ]
                call_datas[call.gt_type]['x'].append(calldata[RC])
                call_datas[call.gt_type]['y'].append(sum(calldata[ACS]))
                call_datas[call.gt_type]['value'].append(gq)
                het_by_sample[sample_name]['num_gt'] += 1
                gqs.append(gq)
                if (call.is_het and gq >= gq_threshold):
                    het_by_sample[sample_name]['num_het'] += 1
                #Should it be == 1 or 2?
                if (gq >= gq_threshold and (call.gt_type == 2 or
                                            call.gt_type == 1)):
                    num_diff_to_ref_gts += 1
                if gq >= gq_threshold:
                    num_called += 1
        missing_calls_prc.append(100 - (num_called * 100 / total_calls))

        if snp.num_called:
            if num_called != 0:
                perc_variable = (num_diff_to_ref_gts / num_called) * 100
#                 if perc_variable == 0.0:
#                     print snp, snp.samples, num_diff_to_ref_gts,  snp.num_called
                data['variable_gt_per_snp'].append(perc_variable)

    return data


def _get_seq_lengths(fhand):
    return {get_name(seq): get_length(seq) for seq in read_seqs([fhand])}


def calc_n_bases_in_chrom_with_snp(counts, ref_fhand):
    n_bases = 0
    ref_lengths = _get_seq_lengths(ref_fhand)
    for ref_name, length in ref_lengths.items():
        if ref_name.strip() in counts:
            n_bases += length
    return n_bases


def calc_density_per_chrom(counts, ref_fhand, size=100):
    densities = {}
    ref_lengths = _get_seq_lengths(ref_fhand)
    for ref_name, length in ref_lengths.items():
        seq_count = counts[ref_name]
        if seq_count == 0:
            densities[ref_name] = 0
        else:
            densities[ref_name] = round((seq_count / length) * 100, 2)
    return densities
