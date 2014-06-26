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
HOM = 3

DP = 'depth'
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
        self._snv_quals = IntCounter()

        self._call_data = {HOM_REF: {'x': array(int_code), 'y':
                                     array(int_code),
                                     'value': array(float_code)},
                          HET: {'x': array(int_code), 'y': array(int_code),
                                'value': array(float_code)},
                          HOM_ALT: {'x': array(int_code), 'y': array(int_code),
                                    'value': array(float_code)}
                         }
        self._missing_calls_prc = []
        self._counts_distribution_in_gt = {}
        self._depths_distribution = {'called': [], 'uncalled': []}
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
        counts_distribution_in_gt = self._counts_distribution_in_gt
        depths_distribution = self._depths_distribution
        snv_quals = self._snv_quals

        for snp in self.reader:
            snv_qual = snp.QUAL
            if snv_qual is not None:
                snv_quals[int(snv_qual)] += 1
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
                calldata = get_call_data(call, vcf_variant)
                gt = calldata[GT]
                depth = calldata[DP]
                if depth is None:
                    depth = 0
                rc = calldata[RC]
                gq = calldata[GQ]
                if depth in counts_distribution_in_gt:
                    if gq < gq_threshold:
                        gt = None
                    if gt in counts_distribution_in_gt[depth]:
                        if rc in counts_distribution_in_gt[depth][gt]:
                            counts_distribution_in_gt[depth][gt][rc] += 1
                        else:
                            counts_distribution_in_gt[depth][gt][rc] = 1
                    else:
                        counts_distribution_in_gt[depth][gt] = \
                                                            IntCounter({rc: 1})
                else:
                    counts_distribution_in_gt[depth] = \
                                                      {gt: IntCounter({rc: 1})}
                if call.called:
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
                        depths_distribution['called'].append(depth)
                    else:
                        depths_distribution['uncalled'].append(depth)
                else:
                    depths_distribution['uncalled'].append(depth)
            missing_calls_prc.append(100 - (num_called * 100 / total_calls))

            if snp.num_called:
                if num_called != 0:
                    perc_variable = (num_diff_to_ref_gts / num_called) * 100
                    variable_gt_per_snp[int(perc_variable)] += 1

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

    @property
    def counts_distribution_in_gt(self):
        return self._counts_distribution_in_gt

    @property
    def depths_distribution(self):
        return self._depths_distribution

    @property
    def snv_quals(self):
        return self._snv_quals


def _get_seq_lengths(fhand):
    return {get_name(seq): get_length(seq) for seq in read_seqs([fhand])}


def calc_n_bases_in_chrom_with_snp(counts, ref_fhand):
    n_bases = 0
    ref_lengths = _get_seq_lengths(ref_fhand)
    for ref_name, length in ref_lengths.items():
        if ref_name.strip() in counts:
            n_bases += length
    return n_bases


def calc_n_bases_per_n_snps_in_chrom(counts, ref_fhand):
    distribution = {}
    ref_lengths = _get_seq_lengths(ref_fhand)
    for ref_name, length in ref_lengths.items():
        if ref_name.strip() in counts:
            n = counts[ref_name.strip()]
            if n in distribution:
                distribution[n] += length
            else:
                distribution[n] = length
    return distribution


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


def create_snp_name(vcf_record):
    return str(vcf_record.CHROM) + '_' + str(vcf_record.POS)


def choose_samples(record, sample_names):
    if sample_names is None:
        chosen_samples = record.samples
    else:
        filter_by_name = lambda x: True if x.sample in sample_names else False
        chosen_samples = filter(filter_by_name, record.samples, )
    return chosen_samples


class VCFcomparisons(object):
    def __init__(self, vcf_path, samples=None):
        reader = Reader(filename=vcf_path)
        self.index = {}
        self.samples = samples
        for vcf_record in reader:
            snp_name = create_snp_name(vcf_record)
            self.index[snp_name] = vcf_record

    def calculate_statistics(self, reader, samples=None):
        n_common_snps = 0
        total_snps = 0
        common_genotypes = 0
        uncalled_genotypes = 0
        different_genotypes = 0
        for vcf_record in reader:
            total_snps += 1
            snp_name = create_snp_name(vcf_record)
            if snp_name not in self.index:
                continue
            n_common_snps += 1
            for call1 in choose_samples(vcf_record, samples):
                for call2 in choose_samples(self.index[snp_name],
                                            self.samples):
                    if call1.gt_type is None or call2.gt_type is None:
                        uncalled_genotypes += 1
                    elif call1.gt_type == call2.gt_type:
                        common_genotypes += 1
                    else:
                        different_genotypes += 1
        common_snps_prc = n_common_snps / float(total_snps) * 100
        statistics = {'common_snps_prc': common_snps_prc,
                      'common': common_genotypes,
                      'uncalled': uncalled_genotypes,
                      'different': different_genotypes}
        return statistics


WINDOWS_SIZE = 50
MAFS = 'mafs'
GT_DEPTHS = 'gt_depths'
GT_QUALS = 'gt_quals'
HETEROZIGOSITY = 'heterozigosity'
GT_TYPES = 'gt_types'
SAMPLE_COUNTERS = (MAFS, GT_DEPTHS, GT_QUALS, HETEROZIGOSITY, GT_TYPES)

SNV_QUALS = 'snv_quals'
SNV_DENSITY = 'snv_density'
HET_IN_SNP = 'heterozigotes_in_snp'


class _AlleleCounts2D(object):
    def __init__(self):
        self._data = {}
        self._genotypes = {HOM_REF: set(['0/0']),
                           HOM_ALT: set(),
                           HET: set()}

    @property
    def genotypes(self):
        return self._genotypes

    def add(self, rc, acs, gt, gq):
        if gt != '0/0':
            if gt[0] == gt[-1]:
                self._genotypes[HOM_ALT].add(gt)
            else:
                self._genotypes[HET].add(gt)

        data = self._data
        if rc not in data:
            data[rc] = {}
        if acs not in data[rc]:
            data[rc][acs] = {}
        if gt not in data[rc][acs]:
            data[rc][acs][gt] = {'num_gt': 0, 'sum_gq': 0}
        data[rc][acs][gt]['num_gt'] += 1
        data[rc][acs][gt]['sum_gq'] += gq

    def _get_data_for_gt_type(self, ref_count, alt_count, gt_type):
        genotypes = self._data.get(ref_count, {}).get(alt_count, None)
        if genotypes is None:
            return None

        gts = self._genotypes[gt_type]
        gt_count = 0
        gq_quals = 0
        for gt in gts:
            if gt in genotypes:
                gt_count += genotypes[gt]['num_gt']
                gq_quals += genotypes[gt]['sum_gq']

        if not gt_count:
            return None

        return gt_count, gq_quals

    def get_gt_count(self, ref_count, alt_count, gt_type):
        count_and_qual = self._get_data_for_gt_type(ref_count, alt_count,
                                                    gt_type)
        if count_and_qual is None:
            return None
        else:
            return count_and_qual[0]

    def get_avg_gt_qual(self, ref_count, alt_count, gt_type):
        counts_and_quals = self._get_data_for_gt_type(ref_count, alt_count,
                                                      gt_type)
        if counts_and_quals is None:
            return None

        gt_count, gq_quals = counts_and_quals

        return gq_quals / gt_count

    def get_avg_gt_quals(self, gt_type, rc_max=None, ac_max=None):
        for rc, allele_counts in self._data.items():
            if rc_max and rc > rc_max:
                continue
            for ac in allele_counts:
                if ac_max and ac > ac_max:
                    continue
                avg_qual = self.get_avg_gt_qual(rc, ac, gt_type)
                if avg_qual is not None:
                    yield rc, ac, avg_qual

    def get_gt_counts(self, gt_type, rc_max=None, ac_max=None):
        for rc, allele_counts in self._data.items():
            if rc_max and rc > rc_max:
                continue
            for ac in allele_counts:
                if ac_max and ac > ac_max:
                    continue
                gt_count = self.get_gt_count(rc, ac, gt_type)
                if gt_count is not None:
                    yield rc, ac, gt_count

    def get_gt_depths_for_coverage(self, coverage):
        for ref_count in range(coverage + 1):
            alt_count = coverage - ref_count

            genotypes = self._data.get(ref_count, {}).get(alt_count, None)
            if genotypes is None:
                continue
            yield ref_count, genotypes


class VcfStats2(object):
    def __init__(self, vcf_path, gq_threshold=None):
        self._reader = Reader(filename=vcf_path)
        self._random_reader = Reader(filename=vcf_path)
        self._vcf_variant = get_snpcaller_name(self._reader)
        self._samples = self._reader.samples
        self._gq_threshold = 0 if gq_threshold is None else gq_threshold
        # sample_counter
        self._sample_counters = {}
        for counter_name in SAMPLE_COUNTERS:
            if counter_name not in self._sample_counters:
                self._sample_counters[counter_name] = {}
            for sample in self._samples:
                if counter_name in (GT_DEPTHS, GT_QUALS):
                    counters = {HOM: IntCounter(), HET: IntCounter()}
                else:
                    counters = IntCounter()
                self._sample_counters[counter_name][sample] = counters

        self._snv_counters = {MAFS: IntCounter(),
                              SNV_QUALS: IntCounter(),
                              HET_IN_SNP: IntCounter(),
                              SNV_DENSITY: IntCounter()}
        self._ac2d = _AlleleCounts2D()
        self._calculate()

    def _add_maf(self, snp):
        mafs = calculate_maf(snp, vcf_variant=self._vcf_variant)
        if mafs:
            for sample, maf in mafs.items():
                if maf:
                    maf = int(maf * 100)
                    if sample == 'all':
                        self._snv_counters[MAFS][maf] += 1
                    else:
                        self._sample_counters[MAFS][sample][maf] += 1

    def _add_snv_qual(self, snp):
        snv_qual = snp.QUAL
        if snv_qual is not None:
            self._snv_counters[SNV_QUALS][int(snv_qual)] += 1

    def _add_snv_density(self, snp):
        windows_size = WINDOWS_SIZE
        pos = snp.POS
        start = pos - windows_size if pos - windows_size > windows_size else 0
        end = pos + windows_size
        chrom = snp.CHROM
        num_snvs = len(list(self._random_reader.fetch(chrom, start, end))) - 1

        self._snv_counters[SNV_DENSITY][num_snvs] += 1

    def _add_snv_het_obs_fraction(self, snp, min_num_samples=6):
        if snp.num_called < min_num_samples:
            return
        het_for_snp = int((snp.num_het / snp.num_called) * 100)
        self._snv_counters[HET_IN_SNP][het_for_snp] += 1

    def _calculate(self):
        snp_counter = 0
        vcf_variant = self._vcf_variant
        for snp in self._reader:
            snp_counter += 1
            self._add_maf(snp)
            self._add_snv_qual(snp)
            self._add_snv_density(snp)
            self._add_snv_het_obs_fraction(snp)

            for call in snp.samples:
                if not call.called:
                    continue
                sample_name = call.sample
                calldata = get_call_data(call, vcf_variant)
                dp = calldata[DP]
                gq = calldata[GQ]
                rc = calldata[RC]
                gt = calldata[GT]
                acs = sum(calldata[ACS])
                gt_type = call.gt_type
                if gq >= self._gq_threshold:
                    gt_broud_type = HET if call.is_het else HOM
                    self._sample_counters[GT_DEPTHS][sample_name][gt_broud_type][dp] += 1
                    self._sample_counters[GT_QUALS][sample_name][gt_broud_type][gq] += 1
                    self._sample_counters[GT_TYPES][sample_name][gt_type] += 1
                self._ac2d.add(rc=rc, acs=acs, gt=gt, gq=gq)

    def _get_sample_counter(self, kind, sample=None, gt_broud_type=None):
        counters = self._sample_counters[kind]
        if sample is not None:
            if gt_broud_type is None:
                return counters[sample]
            else:
                return counters[sample][gt_broud_type]
        all_counters = IntCounter()
        for sample_counter in counters.values():
            if gt_broud_type is None:
                all_counters += sample_counter
            else:
                all_counters += sample_counter[gt_broud_type]
        return all_counters

    def mafs(self, sample=None):
        if sample is None:
            return self._snv_counters[MAFS]
        return self._get_sample_counter(MAFS, sample)

    def gt_depths(self, gt_broud_type, sample=None):
        return self._get_sample_counter(GT_DEPTHS, sample,
                                        gt_broud_type=gt_broud_type)

    def gt_quals(self, gt_broud_type, sample=None):
        return self._get_sample_counter(GT_QUALS, sample,
                                        gt_broud_type=gt_broud_type)

    def heterozigotes_by_sample(self, sample):
        sample_gt_types = self._get_sample_counter(GT_TYPES, sample)
        het_gt = sample_gt_types[HET]
        all_gts = sample_gt_types.count
        print het_gt, all_gts
        return het_gt / all_gts

    def gt_types(self, sample=None):
        return self._get_sample_counter(GT_TYPES, sample)

    @property
    def samples(self):
        return self._samples

    @property
    def snv_density(self):
        return self._snv_counters[SNV_DENSITY]

    @property
    def snv_quals(self):
        return self._snv_counters[SNV_QUALS]

    def het_by_snp(self):
        return

    @property
    def allelecount2d(self):
        return self._ac2d

