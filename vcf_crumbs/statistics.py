from __future__ import division

from collections import Counter

from  vcf import Reader

from crumbs.seq import get_name, get_length
from crumbs.seqio import read_seqs
from crumbs.statistics import IntCounter, IntBoxplot


VARSCAN = 'VarScan'
GATK = 'gatk'
FREEBAYES = 'freebayes'

# taken from pyvcf
HOM_REF = 0
HET = 1
HOM_ALT = 2
HOM = 3

DP = 'depth'
ADS = 'allele_depths'
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


class _AlleleDepths(object):
    def __init__(self, call, snp_caller):
        self._sum_alternatives = None
        self._al_counts = None
        self.has_alternative_counts = None
        self._get_counts(call, snp_caller)

    def _get_counts(self, call, snp_caller):
        if snp_caller == GATK:
            self._get_counts_gatk(call)
        elif snp_caller == VARSCAN:
            self._get_counts_varscan(call)
        elif snp_caller == FREEBAYES:
            self._get_counts_freebayes(call)
        else:
            msg = 'SNP caller not supported yet'
            raise NotImplementedError(msg)

    def _get_counts_gatk(self, call):
        als = [int(a) for a in call.gt_alleles]
        al_counts = {al_: al_count for al_count, al_ in zip(call.data.AD, als)}
        self._al_counts = al_counts
        sum_alt = sum(alc for al, alc in al_counts.items() if al != 0)
        self._sum_alternatives = sum_alt
        self.has_alternative_counts = True

    def _get_counts_varscan(self, call):
        data = call.data
        rd = data.RD
        self._al_counts = {0: rd}
        self._sum_alternatives = data.AD
        self.has_alternative_counts = False

    def _get_counts_freebayes(self, call):
        data = call.data
        rd = data.RO

        ad = data.AO
        if isinstance(ad, int):
            ad = [ad]
        # the number of alternative alleles should be all alleles - 1
        assert len(call.site.alleles) - 1 == len(ad)

        al_counts = {allele + 1: count for allele, count in enumerate(ad)}
        self._sum_alternatives = sum(al_counts.values())
        al_counts[0] = rd
        self._al_counts = al_counts
        self.has_alternative_counts = True

    @property
    def alt_sum_depths(self):
        return self._sum_alternatives

    @property
    def ref_depth(self):
        return self._al_counts.get(0, None)

    @property
    def allele_depths(self):
        return self._al_counts


def get_call_data(call, vcf_variant):
    data = call.data
    gt = data.GT
    try:
        gq = data.GQ
    except AttributeError:
        gq = None
    dp = data.DP

    if call.called:
        depths = _AlleleDepths(call, snp_caller=vcf_variant)
        calldata = {GT: gt, GQ: gq, DP: dp, ADS: depths}
    else:
        calldata = {}
    return calldata


def calculate_maf_dp(snp, vcf_variant):
    total_ad = 0
    total_rd = 0
    mafs = {}
    for call in snp.samples:
        if call.called:
            calldata = get_call_data(call, vcf_variant)
            depths = calldata[ADS]
            rd = depths.ref_depth
            ad = depths.alt_sum_depths
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


def calculate_maf_and_mac(snp):
    n_chroms_sampled = 0
    allele_counts = Counter()
    for call in snp.samples:
        if call.called:
            genotype = call.gt_alleles
            n_chroms_sampled += len(genotype)
            for al in genotype:
                allele_counts[al] += 1
    if not n_chroms_sampled:
        return None
    maf = max(allele_counts.values()) / n_chroms_sampled
    mac = max(allele_counts.values())
    return maf, mac


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
MACS = 'macs'
MAFS_DP = 'mafs_detph'
GT_DEPTHS = 'gt_depths'
GT_QUALS = 'gt_quals'
HETEROZIGOSITY = 'heterozigosity'
GT_TYPES = 'gt_types'
SAMPLE_COUNTERS = (MAFS_DP, GT_DEPTHS, GT_QUALS, HETEROZIGOSITY,
                   GT_TYPES)

SNV_QUALS = 'snv_quals'
SNV_DENSITY = 'snv_density'
HET_IN_SNP = 'heterozigotes_in_snp'


class _AlleleCounts2D(object):
    def __init__(self):
        self._data = {}
        self._genotypes = {HOM_REF: set([(0, 0)]),
                           HOM_ALT: set(),
                           HET: set()}

    @property
    def genotypes(self):
        return self._genotypes

    def add(self, rc, acs, gt, gq):
        gt = tuple(map(int, sorted(gt)))
        if gt != (0, 0):
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


MIN_NUM_SAMPLES = 6
REMARKABLE_DEPTHS = (10, 20, 35, 50, 70)


class VcfStats(object):
    def __init__(self, vcf_path, gq_threshold=None, dp_threshold=100,
                 min_samples_for_heterozigosity=MIN_NUM_SAMPLES,
                 remarkable_coverages=None):
        if remarkable_coverages is None:
            remarkable_depths = REMARKABLE_DEPTHS
        self.remarkable_depths = remarkable_depths

        self._reader = Reader(filename=vcf_path)
        self._random_reader = Reader(filename=vcf_path)

        self._vcf_variant = get_snpcaller_name(self._reader)
        self._samples = self._reader.samples
        self._gq_threshold = 0 if gq_threshold is None else gq_threshold
        self._min_samples_for_heterozigosity = min_samples_for_heterozigosity

        self.dp_threshold = dp_threshold
        self._gt_qual_cov_counter = {HOM: IntBoxplot(), HET: IntBoxplot()}
        self._ac2d = _AlleleCounts2D()

        self.sample_dp_coincidence = {1: IntCounter()}
        for cov in remarkable_depths:
            self.sample_dp_coincidence[cov] = IntCounter()

        self.called_gts = IntCounter()

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
                              MACS: IntCounter(),
                              MAFS_DP: IntCounter(),
                              SNV_QUALS: IntCounter(),
                              HET_IN_SNP: IntCounter(),
                              SNV_DENSITY: IntCounter()}
        self._calculate()

    def _add_maf_and_mac(self, snp):
        maf, mac = calculate_maf_and_mac(snp)
        if maf:
            maf = int(round(maf * 100))
            self._snv_counters[MAFS][maf] += 1
        if mac:
            self._snv_counters[MACS][mac] += 1

    def _add_maf_dp(self, snp):
        mafs = calculate_maf_dp(snp, vcf_variant=self._vcf_variant)
        if mafs:
            for sample, maf in mafs.items():
                if maf:
                    maf = int(round(maf * 100))
                    if sample == 'all':
                        self._snv_counters[MAFS_DP][maf] += 1
                    else:
                        self._sample_counters[MAFS_DP][sample][maf] += 1

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

    def _add_snv_het_obs_fraction(self, snp, min_num_samples):
        if snp.num_called < min_num_samples:
            return
        het_for_snp = int((snp.num_het / snp.num_called) * 100)
        self._snv_counters[HET_IN_SNP][het_for_snp] += 1

    def _num_samples_higher_equal_dp(self, depth, snp):
        vcf_variant = self._vcf_variant
        n_samples = 0
        for call in snp.samples:
            if not call.called:
                continue
            calldata = get_call_data(call, vcf_variant)
            dp = calldata[DP]
            if dp >= depth:
                n_samples += 1
        return n_samples

    def _calculate(self):
        snp_counter = 0
        vcf_variant = self._vcf_variant
        for snp in self._reader:
            snp_counter += 1
            self._add_maf_dp(snp)
            self._add_maf_and_mac(snp)
            self._add_snv_qual(snp)
            self._add_snv_density(snp)
            self._add_snv_het_obs_fraction(snp,
                                          self._min_samples_for_heterozigosity)

            for depth, counter in self.sample_dp_coincidence.viewitems():
                n_samples = self._num_samples_higher_equal_dp(depth, snp)
                counter[n_samples] += 1

            n_gt_called = 0
            for call in snp.samples:
                if not call.called:
                    continue
                n_gt_called += 1
                sample_name = call.sample
                calldata = get_call_data(call, vcf_variant)
                dp = calldata[DP]
                gq = calldata[GQ]
                depths = calldata[ADS]
                rc = depths.ref_depth
                acs = depths.alt_sum_depths
                gt_type = call.gt_type

                gt_broud_type = HET if call.is_het else HOM

                if dp < self.dp_threshold:
                    self._gt_qual_cov_counter[gt_broud_type].append(dp, gq)

                if gq >= self._gq_threshold:
                    self._sample_counters[GT_DEPTHS][sample_name][gt_broud_type][dp] += 1
                    self._sample_counters[GT_QUALS][sample_name][gt_broud_type][gq] += 1
                    self._sample_counters[GT_TYPES][sample_name][gt_type] += 1
                self._ac2d.add(rc=rc, acs=acs, gt=call.gt_alleles, gq=gq)
            self.called_gts[n_gt_called] += 1

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

    def macs(self):
        return self._snv_counters[MACS]

    def mafs(self):
        return self._snv_counters[MAFS]

    def mafs_dp(self, sample=None):
        if sample is None:
            return self._snv_counters[MAFS_DP]
        return self._get_sample_counter(MAFS_DP, sample)

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

    @property
    def het_by_snp(self):
        return self._snv_counters[HET_IN_SNP]

    @property
    def allelecount2d(self):
        return self._ac2d

    @property
    def gt_depths_by_gt_and_qual(self):
        return self._gt_qual_cov_counter
