from __future__ import division
from operator import itemgetter

from vcf import Reader

from crumbs.seq.seq import get_name, get_length
from crumbs.seq.seqio import read_seqs
from crumbs.statistics import IntCounter, IntBoxplot
from crumbs.plot import get_fig_and_canvas, draw_int_boxplot
from crumbs.vcf.snv import (VARSCAN, GATK, FREEBAYES, HOM_REF, HET, HOM_ALT,
                            HOM, DEF_MIN_CALLS_FOR_POP_STATS, VCFReader,
                            pyvcfReader)

# TODO: This must be optional
from crumbs.bam.coord_transforms import ReadRefCoord


# Missing docstring
# pylint: disable=C0111

REMARKABLE_DEPTHS = (5, 10, 20, 35, 50, 70)
DP = 'depth'
ADS = 'allele_depths'
GQ = 'genotype_quality'
GT = 'genotype'


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
INBREED_F_IN_SNP = 'inbreeding_coef_f'
DEPTHS = 'depths'


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
        if rc is None or acs is None:
            return None
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


class VcfStats(object):
    def __init__(self, vcf_fpath, gq_threshold=None, dp_threshold=100,
                 min_calls_for_pop_stats=DEF_MIN_CALLS_FOR_POP_STATS,
                 remarkable_coverages=None, window_size=WINDOWS_SIZE):
        if remarkable_coverages is None:
            remarkable_depths = REMARKABLE_DEPTHS
        self.remarkable_depths = remarkable_depths

        self._reader = VCFReader(open(vcf_fpath),
                               min_calls_for_pop_stats=min_calls_for_pop_stats)

        self._random_reader = pyvcfReader(filename=vcf_fpath)

        self.window_size = window_size
        self._gq_threshold = 0 if gq_threshold is None else gq_threshold

        self.dp_threshold = dp_threshold
        self._gt_qual_depth_counter = {HOM: IntBoxplot(), HET: IntBoxplot()}
        self._ac2d = _AlleleCounts2D()

        self.sample_dp_coincidence = {1: IntCounter()}
        for cov in remarkable_depths:
            self.sample_dp_coincidence[cov] = IntCounter()

        self.called_snvs = 0
        self.called_gts = IntCounter()

        # sample_counter
        self._sample_counters = {}

        for counter_name in SAMPLE_COUNTERS:
            if counter_name not in self._sample_counters:
                self._sample_counters[counter_name] = {}
            for sample in self._reader.samples:
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
                              SNV_DENSITY: IntCounter(),
                              INBREED_F_IN_SNP: IntCounter(),
                              DEPTHS: IntCounter()}
        self._calculate()

    def _add_depth(self, snp):
        depth = snp.depth
        if depth is None:
            depth = 0
        self._snv_counters[DEPTHS][depth] += 1

    def _add_maf_and_mac(self, snp):
        maf = snp.maf
        if maf:
            maf = int(round(maf * 100))
            self._snv_counters[MAFS][maf] += 1
        mac = snp.mac
        if mac:
            self._snv_counters[MACS][mac] += 1

    def _add_maf_dp(self, snp):
        maf_dp = snp.maf_depth
        if maf_dp is not None:
            self._snv_counters[MAFS_DP][int(round(maf_dp * 100))] += 1
        for call in snp.calls:
            maf_dp = call.maf_depth
            if maf_dp is None:
                continue
            sample = call.sample
            maf_depth = int(round(maf_dp * 100))
            self._sample_counters[MAFS_DP][sample][maf_depth] += 1

    def _add_snv_qual(self, snp):
        snv_qual = snp.qual
        if snv_qual is not None:
            self._snv_counters[SNV_QUALS][int(round(snv_qual))] += 1

    def _add_snv_density(self, snp):
        windows_size = self.window_size
        pos = snp.pos
        start = pos - windows_size if pos - windows_size > windows_size else 0
        end = pos + windows_size
        chrom = snp.chrom
        num_snvs = len(list(self._random_reader.fetch(chrom, start, end))) - 1

        self._snv_counters[SNV_DENSITY][num_snvs] += 1

    def _add_snv_het_obs_fraction(self, snp):
        obs_het = snp.obs_het
        if obs_het is None:
            return
        self._snv_counters[HET_IN_SNP][int(round(obs_het * 100))] += 1

        inbreed_coef = snp.inbreed_coef
        if inbreed_coef is None:
            return
        inbreed_coef = int(round(inbreed_coef * 100))
        self._snv_counters[INBREED_F_IN_SNP][inbreed_coef] += 1

    @staticmethod
    def _num_samples_higher_equal_dp(depth, snp):
        n_samples = 0
        for call in snp.calls:
            if not call.called:
                continue
            if call.depth >= depth:
                n_samples += 1
        return n_samples

    def _calculate(self):
        snp_counter = 0
        for snp in self._reader.parse_snvs():
            snp_counter += 1
            self._add_maf_dp(snp)
            self._add_maf_and_mac(snp)
            self._add_snv_qual(snp)
            self._add_snv_density(snp)
            self._add_snv_het_obs_fraction(snp)
            self._add_depth(snp)

            for depth, counter in self.sample_dp_coincidence.viewitems():
                n_samples = self._num_samples_higher_equal_dp(depth, snp)
                counter[n_samples] += 1

            n_gt_called = 0
            for call in snp.calls:
                if not call.called:
                    continue
                n_gt_called += 1
                sample_name = call.sample
                ref_depth = call.ref_depth
                acs = call.alt_sum_depths
                gt_type = call.gt_type

                gt_broud_type = HET if call.is_het else HOM

                depth = call.depth
                gt_qual = call.gt_qual
                if depth is not None and depth < self.dp_threshold:
                    self._gt_qual_depth_counter[gt_broud_type].append(depth,
                                                                      gt_qual)
                # CHECK THIS. This is an special case where the only info we
                # have is the genotype
                if gt_qual is None:
                    self._sample_counters[GT_TYPES][sample_name][gt_type] += 1
                    if depth is not None:
                        self._sample_counters[GT_DEPTHS][sample_name][gt_broud_type][depth] += 1
                elif gt_qual >= self._gq_threshold:
                    self._sample_counters[GT_TYPES][sample_name][gt_type] += 1
                    self._sample_counters[GT_QUALS][sample_name][gt_broud_type][gt_qual] += 1
                    self._sample_counters[GT_DEPTHS][sample_name][gt_broud_type][depth] += 1
                self._ac2d.add(rc=ref_depth, acs=acs, gt=call.int_alleles,
                               gq=gt_qual)
            self.called_gts[n_gt_called] += 1
            self.called_snvs += 1

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

    def heterozigosity_for_sample(self, sample):
        sample_gt_types = self._get_sample_counter(GT_TYPES, sample)

        het_gt = sample_gt_types[HET]
        all_gts = sample_gt_types.count
        try:
            heterozigosity = het_gt / all_gts
        except ZeroDivisionError:
            heterozigosity = 0
        return heterozigosity

    def gt_types(self, sample=None):
        return self._get_sample_counter(GT_TYPES, sample)

    @property
    def samples(self):
        return self._reader.samples

    @property
    def min_calls_for_pop_stats(self):
        return self._reader.min_calls_for_pop_stats

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
    def inbreeding_by_snp(self):
        return self._snv_counters[INBREED_F_IN_SNP]

    @property
    def allelecount2d(self):
        return self._ac2d

    @property
    def gt_depths_by_gt_and_qual(self):
        return self._gt_qual_depth_counter

    @property
    def depths(self):
        return self._snv_counters[DEPTHS]


def calc_snv_read_pos_stats(sam, snvs, max_snps=None, max_pos=None):

    pileup_cols = sam.pileup()
    read_5_pos_cnts_rg = {}
    read_3_pos_cnts_rg = {}
    read_5_pos_box_rg = {}
    read_3_pos_box_rg = {}

    for index, snv in enumerate(snvs):
        if max_snps and index >= max_snps:
            break
        chrom = snv.chrom
        ref_pos = snv.pos
        snv_qual = snv.qual
        snv_col = None
        for col in pileup_cols:
            ref_name = sam.getrname(col.reference_id)
            if ref_name == chrom and col.reference_pos == ref_pos:
                snv_col = col
                break
        if snv_col is None:
            raise RuntimeError('No pileup found for snv {}:{}'.format(chrom,
                                                                      ref_pos))

        for pileup_read in snv_col.pileups:
            try:
                read_group = pileup_read.alignment.opt('RG')
            except KeyError:
                read_group = None
            read_ref_coord = ReadRefCoord(pileup_read.alignment, sam)
            read_pos = read_ref_coord.get_read_pos((chrom, ref_pos))
            read_pos_end = read_ref_coord.get_read_pos_counting_from_end((chrom,
                                                                          ref_pos))
            if read_group not in read_5_pos_cnts_rg:
                read_5_pos_cnts_rg[read_group] = IntCounter()
                read_3_pos_cnts_rg[read_group] = IntCounter()
                read_5_pos_box_rg[read_group] = IntBoxplot()
                read_3_pos_box_rg[read_group] = IntBoxplot()
            read_5_pos_cnts = read_5_pos_cnts_rg[read_group]
            read_3_pos_cnts = read_3_pos_cnts_rg[read_group]
            read_5_pos_box = read_5_pos_box_rg[read_group]
            read_3_pos_box = read_3_pos_box_rg[read_group]

            if (read_pos is not None and
                (not max_pos or read_pos + 1 <= max_pos)):
                    read_5_pos_cnts[read_pos + 1] += 1
                    read_5_pos_box.append(read_pos + 1, snv_qual)
            if (read_pos_end is not None and (not max_pos or
                abs(read_pos_end) <= max_pos)):
                    read_3_pos_cnts[abs(read_pos_end)] += 1
                    read_3_pos_box.append(abs(read_pos_end), snv_qual)

    return {'5_read_pos_counts': read_5_pos_cnts_rg,
            '3_read_pos_counts': read_3_pos_cnts_rg,
            '5_read_pos_boxplot': read_5_pos_box_rg,
            '3_read_pos_boxplot': read_3_pos_box_rg}


def _draw_one_read_pos_stats(stats, axes, box_key, count_key, title):
    axes1 = axes
    xpos = draw_int_boxplot(stats[box_key], axes=axes1)[1]
    axes1.set_ylabel('SNV Qualities')

    axes2 = axes1.twinx()
    ys = zip(*sorted(stats[count_key].items(), key=itemgetter(0)))[1]
    axes2.plot(xpos, ys, c='b')
    axes2.set_ylabel('Num. reads', color='b')
    axes1.set_xlabel('Read position')
    axes1.set_title(title)

    # vertical x labels
    for label in axes1.xaxis.get_ticklabels():
        label.set_rotation(90)


def draw_read_pos_stats(stats, plot_fhand):
    n_read_groups = len(stats['5_read_pos_counts'])
    figure, canvas = get_fig_and_canvas(n_read_groups, 2)
    n_plot = 0
    for read_group in stats['5_read_pos_counts'].keys():
        rg_stat = {'5_read_pos_counts': stats['5_read_pos_counts'][read_group],
                   '3_read_pos_counts': stats['3_read_pos_counts'][read_group],
                   '5_read_pos_boxplot': stats['5_read_pos_boxplot'][read_group],
                   '3_read_pos_boxplot': stats['3_read_pos_boxplot'][read_group]}
        n_plot += 1
        axes1 = figure.add_subplot(n_read_groups, 2, n_plot)
        if read_group is None:
            title = "Pos. counted from 5'"
        else:
            title = "Pos. counted from 5' (rg: {})".format(read_group)
        _draw_one_read_pos_stats(rg_stat, axes1, '5_read_pos_boxplot',
                                 '5_read_pos_counts', title)
        if read_group is None:
            title = "Pos. counted from 3'"
        else:
            title = "Pos. counted from 3' (rg: {})".format(read_group)
        n_plot += 1
        axes2 = figure.add_subplot(n_read_groups, 2, n_plot)
        _draw_one_read_pos_stats(rg_stat, axes2, '3_read_pos_boxplot',
                                 '3_read_pos_counts', title)
    canvas.print_figure(plot_fhand)
    plot_fhand.flush()

# TODO: we should be able to limit the snv qualities
# TODO: Y axes must start in 0
# TODO: implements max snv quality
