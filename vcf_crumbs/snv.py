from __future__ import division

from collections import Counter

from vcf import Reader as pyvcf_Reader

# Missing docstring
# pylint: disable=C0111

VARSCAN = 'VarScan'
GATK = 'gatk'
FREEBAYES = 'freebayes'
HOM_REF = 0
HET = 1
HOM_ALT = 2
HOM = 3

DEF_MIN_CALLS_FOR_POP_STATS = 10
# TODO check if SNV can be coverted in a proxy using the recipes
# http://code.activestate.com/recipes/496741-object-proxying/
# http://code.activestate.com/recipes/496742-shelfproxy/


class VCFReader(object):
    def __init__(self, fhand,
                 min_calls_for_pop_stats=DEF_MIN_CALLS_FOR_POP_STATS):
        self.pyvcf_reader = pyvcf_Reader(fsock=fhand)
        self._min_calls_for_pop_stats = min_calls_for_pop_stats
        self._snpcaller = None

    def parse_snps(self):
        min_calls_for_pop_stats = self._min_calls_for_pop_stats
        for snp in self.pyvcf_reader:
            snp = SNV(snp, reader=self,
                      min_calls_for_pop_stats=min_calls_for_pop_stats)
            yield snp

    @property
    def snpcaller(self):
        if self._snpcaller is not None:
            return self._snpcaller

        metadata = self.pyvcf_reader.metadata
        if 'source' in metadata:
            if 'VarScan2' in metadata['source']:
                snpcaller = VARSCAN
            elif 'freebayes' in metadata['source'][0].lower():
                snpcaller = FREEBAYES
        elif 'UnifiedGenotyper' in metadata:
            snpcaller = GATK
        else:
            raise NotImplementedError('Can not get snp caller of the vcf file')
        self._snpcaller = snpcaller
        return snpcaller


class SNV(object):
    '''A proxy around the pyvcf _Record with some additional functionality'''
    def __init__(self, record, reader,
                 min_calls_for_pop_stats=DEF_MIN_CALLS_FOR_POP_STATS):
        # min_calls_for_pop_stats is defined here because otherwise it would
        # behave like a global variable for all SNPs read from a common
        # reader
        self.record = record
        self.reader = reader
        self.min_calls_for_pop_stats = min_calls_for_pop_stats
        self._mac_analyzed = False
        self._maf = None
        self._mac = None
        self._maf_dp_analyzed = False
        self._allele_depths = None
        self._maf_depth = None
        self._depth = None

    @property
    def obs_het(self):
        snp = self.record
        n_called = snp.num_called
        if n_called >= self.min_calls_for_pop_stats:
            return snp.num_het / n_called
        else:
            return None

    @property
    def exp_het(self):
        snp = self.record
        if snp.num_called >= self.min_calls_for_pop_stats:
            return snp.heterozygosity
        else:
            return None

    @property
    def samples(self):
        return [Call(sample, snv=self) for sample in self.record.samples]

    def _calculate_maf_and_mac(self):
        if self._mac_analyzed:
            return
        self._mac_analyzed = True
        snp = self.record

        n_chroms_sampled = 0
        allele_counts = Counter()
        for call in snp.samples:
            if call.called:
                genotype = call.gt_alleles
                n_chroms_sampled += len(genotype)
                for allele in genotype:
                    allele_counts[allele] += 1
        if not n_chroms_sampled:
            return
        max_allele_count = max(allele_counts.values())
        self._maf = max_allele_count / n_chroms_sampled
        self._mac = n_chroms_sampled - max_allele_count

    @property
    def maf(self):
        'Frequency of the most abundant allele'
        snp = self.record
        if snp.num_called >= self.min_calls_for_pop_stats:
            self._calculate_maf_and_mac()
            return self._maf
        else:
            return None

    @property
    def mac(self):
        'Sum of the allele count of all alleles but the most abundant'
        self._calculate_maf_and_mac()
        return self._mac

    def _calculate_mafs_dp(self):
        if self._maf_dp_analyzed:
            return
        self._maf_dp_analyzed = True

        allele_depths = Counter()
        for call in self.samples:
            if not call.has_alternative_counts:
                continue
            for allele, depth in call.allele_depths.items():
                allele_depths[allele] += depth
        depth = sum(allele_depths.values())
        maf_dp = max(allele_depths.values()) / depth
        self._depth = depth
        self._maf_depth = maf_dp
        self._allele_depths = allele_depths

    @property
    def maf_depth(self):
        self._calculate_mafs_dp()
        return self._maf_depth

    @property
    def allele_depths(self):
        self._calculate_mafs_dp()
        return self._allele_depths

    @property
    def depth(self):
        self._calculate_mafs_dp()
        return self._depth

    @property
    def POS(self):
        return self.record.POS

    @property
    def is_snp(self):
        return self.record.is_snp

    def __str__(self):
        return str(self.record)

    def __unicode__(self):
        return self.record.__unicode__()


class Call(object):
    def __init__(self, call, snv):
        self.call = call
        self.snv = snv
        self._alt_sum_depths = None
        self._allele_depths = {}
        self._has_alternative_counts = None
        self._depths_analyzed = False

    def _get_allele_depths(self):
        if self._depths_analyzed:
            return
        self._depths_analyzed = True
        if not self.call.called:
            return
        snp_caller = self.snv.reader.snpcaller
        if snp_caller is None:
            msg = 'To parse the read depths you have to set the snpcaller'
            raise ValueError(msg)
        elif snp_caller == GATK:
            self._get_depths_gatk()
        elif snp_caller == VARSCAN:
            self._get_depths_varscan()
        elif snp_caller == FREEBAYES:
            self._get_depths_freebayes()
        else:
            msg = 'SNP caller not supported yet'
            raise NotImplementedError(msg)

    def _get_depths_gatk(self):
        call = self.call
        als = [int(a) for a in call.gt_alleles]
        al_counts = {al_: al_count for al_count, al_ in zip(call.data.AD, als)}
        self.allele_depths = al_counts
        sum_alt = sum(alc for al, alc in al_counts.items() if al != 0)
        self._alt_sum_depths = sum_alt
        self._has_alternative_counts = True

    def _get_depths_varscan(self):
        call = self.call
        data = call.data
        ref_depth = data.RD
        self._allele_depths = {0: ref_depth}
        self._alt_sum_depths = data.AD
        self._has_alternative_counts = False

    def _get_depths_freebayes(self):
        call = self.call
        data = call.data
        ref_depth = data.RO

        alt_depths = data.AO
        if isinstance(alt_depths, int):
            alt_depths = [alt_depths]
        # the number of alternative alleles should be all alleles - 1
        assert len(call.site.alleles) - 1 == len(alt_depths)

        al_dps = {allele + 1: count for allele, count in enumerate(alt_depths)}
        self._alt_sum_depths = sum(al_dps.values())
        al_dps[0] = ref_depth
        self._allele_depths = al_dps
        self._has_alternative_counts = True

    @property
    def ref_depth(self):
        self._get_allele_depths()
        return self._allele_depths.get(0, None)

    @property
    def alt_sum_depths(self):
        self._get_allele_depths()
        return self._alt_sum_depths

    @property
    def allele_depths(self):
        self._get_allele_depths()
        return self._allele_depths

    @property
    def has_alternative_counts(self):
        self._get_allele_depths()
        return self._has_alternative_counts

    @property
    def maf_depth(self):
        al_dps = self.allele_depths
        if self.has_alternative_counts:
            return max(al_dps.values()) / sum(al_dps.values())
        else:
            return None

    def __str__(self):
        return str(self.call)

    def __unicode__(self):
        return self.call.__unicode__()
