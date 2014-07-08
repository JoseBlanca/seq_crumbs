
from vcf import Reader as pyvcf_Reader  # Imported here just for convenience
                                        # when importing from other modules

# Missing docstring
# pylint: disable=C0111

VARSCAN = 'VarScan'
GATK = 'gatk'
FREEBAYES = 'freebayes'

# TODO check if SNV can be coverted in a proxy using the recipes
# http://code.activestate.com/recipes/496741-object-proxying/
# http://code.activestate.com/recipes/496742-shelfproxy/


class SNV(object):
    '''A proxy around the pyvcf _Record with some additional functionality'''
    def __init__(self, record, min_calls_for_het=10, snp_caller=None):
        self.record = record
        self.min_calls_for_het = min_calls_for_het
        self.snp_caller = snp_caller

    @property
    def obs_het(self):
        snp = self.record
        n_called = snp.num_called
        if n_called >= self.min_calls_for_het:
            return snp.num_het / n_called
        else:
            return None

    @property
    def exp_het(self):
        snp = self.record
        if snp.num_called >= self.min_calls_for_het:
            return snp.heterozygosity
        else:
            return None

    @property
    def samples(self):
        caller = self.snp_caller
        return [Call(sample, caller) for sample in self.record.samples]


class Call(object):
    def __init__(self, call, snp_caller=None):
        self.call = call
        self.snp_caller = snp_caller
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
        snp_caller = self.snp_caller
        if snp_caller is None:
            msg = 'To parse the read depths you have to set the snp_caller'
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
