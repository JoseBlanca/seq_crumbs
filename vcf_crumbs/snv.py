
from vcf import Reader as pyvcf_Reader  # Imported here just for convenience
                                        # when importing from other modules

# Missing docstring
# pylint: disable=C0111


class SNV(object):
    '''A wrapper around the pyvcf _Record with some additional functionality'''
    def __init__(self, record, min_calls_for_het=10):
        self.record = record
        self.min_calls_for_het = min_calls_for_het

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

    def __getattr__(self, attrname):
        parent_attr = getattr(self.record, attrname)
        return parent_attr
