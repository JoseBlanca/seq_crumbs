from __future__ import division

from vcf_crumbs.snv import VCFReader, VCFWriter

# Missing docstring
# pylint: disable=C0111
# Too few public methods
# pylint: disable=R0903

DEF_PROB_AA_THRESHOLD=0.9999
HW = 'hw'


def run_genotype_filters(in_fhand, out_fhand, gt_filters, template_fhand=None,
                         reader_kwargs=None):
    if reader_kwargs is None:
        reader_kwargs = {}
    reader = VCFReader(in_fhand, **reader_kwargs)

    templa_reader = reader if template_fhand is None else VCFReader(template_fhand)
    writer = VCFWriter(out_fhand, template_reader=templa_reader)

    for snv in reader.parse_snvs():
        for mapper in gt_filters:
            snv = mapper(snv)
        writer.write_snv(snv)


class LowQualityGenotypeFilter(object):

    def __init__(self, min_qual):
        self._min_qual = min_qual

    def __call__(self, snv):
        return snv.remove_gt_from_low_qual_calls(min_qual=self._min_qual)


class HetGenotypeFilter(object):

    def __call__(self, snv):
        return snv.remove_gt_from_het_calls()


class LowEvidenceAlleleFilter(object):
    def __init__(self, prob_aa_threshold=DEF_PROB_AA_THRESHOLD,
                 genotypic_freqs_method=HW):
        self._min_prob = prob_aa_threshold
        self.genotypic_freqs_method = genotypic_freqs_method

    def __call__(self, snv):
        allele_freqs = snv.allele_freqs
        if not allele_freqs:
            def set_all_gt_to_none(call):
                return call.copy_setting_gt(gt=None, return_pyvcf_call=True)
            return snv.copy_mapping_calls(set_all_gt_to_none)

        genotypic_freqs_method = self.genotypic_freqs_method
        calls = []
        for call in snv.calls:
            if not call.called:
                filtered_call = call.call
            else:
                alleles = call.int_alleles
                if len(set(alleles)) > 1:
                    filtered_call = call.call
                else:
                    allele = alleles[0]
                    freq = allele_freqs[allele]
                    allele_depths = call.allele_depths
                    if not allele_depths:
                        msg = 'Allele depths are required for the lowEvidence'
                        msg += 'Allele filter'
                        raise RuntimeError(msg)
                    depth = call.allele_depths[allele]
                    if genotypic_freqs_method == HW:
                        prob = prob_aa_given_n_a_reads(depth, freq)
                    else:
                        msg = 'Method not implemented for genotypic freqs: '
                        msg += genotypic_freqs_method
                        raise NotImplementedError(msg)
                    if prob >= self._min_prob:
                        filtered_call = call.call
                    else:
                        geno = call.call.data.GT[:-1] + '.'
                        filtered_call = call.copy_setting_gt(gt=geno,
                                                        return_pyvcf_call=True)
            calls.append(filtered_call)
        return snv.copy_mapping_calls(calls)


def prob_aa_given_n_a_reads(num_a_reads, freq_a_in_pop):
    'It assumes HW'
    # TODO fix for Backcross
    proba = freq_a_in_pop
    res = proba**2
    res /= proba**2 + 0.5**(num_a_reads - 1) * proba * (1 - proba)
    return res

