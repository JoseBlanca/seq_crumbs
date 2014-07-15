from vcf_crumbs.snv import VCFReader, VCFWriter


def run_genotype_filters(in_fhand, out_fhand, gt_filters, template_fhand=None):
    reader = VCFReader(in_fhand)
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
