from __future__ import division


from collections import Counter
from array import array


from  vcf import Reader

from crumbs.seq import get_name, get_length
from crumbs.seqio import read_seqs


VARSCAN = 'VarScan'
GATK = 'gatk'

# taken from pyvcf
HOM_REF = 0
HET = 1
HOM_ALT = 2


def get_snpcaller_name(reader):
    metadata = reader.metadata
    if 'source' in metadata and 'VarScan2' in metadata['source']:
        return VARSCAN
    if 'UnifiedGenotyper' in metadata:
        return GATK
    return None


def _get_call_data(call, snpcaller):
    data = call.data
    gt = data.GT
    gq = data.GQ
    dp = data.DP
    if snpcaller == GATK:
        rd = data.AD[0]
        ad = sum(data.AD[1:])
    elif snpcaller == VARSCAN:
        rd = data.RD
        ad = data.AD
    return gt, gq, dp, rd, ad


def calculate_maf(snp, snpcaller):
    total_ad = 0
    total_rd = 0
    for call in snp.samples:
        if call.called:
            rd, ad = _get_call_data(call, snpcaller)[3:]
            total_ad += ad
            total_rd += rd
    values = [total_ad, total_rd]
    total = sum(values)
    if not total:
        return None
    maf = max(values) / total
    return maf


def get_data_from_vcf(vcf_path):
    reader = Reader(filename=vcf_path)
    snpcaller = get_snpcaller_name(reader)
    data = {'snps_per_chromo': Counter(),
            'maf_per_snp': [],
            'call_data': {HOM_REF: {'x': array('B'), 'y': array('B'),
                                    'value': array('B')},
                          HET: {'x': array('B'), 'y': array('B'),
                                'value': array('B')},
                          HOM_ALT: {'x': array('B'), 'y': array('B'),
                                    'value': array('B')}
                         }
           }

    chrom_counts = data['snps_per_chromo']
    mafs = data['maf_per_snp']
    call_datas = data['call_data']
    for snp in reader:
        chrom_counts[snp.CHROM] += 1
        mafs.append(calculate_maf(snp, snpcaller))
        #sample_data
        for call in snp.samples:
            if call.called:
                gt, gq, dp, rd, ad = _get_call_data(call, snpcaller)
                call_datas[call.gt_type]['x'].append(rd)
                call_datas[call.gt_type]['y'].append(ad)
                call_datas[call.gt_type]['value'].append(gq)

    return data


def _get_seq_lengths(fhand):
    return {get_name(seq): get_length(seq) for seq in read_seqs([fhand])}


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
