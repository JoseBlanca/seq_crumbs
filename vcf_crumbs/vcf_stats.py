from __future__ import division


from collections import Counter
from array import array


from  vcf import Reader

from crumbs.seq import get_name, get_length
from crumbs.seqio import read_seqs


VARSCAN = 'VarScan'
GATK = 'gatk'
FREEBAYES = 'freebayes'

# taken from pyvcf
HOM_REF = 0
HET = 1
HOM_ALT = 2


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


def _get_call_data(call, vcf_variant):
    data = call.data
    gt = data.GT
    gq = int(data.GQ)
    dp = data.DP
    if vcf_variant == GATK:
        rd = data.AD[0]
        ad = sum(data.AD[1:])
    elif vcf_variant == VARSCAN:
        rd = data.RD
        ad = data.AD
    elif vcf_variant == FREEBAYES:
        rd = data.RO
        ad = data.AO
        if isinstance(ad, list):
            ad = sum(ad)
    else:
        raise NotImplementedError('Not using one of the supported snp callers')
    return gt, gq, dp, rd, ad


def calculate_maf_old(snp, vcf_variant):
    total_ad = 0
    total_rd = 0
    for call in snp.samples:
        if call.called:
            rd, ad = _get_call_data(call, vcf_variant)[3:]
            total_ad += ad
            total_rd += rd
    values = [total_ad, total_rd]
    total = sum(values)
    if not total:
        return None
    maf = max(values) / total
    return maf


def calculate_maf(snp, vcf_variant):
    total_ad = 0
    total_rd = 0
    mafs = {}
    for call in snp.samples:
        if call.called:
            rd, ad = _get_call_data(call, vcf_variant)[3:]
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


def get_data_from_vcf(vcf_path):
    reader = Reader(filename=vcf_path)
    vcf_variant = get_snpcaller_name(reader)
    typecode = 'I'
    data = {'samples': set(),
            'snps_per_chromo': Counter(),
            'maf_per_snp': [],
            'call_data': {HOM_REF: {'x': array(typecode), 'y': array(typecode),
                                    'value': array(typecode)},
                          HET: {'x': array(typecode), 'y': array(typecode),
                                'value': array(typecode)},
                          HOM_ALT: {'x': array(typecode), 'y': array(typecode),
                                    'value': array(typecode)}
                         }
           }

    chrom_counts = data['snps_per_chromo']
    mafs = data['maf_per_snp']
    call_datas = data['call_data']
    samples = data['samples']
    for snp in reader:
        chrom_counts[snp.CHROM] += 1
        maf = calculate_maf(snp, vcf_variant=vcf_variant)
        if maf is not None:
            mafs.append(maf)
        #sample_data
        for call in snp.samples:
            samples.add(call.sample)
            if call.called:
                gt, gq, dp, rd, ad = _get_call_data(call, vcf_variant)
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
