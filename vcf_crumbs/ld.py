from __future__ import division

from vcf_crumbs.statistics import choose_samples
from collections import Counter, namedtuple

# Missing docstring
# pylint: disable=C0111

HaploCount = namedtuple('HaploCount', ['AB', 'ab', 'Ab', 'aB'])


def calculate_r_sqr(snp1, snp2, samples=None):
    calls1 = choose_samples(snp1.record, sample_names=samples)
    calls2 = choose_samples(snp2.record, sample_names=samples)
    haplo_counts = _count_biallelic_haplotypes(calls1, calls2)

    if haplo_counts is None:
        return None
    return _calculate_r_sqr(haplo_counts)


def _calculate_r_sqr(haplo_counts):
    # Invalid name. Due to using uppercases
    # pylint: disable=C0103
    total_samples = sum(haplo_counts)
    cnt_AB = haplo_counts.AB
    cnt_Ab = haplo_counts.Ab
    cnt_aB = haplo_counts.aB
    cnt_ab = haplo_counts.ab

    freqA = (cnt_AB + cnt_Ab) / total_samples
    freqB = (cnt_AB + cnt_aB) / total_samples

    if freqA == 0 or freqB == 0 or freqA == 1 or freqB == 1:
        return None

    rsqr = ((cnt_AB / total_samples) * (cnt_ab / total_samples))
    rsqr -= ((cnt_aB / total_samples) * (cnt_Ab / total_samples))
    rsqr *= rsqr
    rsqr /= ((freqA * (1 - freqA)) * (freqB * (1 - freqB)))
    return rsqr


def _most_freq_alleles(calls):
    gts_counts = Counter()
    gts_counts.update([al for call in calls for al in call.gt_alleles])
    return [al for al, count in gts_counts.most_common(2)]


def _remove_no_calls_and_hets(calls1, calls2):
    filtered_calls1 = []
    filtered_calls2 = []
    for call1, call2 in zip(calls1, calls2):
        if (not call1.called or not call2.called or call1.is_het or
            call2.is_het):
            continue
        filtered_calls1.append(call1)
        filtered_calls2.append(call2)
    return filtered_calls1, filtered_calls2


def _count_biallelic_haplotypes(calls1, calls2):
    # Invalid name. Due to using uppercases
    # pylint: disable=C0103
    calls1, calls2 = _remove_no_calls_and_hets(calls1, calls2)

    freq_alleles = []
    freq_alleles.append(_most_freq_alleles(calls1))
    freq_alleles.append(_most_freq_alleles(calls2))

    haplo_count = Counter()
    for call1, call2 in zip(calls1, calls2):
        al_snp_1 = call1.gt_alleles[0]
        # We're transforming all markers in biallelic
        if al_snp_1 != freq_alleles[0][0]:
            al_snp_1 = freq_alleles[0][1]
        al_snp_2 = call2.gt_alleles[0]
        if al_snp_2 != freq_alleles[1][0]:
            al_snp_2 = freq_alleles[1][1]

        haplo = (al_snp_1, al_snp_2)
        haplo_count[haplo] += 1

    alleles_snp1, alleles_snp2 = zip(*[(haplo[0], haplo[1]) for haplo in haplo_count.keys()])
    alleles_snp1 = set(alleles_snp1)
    alleles_snp2 = set(alleles_snp2)

    haplo_AB = haplo_count.most_common(1)[0][0]
    allele_A = haplo_AB[0]
    allele_B = haplo_AB[1]
    
    non_A_alleles = list(alleles_snp1.difference(allele_A))
    non_B_alleles = list(alleles_snp2.difference(allele_B))

    if not non_A_alleles or not non_B_alleles:
        return None
    allele_a = non_A_alleles[0]
    allele_b = non_B_alleles[0]

    count_AB = haplo_count.get(haplo_AB, 0)
    count_Ab = haplo_count.get((allele_A, allele_b), 0)
    count_aB = haplo_count.get((allele_a, allele_B), 0)
    count_ab = haplo_count.get((allele_a, allele_b), 0)

    return HaploCount(count_AB, count_ab, count_Ab, count_aB)
