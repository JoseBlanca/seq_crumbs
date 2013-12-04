from __future__ import division

import json
from os.path import join
from vcf import Reader

from Bio import SeqIO
from Bio.Restriction.Restriction import CommOnly, RestrictionBatch, Analysis
from vcf_crumbs.prot_change import (get_amino_change, IsIndelError,
                                    BetweenSegments, OutsideAlignment)
from vcf_crumbs.utils import DATA_DIR

COMMON_ENZYMES = ['EcoRI', 'SmaI', 'BamHI', 'AluI', 'BglII', 'SalI', 'BglI',
                  'ClaI', 'TaqI', 'PstI', 'PvuII', 'HindIII', 'EcoRV',
                  'HaeIII', 'KpnI', 'ScaI', 'HinfI', 'DraI', 'ApaI', 'BstEII',
                  'ZraI', 'BanI', 'Asp718I']


def _calculate_maf_for_counts(counts):
    counts = counts.values()
    dp = sum(counts)

    freqs = sorted([count / float(dp) for count in counts])
    return freqs[-1] if freqs else None


def calculate_maf(record):
    counts = count_alleles(record)['all']
    return _calculate_maf_for_counts(counts)


def aggregate_allele_counts(allele_counts):
    all_counts = {}
    for values in allele_counts.values():
        for allele, count in values.items():
            if allele not in all_counts:
                all_counts[allele] = 0
            all_counts[allele] += count
    return {'all': all_counts}


def count_alleles(record, sample_names=None):
    if sample_names is None:
        choosen_samples = record.samples
    else:
        filter_by_name = lambda x: True if x.sample in sample_names else False
        choosen_samples = filter(filter_by_name, record.samples, )

    counts = {}
    alleles = _str_alleles(record)
    for call in choosen_samples:
        if call.gt_bases is None:
            continue
        if call.sample not in counts:
            counts[call.sample] = {}
#         if call.gt_alleles == ['1', '2']:
#             msg = "This code assumes that there aren't genotypes with both"
#             msg += " alleles being alternate alleles"
#             print record.CHROM, record.POS
#             raise NotImplementedError(msg)
        for index, genotype in enumerate(call.gt_alleles):
            genotype = int(genotype)
            allele = alleles[genotype]
            if allele not in counts:
                counts[call.sample][allele] = 0
            if hasattr(call.data, 'RD'):
                #VARSCAN
                allele_counts = call.data.RD if genotype == 0 else call.data.AD
            else:
                #GATK
                allele_counts = call.data.AD[index]
            counts[call.sample][allele] += allele_counts
        if not counts[call.sample]:
            del counts[call.sample]
    if not sample_names:
        return aggregate_allele_counts(counts)

    return counts


def _str_alleles(record):
    alleles = [record.alleles[0]]
    for alt_allele in record.alleles[1:]:
        alleles.append(alt_allele.sequence)
    return alleles


class BaseFilter(object):
    def __call__(self, record):
        raise NotImplementedError()

    @property
    def encode_conf(self):
        encoded_conf = json.dumps(self.conf).replace("\"", "'")
        return encoded_conf

    def _clean_filter(self, record):
        name = self.name
        if record.FILTER and name in record.FILTER:
            record.FILTER.remove(name)

    @property
    def info(self):
        return {}

    @property
    def info_id(self):
        if self.info:
            return self.info['id']
        return None


class CloseToSnvFilter(BaseFilter):
    '''Filter snps with other close snvs.

    Allowed snv_types:  [snp, indel, unknown] '''

    def __init__(self, vcf_fpath, distance=60, max_maf=None, snv_type=None):
        self.distance = distance
        self.reader = Reader(open(vcf_fpath))
        self.maf = max_maf
        self.snv_type = snv_type
        self.conf = {'distance': distance, 'maf': max_maf,
                     'snv_type': snv_type}

    def __call__(self, record):
        self._clean_filter(record)
        chrom = record.CHROM
        pos = record.POS
        start = pos - self.distance if pos - self.distance > 0 else 0
        end = pos + self.distance
        snv_type = self.snv_type
        maf = self.maf
        passed_snvs = 0
        for snv in self.reader.fetch(chrom, start, end):
            if snv.POS == pos:
                continue
            if maf is None and snv_type is None:
                passed_snvs += 1
            elif maf is None and snv_type is not None:
                if snv_type == snv.var_type:
                    passed_snvs += 1
            elif maf is not None  and snv_type is None:
                calculated_maf = calculate_maf(snv)
                if calculated_maf and calculated_maf < maf:
                    passed_snvs += 1
            else:
                calculated_maf = calculate_maf(snv)
                if calculated_maf and calculated_maf < maf and snv_type == snv.var_type:
                    passed_snvs += 1
        if passed_snvs:
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        snv_type = '' if self.snv_type is None else self.snv_type[0]
        maf_str = '' if self.maf is None else '_{:.2f}'.format(self.maf)
        return 'cs{}{}{}'.format(snv_type, self.distance, maf_str)

    @property
    def description(self):
        snv_type = 'snv' if self.snv_type is None else self.snv_type
        maf_str = ', with maf:{:.2f}'.format(self.maf) if self.maf else ''
        desc = 'The snv is closer than {} nucleotides to another {}{} :: {}'
        return desc.format(self.distance, snv_type, maf_str, self.encode_conf)


def _get_lengths(fhand):
    lengths = {}
    for seq in SeqIO.parse(fhand, 'fasta'):
        lengths[seq.id] = len(seq)
    return lengths


class HighVariableRegionFilter(BaseFilter):
    'Filter depending on the variability of the region'

    def __init__(self, max_variability, vcf_fpath, ref_fpath, window=None):
        self.max_variability = max_variability
        self.hg_reader = Reader(open(vcf_fpath))
        self.window = window
        self._lengths = _get_lengths(open(ref_fpath))
        self.conf = {'max_variability': max_variability, 'window': window}

    def __call__(self, record):
        self._clean_filter(record)
        chrom = record.CHROM
        pos = record.POS
        seq_len = self._lengths[chrom]
        if not self.window:
            start = 0
            end = seq_len - 1
            window_len = seq_len
        else:
            w2 = self.window / 2
            start = int(pos - w2)
            if start < 0:
                start = 0
            end = int(pos + w2)
            if end > seq_len - 1:
                end = seq_len - 1
            window_len = self.window

        num_snvs = len(list(self.hg_reader.fetch(chrom, start, end)))

        freq = num_snvs / window_len

        if freq > self.max_variability:
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        return 'hv{}'.format(self.max_variability)

    @property
    def description(self):
        desc = 'The region has more than {} snvs per 100 bases:: {}'
        return desc.format(self.max_variability * 100, self.encode_conf)


class CloseToLimitFilter(BaseFilter):
    'Filter snps close to chromosome limits'

    def __init__(self, distance, ref_fpath):
        self.distance = distance
        self._lengths = _get_lengths(open(ref_fpath))
        self.conf = {'distance': distance}

    def __call__(self, record):
        self._clean_filter(record)
        chrom = record.CHROM
        pos = record.POS
        seq_len = self._lengths[chrom]
        if pos < self.distance or seq_len - pos < self.distance:
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        return 'cl{}'.format(self.distance)

    @property
    def description(self):
        desc = 'The snv is closer than {} nucleotides to the reference edge'
        desc += ':: {}'
        return desc.format(self.distance, self.encode_conf)


class MafLimitFilter(BaseFilter):
    'Filter by maf'

    def __init__(self, max_maf, samples=None):
        self.max_maf = max_maf
        self.samples = samples
        self.conf = {'max_maf': max_maf, 'samples': samples}

    def __call__(self, record):
        self._clean_filter(record)
        allele_counts = count_alleles(record, sample_names=self.samples)
        all_counts = aggregate_allele_counts(allele_counts)
        maf = _calculate_maf_for_counts(all_counts['all'])
        if max and maf > self.max_maf:
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        return 'maf{}'.format(self.max_maf)

    @property
    def description(self):
        samples = 'all' if self.samples is None else ','.join(self.samples)
        desc = 'The most frequent allele in samples: {}.'
        desc += 'frequency greater than {} :: {}'
        return desc.format(samples, self.max_maf, self.encode_conf)


class CapEnzymeFilter(BaseFilter):
    'Filters by detectablity by restriction enzymes'

    def __init__(self, all_enzymes, ref_fpath):
        self.all_enzymes = all_enzymes
        self.ref_index = SeqIO.index(ref_fpath, 'fasta')
        self.conf = {'all_enzymes': all_enzymes}

    def __call__(self, record):
        self._clean_filter(record)
        alleles = _str_alleles(record)
        enzymes = set()
        #we have to make all the posible conbinations
        ref = self.ref_index[record.CHROM]
        used_combinations = []
        for i_index in range(len(alleles)):
            for j_index in range(len(alleles)):
                allelei = alleles[i_index]
                allelej = alleles[j_index]
                if (i_index == j_index or
                    (allelei, allelej) in used_combinations or
                    (allelej, allelei) in used_combinations):
                    continue
                used_combinations.append((allelei, allelej))
                i_j_enzymes = _cap_enzymes_between_alleles(allelei, allelej,
                                                           ref, record.POS,
                                                           record.end,
                                                  all_enzymes=self.all_enzymes)
                enzymes = enzymes.union(i_j_enzymes)
        if not enzymes:
            record.add_filter(self.name)
        else:
            enzymes = [str(e) for e in list(enzymes)]
            record.add_info(info=self.info['id'], value=','.join(enzymes))
        return record

    @property
    def name(self):
        return 'ce{}'.format('t' if self.all_enzymes else 'f')

    @property
    def description(self):
        enzs = 'all' if self.all_enzymes else 'cheap_ones'
        desc = "SNV is not a CAP detectable by the enzymes: {}:: {}"
        return desc.format(enzs, self.encode_conf)

    @property
    def info(self):
        return {'id': 'CAP', 'num': 1, 'type': 'String',
                'desc': 'Enzymes that can be use to differentiate the alleles'}


def _cap_enzymes_between_alleles(allele1, allele2, reference, start, end,
                                 all_enzymes=False):
    '''It looks in the enzymes that differenciate the given alleles.

    It returns a set.
    '''
    # we have to build the two sequences
    if all_enzymes:
        restriction_batch = CommOnly
    else:
        restriction_batch = RestrictionBatch(COMMON_ENZYMES)

    sseq = reference.seq
    post_seq_start = start - 100 if start - 100 > 0 else 0
    prev_seq = sseq[post_seq_start: start - 1]
    post_seq = sseq[end: end + 100]

    seq1 = prev_seq + allele1 + post_seq
    seq2 = prev_seq + allele2 + post_seq
    anal1 = Analysis(restriction_batch, seq1, linear=True)
    enzymes1 = set(anal1.with_sites().keys())
    anal1 = Analysis(restriction_batch, seq2, linear=True)
    enzymes2 = set(anal1.with_sites().keys())
    enzymes = set(enzymes1).symmetric_difference(set(enzymes2))
    return enzymes


class KindFilter(BaseFilter):
    '''Choose Snv by type

    type options = snp, indel, unknown
    '''

    def __init__(self, kind='snp'):
        self.kind = kind
        self.conf = {'kind': kind}

    def __call__(self, record):
        self._clean_filter(record)
        rec_type = record.var_type
        if rec_type != self.kind:
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        return 'vk{}'.format(self.kind[0].lower())

    @property
    def description(self):
        return 'It is not an {} :: {}'.format(self.kind, self.encode_conf)


class IsVariableFilter(BaseFilter):
    'Variable in readgroup'

    def __init__(self, filter_id, samples, max_maf=None, min_reads=6,
                 min_reads_per_allele=3, in_union=True, in_all_groups=True,
                 reference_free=True
                 ):
        self.max_maf = max_maf
        self.samples = samples
        self.min_reads = min_reads
        self.min_reads_per_allele = min_reads_per_allele
        self.in_union = in_union
        self.in_all_groups = in_all_groups
        self.reference_free = reference_free
        self.filter_id = filter_id
        self.conf = {'max_maf': max_maf, 'samples': samples,
                     'min_reads': min_reads, 'in_union': in_union,
                     'min_reads_per_allele': min_reads_per_allele,
                     'in_all_groups': in_all_groups,
                     'reference_free': reference_free}

    def __call__(self, record):
        self._clean_filter(record)
        is_variable = variable_in_samples(record, self.samples, self.in_union,
                                          self.max_maf, self.in_all_groups,
                                          self.reference_free, self.min_reads,
                                          self.min_reads_per_allele)
        if not is_variable:
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        return 'vs{}'.format(self.filter_id)

    @property
    def description(self):
        samples = ','.join(self.samples)
        desc = 'It is not variable, or no data, in the samples : {}. '
        desc += 'All together: {} :: {}'
        return desc.format(samples, str(self.in_union), self.encode_conf)


class IsNotVariableFilter(BaseFilter):
    'Variable in readgroup'

    def __init__(self, filter_id, samples, max_maf=None, min_reads=6,
                 in_union=True, reference_free=True):
        self.max_maf = max_maf
        self.samples = samples
        self.min_reads = min_reads
        self.in_union = in_union
        self.reference_free = reference_free
        self.filter_id = filter_id
        self.conf = {'max_maf': max_maf, 'samples': samples,
                     'min_reads': min_reads, 'in_union': in_union,
                     'reference_free': reference_free}

    def __call__(self, record):
        self._clean_filter(record)
        is_invariable = invariant_in_samples(record, self.samples,
                                    in_union=self.in_union, maf=self.max_maf,
                                    reference_free=self.reference_free,
                                    min_reads=self.min_reads)
        if not is_invariable:
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        return 'nvs{}'.format(self.filter_id)

    @property
    def description(self):
        samples = ','.join(self.samples)
        desc = 'It is variable, or no data, in the samples : {}. '
        desc += 'All together: {} :: {}'
        return desc.format(samples, str(self.in_union), self.encode_conf)


def invariant_in_samples(record, samples, in_union=True,
                           in_all_groups=True, reference_free=True, maf=None,
                           min_reads=None):
    'it check if the given snv is invariant form the given groups'
    allele_counts = count_alleles(record, samples)

    if not allele_counts:
        return False if reference_free else True

    if in_union:
        allele_counts = aggregate_allele_counts(allele_counts)

    invariable_in_samples = []
    alleles = _str_alleles(record)
    for sample_counts in allele_counts.values():
        invariable_in_sample = not _is_variable_in_sample(alleles,
                                                 sample_counts,
                                                 reference_free=reference_free,
                                                 maf=maf,
                                                 min_reads=min_reads,
                                                 min_reads_per_allele=1)
        invariable_in_samples.append(invariable_in_sample)

    if in_all_groups:
        return all(invariable_in_samples)
    else:
        return any(invariable_in_samples)


def variable_in_samples(record, samples, in_union=True, maf=None,
                        in_all_groups=True, reference_free=True,
                        min_reads=6, min_reads_per_allele=3):
    allele_counts = count_alleles(record, samples)

    if in_union:
        allele_counts = aggregate_allele_counts(allele_counts)

    variable_in_samples = []
    alleles = _str_alleles(record)
    for sample_counts in allele_counts.values():
        variable_in_sample = _is_variable_in_sample(alleles, sample_counts,
                                                reference_free=reference_free,
                                                maf=maf, min_reads=min_reads,
                                     min_reads_per_allele=min_reads_per_allele)

        variable_in_samples.append(variable_in_sample)
    if in_all_groups:
        return all(variable_in_samples)
    else:
        return any(variable_in_samples)


def _is_variable_in_sample(alleles, allele_counts, reference_free, maf,
                           min_reads, min_reads_per_allele):
    'It checks if the allele is variable'

    num_alleles = 0
    # first remove  the bad alleles
    for allele in allele_counts:
        if allele_counts[allele] >= min_reads_per_allele:
            num_alleles += 1
    if not num_alleles:
        return False
    ref_allele = alleles[0]

    if (ref_allele in allele_counts and
        allele_counts[ref_allele] >= min_reads_per_allele):
        ref_in_alleles = True
    else:
        ref_in_alleles = False

    if num_alleles == 1:
        if reference_free:
            return False
        elif not reference_free and ref_in_alleles:
            return False

    num_reads = sum(allele_counts.values())
    maf_allele = _calculate_maf_for_counts(allele_counts)
    if ((maf and maf < maf_allele) or
        (min_reads and min_reads > num_reads)):
        return False

    return True


class ChangeAminoFilter(BaseFilter):
    'The aminoacid changes with the alternative allele'

    def __init__(self, ref_fpath, orf_seq_fpath):
        self.orf_suffix = '_orf_seq'
        self.ref_index = SeqIO.index(ref_fpath, 'fasta')
        self.orf_seq_index = SeqIO.index(orf_seq_fpath, 'fasta')

        self.conf = {}

    def __call__(self, record):
        self._clean_filter(record)
        seq_name = record.CHROM
        seq_ref = self.ref_index[seq_name]
        aminos = None
        try:
            seq_estscan = self.orf_seq_index[seq_name + self.orf_suffix]
        except KeyError:
            return record
        try:
            aminos = get_amino_change(seq_ref, seq_estscan, record)
        except IsIndelError:
            record.add_filter(self.name)
            record.add_info(info=self.info_id, value='INDEL')
        except BetweenSegments:
            pass
        except OutsideAlignment:
            pass
        if aminos is None:
            return record

        if set(aminos['ref_amino']) != set(aminos['alt_amino']):
            record.add_filter(self.name)
            info_val = '{}->{}'.format(aminos['ref_amino'],
                                       ','.join(aminos['alt_amino']))
            record.add_info(info=self.info_id, value=info_val)
        return record

    @property
    def name(self):
        return 'caa'

    @property
    def description(self):
        return "Alt alleles change the amino acid"

    @property
    def info(self):
        return {'id': 'AAC', 'num': 1, 'type': 'String',
                'desc': 'Amino acid change'}


def _parse_blossum_matrix(path):
    matrix = {}
    fhand = open(path)
    header = fhand.readline().strip().split('\t')[1:]
    header = [h.strip() for h in header]
    for line in fhand:
        items = line.strip().split('\t')
        aa1 = items[0].strip()
        for aa2, value in zip(header, items[1:]):
            matrix[(aa1, aa2)] = int(value)
    return matrix


class ChangeAminoSeverityFilter(BaseFilter):
    'Check if the change is severe  or not'

    def __init__(self, ref_fpath, orf_seq_fpath):
        self.orf_suffix = '_orf_seq'
        self.ref_index = SeqIO.index(ref_fpath, 'fasta')
        self.orf_seq_index = SeqIO.index(orf_seq_fpath, 'fasta')
        blossum_path = join(DATA_DIR, 'blossum90.csv')
        self.blosum = _parse_blossum_matrix(blossum_path)

        self.conf = {}

    def _is_severe(self, ref_aa, alt_aa):
        return True if self.blosum[(ref_aa, alt_aa)] < 0 else False

    def __call__(self, record):
        self._clean_filter(record)
        seq_name = record.CHROM
        seq_ref = self.ref_index[seq_name]
        aminos = None
        try:
            seq_estscan = self.orf_seq_index[seq_name + self.orf_suffix]
        except KeyError:
            return record
        try:
            aminos = get_amino_change(seq_ref, seq_estscan, record)
        except IsIndelError:
            record.add_filter(self.name)
        except BetweenSegments:
            pass
        except OutsideAlignment:
            pass

        if aminos is None:
            return record

        ref_aa = aminos['ref_amino']
        alt_aas = aminos['alt_amino']
        if set(ref_aa) == set(alt_aas):
            return record
        if any([self._is_severe(ref_aa, alt_aa) for alt_aa in alt_aas]):
            record.add_filter(self.name)
        return record

    @property
    def name(self):
        return 'cas'

    @property
    def description(self):
        return "Alt alleles change the amino acid and the change is severe"
