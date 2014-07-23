from __future__ import division

import sys
from collections import Counter, OrderedDict
import gzip

from vcf import Reader as pyvcfReader
from vcf import Writer as pyvcfWriter
from vcf.model import make_calldata_tuple
# ouch, _Call is a private class, but we don't know how to modify a Call
from vcf.model import _Call as pyvcfCall
from vcf.model import _Record as pyvcfRecord

from vcf_crumbs.iterutils import generate_windows
from crumbs.seqio import read_seqs
from crumbs.seq import get_name, get_length

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
DEF_MIN_NUM_SNPS_IN_WIN = 5

# TODO check if SNV can be converted in a proxy using the recipes
# http://code.activestate.com/recipes/496741-object-proxying/
# http://code.activestate.com/recipes/496742-shelfproxy/


class _SNPQueue(object):
    def __init__(self):
        self.queue = []

    def empty(self):
        self.queue = []

    def pop(self, location):
        queue = self.queue
        index_to_remove = None
        for index, snp in enumerate(queue):
            if snp.pos < location:
                index_to_remove = index
            else:
                break
        if index_to_remove is not None:
            del queue[:index_to_remove + 1]

    def extend(self, snps):
        self.queue.extend(snps)


class _SNPSlidingWindow(object):
    def __init__(self, snp_reader, win_size, win_step, min_num_snps,
                 ref_fhand=None):
        self._snp_queue = _SNPQueue()
        self._reader = snp_reader
        self.win_size = win_size
        self.win_step = win_step
        self.min_num_snps = min_num_snps
        self._ref_fhand = ref_fhand

    def _snp_in_window(self, snp, win):
        if win[0] <= snp.pos < win[1]:
            return True
        else:
            return False

    def _get_chrom_lengths(self):
        chrom_lens = OrderedDict()
        if self._ref_fhand is None:
            vcf_fhand = gzip.open(self._reader.fhand.name)
            for line in vcf_fhand:
                line = line.strip()
                if line.startswith('#'):
                    continue
                items = line.split()
                chrom = items[0]
                loc = int(items[1])
                if chrom not in chrom_lens:
                    chrom_lens[chrom] = loc
                else:
                    if loc > chrom_lens[chrom]:
                        chrom_lens[chrom] = loc

        else:
            for read in read_seqs([self._ref_fhand]):
                chrom_lens[get_name(read)] = get_length(read)
        return chrom_lens

    def windows(self):
        chrom_lengths = self._get_chrom_lengths()
        snp_queue = self._snp_queue
        for chrom, chrom_length in chrom_lengths.items():
            wins = generate_windows(start=0,
                                    size=self.win_size, step=self.win_step,
                                    end=chrom_length + 1)
            snp_queue.empty()
            for win in wins:
                snp_queue.pop(win.start)
                if snp_queue.queue:
                    new_strech_start = snp_queue.queue[-1].pos + 1
                else:
                    new_strech_start = win.start
                new_snps = self._reader.fetch_snvs(chrom, new_strech_start,
                                                   win.end)
                snp_queue.extend(new_snps)
                if len(snp_queue.queue) >= self.min_num_snps:
                    yield {'chrom': chrom, 'start': win.start, 'end': win.end,
                           'snps': snp_queue.queue[:]}


class VCFReader(object):
    def __init__(self, fhand,
                 min_calls_for_pop_stats=DEF_MIN_CALLS_FOR_POP_STATS):
        self.fhand = fhand
        self.pyvcf_reader = pyvcfReader(fsock=fhand)
        self.min_calls_for_pop_stats = min_calls_for_pop_stats
        self._snpcaller = None

    def parse_snvs(self):
        min_calls_for_pop_stats = self.min_calls_for_pop_stats
        last_snp = None
        try:
            for snp in self.pyvcf_reader:
                snp = SNV(snp, reader=self,
                          min_calls_for_pop_stats=min_calls_for_pop_stats)
                last_snp = snp
                yield snp
        except:
            if last_snp is not None:
                chrom = str(last_snp.chrom)
                pos = str(last_snp.pos)
                msg = 'Last parsed SNP was: ' + str(chrom) + ' ' + str(pos + 1)
                msg += '\n'
                sys.stderr.write(msg)
                raise

    def fetch_snvs(self, chrom, start, end=None):
        min_calls_for_pop_stats = self.min_calls_for_pop_stats
        try:
            snvs = self.pyvcf_reader.fetch(chrom, start + 1, end=end)
        except KeyError:
            snvs = []

        for snp in snvs:
            snp = SNV(snp, reader=self,
                      min_calls_for_pop_stats=min_calls_for_pop_stats)
            yield snp

    def sliding_windows(self, size, step=None, ref_fhand=None,
                        min_num_snps=DEF_MIN_NUM_SNPS_IN_WIN):
        random_snp_reader = VCFReader(open(self.fhand.name))
        sliding_window = _SNPSlidingWindow(snp_reader=random_snp_reader,
                                           win_size=size, win_step=step,
                                           min_num_snps=min_num_snps,
                                           ref_fhand=ref_fhand)
        for window in sliding_window.windows():
            yield window

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
            else:
                snpcaller = metadata['source']
        elif 'UnifiedGenotyper' in metadata:
            snpcaller = GATK
        else:
            raise NotImplementedError('Can not get snp caller of the vcf file')
        self._snpcaller = snpcaller
        return snpcaller

    @property
    def samples(self):
        return self.pyvcf_reader.samples

    @property
    def filters(self):
        return self.pyvcf_reader.filters

    @property
    def infos(self):
        return self.pyvcf_reader.infos


class VCFWriter(pyvcfWriter):

    def __init__(self, stream, template_reader, lineterminator="\n"):
        template = template_reader.pyvcf_reader
        super(VCFWriter, self).__init__(stream, template,
                                        lineterminator=lineterminator)

    def write_snv(self, snv):
        super(VCFWriter, self).write_record(snv.record)


class SNV(object):
    '''A proxy around the pyvcf _Record with some additional functionality'''
    def __init__(self, record, reader, ploidy=2,
                 min_calls_for_pop_stats=DEF_MIN_CALLS_FOR_POP_STATS):
        # min_calls_for_pop_stats is defined here because otherwise it would
        # behave like a global variable for all SNPs read from a common
        # reader
        self.record = record
        self.reader = reader
        self.ploidy = ploidy
        self.min_calls_for_pop_stats = min_calls_for_pop_stats
        self._mac_analyzed = False
        self._maf = None
        self._mac = None
        self._maf_dp_analyzed = False
        self._allele_depths = None
        self._maf_depth = None
        self._depth = None

    @property
    def alleles(self):
        alleles = []
        for allele in self.record.alleles:
            try:
                allele = allele.sequence
            except AttributeError:
                pass
            alleles.append(allele)
        return alleles

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
    def calls(self):
        return [Call(sample, snv=self) for sample in self.record.samples]

    def get_call(self, sample):
        call = self.record.genotype(sample)
        return Call(call, self)

    def _calculate_maf_and_mac(self):
        if self._mac_analyzed:
            return
        self._mac_analyzed = True
        ploidy = self.ploidy
        snp = self.record

        n_chroms_sampled = 0
        allele_counts = Counter()
        for call in snp.samples:
            if call.called:
                genotype = call.gt_alleles
                assert len(genotype) == ploidy
                n_chroms_sampled += ploidy
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

    @property
    def inbreed_coef(self):
        obs_het = self.obs_het
        if obs_het is None:
            return None
        try:
            inbreed_f = 1 - (obs_het / self.exp_het)
            return inbreed_f
        except ZeroDivisionError:
            return None

    def _calculate_mafs_dp(self):
        if self._maf_dp_analyzed:
            return None
        self._maf_dp_analyzed = True

        allele_depths = Counter()
        for call in self.calls:
            if not call.has_alternative_counts:
                continue
            for allele, depth in call.allele_depths.items():
                allele_depths[allele] += depth
        depth = sum(allele_depths.values())
        if not depth:
            return None
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
    def chrom(self):
        return self.record.CHROM

    @property
    def pos(self):
        return self.record.POS - 1

    @property
    def end(self):
        return self.record.end

    @property
    def qual(self):
        return self.record.QUAL

    @property
    def id(self):
        return self.record.ID

    @property
    def ref(self):
        return self.record.REF

    @property
    def num_called(self):
        return self.record.num_called

    @property
    def call_rate(self):
        return self.record.num_called / len(self.reader.samples)

    @property
    def is_snp(self):
        return self.record.is_snp

    def __str__(self):
        return str(self.record)

    def __unicode__(self):
        return self.record.__unicode__()

    def __repr__(self):
        return repr(self.record)

    @property
    def filters(self):
        return self.record.FILTER

    @filters.setter
    def filters(self, value):
        self.record.FILTER = value

    def add_filter(self, filter_name):
        self.record.add_filter(filter_name)

    @property
    def infos(self):
        return self.record.INFO

    def add_info(self, *args, **kwargs):
        self.record.add_info(*args, **kwargs)

    @property
    def kind(self):
        # return snv type. [snp, indel, unknown]
        return self.record.var_type

    @property
    def is_indel(self):
        return self.record.is_indel

    @property
    def calldata_class(self):
        return make_calldata_tuple(self.record.FORMAT.split(':'))

    def _filter_out_calls(self, filter_):
        calls = []
        sample_indexes = {}
        for index, call in enumerate(self.calls):
            if filter_(call):
                call = call.copy_setting_gt_to_none(return_pyvcf_call=True)
            else:
                call = call.call
            calls.append(call)
            sample_indexes[call.sample] = index
        record = self.record
        record = pyvcfRecord(record.CHROM, record.POS, record.ID,
                             record.REF, record.ALT, record.QUAL,
                             record.FILTER, record.INFO, record.FORMAT,
                             sample_indexes, calls)
        snv = SNV(record, self.reader,
                  min_calls_for_pop_stats=self.min_calls_for_pop_stats)
        return snv

    def remove_gt_from_het_calls(self):
        def filter_(call):
            if call.is_het:
                return True
            else:
                return False
        return self._filter_out_calls(filter_)

    def remove_gt_from_low_qual_calls(self, min_qual):
        'It returns a new SNV with low qual call set to uncalled'
        def filter_(call):
            if min_qual is not None and call.gt_qual < min_qual:
                return True
            else:
                return False
        return self._filter_out_calls(filter_)

    def filter_calls_by_sample(self, samples, reverse=False, keep_info=False):

        if reverse:
            check = lambda call: call.sample not in samples
        else:
            check = lambda call: call.sample in samples

        calls = [call.call for call in self.calls if check(call)]
        sample_indexes = {call.sample: idx for idx, call in enumerate(calls)}
        int_alleles = set()
        for call in calls:
            if call.called:
                int_alleles.update([int(allele) for allele in call.gt_alleles])

        record = self.record
        alt_alleles = [a for i, a in enumerate(record.ALT) if i + 1 in int_alleles]
        if not alt_alleles:
            alt_alleles.append(None)
        info = record.INFO if keep_info else {}
        record = pyvcfRecord(record.CHROM, record.POS, record.ID,
                             record.REF, alt_alleles, record.QUAL,
                             record.FILTER, info, record.FORMAT,
                             sample_indexes, calls)
        snv = SNV(record, self.reader,
                  min_calls_for_pop_stats=self.min_calls_for_pop_stats)
        return snv


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
        self._allele_depths = al_counts
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
    def depth(self):
        # In freebayes the allele depth (DP) and the sum of allele
        # observations do not match
        if self.snv.reader.snpcaller == FREEBAYES:
            al_dps = self.allele_depths
            depth = sum(al_dps.values())
        else:
            depth = self.call.data.DP
        return depth

    @property
    def gt_qual(self):
        return self.call.data.GQ

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
        depth = sum(al_dps.values())
        if depth:
            major_allele_depth = max(al_dps.values())
            return major_allele_depth / depth
        else:
            return None

    @property
    def sample(self):
        return self.call.sample

    @property
    def called(self):
        return self.call.called

    @property
    def gt_type(self):
        return self.call.gt_type

    @property
    def is_het(self):
        return self.call.is_het

    @property
    def int_alleles(self):
        call = self.call
        if call.called:
            return [int(al) for al in call.gt_alleles]
        else:
            return []

    def __str__(self):
        return str(self.call)

    def __unicode__(self):
        return self.call.__unicode__()

    def __repr__(self):
        return repr(self.call)

    def copy_setting_gt_to_none(self, return_pyvcf_call=False):
        snv = self.snv
        calldata_class = snv.calldata_class
        call = self.call
        # Access to a protected member. In this case namedtuple _fields
        # is not a protected member
        # pylint: disable=W0212
        sampdat = [None if field == 'GT' else getattr(call.data, field) for field in calldata_class._fields]

        pyvcf_call = pyvcfCall(self.snv.record, self.call.sample,
                               calldata_class(*sampdat))
        if return_pyvcf_call:
            call = pyvcf_call
        else:
            call = Call(pyvcf_call, snv)
        return call
