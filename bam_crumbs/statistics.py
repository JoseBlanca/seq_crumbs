from __future__ import division


from numpy import histogram, zeros, median, sum

from crumbs.statistics import draw_histogram, IntCounter, LABELS

# pylint: disable=C0111


DEFAULT_N_BINS = 20


def count_reads(ref_name, bam, start=None, end=None):
    'It returns the count of aligned reads in the region'
    return bam.count(reference=ref_name, start=start, end=end)


class ArrayWrapper(object):
    'A thin wrapper aroung numpy array to have the same interface as IntCounter'
    def __init__(self, array):
        self.array = array
        self.labels = LABELS.copy()

    @property
    def min(self):
        return self.array.min()

    @property
    def max(self):
        return self.array.max()

    @property
    def average(self):
        return self.array.mean()

    @property
    def median(self):
        return median(self.array)

    @property
    def variance(self):
        return self.array.var()

    @property
    def count(self):
        return len(self.array)

    @property
    def sum(self):
        return sum(self.array)

    def calculate_distribution(self, bins=DEFAULT_N_BINS, min_=None,
                               max_=None):
        if min_ is None:
            min_ = self.min
        if max_ is None:
            max_ = self.max

        counts, bins = histogram(self.array, bins=bins, range=(min_, max_))
        return {'bin_limits': bins, 'counts': counts}

    def update_labels(self, labels):
        'It prepares the labels for output files'
        self.labels.update(labels)

    def __str__(self):
        'It writes some basic stats of the values'
        if self.count != 0:
            labels = self.labels
            # now we write some basic stats
            format_num = lambda x: '{:,d}'.format(x) if isinstance(x, int) else '%.2f' % x
            text = '{}: {}\n'.format(labels['minimum'], format_num(self.min))
            text += '{}: {}\n'.format(labels['maximum'], format_num(self.max))
            text += '{}: {}\n'.format(labels['average'],
                                      format_num(self.average))

            if labels['variance'] is not None:
                text += '{}: {}\n'.format(labels['variance'],
                                          format_num(self.variance))
            if labels['sum'] is not None:
                text += '{}: {}\n'.format(labels['sum'],
                                          format_num(self.sum))
            if labels['items'] is not None:
                text += '{}: {}\n'.format(labels['items'], self.count)
            text += '\n'

            distrib = self.calculate_distribution()
            text += draw_histogram(distrib['bin_limits'], distrib['counts'])
            return text
        return ''


class ReferenceStats(object):
    def __init__(self, bams):
        self._bams = bams
        self._rpkms = None
        self._tot_reads = 0
        self._lengths = None
        self._count_reads()

    def _count_reads(self):
        tot_reads = self._tot_reads
        nreferences = self._bams[0].nreferences
        rpks = zeros(nreferences)
        lengths = IntCounter()
        for bam in self._bams:
            if bam.nreferences != nreferences:
                msg = 'BAM files should have the same references'
                raise ValueError(msg)
            # For the references we use the first BAM to make sure that the
            # references are the same in all bams
            for index, ref in enumerate(self._bams[0].header['SQ']):
                length = ref['LN']
                count = count_reads(ref['SN'], bam)
                rpk = count / length
                tot_reads += count
                rpks[index] = rpk
                lengths[length] += 1
            self._lengths = lengths

            # from rpk to rpkms
            million_reads = tot_reads / 1e6
            rpks /= million_reads
            self._rpkms = ArrayWrapper(rpks)

    @property
    def lengths(self):
        return self._lengths

    @property
    def rpkms(self):
        return self._rpkms


class MapqCounter(IntCounter):
    def __init__(self, bams):
        self._bams = bams
        self._count_mapqs()

    def _count_mapqs(self):
        for bam in self._bams:
            for read in bam.fetch():
                self[read.mapq] += 1


class CoverageCounter(IntCounter):
    def __init__(self, bams):
        self._bams = bams
        self._count_cov()

    def _count_cov(self):
        for bam in self._bams:
            for column in bam.pileup():
                self[len(column.pileups)] += 1
