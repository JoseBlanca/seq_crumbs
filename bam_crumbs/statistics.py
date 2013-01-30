from __future__ import division

from array import array

from numpy import array as np_array, histogram

from crumbs.statistics import draw_histogram, IntCounter

# pylint: disable=C0111


def count_reads(ref_name, bam, start=None, end=None):
    'It returns the count of aligned reads in the region'
    return bam.count(reference=ref_name, start=start, end=end)


class RpkmCounter(object):
    def __init__(self, bam):
        self._bam = bam
        self._rpkms = None
        self._count_reads()

    def _count_reads(self):
        tot_reads = 0
        rpks = array('f')
        for ref in self._bam.header['SQ']:
            length = ref['LN']
            count = count_reads(ref['SN'], self._bam)
            rpk = count / length
            tot_reads += count
            rpks.append(rpk)

        # from rpk to rpkms
        million_reads = tot_reads / 1e6
        for index in range(len(rpks)):
            rpks[index] /= million_reads
        rpkms = np_array(rpks)
        self._rpkms = rpkms

    def ascii_histogram(self):
        counts, bins = histogram(self._rpkms)
        return draw_histogram(bins, counts)


class MapqCounter(IntCounter):
    def __init__(self, bam):
        self._bam = bam
        self._count_mapqs()

    def _count_mapqs(self):
        for read in self._bam.fetch():
            self[read.mapq] += 1


class CoverageCounter(IntCounter):
    def __init__(self, bam):
        self._bam = bam
        self._count_cov()

    def _count_cov(self):
        for column in self._bam.pileup():
            self[len(column.pileups)] += 1
