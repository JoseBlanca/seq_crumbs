# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of seq_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with seq_crumbs. If not, see <http://www.gnu.org/licenses/>.

import subprocess
from tempfile import NamedTemporaryFile
from collections import Counter

from crumbs.utils.bin_utils import (get_binary_path, popen,
                                    check_process_finishes)
from crumbs.utils.tags import (PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS)
from crumbs.seqio import write_seqrecords, read_seqrecords
from Bio.SeqFeature import SeqFeature, FeatureLocation


def _run_estscan(seqrecords, pep_out_fpath, dna_out_fpath, matrix_fpath):
    'It runs estscan in the input seqs'
    seq_fhand = write_seqrecords(seqrecords, file_format='fasta')
    seq_fhand.flush()
    binary = get_binary_path('estscan')

    cmd = [binary, '-t', pep_out_fpath, '-o', dna_out_fpath, '-M',
           matrix_fpath, seq_fhand.name]
    process = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])
    seq_fhand.close()


def _read_estcan_result(fhand, result, file_type):
    'It reads a dna or pep ESTscan result file'
    for seq in read_seqrecords([fhand], file_format='fasta'):
        items = [i.strip() for i in seq.description.split(';')]
        strand = -1 if 'minus strand' in items else 1
        start, end = items[0].split(' ', 3)[1:3]
        seqid = seq.id
        try:
            seq_orfs = result[seqid]
        except KeyError:
            seq_orfs = {}
            result[seqid] = seq_orfs
        orf_key = (int(start), int(end), strand)
        if orf_key in seq_orfs:
            orf = seq_orfs[orf_key]
        else:
            orf = {}
            seq_orfs[orf_key] = orf
        orf[file_type] = seq.seq


def _read_estcan_results(pep_fhand, dna_fhand):
    'It reads the fhand result files'
    result = {}
    _read_estcan_result(pep_fhand, result, 'pep')
    _read_estcan_result(dna_fhand, result, 'dna')
    return result


class EstscanOrfAnnotator(object):
    'It annotates the given seqrecords'
    def __init__(self, usage_matrix):
        'Initiator'
        self._usage_matrix = usage_matrix
        self._stats = Counter()

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It runs the actual annotations'
        stats = self._stats
        stats[PROCESSED_PACKETS] += 1
        pep_fhand = NamedTemporaryFile()
        dna_fhand = NamedTemporaryFile()
        _run_estscan(seqrecords, pep_fhand.name, dna_fhand.name,
                     self._usage_matrix)
        # now we read the result files
        estscan_result = _read_estcan_results(open(pep_fhand.name),
                                              open(dna_fhand.name))

        for seq in seqrecords:
            stats[PROCESSED_SEQS] += 1
            seq_name = seq.id
            orfs = estscan_result.get(seq_name, {})
            feats = []
            for (start, end, strand), seqs in orfs.viewitems():
                start -= 1
                # end is fine  -- end[
                feat = SeqFeature(location=FeatureLocation(start, end, strand),
                                  type='ORF', qualifiers=seqs)
                feats.append(feat)
            if feats:
                seq.features.extend(feats)
            stats[YIELDED_SEQS] += 1

        dna_fhand.close()
        pep_fhand.close()
        return seqrecords
