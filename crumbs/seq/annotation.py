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
import os.path
from tempfile import NamedTemporaryFile
from random import randint

from crumbs.utils.optional_modules import SeqFeature, FeatureLocation

from crumbs.utils.bin_utils import (get_binary_path, popen,
                                    check_process_finishes)
from crumbs.utils.tags import FIVE_PRIME, THREE_PRIME
from crumbs.seq.seq import get_description, get_name, get_str_seq
from crumbs.seq.seqio import write_seqs, read_seqs
from crumbs.blast import Blaster
from crumbs.settings import get_setting

# pylint: disable=R0903


def _run_estscan(seqs, pep_out_fpath, dna_out_fpath, matrix_fpath):
    'It runs estscan in the input seqs'
    seq_fhand = write_seqs(seqs, file_format='fasta')
    seq_fhand.flush()
    binary = get_binary_path('estscan')

    cmd = [binary, '-t', pep_out_fpath, '-o', dna_out_fpath, '-M',
           matrix_fpath, seq_fhand.name]
    process = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])
    seq_fhand.close()


def _read_estcan_result(fhand, result, file_type):
    'It reads a dna or pep ESTscan result file'
    for seq in read_seqs([fhand]):
        items = [i.strip() for i in get_description(seq).split(';')]
        strand = -1 if 'minus strand' in items else 1
        start, end = items[0].split(' ', 3)[1:3]
        # estscan changes the name, we have to fix it
        seqid = get_name(seq).strip(';')
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
        orf[file_type] = get_str_seq(seq)


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

    def __call__(self, seqs):
        'It runs the actual annotations'
        if not seqs:
            return seqs
        pep_fhand = NamedTemporaryFile()
        dna_fhand = NamedTemporaryFile()
        _run_estscan(seqs, pep_fhand.name, dna_fhand.name,
                     self._usage_matrix)
        # now we read the result files
        estscan_result = _read_estcan_results(open(pep_fhand.name),
                                              open(dna_fhand.name))
        for seq in seqs:
            seq_name = get_name(seq)
            orfs = estscan_result.get(seq_name, {})
            feats = []
            for (start, end, strand), str_seqs in orfs.viewitems():
                start -= 1
                # end is fine  -- end[
                feat = SeqFeature(location=FeatureLocation(start, end, strand),
                                  type='ORF', qualifiers=str_seqs)
                feats.append(feat)
            if feats:
                seq.object.features.extend(feats)

        dna_fhand.close()
        pep_fhand.close()
        return seqs


def _detect_polya_tail(seq, location, min_len, max_cont_mismatches):
    '''It detects 3' poylA or 5' polyT tails.

    This function is a re-implementation of the EMBOSS's trimest code.
    It will return the position of a poly-A in 3' or a poly-T in 5'.
    It returns the start and end of the tail. The nucleotide in the end
    position won't be included in the poly-A.
    '''
    if location == FIVE_PRIME:
        tail_nucl = 'T'
        inc = 1
        start = 0
        end = len(seq)
    elif location == THREE_PRIME:
        tail_nucl = 'A'
        inc = -1
        start = -1
        end = -len(seq) - 1
    else:
        msg = 'location should be five or three prime'
        raise ValueError(msg)

    mismatch_count = 0
    poly_count = 0
    tail_len = 0
    result = 0

    for index in range(start, end, inc):
        nucl = seq[index].upper()
        if nucl == tail_nucl:
            poly_count += 1
            mismatch_count = 0
        elif nucl == 'N':
            pass
        else:
            poly_count = 0
            mismatch_count += 1
        if poly_count >= min_len:
            result = tail_len + 1
        if mismatch_count > max_cont_mismatches:
            break
        tail_len += 1
    if result and location == FIVE_PRIME:
        start = 0
        end = result
        result = start, end
    elif result and location == THREE_PRIME:
        end = len(seq)
        start = end - result
        result = start, end
    else:
        result = None
    return result


def _annotate_polya(seq, min_len, max_cont_mismatches):
    'It annotates the polyA with the EMBOSS trimest method'
    str_seq = get_str_seq(seq)
    polya = _detect_polya_tail(str_seq, THREE_PRIME, min_len,
                               max_cont_mismatches)
    polyt = _detect_polya_tail(str_seq, FIVE_PRIME, min_len,
                               max_cont_mismatches)
    a_len = polya[1] - polya[0] if polya else 0
    t_len = polyt[1] - polyt[0] if polyt else 0
    chosen_tail = None
    if a_len > t_len:
        chosen_tail = 'A'
    elif t_len > a_len:
        chosen_tail = 'T'
    elif a_len and a_len == t_len:
        if randint(0, 1):
            chosen_tail = 'A'
        else:
            chosen_tail = 'T'
    if chosen_tail:
        strand = 1 if chosen_tail == 'A' else -1
        start, end = polya if chosen_tail == 'A' else polyt
        feat = SeqFeature(location=FeatureLocation(start, end, strand),
                          type='polyA_sequence')
        # We're assuming that the seq has a SeqRecord in it
        seq.object.features.append(feat)


class PolyaAnnotator(object):
    'It annotates the given seqrecords with poly-A or poly-T regions'
    def __init__(self, min_len=get_setting('POLYA_ANNOTATOR_MIN_LEN'),
                max_cont_mismatches=get_setting('POLYA_ANNOTATOR_MISMATCHES')):
        '''It inits the class.

        min_len - minimum number of consecutive As (or Ts) to extend the tail
        max_cont_mismatches - maximum number of consecutive no A (or Ts) to
                              break a tail.
        '''
        self._min_len = min_len
        self._max_cont_mismatches = max_cont_mismatches

    def __call__(self, seqrecords):
        'It runs the actual annotations'
        max_cont_mismatches = self._max_cont_mismatches
        min_len = self._min_len

        for seq in seqrecords:
            _annotate_polya(seq, min_len, max_cont_mismatches)
        return seqrecords


class BlastAnnotator(object):
    'It annotates using blast'
    def __init__(self, blastdb, program, dbtype=None, filters=None,
                 params=None, remote=False):
        'Initializes the class'
        self.blastdb = blastdb
        self._program = program
        self._filters = [] if filters is None else filters
        self._params = params
        self._dbtype = dbtype
        self._remote = remote

    def __call__(self, seqrecords):
        'It does the work'
        if not seqrecords:
            return seqrecords
        matcher = Blaster(seqrecords, self.blastdb, self._program,
                               self._dbtype, filters=self._filters,
                               params=self._params, remote=self._remote)
        blasts = matcher.blasts
        blastdb = os.path.basename(self.blastdb)
        for seqrecord in seqrecords:
            align_result = blasts.get(get_name(seqrecord), None)
            if not align_result:
                continue
            match_counter = 0
            for match in align_result['matches']:
                subject = match['subject']['name']
                match_counter += 1
                for match_part in  match['match_parts']:
                    if match_part['subject_end'] < match_part['subject_start']:
                        strand = -1
                        subject_start = match_part['subject_end']
                        subject_end = match_part['subject_start']
                    else:
                        strand = 1
                        subject_start = match_part['subject_start']
                        subject_end = match_part['subject_end']

                    query_start = match_part['query_start']
                    query_end = match_part['query_end']
                    qualifiers = {}
                    qualifiers['Target'] = {'start': subject_start,
                                            'end': subject_end,
                                            'name': subject}
                    qualifiers['score'] = match_part['scores']['expect']
                    qualifiers['identity'] = match_part['scores']['identity']
                    qualifiers['blastdb'] = blastdb
                    location = FeatureLocation(query_start, query_end, strand)
                    feature = SeqFeature(location=location, type='match_part',
                                         qualifiers=qualifiers,
                                       id='match{0:03d}'.format(match_counter))
                    seqrecord.object.features.append(feature)
        return seqrecords
