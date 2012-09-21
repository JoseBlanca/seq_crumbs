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

import os.path
import subprocess
import tempfile

from crumbs.seqio import seqio, write_seqrecords
from crumbs.utils.bin_utils import (check_process_finishes, popen,
                                    get_binary_path)
from crumbs.alignment_result import (filter_alignments, ELONGATED, QUERY,
                                     covered_segments_from_match_parts,
                                     elongate_match_parts_till_global,
                                     TabularBlastParser)

BLAST_FIELDS = {'query': 'qseqid', 'subject': 'sseqid', 'identity': 'pident',
                'aligment_length': 'length', 'mismatches': 'mismatch',
                'gap_open': 'gapopen', 'query_start': 'qstart',
                'query_end': 'qend', 'subject_start': 'sstart',
                'subject_end': 'send', 'expect': 'evalue', 'score': 'bitscore',
                'subject_length': 'slen', 'query_length': 'qlen'}


def generate_tabblast_format(fmt):
    'Given a list with fields with our names it return one with the blast ones'
    return '6 ' + ' '.join([BLAST_FIELDS[f] for f in fmt])


def _makeblastdb_plus(seq_fpath, dbtype, outputdb=None):
    'It creates the blast db database'
    cmd = [get_binary_path('makeblastdb'), '-in', seq_fpath, '-dbtype', dbtype]
    if outputdb is not None:
        cmd.extend(['-out', outputdb])
    process = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])


def get_blast_db(seq_fpath, dbtype, directory=None):
    '''It returns an fpath to a blast db for the given sequence fpath.

    If the blast database does not exists it will create a new one.
    '''
    assert dbtype in ('nucl', 'prot')

    if directory is None:
        directory = tempfile.gettempdir()

    db_name = os.path.splitext(os.path.basename(seq_fpath))[0]
    db_path = os.path.join(directory, db_name)
    if os.path.exists(db_path):
        return db_path

    seqio([open(seq_fpath)], [open(db_path, 'w')], out_format='fasta',
          copy_if_same_format=False)

    _makeblastdb_plus(db_path, dbtype)

    return db_path


def do_blast(query_fpath, db_fpath, program, out_fpath, params=None):
    'It does a blast'
    if not params:
        params = {}
    evalue = params.get('evalue', 0.001)
    task = params.get('task', 'megablast')
    outfmt = str(params.get('outfmt', 5))
    assert task in ('blastn', 'blastn-short', 'dc-megablast', 'megablast',
                    'rmblastn')

    if program not in ('blastn', 'blastp', 'blastx', 'tblastx', 'tblastn'):
        raise ValueError('The given program is invalid: ' + str(program))
    binary = get_binary_path(program)
    cmd = [binary, '-query', query_fpath, '-db', db_fpath, '-out', out_fpath]
    cmd.extend(['-evalue', str(evalue), '-task', task])
    cmd.extend(['-outfmt', outfmt])
    process = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])


def _do_blast_2(db_fhand, queries, program, blast_format=None, params=None):
    '''It returns an alignment result with the blast.

    It is an alternative interface to the one based on fpaths.
    db_fhand should be a plain sequence file.
    queries should be a SeqRecord list.
    If an alternative blast output format is given it should be tabular, so
    blast_format is a list of fields.
    '''

    query_fhand = write_seqrecords(queries, file_format='fasta')
    query_fhand.flush()

    blastdb = get_blast_db(db_fhand.name, dbtype='nucl')

    if blast_format is None:
        blast_format = ['query', 'subject', 'query_length', 'subject_length',
                        'query_start', 'query_end', 'subject_start',
                        'subject_end', 'expect', 'identity']
    fmt = generate_tabblast_format(blast_format)
    if params is None:
        params = {}
    params['outfmt'] = fmt

    blast_fhand = tempfile.NamedTemporaryFile(suffix='.blast')
    do_blast(query_fhand.name, blastdb, program, blast_fhand.name, params)

    return TabularBlastParser(blast_fhand, blast_format)


class BlastMatcher(object):
    'It matches the given SeqRecords against the reads in the file using Blast'
    def __init__(self, reads_fhand, seqrecords, program, params=None,
                 filters=None, elongate_for_global=False):
        'It inits the class.'
        self._subjects = seqrecords
        self.program = program
        if params is None:
            params = {}
        self.params = params
        if filters is None:
            filters = []
        self.filters = filters
        self.elongate_for_global = elongate_for_global
        self._match_parts = self._look_for_blast_matches(reads_fhand,
                                                         seqrecords)

    def _look_for_blast_matches(self, seq_fhand, oligos):
        'It looks for the oligos in the given sequence files'
        # we need to keep the blast_fhands, because they're temp files and
        # otherwise they might be removed
        blasts = _do_blast_2(seq_fhand, oligos, params=self.params,
                             program=self.program)
        if self.filters is not None:
            blasts = filter_alignments(blasts, config=self.filters)

        # Which are the regions covered in each sequence?
        indexed_match_parts = {}
        one_oligo = True if len(oligos) == 1 else False
        for blast in blasts:
            oligo = blast['query']
            for match in blast['matches']:
                read = match['subject']
                if self.elongate_for_global:
                    elongate_match_parts_till_global(match['match_parts'],
                                                 query_length=oligo['length'],
                                                 subject_length=read['length'],
                                                 align_completely=QUERY)

                #match_parts = [m['match_parts'] for m in blast['matches']]
                match_parts = match['match_parts']
                if one_oligo:
                    indexed_match_parts[read['name']] = match_parts
                else:
                    try:
                        indexed_match_parts[read['name']].extend(match_parts)
                    except KeyError:
                        indexed_match_parts[read['name']] = match_parts
        return indexed_match_parts

    def get_matched_segments_for_read(self, read_name):
        'It returns the matched segments for any oligo'
        ignore_elongation_shorter = 3

        try:
            match_parts = self._match_parts[read_name]
            if read_name == 'seq5':
                pass
        except KeyError:
            # There was no match in the blast
            return None

        # Any of the match_parts has been elongated?
        elongated_match = False
        for m_p in match_parts:
            if ELONGATED in m_p and m_p[ELONGATED] > ignore_elongation_shorter:
                elongated_match = True
        segments = covered_segments_from_match_parts(match_parts,
                                                     in_query=False)
        return segments, elongated_match
