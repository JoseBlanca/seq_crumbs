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

from crumbs.seqio import seqio
from crumbs.utils.bin_utils import (check_process_finishes, popen,
                                    get_binary_path)

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
