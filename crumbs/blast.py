'''
Created on 25/04/2012

@author: jose
'''

import os.path
import subprocess


from crumbs.seqio import seqio
from crumbs.utils import check_process_finishes, popen, get_binary_path


def _makeblastdb_plus(seq_fpath, dbtype, outputdb=None):
    'It creates the blast db database'
    cmd = [get_binary_path('makeblastdb'), '-in', seq_fpath, '-dbtype', dbtype]
    if outputdb is not None:
        cmd.extend(['-out', outputdb])
    process = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])


def get_blast_db(seq_fpath, directory, dbtype):
    '''It returns an fpath to a blast db for the given sequence fpath.

    If the blast database does not exists it will create a new one.
    '''
    assert dbtype in ('nucl', 'prot')

    db_name = os.path.splitext(os.path.basename(seq_fpath))[0]
    db_path = os.path.join(directory, db_name)
    if os.path.exists(db_path):
        return db_path

    seqio([open(seq_fpath)], [open(db_path, 'w')], out_format='fasta',
          copy_if_same_format=False)

    _makeblastdb_plus(db_path, dbtype)

    return db_path


def do_blast(seq_fpath, db_fpath, program, out_fpath, params=None):
    'It does a blast'
    if not params:
        params = {}
    evalue = params.get('evalue', 0.001)
    task = params.get('task', 'megablast')
    assert task in ('blastn', 'blastn-short', 'dc-megablast', 'megablast',
                    'rmblastn')

    if program not in ('blastn', 'blastp', 'blastx', 'tblastx', 'tblastn'):
        raise ValueError('The given program is invalid: ' + program)
    cmd = [program, '-query', seq_fpath, '-db', db_fpath, '-out', out_fpath]
    cmd.extend(['-evalue', str(evalue), '-task', task])
    cmd.extend(['-outfmt', '5'])
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])
