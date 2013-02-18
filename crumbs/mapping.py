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
import shutil
from subprocess import PIPE
from tempfile import NamedTemporaryFile

from crumbs.utils.bin_utils import (check_process_finishes, get_binary_path,
                                    popen, get_num_threads)
from crumbs.settings import get_setting
from crumbs.utils.file_utils import TemporaryDir


def _bwa_index_exists(index_path):
    'It checks if a a bwa index exists giving its path'
    return os.path.exists('{}.bwt'.format(index_path))


def _create_bwa_index(index_fpath):
    binary = get_binary_path('bwa')
    # how many sequences do we have?
    n_seqs = [l for l in open(index_fpath) if l[0] == '>']
    algorithm = 'bwtsw' if n_seqs > 10000 else 'is'
    cmd = [binary, 'index', '-a', algorithm, index_fpath]
    process = popen(cmd, stdout=PIPE, stderr=PIPE)
    check_process_finishes(process, binary=cmd[0])


def get_or_create_bwa_index(fpath, directory=None):
    'It creates the bwa index for the given reference'

    if directory is not None:
        index_fpath = os.path.join(directory, os.path.basename(fpath))
    else:
        index_fpath = fpath

    if not _bwa_index_exists(index_fpath):
        if os.path.exists(index_fpath):
            temp_dir = TemporaryDir()
            tmp_index_fpath = os.path.join(temp_dir.name,
                                           os.path.basename(index_fpath))
            os.symlink(fpath, tmp_index_fpath)
            _create_bwa_index(tmp_index_fpath)
            for file_ in os.listdir(temp_dir.name):
                if file_ == os.path.basename(index_fpath):
                    continue
                shutil.copy(os.path.join(temp_dir.name, file_),
                            os.path.join(os.path.dirname(index_fpath), file_))
        else:
            os.symlink(fpath, index_fpath)
            _create_bwa_index(index_fpath)

    return index_fpath


def map_with_bwa_sw(index_fpath, bam_fpath, unpaired_fpath=None,
                    paired_fpaths=None, readgroup=None, threads=None,
                    log_fpath=None, extra_params=None):
    'It maps with bwa ws algorithm'
    if paired_fpaths is None and unpaired_fpath is None:
        raise RuntimeError('At least one file to map is required')
    elif paired_fpaths is not None and unpaired_fpath is not None:
        msg = 'Bwa can not map unpaired and unpaired reads together'
        raise RuntimeError(msg)

    if readgroup is None:
        readgroup = {}

    if extra_params is None:
        extra_params = []

    binary = get_binary_path('bwa')

    cmd = [binary, 'bwasw', '-t', str(get_num_threads(threads)), index_fpath]
    cmd.extend(extra_params)

    if paired_fpaths is not None:
        cmd.extend(paired_fpaths)
    if unpaired_fpath is not None:
        cmd.append(unpaired_fpath)

    if log_fpath is None:
        stderr = NamedTemporaryFile(suffix='.stderr')
    else:
        stderr = open(log_fpath, 'w')
    # raw_input(' '.join(cmd))
    bwa = popen(cmd, stderr=stderr, stdout=PIPE)

    # add readgroup using picard
    picard_tools = get_setting("PICARD_TOOLS_DIR")
    if readgroup:
        cmd = ['java', '-jar',
           os.path.join(picard_tools, 'AddOrReplaceReadGroups.jar'),
           'INPUT=/dev/stdin', 'OUTPUT={0}'.format(bam_fpath),
           'RGID={0}'.format(readgroup['ID']),
           'RGLB={0}'.format(readgroup['LB']),
           'RGPL={0}'.format(readgroup['PL']),
           'RGSM={0}'.format(readgroup['SM']),
           'VALIDATION_STRINGENCY=LENIENT']
    else:
        cmd = [get_binary_path('samtools'), 'view', '-h', '-b', '-S', '-',
               '-o', bam_fpath]

    samtools = popen(cmd, stdin=bwa.stdout, stderr=stderr)
    bwa.stdout.close()  # Allow p1 to receive a SIGPIPE if samtools exits.
    samtools.communicate()


def _bowtie2_index_exists(index_path):
    'It checks if a a bowtie2 index exists giving its path'
    return os.path.exists('{}.1.bt2'.format(index_path))


def get_or_create_bowtie2_index(fpath, directory=None):
    "it creates the bowtie2 index"
    binary = get_binary_path('bowtie2-build')
    if directory is not None:
        index_fpath = os.path.join(directory, os.path.basename(fpath))
    else:
        index_fpath = fpath

    if not _bowtie2_index_exists(fpath):
        cmd = [binary, '-f', fpath, index_fpath]
        process = popen(cmd, stdout=PIPE, stderr=PIPE)
        check_process_finishes(process, binary=cmd[0])
    return index_fpath


def map_with_bowtie2(index_fpath, bam_fpath, paired_fpaths=None,
                     unpaired_fpaths=None, readgroup=None, threads=None,
                     log_fpath=None, preset='very-sensitive-local',
                     extra_params=None):
    '''It maps with bowtie2.

    paired_seqs is a list of tuples, in which each tuple are paired seqs
    unpaired_seqs is a list of files
    '''
    if readgroup is None:
        readgroup = {}

    if extra_params is None:
        extra_params = []

    if paired_fpaths is None and unpaired_fpaths is None:
        raise RuntimeError('At least one file to map is required')

    binary = get_binary_path('bowtie2')
    cmd = [binary, '-x', index_fpath, '--{0}'.format(preset),
           '-p', str(get_num_threads(threads))]

    cmd.extend(extra_params)
    if unpaired_fpaths:
        cmd.extend(['-U', ','.join(unpaired_fpaths)])
    if paired_fpaths:
        plus = [pairs[0] for pairs in paired_fpaths]
        minus = [pairs[1] for pairs in paired_fpaths]
        cmd.extend(['-1', ','.join(plus), '-2', ','.join(minus)])

    if 'ID' in readgroup.keys():
        for key, value in readgroup.items():
            if key not in ('ID', 'LB', 'SM', 'PL'):
                msg = 'The readgroup header tag is not valid: {}'.format(key)
                raise RuntimeError(msg)
            if key == 'ID':
                cmd.extend(['--rg-id', value])
            else:
                cmd.extend(['--rg', '{0}:{1}'.format(key, value)])

    if log_fpath is None:
        stderr = NamedTemporaryFile(suffix='.stderr')
    else:
        stderr = open(log_fpath, 'w')

#    raw_input(' '.join(cmd))
    bowtie2 = popen(cmd, stderr=stderr, stdout=PIPE)
    # print bowtie2.stdout.read()
    cmd = [get_binary_path('samtools'), 'view', '-h', '-b', '-S', '-', '-o',
           bam_fpath]

    samtools = popen(cmd, stdin=bowtie2.stdout, stderr=stderr)
    bowtie2.stdout.close()  # Allow p1 to receive a SIGPIPE if samtools exits.
    samtools.communicate()
