# Copyright 2013 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of seq_crumbs.
# vcf_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# vcf_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with vcf_crumbs. If not, see <http://www.gnu.org/licenses/>.

import os.path
from subprocess import check_call, Popen, PIPE

TEST_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                             '..', 'test', 'test_data'))

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))


def compress_with_bgzip(in_fhand, compressed_fhand):
    '''It compresses the input fhand.

    The input fhand has to be a real file with a name.
    '''
    cmd = ['bgzip', '-c', in_fhand.name]
    check_call(cmd, stdout=compressed_fhand)


def uncompress_gzip(in_fhand, uncompressed_fhand):
    '''It compresses the input fhand.

    The input fhand has to be a real file with a name.
    '''
    cmd = ['gzip', '-dc', in_fhand.name]
    check_call(cmd, stdout=uncompressed_fhand)


def index_vcf_with_tabix(in_fpath):
    cmd = ['tabix', '-p', 'vcf', in_fpath]
    tabix = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = tabix.communicate()
    if tabix.returncode:
        msg = 'Some problem indexing with tabix.\n'
        msg += 'stdout: ' + stdout
        msg += 'stderr: ' + stderr
        raise RuntimeError(msg)
