
import os.path
from subprocess import check_call
import shutil
from tempfile import NamedTemporaryFile

import pysam

from bam_crumbs.utils.flag import create_flag
from bam_crumbs.settings import get_setting

# pylint: disable=C0111


def filter_bam(in_fpath, out_fpath, min_mapq=0, required_flag_tags=None,
               filtering_flag_tags=None, regions=None):
    cmd = ['-bh']

    # The following line:
    cmd.append('-o' + out_fpath)
    # should be
    # cmd.extend(['-o', out_fpath])
    # but it is a workaround, take a look at:
    # https://groups.google.com/forum/#!msg/pysam-user-group/ooHgIiNVe4c/CcY06d45rzQJ

    if min_mapq:
        cmd.extend(['-q', str(min_mapq)])

    if required_flag_tags:
        flag = create_flag(required_flag_tags)
        cmd.extend(['-f', str(flag)])

    if filtering_flag_tags:
        flag = create_flag(filtering_flag_tags)
        cmd.extend(['-F', str(flag)])

    cmd.extend([in_fpath])

    if regions:
        regions = ['{0}:{1}-{2}'.format(*s) for s in regions.segments]
        cmd.extend(regions)

    pysam.view(*cmd)


def sort_bam(in_bam_fpath, out_bam_fpath=None):

    if out_bam_fpath is None:
        out_bam_fpath = in_bam_fpath

    if out_bam_fpath == in_bam_fpath:
        sorted_fhand = NamedTemporaryFile(suffix='.sorted.bam', delete=False)
        temp_out_fpath = sorted_fhand.name
    else:
        temp_out_fpath = out_bam_fpath

    picard_tools = get_setting("PICARD_TOOLS_DIR")
    cmd = ['java', '-jar', os.path.join(picard_tools, 'SortSam.jar'),
           'INPUT={0}'.format(in_bam_fpath),
           'OUTPUT={0}'.format(temp_out_fpath),
           'SORT_ORDER=coordinate', 'VALIDATION_STRINGENCY=LENIENT']
    stderr = NamedTemporaryFile(suffix='picard.stderr')
    check_call(cmd, stderr=stderr)

    if temp_out_fpath != out_bam_fpath:
        shutil.move(temp_out_fpath, out_bam_fpath)


def index_bam(bam_fpath):
    'It indexes a bam file'
    pysam.index(bam_fpath)
