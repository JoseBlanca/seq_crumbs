
import os.path
from subprocess import check_call
import shutil
from tempfile import NamedTemporaryFile, gettempdir

import pysam

from bam_crumbs.utils.flag import create_flag
from bam_crumbs.settings import get_setting
from crumbs.utils.bin_utils import get_num_threads

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


def _create_sam_reference_index(fpath):
    'It creates a sam index for a reference sequence file'
    index_fpath = fpath + '.fai'
    if os.path.exists(index_fpath):
        return
    pysam.faidx(fpath)


def _create_picard_dict(fpath):
    'It creates a picard dict if if it does not exist'
    dict_path = os.path.splitext(fpath)[0] + '.dict'
    if os.path.exists(dict_path):
        return
    picard_dir = get_setting("PICARD_TOOLS_DIR")
    picard_jar = os.path.join(picard_dir, 'CreateSequenceDictionary.jar')
    cmd = ['java', '-jar', picard_jar,
           'R=%s' % fpath,
           'O=%s' % dict_path]
    stderr = NamedTemporaryFile(suffix='picard.stderr')
    check_call(cmd, stderr=stderr)


def realign_bam(in_bam_fpath, reference_fpath, out_bam_fpath=None,
                threads=False):

    if out_bam_fpath is None:
        out_bam_fpath = in_bam_fpath

    if out_bam_fpath == in_bam_fpath:
        realigned_fhand = NamedTemporaryFile(suffix='.realigned.bam', delete=False)
        temp_out_fpath = realigned_fhand.name
    else:
        temp_out_fpath = out_bam_fpath

    _realign_bam(in_bam_fpath, reference_fpath, temp_out_fpath, threads=False)
    sort_bam(temp_out_fpath)

    if temp_out_fpath != out_bam_fpath:
        shutil.move(temp_out_fpath, out_bam_fpath)


def _realign_bam(bam_fpath, reference_fpath, out_bam_fpath, threads=False):
    'It realigns the bam using GATK Local realignment around indels'
    # reference sam index
    _create_sam_reference_index(reference_fpath)

    # reference picard dict
    _create_picard_dict(reference_fpath)

    # bam index
    index_bam(bam_fpath)

    # the intervals to realign
    gatk_dir = get_setting("GATK_DIR")
    gatk_jar = os.path.join(gatk_dir, 'GenomeAnalysisTK.jar')
    intervals_fhand = NamedTemporaryFile(suffix='.intervals')
    stderr = NamedTemporaryFile(suffix='picard.stderr')
    stdout = NamedTemporaryFile(suffix='picard.stdout')
    cmd = ['java', '-jar', gatk_jar, '-T', 'RealignerTargetCreator',
           '-I', bam_fpath, '-R', reference_fpath, '-o', intervals_fhand.name]
    check_call(cmd, stderr=stderr, stdout=stdout)

    # the realignment itself
    cmd = ['java', '-jar', gatk_jar, '-I', bam_fpath, '-R', reference_fpath,
           '-T', 'IndelRealigner', '-targetIntervals', intervals_fhand.name,
           '-o', out_bam_fpath]

    if threads and threads > 1:
        cmd.extend(['-nt', str(get_num_threads(threads))])
    check_call(cmd, stderr=stderr, stdout=stdout)
    intervals_fhand.close()


def calmd_bam(in_bam_fpath, reference_fpath, out_bam_fpath=None):

    if out_bam_fpath is None:
        out_bam_fpath = in_bam_fpath

    if out_bam_fpath == in_bam_fpath:
        realigned_fhand = NamedTemporaryFile(suffix='.realigned.bam', delete=False)
        temp_out_fpath = realigned_fhand.name
    else:
        temp_out_fpath = out_bam_fpath

    _calmd_bam(in_bam_fpath, reference_fpath, temp_out_fpath)


def _calmd_bam(bam_fpath, reference_fpath, out_bam_fpath):
    out_fhand = open(out_bam_fpath, 'wb')
    for line in pysam.calmd(*["-bAr", bam_fpath, reference_fpath]):
        out_fhand.write(line)
    # out_fhand.write(pysam.calmd(*["-bAr", bam_fpath, reference_fpath]))
    out_fhand.flush()
    out_fhand.close()
