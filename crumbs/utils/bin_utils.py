import sys
import os.path
import platform

from crumbs.settings import get_setting
from crumbs.third_party import cgitb
from crumbs.exceptions import (UnknownFormatError, FileNotFoundError,
                               WrongFormatError, TooManyFiles,
                               MalformedFile, SampleSizeError,
                               ExternalBinaryError, MissingBinaryError,
                               IncompatibleFormatError,
                               UndecidedFastqVersionError, MaxNumReadsInMem,
                               PairDirectionError, InterleaveError,
                               OptionalRequirementError)

from crumbs.utils.tags import (OUTFILE, GUESS_FORMAT, BGZF, GZIP, BZIP2,
                               ERROR_ENVIRON_VARIABLE)

from crumbs import __version__ as version

BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..',
                                       'bin'))


def main(funct):
    'The main function of a script'
    argv = sys.argv
    if '--error_log' in argv:
        error_fpath_index = argv.index('--error_log') + 1
        error_fpath = argv[error_fpath_index]
    else:
        binary = os.path.split(sys.argv[0])[-1]
        error_fpath = binary + '.error'

    stderr = sys.stderr
    try:
        # This code is required to test the error handling
        fail = os.environ.get(ERROR_ENVIRON_VARIABLE, None)
        if fail:
            raise RuntimeError('Generating a test error')
        return(funct())
    except FileNotFoundError, error:
        stderr.write(str(error) + '\n')
        return 3
    except UnknownFormatError, error:
        stderr.write(str(error) + '\n')
        return 4
    except WrongFormatError, error:
        stderr.write(str(error) + '\n')
        return 5
    except TooManyFiles, error:
        stderr.write(str(error) + '\n')
        return 6
    except MalformedFile, error:
        stderr.write(str(error) + '\n')
        return 7
    except SampleSizeError, error:
        stderr.write(str(error) + '\n')
        return 8
    except ExternalBinaryError, error:
        stderr.write(str(error) + '\n')
        return 9
    except MissingBinaryError, error:
        stderr.write(str(error) + '\n')
        return 10
    except IncompatibleFormatError, error:
        stderr.write(str(error) + '\n')
        return 11
    except UndecidedFastqVersionError, error:
        stderr.write(str(error) + '\n')
        return 12
    except MaxNumReadsInMem, error:
        stderr.write(str(error) + '\n')
        return 13
    except PairDirectionError, error:
        stderr.write(str(error) + '\n')
        return 14
    except InterleaveError, error:
        stderr.write(str(error) + '\n')
        return 15
    except KeyboardInterrupt, error:
        stderr.write('Program stopped by user request\n')
        return 16
    except OptionalRequirementError, error:
        stderr.write(str(error) + '\n')
        return 17
    except Exception as error:
        msg = 'An unexpected error happened.\n'
        msg += 'The seq_crumbs developers would appreciate your feedback.\n'
        try:
            fail = os.environ.get(ERROR_ENVIRON_VARIABLE, None)
            if fail:
                # error handling debugging
                fhand = sys.stderr
                error_fpath = None
            else:
                fhand = open(error_fpath, 'a')
        except IOError:
            # the fpath for the error is not writable
            fhand = sys.stderr
            error_fpath = None

        if error_fpath:
            msg += 'Please send them the error log'
            msg += ': ' + error_fpath + '\n\n'
        msg += str(error)
        stderr.write(msg)

        if error_fpath:
            hook = cgitb.Hook(display=0, format='text', logfpath=error_fpath)
            hook.handle()

        fhand.write('\nThe command was:\n' + ' '.join(sys.argv) + '\n')
        fhand.close()
        raise


def build_version_msg():
    'It creates a message with the version.'
    bin_name = os.path.split(sys.argv[0])[-1]
    version_msg = bin_name + ' from seq_crumbs version: ' + version
    return version_msg


# This function has been modified from
# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    'It emulates the which Unix utility'
    def is_exe(fpath):
        'The file is an executable'
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath = os.path.dirname(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def get_binary_path(binary_name):
    '''It return the path to the proper binary. It looks on platform and
    architecture to decide it.

    Fails if there is not binary for that architecture
    '''
    if get_setting('USE_EXTERNAL_BIN_PREFIX'):
        ext_binary_name = get_setting('EXTERNAL_BIN_PREFIX') + binary_name
        if os.path.exists(ext_binary_name):
            return ext_binary_name

    if not get_setting('ADD_PATH_TO_EXT_BIN'):
        # I have to check if the binary is on my current directory.
        # If it is there use it, else assumes that it is on the path
        if os.path.exists(os.path.join(os.getcwd(), ext_binary_name)):
            return os.path.join(os.getcwd(), ext_binary_name)
        # return binary_name

    system = platform.system().lower()
    if system == 'windows':
        binary_name += '.exe'
    arch = platform.architecture()[0]

    join = os.path.join

    third_party_path = get_setting('THIRD_PARTY_DIR')

    binary_path = os.path.abspath(join(third_party_path, system, arch,
                                       binary_name))

    if os.path.exists(binary_path):
        return binary_path
    elif arch == '64bit':
        arch = '32bit'
        binary_path = os.path.abspath(join(third_party_path, system, arch,
                                           binary_name))
        if os.path.exists(binary_path):
            return binary_path

    # At this point there is not available binary for the working platform
    # Is the binary really in the path?
    if which(binary_name):
        return binary_name

    msg = '{} not found in the path. Please install it to use ngs_crumbs'
    raise MissingBinaryError(msg.format(binary_name))


def get_num_threads(threads):
    """It returns num of threads to use in parallel.

    You can pass to the funaction the  memory you want to use each thread.
    It calculates the number of treads
    In megabytes
    """
    phisical_threads = os.sysconf('SC_NPROCESSORS_ONLN')
    if not threads:
        return 1
    elif isinstance(threads, bool):
        return phisical_threads
    else:
        return threads
