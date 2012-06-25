'''
Created on 21/06/2012

@author: jose
'''

import sys
import os

from crumbs.utils import cgitb
from crumbs.utils.exceptions import (UnknownFormatError, FileNotFoundError,
                                     WrongFormatError, TooManyFiles,
                                     MalformedFile)

STDIN = 'stdin'
STDOUT = 'stdout'
INFILES = 'infiles'
OUTFILE = 'output'


def main(funct):
    'The main function of a script'
    if len(sys.argv) == 1:
        sys.argv = sys.argv + ['-h']

    argv = sys.argv
    if '--error_log' in argv:
        error_fpath_index = argv.index('--error_log') + 1
        error_fpath = argv[error_fpath_index]
    else:
        binary = sys.argv[0]
        error_fpath = binary + '.error'

    stderr = sys.stderr
    try:
        return(funct())
    except FileNotFoundError, error:
        stderr.write(str(error) + '\n')
        return 2
    except UnknownFormatError, error:
        stderr.write(str(error) + '\n')
        return 3
    except WrongFormatError, error:
        stderr.write(str(error) + '\n')
        return 4
    except TooManyFiles, error:
        stderr.write(str(error) + '\n')
        return 5
    except MalformedFile, error:
        stderr.write(str(error) + '\n')
        return 6
    except Exception as error:
        msg = 'An unexpected error happened.\n'
        msg += 'The seq crumbs developers would appreciate your feedback\n'
        msg += 'Please send them the error log: '
        msg += error_fpath + '\n\n'
        msg += str(error)
        stderr.write(msg)
        hook = cgitb.Hook(display=0, format='text', logfpath=error_fpath)
        hook.handle()
        fhand = open(error_fpath, 'a')
        fhand.write('\nThe command was:\n' + ' '.join(sys.argv) + '\n')
        fhand.close()
        raise


def get_inputs_from_args(parsed_args):
    'It returns the input fhand'
    in_fpaths = getattr(parsed_args, INFILES)
    if in_fpaths == STDIN:
        in_fhands = [sys.stdin]
    else:
        in_fhands = []
        for in_fpath in in_fpaths:
            if os.path.exists(in_fpath):
                in_fhand = open(in_fpath, 'rt')
            else:
                raise FileNotFoundError('A file was not found: ' + in_fpath)
            in_fhands.append(in_fhand)
    return in_fhands


def get_output_from_args(parsed_args):
    'It returns the out_fhand'
    out_fpath = getattr(parsed_args, OUTFILE)
    if out_fpath == STDOUT:
        out_fhand = sys.stdout
    else:
        out_fhand = open(out_fpath, 'w')
    return out_fhand
