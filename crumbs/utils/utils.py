'''
Created on 21/06/2012

@author: jose
'''

import sys

from crumbs.utils import cgitb
from crumbs.utils.exceptions import UnknownFormatError, FileNotFoundError


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
