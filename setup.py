'''
Created on 20/06/2012

@author: jose
'''
import sys

try:
    from setuptools import setup
    _SETUPTOOLS = True
except ImportError:
    from distutils.core import setup
    _SETUPTOOLS = False

# The next three lines are modified from Biopython
__version__ = "Undefined"
for line in open('crumbs/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())
    break


def check_dependencies():
    'If a dependency is not met it stops the installation'
    if not _SETUPTOOLS:
        return  # No need to check anything because it will get installed
    dependencies = []
    try:
        import Bio
    except ImportError:
        dependencies.append('biopython')
    try:
        import argparse
    except ImportError:
        dependencies.append('argparse')

    if not dependencies:
        return

    if len(dependencies) > 1:
        msg = 'These python packages are required dependencies of seq_crumbs, '
        msg += 'please install them: '
    else:
        msg = 'seq_crumbs depends on this package, please install it: '
    msg += ' '.join(dependencies)
    print msg
    sys.exit(-1)


SCRIPTS = ['bin/sff_extract', 'bin/seq_head', 'bin/guess_seq_format']


setup_args = {
              'name': 'seq_crumbs',
              'version': __version__,
              'description': 'Small utilities for sequence files manipulation',
              'author': 'Jose Blanca & Peio Ziarsolo',
              'author_email': 'jblanca@upv.es',
              'url': 'http://bioinf.comav.upv.es/seq_crumbs/',
              'packages': ['crumbs'],
              'scripts': SCRIPTS,
              'license': 'AGPL'
              }

if _SETUPTOOLS:
    requires = ['biopython']
    try:
        import argparse
    except ImportError:
        requires.append(argparse)
    setup_args['install_requires'] = requires

setup(**setup_args)
