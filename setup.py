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

'''
Created on 20/06/2012

@author: jose
'''
import sys
import os

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


def get_scripts():
    scripts = []
    for file_ in os.listdir('bin'):
        if not file_.endswith('.error'):
            scripts.append(os.path.join('bin', file_))
    return scripts


setup_args = {
              'name': 'seq_crumbs',
              'version': __version__,
              'description': 'Small utilities for sequence files manipulation',
              'author': 'Jose Blanca & Peio Ziarsolo',
              'author_email': 'jblanca@upv.es',
              'url': 'http://bioinf.comav.upv.es/seq_crumbs/',
              'packages': ['crumbs', 'crumbs.third_party', 'crumbs.utils'],
              'scripts': get_scripts(),
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
