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

import sys
from sys import version_info
import os
import fnmatch
import glob
import platform
import subprocess

import distutils.command.install_data

try:
    from setuptools import setup
    from setuptools.command import install
    _SETUPTOOLS = True

except ImportError:
    from distutils.core import setup
    from distutils.command import install
    _SETUPTOOLS = False

#print "using_setuptools", _SETUPTOOLS

# The next three lines are modified from Biopython
__version__ = "Undefined"
for line in open('crumbs/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())
        break


def check_dependencies():
    'If a dependency is not met it stops the installation'
    if _SETUPTOOLS:
        return  # No need to check anything because it will get installed
    msg = None
    try:
        import Bio
    except ImportError:
        msg = 'You have to install Bioypython >= 1.60\n'

    if not msg:
        try:
            from Bio.bgzf import BgzfWriter
        except ImportError:
            msg = 'You have an old version of Biopython installed, '
            msg += 'please update to >= 1.60\n'

    if not msg:
        return
    sys.stderr.write(msg)
    sys.exit(-1)

check_dependencies()


def opj(*args):
    path = os.path.join(*args)
    return os.path.normpath(path)


def find_data_file(srcdir, *wildcards, **kw):
    # get a list of all files under the srcdir matching wildcards,
    # returned in a format to be used for install_data
    def walk_helper(arg, dirname, files):
        if '.svn' in dirname or '.git' in dirname:
            return
        names = []
        lst, wildcards = arg
        for wildcard in wildcards:
            wc_name = opj(dirname, wildcard)
            for fpath in files:
                filename = opj(dirname, fpath)
                if (fnmatch.fnmatch(filename, wc_name) and
                    not os.path.isdir(filename)):
                    names.append(filename)
        if names:
            lst.append((dirname, names))

    file_list = []
    recursive = kw.get('recursive', True)
    if recursive:
        os.path.walk(srcdir, walk_helper, (file_list, wildcards))
    else:
        walk_helper((file_list, wildcards),
                    srcdir,
                    [os.path.basename(f) for f in glob.glob(opj(srcdir, '*'))])
    return file_list


def get_platform_bin_dir():
    '''It returns the platform specific bindir. It returns the relative path
    from the source code root dir'''
    system = platform.system().lower()
    arch = platform.architecture()[0]

    return os.path.join('crumbs', 'third_party', 'bin', system, arch)


platform_bin_dir = get_platform_bin_dir()
external_executables = find_data_file(platform_bin_dir, '*')


def get_scripts():
    scripts = []
    for file_ in os.listdir('bin'):
        if not file_.endswith('.error'):
            scripts.append(os.path.join('bin', file_))
    return scripts


class SmartInstall(install.install):
    def run(self):
        result = install.install.run(self)
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')

        # install the manpages
        #check we have rst2man
        try:
            import docutils
            have_rst2man = True
        except ImportError:
            have_rst2man = False

        man_dir = os.path.join(sys.prefix, 'share/man/man1')
        if have_rst2man:
            if not os.path.exists(man_dir):
                os.makedirs(man_dir)

            for fpath in os.listdir('doc'):
                if not fpath.endswith('.rst'):
                    continue
                rst_fpath = os.path.join('doc', fpath)
                man_fpath = os.path.join(man_dir,
                                         os.path.splitext(fpath)[0] + '.1')
                #print 'generating manpage: ', man_fpath
                subprocess.call(['rst2man.py', rst_fpath, man_fpath])

        return result


class InstallData(distutils.command.install_data.install_data):
    """need to change self.install_dir to the actual library dir"""
    def run(self):
        '''It modifies the place in which the thrid_party_binaries will be
        installed.'''
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return distutils.command.install_data.install_data.run(self)


setup_args = {
              'name': 'seq_crumbs',
              'version': __version__,
              'description': 'Small utilities for sequence files manipulation',
              'author': 'Jose Blanca & Peio Ziarsolo',
              'author_email': 'jblanca@upv.es',
              'url': 'http://bioinf.comav.upv.es/seq_crumbs/',
              'packages': ['crumbs', 'crumbs.third_party', 'crumbs.utils'],
              'include_package_data': True,
              'data_files': external_executables,
              'scripts': get_scripts(),
              'license': 'AGPL',
              'cmdclass': {'install': SmartInstall,
                           'install_data': InstallData}
              }

if _SETUPTOOLS:
    setup_args['install_requires'] = ['biopython >= 1.60']
    if version_info[0] < 3 or (version_info[0] == 3 and version_info[1] < 3):
        # until python 3.3 the standard file module has no support for 
        # wrapping file object and required to open a new file
        # bz2file is a backport of the python 3.3 std library module
        setup_args['install_requires'].append('bz2file')

setup(**setup_args)
