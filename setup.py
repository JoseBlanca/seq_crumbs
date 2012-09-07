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
import fnmatch
import glob
import platform
import subprocess

import distutils.command.install_data

try:
    from setuptools import setu
    from setuptools.command import install
    _SETUPTOOLS = True

except ImportError:
    from distutils.core import setup
    from distutils.command import install
    _SETUPTOOLS = False

print "using_setuptools", _SETUPTOOLS

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
        for wc in wildcards:
            wc_name = opj(dirname, wc)
            for f in files:
                filename = opj(dirname, f)
                if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
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


class smart_install(install.install):
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
                print 'generating manpage: ', man_fpath
                subprocess.call(['rst2man.py', rst_fpath, man_fpath])

        return result


class install_data(distutils.command.install_data.install_data):
    """need to change self.install_dir to the actual library dir"""
    def run(self):
        # cambiar el sitio donde se van a instalar los third_party_binaries
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
              'data_files': external_executables,
              'scripts': get_scripts(),
              'license': 'AGPL',
              'cmdclass': {'install': smart_install,
                           'install_data': install_data}
              }

if _SETUPTOOLS:
    requires = ['biopython']
    try:
        import argparse
    except ImportError:
        requires.append(argparse)
    setup_args['install_requires'] = requires

setup(**setup_args)
