#!/usr/bin/env python

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

import argparse
import os
import platform
import shutil
import subprocess
import tarfile
import sys

from os.path import join

from crumbs import settings


def _setup_argparse():
    'It prepares the command line argument parsing.'
    description = 'creates the binary distribution'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-p', '--pyinstaller_dir', dest='pyinstaller_dir',
                        help='Sequence input files to process')
    parser.add_argument('-a', '--app_dir', dest='app_dir',
                        help='app source code root dir')
    return parser


def _parse_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    parsed_args = parser.parse_args()

    if parsed_args.pyinstaller_dir is None:
        parser.error('pyinstaller dir is mandatory')
    else:
        pyinstaller_dir = os.path.abspath(parsed_args.pyinstaller_dir)

    if parsed_args.app_dir is None:
        parser.error('app dir is mandatory')
    else:
        app_dir = os.path.abspath(parsed_args.app_dir)

    return pyinstaller_dir, app_dir


def get_version(app_dir):
    'it gets the version of the app'
    version = "Undefined"
    for line in open(join(app_dir, 'crumbs/__init__.py')):
        if (line.startswith('__version__')):
            version = line.split('=')[1].strip()
    version = version.replace("'", "")
    return version


def check_pyinstaller_sanity(dir_):
    'It checks that the pyinstaller dir is a valid pyinstaller directory'
    files_ = os.listdir(dir_)
    assert 'PyInstaller' in files_
    assert 'pyinstaller.py' in files_


def get_platform_bin_dir(app_dir):
    '''It returns the platform specific bindir. It returns the relative path
    from the source code root dir'''
    system = platform.system().lower()
    arch = platform.architecture()[0]

    return os.path.join(app_dir, 'crumbs', 'third_party', 'bin', system, arch)


def get_bin_dist_directory(app_dir):
    'It creates the bin dist name'
    # get version
    version = get_version(app_dir)
    arch = platform.architecture()[0]
    arch = 'i686' if arch == '32bits' else 'x64'

    bin_dist_name = 'seq_crumbs-{0:s}-{1:s}-linux'.format(version, arch)
    return join(app_dir, 'dist', bin_dist_name)


def copy_and_rename_ext_bin(app_dir, bin_dist_dir):
    ''' It copies and  changes the name of the platform specific external
     binaries.
     It copies them to the binary dist directory'''

    external_bin_dir = get_platform_bin_dir(app_dir)
    for ext_bin in os.listdir(external_bin_dir):
        shutil.copy(join(external_bin_dir, ext_bin),
                join(bin_dist_dir, settings.EXTERNAL_BIN_PREFIX + ext_bin))


def write_lines(lines, fpath):
    'it writes a list of lines to a fhand'
    settings_fhand = open(fpath, 'w')
    for line in lines:
        settings_fhand.write(line)
    settings_fhand.flush()
    settings_fhand.close()


def modify_settings_content(lines):
    'it modifies the settings content to use prefix and remove path'
    new_lines = []
    for line in lines:
        if line.startswith('USE_EXTERNAL_BIN_PREFIX'):
            line = 'USE_EXTERNAL_BIN_PREFIX = True\n'
        if line.startswith('ADD_PATH_TO_EXT_BIN'):
            line = 'ADD_PATH_TO_EXT_BIN = False\n'
        new_lines.append(line)
    return new_lines


def make_binary(script_path, pyinstaller_dir, bin_dist_dir):
    'It converts a script into a binary file using pyinstaller'
    initial_wd = os.getcwd()
    script_name = os.path.basename(script_path)
    os.chdir(pyinstaller_dir)
    subprocess.check_call([sys.executable, 'pyinstaller.py', '--onefile',
                           script_path])

    subprocess.check_call([sys.executable, 'pyinstaller.py', join(script_name,
                                                      script_name + '.spec')])

    shutil.copy(join(script_name, 'dist', script_name),
                join(bin_dist_dir, script_name))

    shutil.rmtree(script_name)
    os.chdir(initial_wd)


def main():
    'the main part'
    parser = _setup_argparse()
    pyinstaller_dir, app_dir = _parse_args(parser)

    check_pyinstaller_sanity(pyinstaller_dir)

    # make directory where all bynaries will go
    bin_dist_dir = get_bin_dist_directory(app_dir)
    if os.path.exists(bin_dist_dir):
        shutil.rmtree(bin_dist_dir)
    os.makedirs(bin_dist_dir)
    try:
        # copy the external executables changing the name
        copy_and_rename_ext_bin(app_dir, bin_dist_dir)

        # modify the conf to compile the binaries with the get_binary()
        # modified
        settings_fpath = join(app_dir, 'crumbs', 'settings.py')
        settings_content = open(settings_fpath).readlines()
        modified_settings = modify_settings_content(settings_content)
        write_lines(modified_settings, settings_fpath)

        # make binary for each executable
        for command in os.listdir(join(app_dir, 'bin')):
            script_path = join(app_dir, 'bin', command)
            first_line = open(script_path).readline()
            if first_line.startswith('#!/usr/bin/env'):
                make_binary(script_path, pyinstaller_dir, bin_dist_dir)
    except Exception:
        # return the settings file to its origin
        write_lines(settings_content, settings_fpath)

    write_lines(settings_content, settings_fpath)
    #make targz
    tar_fpath = bin_dist_dir + '.tar.gz'
    tar = tarfile.open(tar_fpath, 'w:gz')
    tar.add(bin_dist_dir, os.path.basename(bin_dist_dir))
    tar.close()
    shutil.rmtree(bin_dist_dir)

if __name__ == '__main__':
    main()
