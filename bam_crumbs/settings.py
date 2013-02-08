# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of bam_crumbs.
# bam_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# bam_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with bam_crumbs. If not, see <http://www.gnu.org/licenses/>.
import os
import tempfile


# # Use this to modify how get_binary path works
# if need to modify the binary's name
_USE_EXTERNAL_BIN_PREFIX = False
# prefix to add to the binary name
_EXTERNAL_BIN_PREFIX = 'bams_crumbs'
# mark True if need the path or assumes that is on the path
_ADD_PATH_TO_EXT_BIN = True


_THIRD_PART_JAVA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                   'third_party', 'java')
_PICARD_TOOLS_DIR = os.path.join(_THIRD_PART_JAVA_DIR, 'picard-tools')

_TEMP_DIR = None

_DEFAULT_N_BINS = 80
_DEFAULT_N_MOST_ABUNDANT_REFERENCES = 40


class _Settings(dict):
    '''A class that stores the seq_crumbs settings.'''
    def __init__(self):
        'It inits the class'
        super(_Settings, self).__init__()
        self.load_settings()
        tempfile.tempdir = self.__getitem__('TEMP_DIR')

    def load_settings(self):
        'It loads the settings defined in this module'
        for key, val in globals().viewitems():
            if not key.isupper():
                continue
            key = key[1:]  # strip the underscore
            super(_Settings, self).__setitem__(key, val)

        # Are there any environment variable to update the settings?
        for key, value in os.environ.items():
            if key.startswith('SEQ_CRUMBS_'):
                key = key[11:]
                if key in self.viewkeys():
                    value = type(key)(value)
                    super(_Settings, self).__setitem__(key, value)

_settings = _Settings()


def get_settings():
    'It returns the settings'
    # load the settings defined in this module
    return _settings


def get_setting(key):
    'It returns the value for one setting'
    return _settings[key]
