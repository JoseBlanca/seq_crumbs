'''
Created on 21/06/2012

@author: jose
'''


class UnknownFormatError(Exception):
    'Raised when the format of a sequence file cannot be guessed'
    pass


class FileNotFoundError(Exception):
    'The file does not exists'
    pass
