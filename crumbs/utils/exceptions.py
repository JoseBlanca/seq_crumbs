'''
Created on 21/06/2012

@author: jose
'''


class UnknownFormatError(Exception):
    'Raised when the format of a sequence file cannot be guessed'
    pass


class WrongFormatError(Exception):
    'Raised when the file have an unexpected format'
    pass


class FileNotFoundError(Exception):
    'The file does not exists'
    pass


class TooManyFiles(Exception):
    pass


class MalformedFile(Exception):
    'The input file is malformed'
    pass
