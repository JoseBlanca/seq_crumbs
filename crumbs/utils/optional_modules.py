# Copyright 2012-2013 Jose Blanca, Peio Ziarsolo,
# COMAV-Univ. Politecnica Valencia
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
In this module all the imports of the optional modules are handled
'''

from crumbs.exceptions import OptionalRequirementError

MSG = 'A python package to run this executable is required,'
MSG += ' but it is not installed: '
BIO = 'biopython'
BIO_BGZF = 'biopython with Bgzf support'


def create_fake_class(msg):
    class FakePythonRequiredClass(object):
        def __init__(self, *args, **kwargs):
            raise OptionalRequirementError(msg)
    return FakePythonRequiredClass


def create_fake_funct(msg):
    def FakeRequiredfunct(*args, **kwargs):
        raise OptionalRequirementError(msg)
    return FakeRequiredfunct

try:
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
except ImportError:
    SeqRecord = create_fake_class(MSG + BIO)
    Seq = create_fake_class(MSG + BIO)
    SeqFeature = create_fake_class(MSG + BIO)
    FeatureLocation = create_fake_class(MSG + BIO)

try:
    from Bio.SeqIO.SffIO import SffIterator
except ImportError:
    SffIterator = create_fake_class(MSG + BIO)

try:
    from Bio.bgzf import BgzfWriter
except ImportError:
    BgzfWriter = create_fake_class(MSG + BIO_BGZF)


try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError:
    from crumbs.utils.biopython_code import FastqGeneralIterator

try:
    from Bio.SeqIO.FastaIO import FastaIterator
except ImportError:
    FastaIterator = create_fake_class(MSG + BIO)

try:
    from Bio.SeqIO.QualityIO import QualPhredIterator
    from Bio.SeqIO.QualityIO import PairedFastaQualIterator
    from Bio.SeqIO.QualityIO import FastqPhredIterator
    from Bio.SeqIO.QualityIO import FastqSolexaIterator
    from Bio.SeqIO.QualityIO import FastqIlluminaIterator
    from Bio.SeqIO import parse as parse_into_seqrecs
    from Bio.SeqIO import write as write_seqrecs
except ImportError:
    QualPhredIterator = create_fake_class(MSG + BIO)
    PairedFastaQualIterator = create_fake_class(MSG + BIO)
    FastqPhredIterator = create_fake_class(MSG + BIO)
    FastqSolexaIterator = create_fake_class(MSG + BIO)
    FastqIlluminaIterator = create_fake_class(MSG + BIO)
    parse_into_seqrecs = create_fake_funct(MSG + BIO)
    write_seqrecs = create_fake_funct(MSG + BIO)

try:
    from Bio.Blast import NCBIXML, NCBIWWW
except ImportError:
    NCBIXML = create_fake_class(MSG + BIO)
    NCBIWWW = create_fake_class(MSG + BIO)

try:
    from Bio._py3k import _bytes_to_string, _as_bytes
except ImportError:
    _bytes_to_string = create_fake_funct(MSG + BIO)
    _as_bytes = create_fake_funct(MSG + BIO)

try:
    from Bio.SeqIO._index import SeqFileRandomAccess
    from Bio.SeqIO._index import _FormatToRandomAccess
    from Bio.SeqIO import _index
except ImportError:
    SeqFileRandomAccess = create_fake_class(MSG + BIO)
    _FormatToRandomAccess = create_fake_class(MSG + BIO)
    index = create_fake_class(MSG + BIO)

try:
    from Bio.Alphabet import Alphabet, AlphabetEncoder
except ImportError:
    Alphabet = create_fake_class(MSG + BIO)
    AlphabetEncoder = create_fake_class(MSG + BIO)

try:
    from toolz.itertoolz.core import merge_sorted, first
except ImportError:
    merge_sorted = create_fake_funct(MSG + 'toolz')
    first = create_fake_funct(MSG + 'toolz')

try:
    from pysam import Samfile
except ImportError:
    Samfile = create_fake_funct(MSG + 'pysam')
