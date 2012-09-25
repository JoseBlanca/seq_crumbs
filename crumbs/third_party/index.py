# Copyright 2009-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

try:
    from collections import UserDict as _dict_base
except ImportError:
    from UserDict import DictMixin as _dict_base

from Bio._py3k import _bytes_to_string, _as_bytes
from Bio.SeqIO._index import SeqFileRandomAccess
from Bio.Alphabet import Alphabet, AlphabetEncoder
from Bio.SeqIO._index import _FormatToRandomAccess


def index(filename, format, alphabet=None, key_function=None):
    """Indexes a sequence file and returns a dictionary like object.

     - filename - string giving name of file to be indexed
     - format   - lower case string describing the file format
     - alphabet - optional Alphabet object, useful when the sequence type
                  cannot be automatically inferred from the file itself
                  (e.g. format="fasta" or "tab")
     - key_function - Optional callback function which when given a
                  SeqRecord identifier string should return a unique
                  key for the dictionary.

    This indexing function will return a dictionary like object, giving the
    SeqRecord objects as values:

    >>> from Bio import SeqIO
    >>> records = SeqIO.index("Quality/example.fastq", "fastq")
    >>> len(records)
    3
    >>> sorted(records)
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print records["EAS54_6_R1_2_1_540_792"].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> "EAS54_6_R1_2_1_540_792" in records
    True
    >>> print records.get("Missing", None)
    None

    If the file is BGZF compressed, this is detected automatically. Ordinary
    GZIP files are not supported:

    >>> from Bio import SeqIO
    >>> records = SeqIO.index("Quality/example.fastq.bgz", "fastq")
    >>> len(records)
    3
    >>> print records["EAS54_6_R1_2_1_540_792"].seq
    TTGGCAGGCCAAGGCCGATGGATCA

    Note that this psuedo dictionary will not support all the methods of a
    true Python dictionary, for example values() is not defined since this
    would require loading all of the records into memory at once.

    When you call the index function, it will scan through the file, noting
    the location of each record. When you access a particular record via the
    dictionary methods, the code will jump to the appropriate part of the
    file and then parse that section into a SeqRecord.

    Note that not all the input formats supported by Bio.SeqIO can be used
    with this index function. It is designed to work only with sequential
    file formats (e.g. "fasta", "gb", "fastq") and is not suitable for any
    interlaced file format (e.g. alignment formats such as "clustal").

    For small files, it may be more efficient to use an in memory Python
    dictionary, e.g.

    >>> from Bio import SeqIO
    >>> records = SeqIO.to_dict(SeqIO.parse(open("Quality/example.fastq"), "fastq"))
    >>> len(records)
    3
    >>> sorted(records)
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print records["EAS54_6_R1_2_1_540_792"].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>

    As with the to_dict() function, by default the id string of each record
    is used as the key. You can specify a callback function to transform
    this (the record identifier string) into your prefered key. For example:

    >>> from Bio import SeqIO
    >>> def make_tuple(identifier):
    ...     parts = identifier.split("_")
    ...     return int(parts[-2]), int(parts[-1])
    >>> records = SeqIO.index("Quality/example.fastq", "fastq",
    ...                       key_function=make_tuple)
    >>> len(records)
    3
    >>> sorted(records)
    [(413, 324), (443, 348), (540, 792)]
    >>> print records[(540, 792)].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> (540, 792) in records
    True
    >>> "EAS54_6_R1_2_1_540_792" in records
    False
    >>> print records.get("Missing", None)
    None

    Another common use case would be indexing an NCBI style FASTA file,
    where you might want to extract the GI number from the FASTA identifer
    to use as the dictionary key.

    Notice that unlike the to_dict() function, here the key_function does
    not get given the full SeqRecord to use to generate the key. Doing so
    would impose a severe performance penalty as it would require the file
    to be completely parsed while building the index. Right now this is
    usually avoided.

    See also: Bio.SeqIO.index_db() and Bio.SeqIO.to_dict()
    """
    #Try and give helpful error messages:
    if not isinstance(filename, basestring):
        raise TypeError("Need a filename (not a handle)")
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or \
                                     isinstance(alphabet, AlphabetEncoder)):
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    #Map the file format to a sequence iterator:
    return _IndexedSeqFileDict(filename, format, alphabet, key_function)


class _IndexedSeqFileDict(_dict_base):
    """Read only dictionary interface to a sequential sequence file.

    Keeps the keys and associated file offsets in memory, reads the file to
    access entries as SeqRecord objects using Bio.SeqIO for parsing them.
    This approach is memory limited, but will work even with millions of
    sequences.

    Note - as with the Bio.SeqIO.to_dict() function, duplicate keys
    (record identifiers by default) are not allowed. If this happens,
    a ValueError exception is raised.

    By default the SeqRecord's id string is used as the dictionary
    key. This can be changed by suppling an optional key_function,
    a callback function which will be given the record id and must
    return the desired key. For example, this allows you to parse
    NCBI style FASTA identifiers, and extract the GI number to use
    as the dictionary key.

    Note that this dictionary is essentially read only. You cannot
    add or change values, pop values, nor clear the dictionary.
    """
    def __init__(self, filename, format, alphabet, key_function):
        #Use key_function=None for default value
        try:
            proxy_class = _FormatToRandomAccess[format]
        except KeyError:
            raise ValueError("Unsupported format '%s'" % format)
        random_access_proxy = proxy_class(filename, format, alphabet)
        self._proxy = random_access_proxy
        self._key_function = key_function
        if key_function:
            offset_iter = ((key_function(k), o, l) for (k, o, l) in random_access_proxy)
        else:
            offset_iter = random_access_proxy
        offsets = {}
        for key, offset, length in offset_iter:
            #Note - we don't store the length because I want to minimise the
            #memory requirements. With the SQLite backend the length is kept
            #and is used to speed up the get_raw method (by about 3 times).
            #The length should be provided by all the current backends except
            #SFF where there is an existing Roche index we can reuse (very fast
            #but lacks the record lengths)
            #assert length or format in ["sff", "sff-trim"], \
            #       "%s at offset %i given length %r (%s format %s)" \
            #       % (key, offset, length, filename, format)
            if key in offsets:
                self._proxy._handle.close()
                raise ValueError("Duplicate key '%s'" % key)
            else:
                offsets[key] = offset
        self._offsets = offsets

    def __repr__(self):
        return "SeqIO.index(%r, %r, alphabet=%r, key_function=%r)" \
               % (self._proxy._handle.name, self._proxy._format,
                  self._proxy._alphabet, self._key_function)

    def __str__(self):
        if self:
            return "{%s : SeqRecord(...), ...}" % repr(self.keys()[0])
        else:
            return "{}"

    def __contains__(self, key) :
        return key in self._offsets

    def __len__(self):
        """How many records are there?"""
        return len(self._offsets)

    if hasattr(dict, "iteritems"):
        #Python 2, use iteritems but not items etc
        def values(self):
            """Would be a list of the SeqRecord objects, but not implemented.

            In general you can be indexing very very large files, with millions
            of sequences. Loading all these into memory at once as SeqRecord
            objects would (probably) use up all the RAM. Therefore we simply
            don't support this dictionary method.
            """
            raise NotImplementedError("Due to memory concerns, when indexing a "
                                      "sequence file you cannot access all the "
                                      "records at once.")

        def items(self):
            """Would be a list of the (key, SeqRecord) tuples, but not implemented.

            In general you can be indexing very very large files, with millions
            of sequences. Loading all these into memory at once as SeqRecord
            objects would (probably) use up all the RAM. Therefore we simply
            don't support this dictionary method.
            """
            raise NotImplementedError("Due to memory concerns, when indexing a "
                                      "sequence file you cannot access all the "
                                      "records at once.")

        def keys(self) :
            """Return a list of all the keys (SeqRecord identifiers)."""
            #TODO - Stick a warning in here for large lists? Or just refuse?
            return self._offsets.keys()

        def itervalues(self):
            """Iterate over the SeqRecord) items."""
            for key in self.__iter__():
                yield self.__getitem__(key)

        def iteritems(self):
            """Iterate over the (key, SeqRecord) items."""
            for key in self.__iter__():
                yield key, self.__getitem__(key)

        def iterkeys(self):
            """Iterate over the keys."""
            return self.__iter__()

    else:
        #Python 3 - define items and values as iterators
        def items(self):
            """Iterate over the (key, SeqRecord) items."""
            for key in self.__iter__():
                yield key, self.__getitem__(key)

        def values(self):
            """Iterate over the SeqRecord items."""
            for key in self.__iter__():
                yield self.__getitem__(key)

        def keys(self):
            """Iterate over the keys."""
            return self.__iter__()

    def __iter__(self):
        """Iterate over the keys."""
        return iter(self._offsets)

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        #Pass the offset to the proxy
        record = self._proxy.get(self._offsets[key])
        if self._key_function:
            key2 = self._key_function(record.id + ' ' + record.description)
        else:
            #key2 = record.id + ' ' + record.description
            key2 = record.description
        if key != key2:
            raise ValueError("Key did not match (%s vs %s)" % (key, key2))
        return record

    def get(self, k, d=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self.__getitem__(k)
        except KeyError:
            return d

    def get_raw(self, key):
        """Similar to the get method, but returns the record as a raw string.

        If the key is not found, a KeyError exception is raised.

        Note that on Python 3 a bytes string is returned, not a typical
        unicode string.

        NOTE - This functionality is not supported for every file format.
        """
        #Pass the offset to the proxy
        return self._proxy.get_raw(self._offsets[key])

    def __setitem__(self, key, value):
        """Would allow setting or replacing records, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    def update(self, *args, **kwargs):
        """Would allow adding more values, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    def pop(self, key, default=None):
        """Would remove specified record, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    def popitem(self):
        """Would remove and return a SeqRecord, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    def clear(self):
        """Would clear dictionary, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    def fromkeys(self, keys, value=None):
        """A dictionary method which we don't implement."""
        raise NotImplementedError("An indexed a sequence file doesn't "
                                  "support this.")

    def copy(self):
        """A dictionary method which we don't implement."""
        raise NotImplementedError("An indexed a sequence file doesn't "
                                  "support this.")


class FastqRandomAccess(SeqFileRandomAccess):
    """Random access to a FASTQ file (any supported variant).

    With FASTQ the records all start with a "@" line, but so can quality lines.
    Note this will cope with line-wrapped FASTQ files.
    """
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        id = None
        start_offset = handle.tell()
        line = handle.readline()
        if not line:
            #Empty file!
            return
        at_char = _as_bytes("@")
        plus_char = _as_bytes("+")
        if line[0:1] != at_char:
            raise ValueError("Problem with FASTQ @ line:\n%s" % repr(line))
        while line:
            #assert line[0]=="@"
            #This record seems OK (so far)
            id = line[1:].rstrip()
            #Find the seq line(s)
            seq_len = 0
            length = len(line)
            while line:
                line = handle.readline()
                length += len(line)
                if line.startswith(plus_char) : break
                seq_len += len(line.strip())
            if not line:
                raise ValueError("Premature end of file in seq section")
            #assert line[0]=="+"
            #Find the qual line(s)
            qual_len = 0
            while line:
                if seq_len == qual_len:
                    #Should be end of record...
                    end_offset = handle.tell()
                    line = handle.readline()
                    if line and line[0:1] != at_char:
                        ValueError("Problem with line %s" % repr(line))
                    break
                else:
                    line = handle.readline()
                    qual_len += len(line.strip())
                    length += len(line)
            if seq_len != qual_len:
                raise ValueError("Problem with quality section")
            yield _bytes_to_string(id), start_offset, length
            start_offset = end_offset
        #print "EOF"

    def get_raw(self, offset):
        """Similar to the get method, but returns the record as a raw string."""
        #TODO - Refactor this and the __init__ method to reduce code duplication?
        handle = self._handle
        handle.seek(offset)
        line = handle.readline()
        data = line
        at_char = _as_bytes("@")
        plus_char = _as_bytes("+")
        if line[0:1] != at_char:
            raise ValueError("Problem with FASTQ @ line:\n%s" % repr(line))
        identifier = line[1:].rstrip()
        #Find the seq line(s)
        seq_len = 0
        while line:
            line = handle.readline()
            data += line
            if line.startswith(plus_char) : break
            seq_len += len(line.strip())
        if not line:
            raise ValueError("Premature end of file in seq section")
        assert line[0:1] == plus_char
        #Find the qual line(s)
        qual_len = 0
        while line:
            if seq_len == qual_len:
                #Should be end of record...
                pos = handle.tell()
                line = handle.readline()
                if line and line[0:1] != at_char:
                    ValueError("Problem with line %s" % repr(line))
                break
            else:
                line = handle.readline()
                data += line
                qual_len += len(line.strip())
        if seq_len != qual_len:
            raise ValueError("Problem with quality section")
        return data
