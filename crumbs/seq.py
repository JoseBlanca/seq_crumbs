# Copyright 2013 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
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

import itertools
from copy import deepcopy
from collections import namedtuple

from crumbs.utils.optional_modules import SeqRecord
from crumbs.utils.tags import (SEQITEM, SEQRECORD, ILLUMINA_QUALITY,
                               SANGER_QUALITY, SANGER_FASTQ_FORMATS,
                               ILLUMINA_FASTQ_FORMATS)

# pylint: disable=C0111


SeqWrapper = namedtuple('SeqWrapper', ['kind', 'object', 'file_format'])
_SeqItem = namedtuple('SeqItem', ['name', 'lines', 'annotations'])


class SeqItem(_SeqItem):
    def __new__(cls, name, lines, annotations=None):
        # This subclass is required to have a default value in a namedtuple
        if annotations is None:
            annotations = {}
        # add default values
        return super(SeqItem, cls).__new__(cls, name, lines, annotations)


def get_title(seq):
    'Given a seq it returns the title'
    seq_class = seq.kind
    seq = seq.object

    if seq_class == SEQITEM:
        title = seq.lines[0][1:].rstrip()
    elif seq_class == SEQRECORD:
        title = seq.id + ' ' + seq.description
    else:
        msg = 'Do not know how to guess title form this seq class'
        raise NotImplementedError(msg)
    return title


def get_description(seq):
    seq_class = seq.kind
    seq = seq.object
    if seq_class == SEQITEM:
        title_items = seq.lines[0].split(' ', 1)
        desc = title_items[1] if len(title_items) == 2 else None
    elif seq_class == SEQRECORD:
        desc = seq.description
        if desc == '<unknown description>':  # BioPython default
            return None
    return desc


def get_name(seq):
    if 'SeqRecord' in seq.__class__.__name__:
        seq_class = SEQRECORD
    else:
        seq_class = seq.kind
        seq = seq.object
    if seq_class == SEQITEM:
        name = seq.name
    elif seq_class == SEQRECORD:
        name = seq.id
    return name


def get_file_format(seq):
    seq_class = seq.kind

    if seq_class == SEQITEM:
        fmt = seq.file_format
    elif seq_class == SEQRECORD:
        fmt = None
    return fmt


def _break():
    raise StopIteration


def _is_fastq_plus_line(line, seq_name):
    if line == '+\n' or line.startswith('+') and seq_name in line:
        return True
    else:
        return False


def _get_seqitem_quals(seq):
    fmt = seq.file_format
    sitem = seq.object
    if 'fastq' in fmt:
        quals = sitem.lines[3].rstrip()
    else:
        quals = None
    return quals


def get_str_seq(seq):
    seq_class = seq.kind
    if seq_class == SEQITEM:
        seq = seq.object.lines[1].rstrip()
    elif seq_class == SEQRECORD:
        seq = str(seq.object.seq)
    return seq


def get_length(seq):
    seq_class = seq.kind
    if seq_class == SEQITEM:
        length = len(seq.object.lines[1]) -1
    elif seq_class == SEQRECORD:
        length = len(seq.object)
    return length


SANGER_QUALS = {chr(i): i - 33 for i in range(33, 127)}
ILLUMINA_QUALS = {chr(i): i - 64 for i in range(64, 127)}


def _get_seqitem_qualities(seqwrap):
    fmt = seqwrap.file_format.lower()
    if 'fasta' in fmt:
        raise AttributeError('A fasta file has no qualities')
    elif 'fastq' in fmt:
        if 'illumina' in fmt:
            quals_map = ILLUMINA_QUALS
        else:
            quals_map = SANGER_QUALS
        encoded_quals = seqwrap.object.lines[3].rstrip()
        quals = [quals_map[qual] for qual in encoded_quals]
    else:
        raise RuntimeError('Qualities requested for an unknown SeqItem format')
    return quals


def get_int_qualities(seq):
    seq_class = seq.kind
    if seq_class == SEQITEM:
        return _get_seqitem_qualities(seq)
    elif seq_class == SEQRECORD:
        try:
            quals = seq.object.letter_annotations['phred_quality']
        except KeyError:
            msg = 'The given SeqRecord has no phred_quality'
            raise AttributeError(msg)
        return quals


SANGER_STRS = {i - 33: chr(i)  for i in range(33, 127)}
ILLUMINA_STRS = {i - 64: chr(i) for i in range(64, 127)}


def _int_quals_to_str_quals(int_quals, out_format):
    if out_format == SANGER_QUALITY:
        quals_map = SANGER_STRS
    elif out_format == ILLUMINA_QUALITY:
        quals_map = ILLUMINA_STRS
    else:
        msg = 'Unknown or not supported quality format'
        raise ValueError(msg)
    return (quals_map[int_quality] for int_quality in int_quals)


def get_str_qualities(seq, out_format=None):
    if out_format is None:
        out_format = seq.file_format
    if out_format in SANGER_FASTQ_FORMATS:
        out_format = SANGER_QUALITY
    elif out_format in ILLUMINA_FASTQ_FORMATS:
        out_format = ILLUMINA_QUALITY

    seq_class = seq.kind
    if seq_class == SEQITEM:
        in_format = seq.file_format
        if 'fasta' in in_format:
            raise ValueError('A fasta file has no qualities')
        if in_format in SANGER_FASTQ_FORMATS:
            in_format = SANGER_QUALITY
        elif in_format in ILLUMINA_FASTQ_FORMATS:
            in_format = ILLUMINA_QUALITY
        else:
            msg = 'Unknown or not supported quality format: '
            msg += in_format
            raise ValueError(msg)
        if in_format == out_format:
            quals = seq.object.lines[3].rstrip()
        else:
            int_quals = get_int_qualities(seq)
            quals = ''.join(_int_quals_to_str_quals(int_quals, out_format))
    elif seq_class == SEQRECORD:
        int_quals = get_int_qualities(seq)
        quals = ''.join(_int_quals_to_str_quals(int_quals, out_format))
    return quals


def get_annotations(seq):
    return seq.object.annotations


def _copy_seqrecord(seqrec, seq=None, name=None, id_=None):
    'Given a seqrecord it returns a new seqrecord with seq or qual changed.'
    if seq is None:
        seq = seqrec.seq
    if id_ is  None:
        id_ = seqrec.id
    if name is None:
        name = seqrec.name

    # the letter annotations
    let_annot = {annot: v for annot, v in seqrec.letter_annotations.items()}

    # the rest of parameters
    description = seqrec.description
    dbxrefs = seqrec.dbxrefs[:]
    features = seqrec.features[:]  # the features are not copied
    annotations = deepcopy(seqrec.annotations)

    # the new sequence
    new_seq = SeqRecord(seq=seq, id=id_, name=name, description=description,
                        dbxrefs=dbxrefs, features=features,
                        annotations=annotations, letter_annotations=let_annot)

    return new_seq


def _copy_seqitem(seqwrapper, seq=None, name=None):
    seq_item = seqwrapper.object
    lines = seq_item.lines
    fmt = seqwrapper.file_format
    if seq is None:
        lines = lines[:]
    else:
        if 'fasta' in fmt:
            lines = [lines[0], seq + '\n']
        elif 'fastq' in fmt:
            lines = [lines[0], seq + '\n', lines[2], lines[3]]
            if len(lines[1]) != len(lines[3]):
                msg = 'Sequence and quality line length do not match'
                raise ValueError(msg)
        else:
            raise RuntimeError('Unknown format for a SequenceItem')

    if name:
        # name title line
        lines[0] = lines[0][0] + name + '\n'
        # change + line in case has the name in it.
        if 'fastq' in fmt:
            lines[2] = '+\n'
    name = seq_item.name if name is None else name

    annotations = seq_item.annotations
    if annotations is not None:
        annotations = annotations.copy()
    seq = SeqWrapper(kind=seqwrapper.kind,
                     object=SeqItem(name, lines, annotations),
                     file_format=fmt)
    return seq


def copy_seq(seqwrapper, seq=None, name=None):
    seq_class = seqwrapper.kind
    seq_obj = seqwrapper.object
    if seq_class == SEQITEM:
        seq = _copy_seqitem(seqwrapper, seq=seq, name=name)
    elif seq_class == SEQRECORD:
        seq_obj = _copy_seqrecord(seq_obj, seq=seq, name=name, id_=name)
        seq = SeqWrapper(kind=seqwrapper.kind, object=seq_obj,
                         file_format=seqwrapper.file_format)
    return seq


def _slice_seqitem(seqwrap, start, stop):
    fmt = seqwrap.file_format
    seq_obj = seqwrap.object
    lines = seq_obj.lines
    seq_str = get_str_seq(seqwrap)
    seq_str = seq_str[start: stop] + '\n'
    if 'fasta' in fmt:
        lines = [lines[0], seq_str]
    elif 'fastq' in fmt:
        qual_str = get_str_qualities(seqwrap)
        qual_str = qual_str[start: stop]
        qual_str += '\n'
        lines = [lines[0], seq_str, '+\n', qual_str]
    else:
        raise ValueError('Unknown SeqItem type')
    seq_obj = SeqItem(name=seq_obj.name, lines=lines,
                      annotations=seq_obj.annotations)
    return seq_obj


def slice_seq(seq, start=None, stop=None):
    seq_class = seq.kind
    if seq_class == SEQITEM:
        seq_obj = _slice_seqitem(seq, start, stop)
    elif seq_class == SEQRECORD:
        seq_obj = seq.object[start:stop]
    return SeqWrapper(seq.kind, object=seq_obj, file_format=seq.file_format)


def assing_kind_to_seqs(kind, seqs, file_format):
    'It puts each seq into a NamedTuple named Seq'
    return (SeqWrapper(kind, seq, file_format) for seq in seqs)
