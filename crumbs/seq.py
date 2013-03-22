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

from Bio.SeqRecord import SeqRecord

from crumbs.utils.file_formats import _remove_multiline
from crumbs.utils.tags import SEQITEM, SEQRECORD

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
        title = seq.lines[0][1:]
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


def _sitem_fastq_plus_line_index(lines, seq_name):
    for index, line in enumerate(lines):
        if _is_fastq_plus_line(line, seq_name):
            return index
    raise ValueError('No fastq plus line in the given lines')


def _get_seqitem_str_lines(seq):
    fmt = seq.file_format
    sitem = seq.object
    sname = sitem.name
    if 'fastq' in fmt and not 'multiline' in fmt:
        lines = sitem.lines[1:2]
    else:
        lines = (_break() if _is_fastq_plus_line(l, sname) else l for l in sitem.lines[1:])
    return lines


def _get_seqitem_qual_lines(seq):
    fmt = seq.file_format
    sitem = seq.object
    if 'fastq' in fmt:
        if 'multiline' in fmt:
            lines = sitem.lines[_sitem_fastq_plus_line_index(sitem.lines, get_name(seq)) + 1:]
        else:
            lines = sitem.lines[3:4]
    return lines


def get_str_seq(seq):
    seq_class = seq.kind
    if seq_class == SEQITEM:
        seq = ''.join((line.rstrip() for line in _get_seqitem_str_lines(seq)))
    elif seq_class == SEQRECORD:
        seq = str(seq.object.seq)
    return seq


def get_length(seq):
    seq_class = seq.kind
    if seq_class == SEQITEM:
        length = lambda l: len(l) - 1   # It assumes line break and no spaces
        length = sum(map(length, _get_seqitem_str_lines(seq)))
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
        if 'multiline' in fmt:
            lines = (l.rstrip() for l in _get_seqitem_qual_lines(seqwrap))
        else:
            lines = seqwrap.object.lines[-1:]
        lines = [line.strip() for line in lines]
        quals = itertools.chain(quals_map[char] for l in lines for char in l)
    else:
        raise RuntimeError('Qualities requested for an unknown SeqItem format')
    return quals


def get_qualities(seq):
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
        elif 'multiline' in fmt and 'fastq' in fmt:
            qline = ''.join([qline.strip() for qline in _get_seqitem_qual_lines(seqwrapper)])
            lines = [lines[0], seq + '\n', '+\n', qline + '\n']
            fmt = _remove_multiline(fmt)
            if len(lines[1]) != len(lines[3]):
                msg = 'Sequence and quality line length do not match'
                raise ValueError(msg)
        elif 'fastq' in fmt:
            lines = [lines[0], seq + '\n', lines[2], lines[3]]
            if len(lines[1]) != len(lines[3]):
                msg = 'Sequence and quality line length do not match'
                raise ValueError(msg)
        else:
            raise RuntimeError('Unknown format for a SequenceItem')

    if name:
        lines[0] = lines[0][0] + name + '\n'
        if 'fastq' in fmt:
            if 'multiline' in fmt:
                line_plus_index = _sitem_fastq_plus_line_index
            else:
                line_plus_index = 2
            lines[line_plus_index] = '+\n'
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
        qual_str = ''.join(line.rstrip() for line in _get_seqitem_qual_lines(seqwrap))
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
    return SeqWrapper(seq.kind, object=seq_obj,
                      file_format=_remove_multiline(seq.file_format))


def assing_kind_to_seqs(kind, seqs, file_format):
    'It puts each seq into a NamedTuple named Seq'
    return (SeqWrapper(kind, seq, file_format) for seq in seqs)
