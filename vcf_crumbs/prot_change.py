'''
Created on 2013 mai 24

@author: peio
'''
import subprocess
import os
from tempfile import NamedTemporaryFile

from Bio import SeqIO
from Bio.Seq import Seq


class IsIndelError(Exception):
    pass


class OutsideAlignment(Exception):
    pass


class BetweenSegments(Exception):
    pass


class SeqCoords(object):
    '''This class creates a coord system relation between tow secuences.

    First it alignes the sequences and creates de coord system.

    It uses internally 0 coord system.

    The coord system is a list of start, stop segments in both secuences
    ej:

              1         2         3         4         5         6
    01234567890123456789012345678901234567890123456789012345678901234
    ATCTAGGCTGCTACGATTAGCTGATCGATGTTATCGTAGATCTAGCTGATCATGCTAGCTGATCG
    ATCTAGGCTGCTACGA-TAGCTGATCGATGTTATCGTAGATCTAGCTGATCATGC-AGCTGATCG
    0-15, 17-54, 56-64
    0-15, 16-53, 54-62


              1         2         3         4         5         6
    01234567890123456789012345678901234567890123456789012345678901234
    012345678901234567890123 4567890123456789012345678901 2345678901234
    ATCTAGGCTGCTACGATTAGCTGA-CGATGTTATCGTAGATCTAGCTGATCAT-CTAGCTGATCG
    ATCTAGGCT-CTACGATTAGCTGATCGATGTTATCGTAGATC-AGCTGATCATGCTAGCTGATCG
    012345678 90123456789012345678901234567890 123456789012345678901234
    0-8, 10-23, 24-40, 42-51, 52-62
    0-8, 9-22,  24-40, 41-50, 52-62

    '''

    def __init__(self, seq1, seq2):
        "Both secuences are biopython secuences"
        self.coord_system = self._get_coord_system(seq1, seq2)
        self.seq1_name = seq1.id
        self.seq2_name = seq2.id
        self._seq2_len = len(seq2)

    def _get_coord_system(self, seq1, seq2):
        out_fhand, reverse = get_water_alignment(seq1, seq2)
        self.reverse = reverse
        coord_system = build_relations_from_aligment(out_fhand)
        out_fhand.close()
        return coord_system

    def _reverse_pos(self, pos):
        reverse = self.reverse
        if reverse:
            return self._seq2_len - pos - 1
        else:
            return pos

    def _get_segment(self, pos, seq_name):
        'returns the segment index of the given position'
        segments = self.coord_system[seq_name]
        for index, (start, stop) in enumerate(segments):
            if pos >= start and  pos <= stop:
                return index, (start, stop)
        if pos < segments[0][0] or pos > segments[-1][1]:
            raise OutsideAlignment
        else:
            raise BetweenSegments

    def _to_seq_pos(self, pos, to_seq1=True):
        if to_seq1:
            seq1_name = self.seq1_name
            seq2_name = self.seq2_name
        else:
            seq2_name = self.seq1_name
            seq1_name = self.seq2_name

        segment2 = self._get_segment(pos, seq2_name)
        segment2_index, segment2 = segment2
        segment1 = self.coord_system[seq1_name][segment2_index]
        return segment1[0] + pos - segment2[0]

    def to_seq1_pos(self, seq2_pos):
        seq2_pos = self._reverse_pos(seq2_pos)
        return self._to_seq_pos(seq2_pos, to_seq1=True)

    def to_seq2_pos(self, seq1_pos):
        seq2_pos = self._to_seq_pos(seq1_pos, to_seq1=False)
        return self._reverse_pos(seq2_pos)

    def _to_seq_slice(self, start, end, to_seq1=True):
        if to_seq1:
            seq1_name = self.seq1_name
            seq2_name = self.seq2_name
        else:
            seq2_name = self.seq1_name
            seq1_name = self.seq2_name

        stop = end - 1
        segment2_start = self._get_segment(start, seq2_name)
        segment2_stop = self._get_segment(stop, seq2_name)

        segment2_index_start, segment2_start = segment2_start
        segment2_index_stop, segment2_stop = segment2_stop

        if segment2_index_start != segment2_index_stop:
            raise BetweenSegments

        segment1 = self.coord_system[seq1_name][segment2_index_start]
        start = segment1[0] + start - segment2_start[0]
        stop = segment1[0] + stop - segment2_stop[0]
        return (start, stop + 1)

    def to_seq1_slice(self, start, end):
        if self.reverse:
            start = self._reverse_pos(start)
            end = self._reverse_pos(end)
        slice2 = self._to_seq_slice(start, end, to_seq1=True)
        if self.reverse:
            return slice2[1], slice2[0]
        return slice2

    def to_seq2_slice(self, start, end):
        slice1 = self._to_seq_slice(start, end, to_seq1=False)
        if self.reverse:
            start = self._reverse_pos(slice1[1])
            end = self._reverse_pos(slice1[0])
        else:
            start, end = slice1
        return (start, end)


def build_relations_from_aligment(fhand):
    'It returns a relations dict given an alignment in markx10 format'
    #print open(fhand.name).read()
    #we parse the aligment
    in_seq_section = 0
    seq, al_start, seq_len = None, None, None
    seq0_name = None
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>' and line[1] != '>':
            if seq0_name is None:
                seq0_name = line.split()[0][1:]
            else:
                seq1_name = line.split()[0][1:]
            if in_seq_section:
                seq0 = {'seq': seq,
                        'length': seq_len,
                        'al_start': al_start - 1,
                        'name': seq0_name}
            in_seq_section += 1
            seq = ''
            continue
        if not in_seq_section:
            continue
        if '; sq_len:' in line:
            seq_len = int(line.split(':')[-1])
        if '; al_display_start:' in line:
            al_start = int(line.split(':')[-1])
        if line[0] not in (';', '#'):
            seq += line

    seq1 = {'seq': seq,
            'length': seq_len,
            'al_start': al_start - 1,
            'name': seq1_name}

    #now we get the segments
    gap = '-'
    segments = []
    segment0, segment1 = None, None

    seq0_start, seq1_start = seq0['al_start'], seq1['al_start']
    seq0_start_delta, seq1_start_delta = seq0_start, seq1_start
    seq0_delta, seq1_delta = 0, 0
    for index, (nucl0, nucl1) in enumerate(zip(seq0['seq'], seq1['seq'])):
        seq0_index = seq0_start_delta + index - seq0_delta
        seq1_index = seq1_start_delta + index - seq1_delta
        if nucl0 == gap:
            segment0 = seq0_start, seq0_index - 1
            segment1 = seq1_start, seq1_index - 1
            seq0_start = seq0_index
            seq1_start = seq1_index + 1
            seq0_delta += 1

        elif nucl1 == gap:
            segment0 = seq0_start, seq0_index - 1
            segment1 = seq1_start, seq1_index - 1
            seq1_start = seq1_index
            seq0_start = seq0_index + 1
            seq1_delta += 1

        if segment0 and segment1:
            segment = {seq0['name']: segment0, seq1['name']: segment1}
            segments.append(segment)
            segment0, segment1 = None, None
    else:
        segment0 = seq0_start, seq0_index
        segment1 = seq1_start, seq1_index
        segment = {seq0['name']: segment0, seq1['name']: segment1}
        segments.append(segment)

    relations = {}
    for seg in segments:
        for seq_name, limits in seg.items():
            if seq_name not in relations:
                relations[seq_name] = []
            relations[seq_name].append(limits)
    return relations


def _get_water_score(fhand):
    for line in fhand:
        if line.startswith('# Score:'):
            return float(line.split(':')[1].strip())
    return None


def get_water_alignment(seq1, seq2, gap_open=10.0, gap_extend=0.5,
                        out_fmt='markx10'):
    out_fhand = NamedTemporaryFile()
    _do_water_alignment(seq1, seq2,  out_fhand, gap_open=10.0, gap_extend=0.5,
                       out_fmt='markx10', reverse2=False)
    out_fhand2 = NamedTemporaryFile()
    _do_water_alignment(seq1, seq2,  out_fhand2, gap_open=10.0, gap_extend=0.5,
                       out_fmt='markx10', reverse2=True)
    forw_score = _get_water_score(out_fhand)
    rev_score = _get_water_score(out_fhand2)
    if forw_score > rev_score:
        out_fhand.seek(0)
        return out_fhand, False
    else:
        out_fhand2.seek(0)
        return out_fhand2, True


def _do_water_alignment(seq1, seq2,  out_fhand, gap_open=10.0, gap_extend=0.5,
                       out_fmt='markx10', reverse2=False):
    seq1_fhand = NamedTemporaryFile()
    seq2_fhand = NamedTemporaryFile()

    SeqIO.write(seq1, seq1_fhand, 'fasta')
    SeqIO.write(seq2, seq2_fhand, 'fasta')
    seq1_fhand.flush()
    seq2_fhand.flush()
    cmd = ['water', '-asequence', seq1_fhand.name, '-bsequence',
           seq2_fhand.name, '-outfile', out_fhand.name, '-gapopen',
           str(gap_open), '-gapextend', str(gap_extend), '-aformat3', out_fmt]
    if reverse2:
        cmd.append('-sreverse2')
    stdout = open(os.devnull, 'w')
    stderr = open(os.devnull, 'w')
    subprocess.check_call(cmd, stdout=stdout, stderr=stderr)


def get_amino_change(seq_ref, seq_estscan, record):
    if record.is_indel:
        raise IsIndelError()
    position = record.POS - 1
    alt_allele = record.alleles[1].sequence
    seq_coord = SeqCoords(seq_ref, seq_estscan)

    estscan_pos = seq_coord.to_seq2_pos(position)
    if estscan_pos is None:
        return None

    estscan_frame = (estscan_pos % 3) + 1
    estscan_start = estscan_pos + estscan_frame - 1
    estscan_stop = estscan_start + 2

    # check if there is a frameshift in the ref_seq
    ref_slice = seq_coord.to_seq1_slice(estscan_start, estscan_stop)
    if ref_slice is None:
        return None

    ref_seq_aa = seq_ref[ref_slice[0]: ref_slice[1] + 1].seq[:3].translate()
    estscan_seq_aa = seq_estscan[estscan_start: estscan_stop + 1].seq[:3]

    ref_aa = str(estscan_seq_aa.translate())

    if str(ref_seq_aa) != str(ref_aa):
        return None
    aminos = {'ref_amino': ref_aa, 'alt_amino': []}

    for alt_allele in record.alleles[1:]:
        alt_allele = alt_allele.sequence

        alt_seq = [nucl for nucl in (estscan_seq_aa)]
        alt_seq[estscan_frame - 1] = alt_allele
        alt_seq = Seq("".join(alt_seq))
        alt_aa = str(alt_seq.translate())
        aminos['alt_amino'].append(alt_aa)

    return aminos
