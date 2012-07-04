
# similar software SffToCA
from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.blast import get_blast_db, do_blast, generate_tabblast_format
from crumbs.seqio import write_seqrecords, read_seqrecords, guess_format
from crumbs.alignment_result import (TabularBlastParser, filter_alignments,
                                     covered_segments, ELONGATED,
                                     elongate_match_parts_till_global)


FLX_LINKER = 'GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC'
# The titanium linker can be forward or reverse
TITANIUM_LINKER = 'TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG'
TITANIUM_LINKER_REV = 'CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA'

LINKERS = [SeqRecord(Seq(FLX_LINKER), id='flx_linker'),
           SeqRecord(Seq(TITANIUM_LINKER), id='titanium_linker')]


def _do_blast(seq_fhand, oligos_fhand):
    'It returns an alignment result with the blast.'
    # TODO check if it matters if the oligo or the sequences are used to build
    # the database. Usually there will be only one linker.
    blastdb = get_blast_db(seq_fhand.name, dbtype='nucl')

    blast_format = ['query', 'query_length', 'subject', 'subject_length',
                    'query_start', 'query_end', 'subject_start',
                    'subject_end', 'expect', 'identity']
    fmt = generate_tabblast_format(blast_format)
    params = {'task': 'blastn-short', 'outfmt': fmt}

    blast_fhand = NamedTemporaryFile(prefix='oligos_vs_seqs.', suffix='.blast')
    do_blast(oligos_fhand.name, blastdb, 'blastn', blast_fhand.name, params)
    return TabularBlastParser(open(blast_fhand.name), blast_format)


def look_for_matching_segments(seq_fhand, oligos_fhand):
    'It looks for the oligos in the given sequence files'
    blasts = _do_blast(oligos_fhand, seq_fhand)
    min_identity = 87.0
    min_len = 17
    ignore_elongation_shorter = 3
    filters = [{'kind': 'min_length', 'min_num_residues': min_len,
                'length_in_query': True, 'filter_match_parts': True},
               {'kind': 'score_threshold', 'score_key': 'identity',
               'min_score': min_identity}]
    blasts = filter_alignments(blasts, config=filters)

    # Which are the regions covered in each sequence?
    for blast in blasts:
        read = blast['query']
        for match in blast['matches']:
            oligo = match['subject']
            elongate_match_parts_till_global(match['match_parts'],
                                             query_length=read['length'],
                                             subject_length=oligo['length'])

        match_parts = [m['match_parts'] for m in blast['matches']]
        match_parts = [item for sublist in match_parts for item in sublist]
        # Any of the match_parts has been elongated?
        elongated_match = False
        for m_p in match_parts:
            if ELONGATED in m_p and m_p[ELONGATED] > ignore_elongation_shorter:
                elongated_match = True

        segments = covered_segments(match_parts)
        yield read, segments, elongated_match


def split_by_mate_linker(seqrec, (segments, is_partial)):
    'It splits the seqs using segments'

    if not segments:
        return [seqrec]

    elongated_match = is_partial
    if len(segments) == 1:
        segment_start = segments[0][0]
        segment_end = segments[0][1]
        if segment_start == 0:
            new_seqrec = seqrec[segment_end + 1:]
            new_seqrec.id = seqrec.id + '.fn'
            return [new_seqrec]
        elif segment_end == len(seqrec) - 1:
            new_seqrec = seqrec[:segment_start]
            new_seqrec.id = seqrec.id + '.fn'
            return [new_seqrec]
        else:
            new_seqrec1 = seqrec[:segment_start]
            new_seqrec2 = seqrec[segment_end + 1:]
            id_ = seqrec.id
            if elongated_match:
                id_ = seqrec.id + '_pl'
                new_seqrec1.id = id_ + '.part1'
                new_seqrec2.id = id_ + '.part2'
            else:
                new_seqrec1.id = id_ + '.f'
                new_seqrec2.id = id_ + '.r'
            return [new_seqrec1, new_seqrec2]
    else:
        seqrecords = []
        counter = 1
        seq_start = 0
        for segment_start, segment_end in segments:
            if segment_start == 0:
                continue
            seqrecord = seqrec[seq_start:segment_start]
            seqrecord.id = seqrec.id + '_mlc.part{0:d}'.format(counter)
            seqrecords.append(seqrecord)
            counter += 1
            seq_start = segment_end + 1
        else:
            if segment_end != len(seqrec) + 1:
                seqrecord = seqrec[segment_end + 1:]
                seqrecord.id = seqrec.id + '_mlc.part{0:d}'.format(counter)
                seqrecords.append(seqrecord)
        return seqrecords


def split_mates(seq_fhands, out_fhand, linkers=None):
    'It splits the input sequences with the provided linkers.'

    out_seq_fmt = guess_format(seq_fhands[0])

    if linkers is None:
        linkers = LINKERS

    def get_next_segment(linker_segments):
        'It returns the next segments available'
        try:
            segments = linker_segments.next()
        except StopIteration:
            segments = None
        return segments

    linkers_fhand = NamedTemporaryFile(prefix='linkers.', suffix='.fasta')
    write_seqrecords(linkers_fhand, LINKERS, file_format='fasta')

    for seq_fhand in seq_fhands:
        linker_segments = look_for_matching_segments(seq_fhand, linkers_fhand)
        segments = get_next_segment(linker_segments)
        for seqrec in read_seqrecords([seq_fhand]):
            if segments is not None and segments[0]['name'] == seqrec.id:
                split_seqs = split_by_mate_linker(seqrec, segments[1:])
                segments = get_next_segment(linker_segments)
            else:
                split_seqs = [seqrec]
            write_seqrecords(out_fhand, split_seqs, out_seq_fmt)
