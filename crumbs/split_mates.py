
# similar software SffToCA
from tempfile import NamedTemporaryFile


from crumbs.blast import get_blast_db, do_blast, generate_tabblast_format
from crumbs.seqio import write_seqrecords
from crumbs.alignment_result import (TabularBlastParser, filter_alignments,
                                     covered_segments, ELONGATED, QUERY,
                                     elongate_match_parts_till_global)
from crumbs.settings import (PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS,
                             LINKERS)


def _do_blast(seq_fhand, oligos):
    'It returns an alignment result with the blast.'

    oligos_fhand = NamedTemporaryFile(prefix='oligo.', suffix='.fasta')
    write_seqrecords(oligos_fhand, oligos, file_format='fasta')

    blastdb = get_blast_db(seq_fhand.name, dbtype='nucl')

    blast_format = ['query', 'subject', 'query_length', 'subject_length',
                    'query_start', 'query_end', 'subject_start',
                    'subject_end', 'expect', 'identity']
    fmt = generate_tabblast_format(blast_format)
    params = {'task': 'blastn-short', 'outfmt': fmt}

    blast_fhand = NamedTemporaryFile(prefix='oligos_vs_reads.',
                                     suffix='.blast')
    do_blast(oligos_fhand.name, blastdb, 'blastn', blast_fhand.name, params)

    return TabularBlastParser(blast_fhand, blast_format)


class _BlastMatcher(object):
    'It matches the given oligos against the reads using Blast'
    def __init__(self, reads_fhand, oligos):
        'It inits the class.'
        self._match_parts = self._look_for_blast_matches(reads_fhand, oligos)
        self._oligos = oligos

    @staticmethod
    def _look_for_blast_matches(seq_fhand, oligos):
        'It looks for the oligos in the given sequence files'
        # we need to keep the blast_fhands, because they're temp files and
        # otherwise they might be removed
        blasts = _do_blast(seq_fhand, oligos)
        min_identity = 87.0
        min_len = 17
        filters = [{'kind': 'min_length', 'min_num_residues': min_len,
                    'length_in_query': False, 'filter_match_parts': True},
                   {'kind': 'score_threshold', 'score_key': 'identity',
                   'min_score': min_identity}]
        blasts = filter_alignments(blasts, config=filters)

        # Which are the regions covered in each sequence?
        indexed_match_parts = {}
        one_oligo = True if len(oligos) == 1 else False
        for blast in blasts:
            oligo = blast['query']
            for match in blast['matches']:
                read = match['subject']
                elongate_match_parts_till_global(match['match_parts'],
                                                 query_length=oligo['length'],
                                                 subject_length=read['length'],
                                                 align_completely=QUERY)

                #match_parts = [m['match_parts'] for m in blast['matches']]
                match_parts = match['match_parts']
                if one_oligo:
                    indexed_match_parts[read['name']] = match_parts
                else:
                    try:
                        indexed_match_parts[read['name']].extend(match_parts)
                    except KeyError:
                        indexed_match_parts[read['name']] = match_parts
        return indexed_match_parts

    def get_matched_segments_for_read(self, read_name):
        'It returns the matched segments for any oligo'
        ignore_elongation_shorter = 3

        try:
            match_parts = self._match_parts[read_name]
            if read_name == 'seq5':
                pass
        except KeyError:
            # There was no match in the blast
            return None

        # Any of the match_parts has been elongated?
        elongated_match = False
        for m_p in match_parts:
            if ELONGATED in m_p and m_p[ELONGATED] > ignore_elongation_shorter:
                elongated_match = True
        segments = covered_segments(match_parts, in_query=False)
        return segments, elongated_match


class MatePairSplitter(object):
    'It splits the input sequences with the provided linkers.'

    def __init__(self, linkers=None):
        'The initiator'
        self.linkers = LINKERS if linkers is None else linkers
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqs):
        'It split a list of sequences with the provided linkers'
        stats = self._stats
        stats[PROCESSED_PACKETS] += 1
        seq_fhand = NamedTemporaryFile(suffix='.fasta')
        write_seqrecords(seq_fhand, seqs, 'fasta')
        matcher = _BlastMatcher(open(seq_fhand.name), self.linkers)
        new_seqs = []
        for seqrec in seqs:
            stats[PROCESSED_SEQS] += 1
            segments = matcher.get_matched_segments_for_read(seqrec.id)
            if segments is not None:
                split_seqs = self._split_by_mate_linker(seqrec, segments)
            else:
                split_seqs = [seqrec]
            for seq in split_seqs:
                new_seqs.append(seq)
                stats[YIELDED_SEQS] += 1
        return new_seqs

    def _split_by_mate_linker(self, seqrec, (segments, is_partial)):
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
