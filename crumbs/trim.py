'''
Created on 2012 eka 26

@author: peio
'''
import random
import itertools


def _get_uppercase_segments(string):
    '''It detects the unmasked regions of a sequence

    It returns a list of (start, end) tuples'''
    start = 0
    for is_upper, group in itertools.groupby(string, lambda x: x.isupper()):
        group = list(group)
        end = start + len(group) - 1
        if is_upper:
            yield start, end
        start = end + 1


def _get_longest_segment(segments):
    'it returns the longest segment'
    longest = None
    longest_size = None
    for segment in segments:
        size = segment[1] - segment[0]
        if longest is None:
            longest = [segment]
            longest_size = size
        elif size > longest_size:
            longest = [segment]
            longest_size = size
        elif size == longest_size:
            longest.append(segment)

    if longest is None:
        return None
    elif len(longest) == 1:
        return longest[0]
    else:
        return random.choice(longest)


def trim_lowercased_seqs(seqrecords):
    'It trims the masked seqments of the seqrecords'
    for seqrecord in seqrecords:
        seq = str(seqrecord.seq)
        unmasked_segments = _get_uppercase_segments(seq)
        segment = _get_longest_segment(unmasked_segments)
        if segment is not None:
            yield seqrecord[segment[0]:segment[1] + 1]
