# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
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

from operator import itemgetter
import random
import copy

from crumbs.utils.tags import START, END


def get_all_segments(segments, seq_len):
    '''Given a set of some non overlapping regions it returns all regions.

    input:  --- ----       ----
    output: ---+----+++++++----

    '''
    segments_ = copy.deepcopy(segments)
    if not segments_:
        return [((0, seq_len - 1), False)]

    all_segments = []
    #If the first segment doesn't start at zero we create a new one
    if segments_[0][0] == 0:
        start = segments_[0][1] + 1
        all_segments.append((segments_.pop(0), True))
    else:
        start = 0

    for loc in segments_:
        all_segments.append(((start, loc[0] - 1), False))
        all_segments.append((loc, True))
        start = loc[1] + 1
    else:
        #if the last segment does not ends at the end of the sequence we add
        #an extra one
        end = seq_len - 1
        if start <= end:
            all_segments.append(((start, end), False))
    return all_segments


def get_complementary_segments(segments, seq_len):
    'Given a set of regions in a seq it returns the complementary ones'
    non_matched = []
    for segment in get_all_segments(segments, seq_len):
        if not segment[1]:
            non_matched.append(segment[0])
    return non_matched


def get_longest_segment(segments):
    'It returns the longest segment'
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


def get_longest_complementary_segment(segments, seq_len):
    'Given a seq and the locations it returns the longest region not covered'

    segments = merge_overlaping_segments(segments)

    #we want the non-matched locations
    segments = get_complementary_segments(segments, seq_len)

    #now we return the longest one
    return get_longest_segment(segments)


def merge_overlaping_segments(segments, merge_segments_closer=1):
    '''Given a list of segments it returns the overlapped segments.

       segments 1  -------        ----->    -----------
       segments 2       ------
       returns     -----------    ------    -----------
       merge_segments_closer is an integer. Segments closer than the given
       number of residues will be merged.
    '''

    # we collect all start and ends
    limits = []
    for segment in segments:
        start = segment[0]
        end = segment[1]
        if start > end:     # a reversed item
            start, end = end, start
        limit_1 = (START, start)
        limit_2 = (END, end)
        limits.append(limit_1)
        limits.append(limit_2)

    # sort by secondary key: start before end
    limits.sort(key=itemgetter(0))
    #sort by location (primary key)
    limits.sort(key=itemgetter(1))

    # merge the ends and start that differ in only one base
    filtered_limits = []
    previous_limit = None
    for limit in limits:
        if previous_limit is None:
            previous_limit = limit
            continue
        if (previous_limit[0] == END and limit[0] == START and
            previous_limit[1] >= limit[1] - merge_segments_closer):
            #These limits cancelled each other
            previous_limit = None
            continue
        filtered_limits.append(previous_limit)
        previous_limit = limit
    else:
        filtered_limits.append(limit)
    limits = filtered_limits

    # now we create the merged hsps
    starts = 0
    segments = []
    for limit in limits:
        if limit[0] == START:
            starts += 1
            if starts == 1:
                segment_start = limit[1]
        elif limit[0] == END:
            starts -= 1
            if starts == 0:
                segment = (segment_start, limit[1])
                segments.append(segment)
    return segments
