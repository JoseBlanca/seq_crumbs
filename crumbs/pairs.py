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

import re

from Bio import SeqIO

from crumbs.exceptions import MaxNumReadsInMem, PairDirectionError
from crumbs.utils.tags import FWD, REV


def _parse_pair_direction_and_name(seq):
    'It parses the description field to get the name and the pair direction'
    title = seq.id + ' ' + seq.description
    reg_exps = ['(.+)/(\d+)', '(.+)\s(\d+):.+', '(.+)\.(\w)\s?']

    for reg_exp in reg_exps:
        match = re.match(reg_exp, title)
        if match:
            name, direction = match.groups()
            if direction in ('1', 'f'):
                direction = FWD
            elif direction in ('2', '.r'):
                direction = REV
            else:
                raise PairDirectionError('unknown direction descriptor')
            return name, direction
    raise PairDirectionError('Unable to detect the direction of the seq')


def match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                memory_limit=500000):
    'It matches the seq pairs in an iterator and splits the orphan seqs'
    buf_fwd = {'index': {}, 'items': []}
    buf_rev = {'index': {}, 'items': []}
    for seq in seqs:
        seq_name, direction = _parse_pair_direction_and_name(seq)
        if direction == FWD:
            buf1 = buf_rev
            buf2 = buf_fwd
        else:
            buf1 = buf_fwd
            buf2 = buf_rev

        try:
            matching_seq_index = buf1['index'][seq_name]
        except KeyError:
            matching_seq_index = None

        if matching_seq_index is None:
            #add to buff
            buf2['items'].append(seq)
            buf2['index'][seq_name] = len(buf2['items']) - 1
            #check mem limit
            sum_items = len(buf1['items'] + buf2['items'])
            if memory_limit is not None and sum_items >= memory_limit:
                error_msg = 'There are too many consecutive non matching seqs'
                error_msg += ' in your input. We have reached the memory limit'
                raise MaxNumReadsInMem(error_msg)
        else:
            orphan_seqs = buf1['items'][:matching_seq_index]
            matching_seq = buf1['items'][matching_seq_index]

            # fix buffers
            buf1['items'] = buf1['items'][matching_seq_index + 1:]
            buf1['index'] = {s: i for i, s in enumerate(buf1['items'])}

            #write seqs in file
            SeqIO.write(orphan_seqs, orphan_out_fhand, out_format)
            SeqIO.write([matching_seq, seq], out_fhand, out_format)

    else:
        orphan_seqs = buf1['items'] + buf2['items']
        SeqIO.write(orphan_seqs, orphan_out_fhand, out_format)
