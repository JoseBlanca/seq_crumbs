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

import random


def sample(iterator, length, sample_size):
    'It makes a sample for the given iterator'

    if not 0 <= sample_size <= length:
        raise ValueError("sample larger than population")

    if sample_size > length / 2:
        num_items_to_select = length - sample_size
        invert = True
    else:
        num_items_to_select = sample_size
        invert = False

    selected = set(random.randint(0, length - 1) for num in range(num_items_to_select))
    selected_add = selected.add
    while len(selected) < num_items_to_select:
        selected_add(random.randint(0, length - 1))

    for index, item in enumerate(iterator):
        item_in_selected = index in selected
        if item_in_selected and not invert or not item_in_selected and invert:
            yield item


def length(iterator):
    'It counts the items in an iterator. It consumes the iterator'
    count = 0
    for item in iterator:
        count += 1
    return count


class group_in_packets(object):
    '''It groups an iterable into packets of equal number of elements

    [1, 2, 3, 4, 5] -> [[1,2], [3, 4] [5]]
    '''
    def __init__(self, iterable, packet_size):
        'It inits the class'
        self._packet_size = packet_size - 1
        self._iterable = iter(iterable)
        self._current_count = 0
        self._current_item = None

    def __iter__(self):
        'Part of the iterator interface'
        return self

    def next(self):
        'It returns a packet'
        self._current_item = next(self._iterable)    # Exit on StopIteration
        self._current_count = 0
        return self._grouper()

    def _grouper(self):
        'It yields a packet_size number of items'
        while True:
            if self._current_count < self._packet_size:
                yield self._current_item    # Exit on StopIteration
                self._current_item = next(self._iterable)
                self._current_count += 1
            elif self._current_count == self._packet_size:
                yield self._current_item    # Exit on StopIteration
                self._current_count += 1
            else:
                raise StopIteration
