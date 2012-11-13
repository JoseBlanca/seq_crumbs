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
from itertools import izip_longest

from crumbs.exceptions import error_quality_disagree, MalformedFile


def sample(iterator, sample_size):
    'It makes a sample from the given iterator'
    # This implementation holds the sampled items in memory
    # Example of the algorithm seen in:
    # http://nedbatchelder.com/blog/201208/selecting_randomly_from_an_unknown_
    # sequence.html
    # See also:
    # http://stackoverflow.com/questions/12128948/python-random-lines-from-
    # subfolders/

    sample_ = []
    for index, elem in enumerate(iterator):
        if len(sample_) < sample_size:
            sample_.append(elem)
        else:
            if random.randint(0, index) < sample_size:
                # Ned Bactchelder proposed an algorithm with random order
                # sample_[random.randint(0, sample_size - 1)] = elem
                # we prefer to keep the order, so we always add at the end
                sample_.pop(random.randint(0, sample_size - 1))
                sample_.append(elem)
    return sample_


def sample_2(iterator, iter_length, sample_size):
    'It makes a sample from the given iterator'
    # This implementation will use less memory when the number of sampled items
    # is quite high.
    # It requires to know the number of items beforehand

    if not 0 <= sample_size <= iter_length:
        raise ValueError("sample larger than population")

    if sample_size > iter_length / 2:
        num_items_to_select = iter_length - sample_size
        invert = True
    else:
        num_items_to_select = sample_size
        invert = False

    selected = set(random.randint(0, iter_length - 1)
                                           for n in range(num_items_to_select))
    selected_add = selected.add
    while len(selected) < num_items_to_select:
        selected_add(random.randint(0, iter_length - 1))

    for index, item in enumerate(iterator):
        item_in_selected = index in selected
        if item_in_selected and not invert or not item_in_selected and invert:
            yield item


def length(iterator):
    'It counts the items in an iterator. It consumes the iterator'
    count = 0
    # pylint: disable=W0612
    for item in iterator:
        count += 1
    return count


class group_in_packets(object):
    '''It groups an iterable into packets of equal number of elements

    [1, 2, 3, 4, 5] -> [[1,2], [3, 4] [5]]
    '''
    def __init__(self, iterable, packet_size):
        'It inits the class'
        self._packet_size = packet_size
        self._iterable = iter(iterable)

    def __iter__(self):
        'Part of the iterator interface'
        return self

    def next(self):
        'It returns a packet'
        # We are returning lists and not generators (as before) because
        # otherwise we would not be able to use multiproces.Pool
        packet = []
        count = 0
        size = self._packet_size
        try:
            for item in self._iterable:
                packet.append(item)
                count += 1
                if count >= size:
                    break
        except ValueError as error:
            if error_quality_disagree(error):
                raise MalformedFile(str(error))
            raise
        if not packet:
            raise StopIteration
        return packet


def flat_zip_longest(iter1, iter2, fillvalue=None):
    '''It yields items alternatively from both iterators.

    It won't return the items equal to the fillvalue.
    '''
    for item1, item2 in izip_longest(iter1, iter2, fillvalue=fillvalue):
        if item1 != fillvalue:
            yield item1
        if item2 != fillvalue:
            yield item2


def rolling_window(iterator, window):
    'It yields lists of items with a window number of elements'
    try:
        length_ = len(iterator)
    except TypeError:
        length_ = None
    if length_ is None:
        return _rolling_window_iter(iterator, window)
    else:
        return _rolling_window_serie(iterator, window, length_)


def _rolling_window_serie(serie, window, length_):
    '''It yields lists of items with a window number of elements'''
    return (serie[i:i + window] for i in range(length_ - window + 1))


def _rolling_window_iter(iterator, window):
    '''It yields lists of items with a window number of elements giving
     an iterator'''
    items = []
    for item in iterator:
        if len(items) >= window:
            yield items
            items = items[1:]
        items.append(item)
    else:
        if len(items) >= window:
            yield items
