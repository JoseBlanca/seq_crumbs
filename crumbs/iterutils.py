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
from itertools import izip_longest, islice
import cPickle as pickle
from tempfile import NamedTemporaryFile
from multiprocessing import Pool as pool

from toolz.itertoolz.core import merge_sorted


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


def group_in_packets_fill_last(iterable, packet_size, fillvalue=None):
    'ABCDE -> (A, B), (C, D), (E, None)'
    # It is faster than group_in_packets
    iterables = [iter(iterable)] * packet_size
    kwargs = {'fillvalue': fillvalue}
    return izip_longest(*iterables, **kwargs)


def group_in_packets(iterable, packet_size):
    'ABCDE -> (A, B), (C, D), (E,)'
    iterable = iter(iterable)
    while True:
        chunk = tuple(islice(iterable, packet_size))
        if not chunk:
            break
        yield chunk


def flat_zip_longest(iter1, iter2, fillvalue=None):
    '''It yields items alternatively from both iterators.

    It won't return the items equal to the fillvalue.
    '''
    for item1, item2 in izip_longest(iter1, iter2, fillvalue=fillvalue):
        if item1 != fillvalue:
            yield item1
        if item2 != fillvalue:
            yield item2


def rolling_window(iterator, window, step=1):
    'It yields lists of items with a window number of elements'
    try:
        length_ = len(iterator)
    except TypeError:
        length_ = None
    if length_ is None:
        return _rolling_window_iter(iterator, window, step)
    else:
        return _rolling_window_serie(iterator, window, length_, step)


def _rolling_window_serie(serie, window, length_, step):
    '''It yields lists of items with a window number of elements'''
    return (serie[i:i + window] for i in range(0, length_ - window + 1, step))


def _rolling_window_iter(iterator, window, step):
    '''It yields lists of items with a window number of elements giving
     an iterator'''
    items = []
    for item in iterator:
        if len(items) >= window:
            yield items
            items = items[step:]
        items.append(item)
    else:
        if len(items) >= window:
            yield items


def _pickle_items(items, temp_dir):
    fhand = NamedTemporaryFile(suffix='.pickle', dir=temp_dir)
    for item in items:
        fhand.write(pickle.dumps(item))
        fhand.write('\n\n')
    fhand.flush()
    return fhand


def _unpickle_items(fhand):
    str_item = ''
    for line in fhand:
        if line == '\n':
            yield pickle.loads(str_item)
            str_item = ''
        else:
            str_item += line
    if str_item:
        yield pickle.loads(str_item)
    fhand.close()


def unique(items, key=None):
    '''It yields the unique items.

    The items must be sorted. It only compares contiguous items.
    '''
    prev_item = None
    for item in items:
        if prev_item == None:
            duplicated = False
        else:
            if key:
                duplicated = key(item) == key(prev_item)
            else:
                duplicated = item == prev_item
        if not duplicated:
            yield item
        prev_item = item


def sorted_items(items, key=None, max_items_in_memory=None,
                        temp_dir=None):
    if max_items_in_memory:
        grouped_items = group_in_packets(items, max_items_in_memory)
    else:
        grouped_items = [items]

    # sort and write to disk all groups but the last one
    write_to_disk = True if max_items_in_memory else False
    sorted_groups = []
    for group in grouped_items:
        sorted_group = sorted(group, key=key)
        del group
        if write_to_disk:
            group_fhand = _pickle_items(sorted_group, temp_dir=temp_dir)
            sorted_groups.append(_unpickle_items(open(group_fhand.name)))
        else:
            sorted_groups.append(sorted_group)
        del sorted_group
    if len(sorted_groups) > 1:
        if key is None:
            sorted_items = merge_sorted(*sorted_groups)
        else:
            sorted_items = merge_sorted(*sorted_groups, key=key)
    else:
        sorted_items = sorted_groups[0]
    return sorted_items
