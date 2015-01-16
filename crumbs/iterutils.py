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
from itertools import izip_longest, islice, tee, izip
import cPickle as pickle
from tempfile import NamedTemporaryFile
import sqlite3
from collections import namedtuple

from crumbs.utils.optional_modules import merge_sorted
from crumbs.exceptions import SampleSizeError


class _ListLikeDb(object):
    def __init__(self):
        self._db_fhand = NamedTemporaryFile(suffix='.sqlite.db')
        self._conn = sqlite3.connect(self._db_fhand.name)
        create = 'CREATE TABLE items(idx INTEGER PRIMARY KEY, item BLOB);'
        self._conn.execute(create)
        self._conn.commit()
        self._last_item_returned = None

    def __len__(self):
        conn = self._conn
        cursor = conn.execute('SELECT COUNT(*) FROM items;')
        return cursor.fetchone()[0]

    def append(self, item):
        conn = self._conn
        insert = "INSERT INTO items(item) VALUES (?)"
        conn.execute(insert, (pickle.dumps(item),))
        conn.commit()

    def __setitem__(self, key, item):
        update = "UPDATE items SET item=? WHERE idx=?"
        conn = self._conn
        conn.execute(update, (pickle.dumps(item), key))
        conn.commit()

    def __iter__(self):
        return self

    def next(self):
        if self._last_item_returned is None:
            idx = 1
        else:
            idx = self._last_item_returned + 1
        select = 'SELECT item FROM items WHERE idx=?'
        conn = self._conn
        cursor = conn.execute(select, (idx, ))
        item = cursor.fetchone()
        if item is None:
            raise StopIteration()

        item = str(item[0])

        item = pickle.loads(str(item))
        self._last_item_returned = idx
        return item


def sample(iterator, sample_size, in_disk=False):
    '''It makes a sample from the given iterator.

    It does not keep the order.
    Since it does not know before hand the size of the iterator it has to
    keep a buffer as large as the sample size in memory (default) or in disk.
    '''
    # This implementation holds the sampled items in memory
    # Example of the algorithm seen in:
    # http://nedbatchelder.com/blog/201208/selecting_randomly_from_an_unknown_
    # sequence.html
    # See also:
    # http://stackoverflow.com/questions/12128948/python-random-lines-from-
    # subfolders/

    if in_disk:
        sample_ = _ListLikeDb()
    else:
        sample_ = []
    too_big_sample = True
    for index, elem in enumerate(iterator):
        if len(sample_) < sample_size:
            sample_.append(elem)
        else:
            too_big_sample = False
            if random.randint(0, index) < sample_size:
                sample_[random.randint(0, sample_size - 1)] = elem
    if too_big_sample:
        raise SampleSizeError('Sample larger than population')
    return iter(sample_)


def sample_low_mem(iterator, iter_length, sample_size):
    'It makes a sample from the given iterator'
    # This implementation will use less memory when the number of sampled items
    # is quite high.
    # It requires to know the number of items beforehand

    if sample_size <= 0:
        raise SampleSizeError('No items to sample')
    elif sample_size > iter_length:
        raise SampleSizeError('Sample larger than population')

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


def _pickle_items(items, tempdir):
    fhand = NamedTemporaryFile(suffix='.pickle', dir=tempdir)
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
    prev_key = None
    for item in items:
        key_for_item = item if key is None else key(item)
        if prev_key is None:
            duplicated = False
        else:
            duplicated = key_for_item == prev_key
        if not duplicated:
            yield item
        prev_key = key_for_item
    # An alternative implementation would be:
    # return (first(groups[1]) for groups in groupby(items, key))
    # But it is a little bit slower


def sorted_items(items, key=None, max_items_in_memory=None, tempdir=None):
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
            group_fhand = _pickle_items(sorted_group, tempdir=tempdir)
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


def unique_unordered(items, key=None):
    '''It yields the unique items.'''
    unique_keys = set()
    prev_size = len(unique_keys)
    for item in items:
        key_for_item = item if key is None else key(item)
        unique_keys.add(key_for_item)
        current_size = len(unique_keys)
        if current_size > prev_size:
            yield item
            prev_size = current_size


# Taken from a itertools stdlib recipe
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


Window = namedtuple('Window', ['start', 'end'])


def generate_windows(size, step=None, start=0, end=None):
    if step is None:
        step = size

    win_start = None
    while True:
        if win_start is None:
            win_start = start
        else:
            win_start += step
        win_end = win_start + size
        if end:
            if win_end > end:
                break
        yield Window(win_start, win_end)


class PeekableIterator(object):
    def __init__(self, iterable):
        self._stream = iterable
        self._buffer = []

    def __iter__(self):
        return self

    def next(self):
        if self._buffer:
            return self._buffer.pop()
        return self._stream.next()

    def peek(self):
        try:
            item = self._stream.next()
        except StopIteration:
            raise
        self._buffer.append(item)
        return item


class RandomAccessIterator(object):
    def __init__(self, iterable, rnd_access_win):
        self._stream = iterable
        if not rnd_access_win % 2:
            msg = 'rnd_access_win should be odd'
            raise ValueError(msg)
        self._rnd_access_win = rnd_access_win
        self._buff = []
        self._curr_item_in_buff = None
        self._buffer_pos = 0
        self._stream_consumed = False
        self._fill_buffer()

    def _fill_buffer(self):
        half_win = (self._rnd_access_win - 1) // 2
        self._buff = list(islice(self._stream, half_win))
        if len(self._buff) < half_win:
            self._stream_consumed = True

    def __iter__(self):
        return self

    def next(self):
        in_buffer = self._buff
        try:
            stream_next = self._stream.next()
        except StopIteration:
            self._stream_consumed = True

        if not self._stream_consumed:
            in_buffer.append(stream_next)

        if self._curr_item_in_buff is None:
            self._curr_item_in_buff = 0
        else:
            self._curr_item_in_buff += 1

        try:
            item_to_yield = in_buffer[self._curr_item_in_buff]
        except IndexError:
            raise StopIteration

        if len(in_buffer) > self._rnd_access_win:
            in_buffer.pop(0)
            self._curr_item_in_buff -= 1
            self._buffer_pos += 1

        return item_to_yield

    def _getitem(self, start):
        start -= self._buffer_pos
        if start < 0:
            raise IndexError('Negative indexes not supported')
        return self._buff[start]

    def __getitem__(self, index):
        if isinstance(index, int):
            start = index
            return self._getitem(start)
        else:
            start = index.start
            stop = index.stop
            step = index.step
        if start is None:
            start = 0

        if start < 0 or stop < 0:
            raise IndexError('Negative indexes not supported')

        buff_pos = self._buffer_pos
        start -= buff_pos
        stop -= buff_pos

        buff = self._buff
        if start < 0:
            raise IndexError('Index outside buffered window')
        if (not self._stream_consumed and
           (stop is None or stop > len(buff))):
            raise IndexError('Index outside buffered window')

        return buff[start:stop:step]
