
from collections import namedtuple
from itertools import islice

# Missing docstring
# pylint: disable=C0111

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
