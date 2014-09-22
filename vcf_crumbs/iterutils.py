
from collections import namedtuple

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
        self._rnd_access_win = (rnd_access_win  - 1) / 2
        self._next_buff = []
        self._prev_buff = []
        self._fill_buffers()

    def _fill_buffers(self):
        rnd_access_win = self._rnd_access_win
        next_buff = self._next_buff
        while True:
            if len(next_buff) >= rnd_access_win:
                break
            try:
                nxt_item = self._stream.next()
            except StopIteration:
                break
            next_buff.append(nxt_item)

    def __iter__(self):
        return self

    def next(self):
        next_buffer = self._next_buff
        prev_buff = self._prev_buff
        if not next_buffer:
            raise StopIteration()
        nxt_item = next_buffer.pop(0)
        # fill the buffers
        stream_empty = False
        try:
            nxt_nxt_item = self._stream.next()
        except StopIteration:
            stream_empty = True

        # next buffer is refilled
        if not stream_empty:
            next_buffer.append(nxt_nxt_item)

        # prev buffer is updated
        if len(prev_buff) >= self._rnd_access_win:
            prev_buff.pop(0)
        prev_buff.append(nxt_item)
        return nxt_item

    def next_items(self, start, end=0):
        next_buffer = self._next_buff
        if start == 0:
            msg = 'Returning the current item is not implemented'
            raise NotImplementedError(msg)

        start -= 1
        if start >= len(next_buffer):
            msg = 'Index outside the buffered next window'
            raise IndexError(msg)
        if end == 0:
            return next_buffer[start]
        else:
            if end is not None:
                end -= 1
            return next_buffer[start: end]


    def prev_items(self, start, end=0):
        prev_buff = self._prev_buff
        if start > 0 or end > 0:
            raise ValueError('Only negative indexes allowed')
        if end == 0:
            return prev_buff[start]
        else:
            start += 1
            if start == 0:
                start = None
            return prev_buff[end: start]
