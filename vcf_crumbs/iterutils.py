
def generate_windows(size, step=None, start=0):
    if step is None:
        step = size

    win_start = None
    while True:
        if win_start is None:
            win_start = start
        else:
            win_start += step
        win_end = win_start + size
        yield win_start, win_end


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
