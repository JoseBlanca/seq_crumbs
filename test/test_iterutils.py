
import unittest

from vcf_crumbs.iterutils import generate_windows, PeekableIterator

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111


class GenerateWindowsTests(unittest.TestCase):
    def generate_wins(self, size, step, number):
        windows = generate_windows(size=size, step=step)
        wins = []
        for index, win in enumerate(windows):
            if index >= number:
                break
            wins.append(win)
        return wins

    def test_generate_windows(self):
        assert self.generate_wins(10, 10, 3) == [(0, 10), (10, 20), (20, 30)]
        assert self.generate_wins(10, None, 3) == [(0, 10), (10, 20), (20, 30)]
        assert self.generate_wins(10, 5, 3) == [(0, 10), (5, 15), (10, 20)]

    def test_peek(self):
        stream = iter(range(5))
        stream = PeekableIterator(stream)
        assert stream.peek() == 0
        assert stream.next() == 0
        assert stream.next() == 1
        assert stream.peek() == 2
        assert stream.next() == 2
        assert stream.next() == 3
        assert stream.next() == 4
        try:
            stream.peek()
            self.fail('StopIteration expected')
        except StopIteration:
            pass

if __name__ == '__main__':
    # import sys; sys.argv = ['', 'TrimChimericRegions']
    unittest.main()
