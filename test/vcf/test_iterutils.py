
import unittest

from crumbs.vcf.iterutils import (generate_windows, PeekableIterator,
                                  RandomAccessIterator)

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

        res = [(0, 10), (10, 20), (20, 30)]
        assert list(generate_windows(size=10, step=10, end=30)) == res


class PeekableIterTest(unittest.TestCase):
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


class RandomAccessTest(unittest.TestCase):
    def test_iter_access(self):
        seq1 = range(100)
        seq2 = RandomAccessIterator(iter(seq1), 11)
        assert list(seq1) == list(seq2)

    def test_next_items(self):
        seq1 = range(10)
        seq2 = RandomAccessIterator(iter(seq1), 7)
        assert seq1[0] == seq1[0]
        assert seq2[0:3] == seq1[0:3]
        assert seq2[:3] == seq1[:3]
        assert seq2[2] == seq1[2]

        first = seq2.next()
        assert first == 0
        assert seq2[0:4] == seq1[0:4]

        item = seq2.next()
        assert item == 1
        assert seq2[0:5] == seq1[0:5]

        item = seq2.next()
        assert item == 2
        assert seq2[0:6] == seq1[0:6]

        item = seq2.next()
        assert item == 3
        assert seq2[0:7] == seq1[0:7]

        try:
            seq2[0:8]
            self.fail('IndexError expexted')
        except IndexError:
            pass

        item = seq2.next()
        assert item == 4
        assert seq2[1:8] == seq1[1:8]
        assert seq2[1] == seq1[1]
        try:
            seq2[0]
            self.fail('IndexError expexted')
        except IndexError:
            pass

        item = seq2.next()
        assert item == 5
        assert seq2[2:9] == seq1[2:9]

        try:
            seq2[1:8]
            self.fail('IndexError expexted')
        except IndexError:
            pass

        item = seq2.next()
        assert item == 6
        assert seq2[3:10] == seq1[3:10]

        item = seq2.next()
        assert item == 7
        assert seq2[3:10] == seq1[3:10]
        item = seq2.next()
        assert item == 8
        item = seq2.next()
        assert item == 9
        assert seq2[3:10] == seq1[3:10]

    def test_short_iter(self):
        seq1 = range(3)
        seq2 = RandomAccessIterator(iter(seq1), 7)
        assert seq1[0] == seq1[0]
        assert seq2[0:3] == seq1[0:3]
        assert seq2[:3] == seq1[:3]

        assert seq1[0] == seq1[0]
        assert seq2[0:3] == seq1[0:3]
        assert seq2[:3] == seq1[:3]

        first = seq2.next()
        assert first == 0
        assert seq2[0:4] == seq1[0:4]

        item = seq2.next()
        assert item == 1
        assert seq2[0:4] == seq1[0:4]

        item = seq2.next()
        assert item == 2
        assert seq2[0:4] == seq1[0:4]

        try:
            item = seq2.next()
            self.fail('StopIteration expected')
        except StopIteration:
            pass

        seq1 = range(2)
        seq2 = RandomAccessIterator(iter(seq1), 7)
        assert seq1[0] == seq1[0]
        assert seq2[0:3] == seq1[0:3]
        assert seq2[:3] == seq1[:3]

    def test_empty_iter(self):
        seq1 = []
        seq2 = RandomAccessIterator(iter(seq1), 7)
        assert seq2[0:3] == seq1[0:3]

        try:
            seq2.next()
            self.fail('StopIteration expected')
        except StopIteration:
            pass

if __name__ == '__main__':
    # import sys; sys.argv = ['', 'RandomAccessTest.test_prev_items']
    unittest.main()
