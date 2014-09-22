
import unittest

from vcf_crumbs.iterutils import (generate_windows, PeekableIterator,
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
        seq1 = range(100)
        seq2 = RandomAccessIterator(iter(seq1), 11)
        first = seq2.next()
        assert first == 0
        assert seq2.next_items(1) == 1
        assert seq2.next_items(1, 2) == [1]

        try:
            seq2.next_items(6)
            self.fail('IndexError expected')
        except:
            pass
        assert seq2.next_items(1, 6) == [1, 2, 3, 4, 5]
        assert seq2.next_items(1, 6) == seq2.next_items(1, 7)
        assert seq2.next_items(1, 6) == seq2.next_items(1, None)

        assert list(seq2) == list(range(1, 100))

        try:
            seq2.next_items(1)
            self.fail('IndexError expected')
        except:
            pass

    def test_prev_items(self):
        seq1 = range(100)
        seq2 = RandomAccessIterator(iter(seq1), 11)
        
        assert not seq2.prev_items(-1, None)
        
        first = seq2.next()
        assert first == 0
        assert seq2.prev_items(-1) == 0
        assert seq2.prev_items(-1, -2) == [0]

        try:
            seq2.prev_items(-2)
            self.fail('IndexError expected')
        except:
            pass

        seq2.next()
        seq2.next()
        assert seq2.prev_items(-1, -4) == [0, 1, 2]
        assert seq2.prev_items(-1, -4) == seq2.prev_items(-1, -5)
        assert seq2.prev_items(-1, -4) == seq2.prev_items(-1, None)

        assert list(seq2) == list(range(3, 100))

        assert seq2.prev_items(-1, None) == [95, 96, 97, 98, 99]

    # def test_buffered_items(-5, 5)
    # buffered_items(None, None)

if __name__ == '__main__':
    # import sys; sys.argv = ['', 'RandomAccessTest.test_prev_items']
    unittest.main()
