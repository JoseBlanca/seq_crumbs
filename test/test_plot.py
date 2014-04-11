'''
Created on 2013 urr 10

@author: peio
'''

import unittest
from tempfile import NamedTemporaryFile
from array import array

from bam_crumbs.plot import draw_histogram, draw_scatter, draw_density_plot
import numpy


class PlotTests(unittest.TestCase):
    def test_do_histogram(self):
        values = [1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 1, 4]
        fhand = NamedTemporaryFile(suffix='.png')
        draw_histogram(values, fhand, bins=100, title='hola')
        #raw_input(fhand.name)

    def test_draw_scatter(self):
        fhand = NamedTemporaryFile(suffix='.png')
        data = [{'y': array('B', [40, 60, 13]),
                 'x': array('B', [20, 45, 45]),
                 'value': array('B', [9, 90, 90])}]

        draw_scatter(data, fhand, xlim=0, ylim=0)
        #raw_input(fhand.name)

    def test_draw_scatter_line(self):
        fhand = NamedTemporaryFile(suffix='.png')
        data = [{'y': array('B', [40, 60, 13]),
                 'x': array('B', [20, 45, 45])}]
        draw_scatter(data, fhand, plot_lines=True)
        #raw_input(fhand.name)

    def test_draw_density_plot(self):
        fhand = NamedTemporaryFile(suffix='.png')
        data = {'y': array('f', numpy.random.normal(40, 5, 40000)),
                 'x': array('f', numpy.random.normal(40, 5, 40000))}
        draw_density_plot(data['x'], data['y'], fhand, n_bins=200)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
