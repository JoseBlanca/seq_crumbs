'''
Created on 2013 urr 8

@author: peio
'''
import os

from matplotlib.figure import Figure
from matplotlib import colors, cm
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

FIGURE_SIZE = (15.0, 11.0)  # inche
COLORMAPS = ['Blues', 'Reds', 'Greens', 'Accent', 'Winter', 'Bone', 'Binary']
COLORS = ['blue', 'red', 'green', 'gray']
MARKERGLHYPS = ['^', 's', 'o', '.', 'D', 'h', 's']
MARKER_SIZE2 = 100
MARKER_SIZE3 = 5.0


def _get_canvas_and_axes(figure_size=FIGURE_SIZE, left=0.1, right=0.9, top=0.9,
                         bottom=0.1):
    'It returns a matplotlib canvas and axes instance'
    fig = Figure(figsize=figure_size)
    canvas = FigureCanvas(fig)
    axes = fig.add_subplot(111)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)

    return canvas, axes


def _get_format_from_fname(fname):
    'It returns the extension as the format for the file'
    return os.path.splitext(fname)[1][1:]


def draw_histogram(values, fhand, bins=10, range_=None, **kwargs):
    'It draws a histogram of a pandas Series into a file'
    canvas, axes = _get_canvas_and_axes()
    axes.hist(values, bins=bins, range=range_)
    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)

    canvas.print_figure(fhand, format=_get_format_from_fname(fhand.name))
    fhand.flush()


def draw_scatter(groups, fhand, plot_lines=False, **kwargs):
    # groups is a list of x,y and genotype values
    canvas, axes = _get_canvas_and_axes()
    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)

    for index, group in enumerate(groups):
        x_vals = group['x']
        y_vals = group['y']
        sct_kwargs = {}
        sct_kwargs['marker'] = MARKERGLHYPS[index]
        # TODO. Review the API for the colors.
        # The following cases should be posible for color/value:
        #    - a tuple RGB
        #    - nothing. It has to chose one color by default
        #    - a string: 'blue', 'green', etc.
        #    - a list of numbers (intensities) and ColorMap
        # A possible option could be:
        # group['color] for anything that matplotlib could digest for a single
        # color: RGB tuple or string.
        # group['color_intensities'] for the values of the color. In this case
        # a color map should be given or a default one should be used.
        # In this is used 'color' and 'color_intensities' should be
        # incompatible
        # What about plot_lines? That should be incompatible with
        # color_intensities
        color = group.get('color', COLORS[index])
        if 'value' in group:
            # value is a list with color intensities
            sct_kwargs['c'] = group['value']
            sct_kwargs['norm'] = colors.Normalize()
            sct_kwargs['cmap'] = cm.get_cmap(COLORMAPS[index])
            sct_kwargs['s'] = 50
            sct_kwargs['alpha'] = 0.5
            sct_kwargs['edgecolor'] = 'white'
        else:
            sct_kwargs['c'] = color

        if plot_lines:
            sct_kwargs['ms'] = MARKER_SIZE3
            sct_kwargs['mec'] = color
            sct_kwargs['mfc'] = color
            axes.plot(x_vals, y_vals, **sct_kwargs)
        else:
            sct_kwargs['s'] = MARKER_SIZE2
            axes.scatter(x_vals, y_vals, **sct_kwargs)

    canvas.print_figure(fhand, format=_get_format_from_fname(fhand.name))
    fhand.flush()
