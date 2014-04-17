'''
Created on 2014 api 17

@author: peio
'''

from os.path import splitext

try:
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib import colors, cm
except ImportError:
    pass

from crumbs.exceptions import OptionalRequirementError


FIGURE_SIZE = (15.0, 11.0)  # inche
COLORMAPS = ['Blues', 'Reds', 'Greens', 'Accent', 'Winter', 'Bone', 'Binary']
COLORS = ['blue', 'red', 'green', 'gray']
MARKERGLHYPS = ['^', 's', 'o', '.', 'D', 'h', 's']
MARKER_SIZE2 = 100
MARKER_SIZE3 = 5.0


def _guess_output_for_matplotlib(fhand):
    'Given an fhand it guesses if we need png or svg'
    output = None
    if fhand is not None:
        output = splitext(fhand.name)[-1].strip('.')
    if not output:
        output = 'png'
    return output


def _get_canvas_and_axes(figure_size=FIGURE_SIZE, left=0.1, right=0.9, top=0.9,
                         bottom=0.1):
    'It returns a matplotlib canvas and axes instance'
    fig = Figure(figsize=figure_size)
    canvas = FigureCanvas(fig)
    axes = fig.add_subplot(111)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)

    return canvas, axes


def draw_histogram(counts, bin_limits, title=None, xlabel=None,
                   ylabel=None, fhand=None):
    'It draws an histogram and if the fhand is given it saves it'
    plot_format = _guess_output_for_matplotlib(fhand)
    canvas, axes = _get_canvas_and_axes(figure_size=FIGURE_SIZE, left=0.1,
                                        right=0.9, top=0.9, bottom=0.1)

    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)
    if title:
        axes.set_title(title)

    xvalues = range(len(counts))

    axes.bar(xvalues, counts)

    #the x axis label
    xticks_pos = [value + 0.5 for value in xvalues]

    left_val = None
    right_val = None
    xticks_labels = []
    for value in bin_limits:
        right_val = value
        if left_val:
            xticks_label = (left_val + right_val) / 2.0
            if xticks_label >= 10:
                fmt = '%d'
            elif xticks_label >= 0.1 and xticks_label < 10:
                fmt = '%.1f'
            elif xticks_label < 0.1:
                fmt = '%.1e'
            xticks_label = fmt % xticks_label
            xticks_labels.append(xticks_label)
        left_val = right_val

    #we don't want to clutter the plot
    num_of_xlabels = 15
    step = int(len(counts) / float(num_of_xlabels))
    step = 1 if step == 0 else step
    xticks_pos = xticks_pos[::step]
    xticks_labels = xticks_labels[::step]
    axes.set_xticks(xticks_pos)
    axes.set_xticklabels(xticks_labels)

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def build_histogram(values, fhand, bins=10, range_=None, stacked=False,
                   color=None, label=None, log=False, **kwargs):
    'It draws a histogram of a pandas Series into a file'
    canvas, axes = _get_canvas_and_axes()
    plot_format = _guess_output_for_matplotlib(fhand)
    if color is None:
        if label is None:
            axes.hist(values, bins=bins, range=range_, stacked=stacked,
                      log=log)
        else:
            axes.hist(values, bins=bins, range=range_, stacked=stacked,
                      label=label, log=log)
    else:
        axes.hist(values, bins=bins, range=range_, stacked=stacked,
                  label=label, color=color, log=log)
    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)
    if label is not None:
        axes.legend()

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def draw_scatter(groups, fhand, plot_lines=False, **kwargs):
    # groups is a list of x,y and genotype values
    canvas, axes = _get_canvas_and_axes()
    plot_format = _guess_output_for_matplotlib(fhand)

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

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def draw_density_plot(xs, ys, fhand, n_bins=40, canvas=None, axes=None,
                      range_=None, **kwargs):
    plot_format = _guess_output_for_matplotlib(fhand)
    if canvas is None and axes is None:
        canvas, axes = _get_canvas_and_axes()
    elif canvas is not None and axes is not None:
        pass
    else:
        msg = 'If an axes is given the canvas is also required'
        raise NotImplementedError(msg)

    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)

    # TODO check with norm=LogNorm() and/or normed=True
    if range_ is None:
        axes.hist2d(xs, ys, bins=n_bins)
    else:
        axes.hist2d(xs, ys, bins=n_bins, range=range_)

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()
