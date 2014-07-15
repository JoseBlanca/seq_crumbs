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

BAR = 'bar'
LINE = 'line'


def _guess_output_for_matplotlib(fhand):
    'Given an fhand it guesses if we need png or svg'
    output = None
    if fhand is not None:
        output = splitext(fhand.name)[-1].strip('.')
    if not output:
        output = 'png'
    return output


def get_fig_and_canvas(num_rows=1, num_cols=1, figsize=None):
    if figsize is None:
        height = 5.0 * num_rows
        width = 7.5 * num_cols
        if height > 320.0:
            height = 320.0
        figsize = (width, height)
    try:
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    except NameError:
        msg = 'Matplotlib module is required to draw graphical histograms'
        raise OptionalRequirementError(msg)
    return fig, canvas


def get_canvas_and_axes(figure_size=FIGURE_SIZE, left=0.1, right=0.9, top=0.9,
                         bottom=0.1, plot_type=111):
    'It returns a matplotlib canvas and axes instance'
    try:
        fig = Figure(figsize=FIGURE_SIZE)
        canvas = FigureCanvas(fig)
    except NameError:
        msg = 'Matplotlib module is required to draw graphical histograms'
        raise OptionalRequirementError(msg)

    axes = fig.add_subplot(plot_type)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)

    return canvas, axes


def draw_histogram_in_axes(counts, bin_limits, kind=BAR, axes=None, title=None,
                           xlabel=None, ylabel=None, distrib_label=None,
                           linestyle=None):

    if axes is None:
        canvas, axes = get_canvas_and_axes(figure_size=FIGURE_SIZE, left=0.1,
                                           right=0.9, top=0.9, bottom=0.1)
    else:
        canvas = None

    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)
    if title:
        axes.set_title(title)

    if kind == BAR:
        xvalues = range(len(counts))
        axes.bar(xvalues, counts)

        # the x axis label
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

        # we don't want to clutter the plot
        num_of_xlabels = 15
        step = int(len(counts) / float(num_of_xlabels))
        step = 1 if step == 0 else step
        xticks_pos = xticks_pos[::step]
        xticks_labels = xticks_labels[::step]
        axes.set_xticks(xticks_pos)
        axes.set_xticklabels(xticks_labels)
    elif LINE:
        kwargs = {}
        if distrib_label is not None:
            kwargs['label'] = distrib_label
        if linestyle is not None:
            kwargs['linestyle'] = linestyle

        x_values = []
        for index, i in enumerate(bin_limits):
            try:
                i2 = bin_limits[index + 1]
            except IndexError:
                break
            x_values.append((i + i2) / 2.0)
        y_values = counts
        axes.plot(x_values, y_values, **kwargs)
    return axes, canvas


def draw_histogram_in_fhand(counts, bin_limits, title=None, xlabel=None,
                            ylabel=None, fhand=None, kind=BAR):
    'It draws an histogram and if the fhand is given it saves it'
    plot_format = _guess_output_for_matplotlib(fhand)
    canvas, axes = get_canvas_and_axes(figure_size=FIGURE_SIZE, left=0.1,
                                       right=0.9, top=0.9, bottom=0.1)
    draw_histogram_in_axes(counts, bin_limits, axes=axes, title=title,
                           xlabel=xlabel, ylabel=ylabel, kind=kind)

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def draw_histograms(counters, fhand, distrib_labels=None, num_cols=2,
                    plots_per_chart=3, xlabel=None, ylabel=None, titles=None,
                    kind=LINE, xmax=None, xmin=None, linestyles=None):
    if plots_per_chart > 1 and kind == BAR:
        raise ValueError('if kind is BAR only one plot per chart is allowed')

    plot_format = _guess_output_for_matplotlib(fhand)
    num_plots, mod = divmod(len(counters), plots_per_chart)
    if mod != 0:
        num_plots += 1

    num_rows, mod = divmod(num_plots, num_cols)
    if mod != 0:
        num_rows += 1
    fig, canvas = get_fig_and_canvas(num_rows=num_rows, num_cols=num_cols)

    counter_index = 0
    for plot_num in range(1, num_plots + 1):
        axes = fig.add_subplot(num_rows, num_cols, plot_num)
        for i in range(plots_per_chart):
            try:
                counter = counters[counter_index]
                if distrib_labels is None:
                    distrib_label = None
                else:
                    distrib_label = distrib_labels[counter_index]
                if linestyles is None:
                    linestyle = None
                else:
                    linestyle = linestyles[counter_index]
            except IndexError:
                break
            title = titles[counter_index] if titles else None
            try:
                distrib = counter.calculate_distribution(max_=xmax,
                                                         min_=xmin)
            except RuntimeError:
                axes.set_title(title + ' (NO DATA)')
                counter_index += 1
                continue
            except AttributeError as error:
                # if distributions is None
                err_msg = "'NoneType' object has no attribute "
                err_msg += "'calculate_distribution'"
                if err_msg in error:
                    axes.set_title(title + ' (NO DATA)')
                    counter_index += 1
                    continue
                raise
            print distrib_label
            title = titles[counter_index] if titles else None
            draw_histogram_in_axes(distrib['counts'], distrib['bin_limits'],
                                   kind=kind, axes=axes, ylabel=ylabel,
                                   distrib_label=distrib_label, xlabel=xlabel,
                                   title=title, linestyle=linestyle)
            counter_index += 1

        if distrib_labels is not None:
            axes.legend()

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def build_histogram(values, fhand, bins=10, range_=None, stacked=False,
                    color=None, label=None, log=False, **kwargs):
    'It draws a histogram of a pandas Series into a file'
    canvas, axes = get_canvas_and_axes()
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
    # groups is a list of x,y and color_intensity
    canvas, axes = get_canvas_and_axes()
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
        if 'color_intensity' in group:
            # value is a list with color intensities
            sct_kwargs['c'] = group['color_intensity']
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


def lines_to_plot(series, axes):
    for line in series:
        kwargs = {}
        if 'label' in line:
            kwargs['label'] = line['label']
        axes.plot(line['x'], line['y'], **kwargs)
        axes.set_title('title')
    axes.legend()


def draw_lines(series, fhand):
    ''''It draws a line plot with the given series.
    Series is a list of dictionaries, each dictionary at lest needs to contain
    x and y as lists'''
    canvas, axes = get_canvas_and_axes()
    plot_format = _guess_output_for_matplotlib(fhand)
    lines_to_plot(series, axes)

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def draw_density_plot(xs, ys, fhand=None, n_bins=40, canvas=None, axes=None,
                      range_=None, **kwargs):
    plot_format = _guess_output_for_matplotlib(fhand)
    if canvas is None and axes is None:
        canvas, axes = get_canvas_and_axes()
        using_given_axes = False
    elif canvas is not None and axes is not None and fhand is None:
        using_given_axes = True
    else:
        msg = 'If an axes is given the canvas is also required and no fhand'
        raise NotImplementedError(msg)

    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)

    # TODO check with norm=LogNorm() and/or normed=True
    if range_ is None:
        axes.hist2d(xs, ys, bins=n_bins)
    else:
        axes.hist2d(xs, ys, bins=n_bins, range=range_)

    if not using_given_axes:
        canvas.print_figure(fhand, format=plot_format)
        fhand.flush()


def draw_int_boxplot(boxplot, fhand=None, axes=None, title=None,
                     xlabel=None, ylabel=None):

    if axes is None:
        canvas, axes = get_canvas_and_axes()
        using_given_axes = False
    elif fhand is None:
        using_given_axes = True
    else:
        msg = 'If an axes is not given the fhand also required'
        raise NotImplementedError(msg)

    bar_width = 0.8

    x_vals = sorted(boxplot.counts.keys())

    # for the following line to work we would have to create all points in
    # memory, this is why we have reimplemented the boxplot
    # axes.boxplot([list(boxplot.counts[x_val].elements()) for x_val in x_vals]

    # we don't want more than 50 bars
    max_n_boxes = 40
    n_boxes = len(x_vals) - 1
    xpurge = n_boxes // max_n_boxes + 1

    xticks_lables = []
    xticks = []
    for index, x_val in enumerate(x_vals):
        if index % xpurge != 0:
            continue

        xticks_lables.append(str(x_val))
        # we add 0.5 to have the x ticks in  the
        # middle of the bar
        xticks.append(index + 0.5)

        intcounter = boxplot.counts[x_val]
        if intcounter.count < 5:
            # at least for values are required to calculate the qualties
            continue
        quart = intcounter.quartiles
        min_ = intcounter.min
        max_ = intcounter.max
        axes.bar(index + 0.5, quart[2] - quart[0], bottom=quart[0],
                 width=bar_width, facecolor='none', align='center')
        # median
        axes.plot([index + 0.5 - bar_width / 2, index + 0.5 + bar_width / 2],
                  [quart[1], quart[1]], color='red')
        # max
        axes.plot([index + 0.5 - bar_width / 4, index + 0.5 + bar_width / 4],
                  [max_, max_], color='black')
        axes.plot([index + 0.5, index + 0.5], [quart[0], min_], '--',
                  color='black')
        axes.plot([index + 0.5 - bar_width / 4, index + 0.5 + bar_width / 4],
                  [min_, min_], color='black')
        axes.plot([index + 0.5, index + 0.5], [quart[2], max_], '--',
                  color='black')

    max_x = len(x_vals) - 1
    axes.set_xlim(0, max_x + 1)

    axes.set_xticks(xticks)
    axes.set_xticklabels(xticks_lables)

    if title:
        axes.set_title(title)
    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)

    if not using_given_axes:
        canvas.print_figure(fhand)
        fhand.flush()
