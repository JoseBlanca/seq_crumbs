
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from crumbs.plot import get_fig_and_canvas
from crumbs.vcf.snv import VCFReader

FILTER_ALLELES_GT = None    # it could be an integer, i.e. 1 to keep only the
                            # 2 parental alleles in a F2

REFERENCE = 'Reference'
ALLELE_FREQS = 'Allele freqs'
REF_SAMPLES = 'Ref samples'


def _get_alleles(call, filter_alleles_gt):
    alleles = call.int_alleles

    if filter_alleles_gt is not None:
        alleles = [allele for allele in alleles if allele <= filter_alleles_gt]

    return alleles


def _flatten_data(x_data, y_data):
    new_x_data = []
    new_y_data = []
    for x, ys in zip(x_data, y_data):
        for y in ys:
            new_x_data.append(x)
            new_y_data.append(y)
    return new_x_data, new_y_data


def plot_haplotypes(vcf_fhand, plot_fhand, genotype_mode=REFERENCE,
                    filter_alleles_gt=FILTER_ALLELES_GT):
    reader = VCFReader(vcf_fhand)

    # collect data
    genotypes = None
    samples = []
    for snv in reader.parse_snvs():
        if genotypes is None:
            genotypes = {}
            for call in snv.calls:
                sample = call.sample
                genotypes[sample] = []
                samples.append(sample)

        for call in snv.calls:
            alleles = _get_alleles(call, filter_alleles_gt=filter_alleles_gt)
            genotypes[call.sample].append(alleles)

    # draw
    n_samples = len(samples)
    xsize = len(genotypes[sample]) / 100
    if xsize >= 100:
        xsize = 100
    if xsize <= 8:
        xsize = 8
    ysize = n_samples * 2
    if ysize >= 100:
        ysize = 100
    # print xsize, ysize
    figure_size = (xsize, ysize)

    fig = Figure(figsize=figure_size)

    for index, sample in enumerate(samples):
        axes = fig.add_subplot(n_samples, 1, index)
        axes.set_title(sample)
        y_data = genotypes[sample]
        x_data = [i + 1 for i in range(len(y_data))]
        x_data, y_data = _flatten_data(x_data, y_data)

        axes.plot(x_data, y_data, marker='o',
                  linestyle='None', markersize=3.0, markeredgewidth=0,
                  markerfacecolor='red')
        ylim = axes.get_ylim()
        ylim = ylim[0] - 0.1, ylim[1] + 0.1
        axes.set_ylim(ylim)
        axes.tick_params(axis='x', bottom='off', top='off', which='both',
                         labelbottom='off')
        axes.tick_params(axis='y', left='on', right='off', labelleft='off')
        axes.set_ylabel(sample)

    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(plot_fhand, dpi=300)
    plot_fhand.flush()


def window_data_by_chrom(iterator):

    x_vals = []
    y_vals = {}
    prev_chrom = None
    for item in iterator:
        chrom = item['chrom']
        start = item['start']
        end = item['end']
        item['values']

        if prev_chrom is None:
            x_vals = []
            y_vals = {}
        elif prev_chrom != chrom:
            yield prev_chrom, x_vals, y_vals
            x_vals = []
            y_vals = {}

        x_vals.append((start + end) / 2)
        for key, value in item['values'].items():
            if key not in y_vals:
                y_vals[key] = []
            y_vals[key].append(value)
        prev_chrom = chrom
    else:
        if x_vals:
            yield chrom, x_vals, y_vals


def plot_in_genome(iterator, out_base, labels):
    xsize = 100
    ysize = 10

    for chrom, x_vals, y_vals in window_data_by_chrom(iterator):
        fig, canvas = get_fig_and_canvas(figsize=(xsize, ysize))

        num_axes = len(y_vals)
        for index, plotname in enumerate(labels.keys()):
            plot_label = labels[plotname]
            axes = fig.add_subplot(num_axes, 1, index + 1)
            axes.plot(x_vals, y_vals[plotname])
            title = plot_label['title'] + ' for chromosome {}'.format(chrom)
            axes.set_title(title)
            axes.set_ylabel(plot_label['ylabel'])

        axes.set_xlabel('Position in chromosome')

        fhand = open(out_base + '.' + str(chrom) + '.png', 'w')

        canvas.print_figure(fhand, dpi=300)
        fhand.flush()
