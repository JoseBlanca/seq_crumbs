
import random
from operator import itemgetter
from collections import Counter

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from  vcf import Reader
from vcf_crumbs.statistics import get_snpcaller_name, get_call_data, ADS

FILTER_ALLELES_GT = None    # it could be an integer, i.e. 1 to keep only the
                            # 2 parental alleles in a F2

REFERENCE = 'Reference'
ALLELE_FREQS = 'Allele freqs'
REF_SAMPLES = 'Ref samples'


def _get_alleles(call, filter_alleles_gt):
    if not call.called:
        return []

    alleles = list((map(int, call.gt_alleles)))

    if filter_alleles_gt is not None:
        alleles = [allele for allele in alleles if allele <= filter_alleles_gt]

    return alleles


def _get_genotype(call, filter_alleles_gt, snp_caller):
    alleles = _get_sorted_alleles(call, filter_alleles_gt)
    if not alleles:
        return None
    n_alleles = set(alleles)

    # how do we resume both alleles in one number?
    if n_alleles == 1:
        genotype = alleles[0]
    allele_depths = get_call_data(call, snp_caller)[ADS].allele_depths

    al0_depth = allele_depths.get(alleles[0], 0)
    al1_depth = allele_depths.get(alleles[1], 0)
    if alleles > 2 or abs(alleles[1] - alleles[0]) != 1:
        # most abundant allele
        if al0_depth > al1_depth:
            genotype = alleles[0]
        elif al1_depth > al0_depth:
            genotype = alleles[1]
        else:
            genotype = random.choice(alleles)
    else:
        genotype = (max(alleles) - min(alleles)) * (al0_depth + al1_depth) / 2 + min(alleles)

    return genotype


def plot_haplotypes(vcf_fpath, plot_fpath, genotype_mode=REFERENCE,
                    filter_alleles_gt=FILTER_ALLELES_GT):
    reader = Reader(filename=vcf_fpath)

    snp_caller = get_snpcaller_name(reader)
    limite marcadores o skip
    no se ven los heterocigotos

    # collect data
    nan = float('nan')
    genotypes = None
    samples = []
    for snp in reader:
        if genotypes is None:
            genotypes = {}
            for call in snp.samples:
                sample = call.sample
                genotypes[sample] = []
                samples.append(sample)

        snp_alleles = {}
        al_depths = Counter()
        for call in snp.samples:
            alleles = _get_alleles(call, filter_alleles_gt=filter_alleles_gt)
            if alleles:
                alleles = Counter(alleles)
                call_depth = get_call_data(call, snp_caller)[ADS].allele_depths
                for allele in alleles:
                    al_depths[allele] += call_depth[allele]
            else:
                alleles = {}
            sample = call.sample
            snp_alleles[sample] = alleles

        # we map the alleles to the new numbers according to the genetic mode
        if genotype_mode == ALLELE_FREQS:
            depths = reversed(sorted(al_depths.iteritems(), key=itemgetter(1)))
            al_map = {depth[0]: i for i, depth in enumerate(depths)}
        else:
            al_map = None
        if al_map:
            if al_depths:
                al_depths = {al_map[al]: depth for al, depth in al_depths.iteritems()}
            if alleles:
                new_alleles = {}
                for sample, geno in snp_alleles.iteritems():
                    geno = {al_map[al]: count for al, count in geno.iteritems()}
                    new_alleles[sample] = geno
                snp_alleles = new_alleles

        # We resume both alleles in one number
        for sample, geno in snp_alleles.iteritems():
            n_alleles = len(geno)
            if n_alleles == 0:
                geno = nan
            elif n_alleles == 1:
                geno = list(geno.keys())[0]
            elif n_alleles == 2:
                alleles = geno.keys()
                if abs(alleles[0] - alleles[1]) == 1:
                    geno = min(alleles) + (geno[alleles[0]] + geno[alleles[0]]) / 2
                else:
                    geno = random.choice(alleles)
            else:
                geno = random.choice(alleles)
            genotypes[sample].append(geno)

    # draw
    n_samples = len(samples)
    xsize = len(genotypes[sample]) / 100
    if xsize >= 100:
        xsize = 100
    ysize = n_samples * 2
    if ysize >= 100:
        ysize = 100
    print xsize, ysize
    figure_size = (xsize, ysize)

    fig = Figure(figsize=figure_size)

    for index, sample in enumerate(samples):
        axes = fig.add_subplot(n_samples, 1, index)
        axes.set_title(sample)
        y_data = genotypes[sample]
        axes.plot([i + 1 for i in range(len(y_data))], y_data, marker='o',
                  linestyle='None', markersize=1.0, markeredgewidth=0,
                  markerfacecolor='red')
        ylim = axes.get_ylim()
        ylim = ylim[0] - 0.1, ylim[1] + 0.1
        axes.set_ylim(ylim)
        axes.tick_params(axis='x', bottom='off', top='off', which='both',
                         labelbottom='off')
        axes.tick_params(axis='y', left='on', right='off', labelleft='off')
        axes.set_ylabel(sample)

    canvas = FigureCanvasAgg(fig)
    fhand = open(plot_fpath, 'w')
    canvas.print_figure(fhand, dpi=300)
    fhand.flush()
