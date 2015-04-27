#!/usr/bin/env python
import argparse

import glob

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 10})
import matplotlib.pyplot as plt

def mad(arr):
    """
    Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

def nmad(arr):
    return 1.4826*mad(arr)

def add_stat_legend(x):
    textstr = ''
    textstr += '$\mathrm{N}=%d$\n' % len(x)
    textstr += '$\mathrm{mean}=%.2f$\n' % np.nanmean(x)
    textstr += '$\mathrm{median}=%.2f$\n' % np.nanmedian(x)
    textstr += '$\mathrm{nmad}=%.2f$' % nmad(x)
    props = dict(boxstyle='round', facecolor='white')
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, 
        va='top', ha='right', bbox=props)

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", action="store_true",
        help="print verbose output")
    parser.add_argument("-o", "--output", type=str, default=None,
        help="output file")
    parser.add_argument("-i", "--input", type=str, default=None,
        help="required input file")
    parser.add_argument("--mask-outliers", action="store_true",
        help="mask outliers in samples")
    args = parser.parse_args()

    filenames = glob.glob(args.input)

    nfiles = len(filenames)
    nparams = 3
    summary = np.empty((nfiles, nparams))
    for i,filename in enumerate(filenames):
        # data columns:
        # x0 a1 a2
        data = np.loadtxt(filename, ndmin=2)
        summary[i] = np.median(data, axis=0)

    outlier_thres = 20

    # save results summary plot
    fig = plt.figure(figsize=(12,12))

    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (1,1))
    ax3 = plt.subplot2grid((3,3), (2,2))

    ax4 = plt.subplot2grid((3,3), (1,0))
    ax5 = plt.subplot2grid((3,3), (2,0))
    ax6 = plt.subplot2grid((3,3), (2,1))

    # plot the fit parameter distributions
    def plot_param_dist(x, binspec, color, xlabel):
        xmin, xmax, nxbins = [element for tupl in binspec for element in tupl]

        print 'Stats for %s' % xlabel
        print 'min, max: %.4g, %.4g' % (xmin, xmax)

        if args.mask_outliers:
            mask = ~is_outlier(x, thresh=outlier_thres)
            print 'number of masked entries: %d' % len(np.flatnonzero(~mask))
            print 'list of masked entries and values: '
            for index in np.flatnonzero(~mask):
                filename = filenames[index]
                print '-'.join(filename.split('.')[0].split('-')[-2:]), x[index]
            x = x[mask]
            print 'min, max (masked): %f, %f' % (np.min(x), np.max(x))
            print

        plt.hist(x, bins=np.linspace(xmin, xmax, nxbins+1), facecolor=color, alpha=.5, histtype='stepfilled')
        plt.xlabel(xlabel)
        plt.xlim([xmin, xmax])
        plt.grid()
        add_stat_legend(x)

    plimits = ((3000, 6500), (-.3,3.2), (-.8,.6))
    pbins = ((50,), (50,), (50,))

    plt.sca(ax1)
    plot_param_dist(summary[:,0], (plimits[0],pbins[0]), 'blue', r'$x_0$')
    plt.sca(ax2)
    plot_param_dist(summary[:,1], (plimits[1],pbins[1]), 'green', r'$a_1$')
    plt.sca(ax3)
    plot_param_dist(summary[:,2], (plimits[2],pbins[2]), 'red', r'$a_2$')

    # plot the fit parameter distributions
    def plot_param_scatter(x, y, xlim, ylim, xlabel, ylabel):
        if args.mask_outliers:
            mask = ~is_outlier(x, thresh=outlier_thres) & ~is_outlier(y, thresh=outlier_thres) 
            x = x[mask]
            y = y[mask]
        plt.plot(x, y, '+', ms=5)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.grid()
        # calculate correlation coefficient
        corr =  np.corrcoef(x, y)
        rho = corr[0,1]
        # add text box
        textstr = ''
        textstr += '$\mathrm{N}=%d$\n' % len(x)
        textstr += r'$\rho=%.2f$' % rho
        props = dict(boxstyle='round', facecolor='white')
        plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, va='top', ha='right', bbox=props)

    plt.sca(ax4)
    plot_param_scatter(summary[:,0], summary[:,1], plimits[0], plimits[1], r'$x_0$', r'$a_1$')
    plt.sca(ax5)
    plot_param_scatter(summary[:,0], summary[:,2], plimits[0], plimits[2], r'$x_0$', r'$a_2$')
    plt.sca(ax6)
    plot_param_scatter(summary[:,1], summary[:,2], plimits[1], plimits[2], r'$a_1$', r'$a_2$')

    plt.tight_layout()
    fig.savefig(args.output)

if __name__ == '__main__':
    main()
