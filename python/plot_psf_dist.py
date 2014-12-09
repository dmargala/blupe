#!/usr/bin/env python
import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
#mpl.rcParams.update({'font.size': 8})
import matplotlib.pyplot as plt

def add_stat_legend(x):
    textstr = '$\mathrm{N}=%d$\n$\mathrm{mean}=%.2f$\n$\mathrm{median}=%.2f$\n$\mathrm{std}=%.2f$' % (
        len(x), np.nanmean(x), np.nanmedian(x), np.nanstd(x))
    props = dict(boxstyle='round', facecolor='white')
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, va='top', ha='right', bbox=props)

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
    args = parser.parse_args()

    if args.verbose:
        print 'Reading file: %s' % args.input

    # data columns:
    # plate mjd mapmjd mapname mean_ha mean_psf_fwhm
    data = np.loadtxt(args.input, usecols=(0,1,2,4,5), unpack=True)

    print data.shape

    psf = data[4][data[4] > 0]

    print len(psf), np.nanmean(psf), np.nanmedian(psf), np.nanstd(psf)

    fig = plt.figure(figsize=(8,6))
    plt.hist(psf, bins=50, histtype='stepfilled', alpha=0.5)
    plt.xlabel('SEEING50 (arcseconds)')
    plt.ylabel('Counts')

    add_stat_legend(psf)

    plt.grid(True)

    fig.savefig(args.output, bbox_inches='tight')

if __name__ == '__main__':
    main()
