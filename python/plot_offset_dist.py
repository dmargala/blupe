#!/usr/bin/env python
import argparse

import numpy as np

import glob

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 10})
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
        help="output file base name")
    parser.add_argument("-i", "--input", type=str, default=None,
        help="required input file")
    args = parser.parse_args()


    # the first is the fiberid and the next two columns are xfocal and yfocal positions of the target
    nidtokens = 3
    # the rest are the tabulated throughput correction values
    npoints = 71

    # the throughput correction vectors span the range 3500A to 10500A
    xvalues = np.linspace(3500, 10500, npoints, endpoint=True)

    data = np.loadtxt(args.input, ndmin=2)

    offsets = data[:,nidtokens+0::3]

    print 'Read offsets for %d targets ' % len(offsets)

    print offsets.shape

    for i,x in enumerate(xvalues):
        offsets_wave_slice = offsets[:,i]

        fig = plt.figure(figsize=(8,6))
        plt.hist(offsets_wave_slice, bins=50, histtype='stepfilled', alpha=0.5)
        plt.xlabel('Centroid offset (arcseconds)')
        plt.ylabel('Counts')
        plt.title(r'$\lambda = %s$' % x)
        plt.xlim([0, 2])

        add_stat_legend(offsets_wave_slice)

        plt.grid(True)

        fig.savefig(args.output+'-%s.png'%x, bbox_inches='tight')
        fig.close()

if __name__ == '__main__':
    main()
