#!/usr/bin/env python
import argparse

import numpy as np

import glob

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 8})
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

    nfields = 4

    # the throughput correction vectors span the range 3500A to 10500A
    xvalues = np.linspace(3500, 10500, npoints, endpoint=True)

    data = np.loadtxt(args.input, ndmin=2)

    dx = data[:, nidtokens+0::nfields]
    dy = data[:, nidtokens+1::nfields]

    print 'Read offsets for %d targets ' % len(dx)
    print 'Num entries per target: %d' % dx.shape[1]

    edges = np.linspace(-1.5, 1.5,101)

    mywaves = [3500, 4000, 4500, 5500, 7000, 10500]

    fig = plt.figure(figsize=(8,9))

    ax1 = plt.subplot2grid((3,2), (0,0))
    ax2 = plt.subplot2grid((3,2), (1,0))
    ax3 = plt.subplot2grid((3,2), (2,0))

    ax4 = plt.subplot2grid((3,2), (0,1))
    ax5 = plt.subplot2grid((3,2), (1,1))
    ax6 = plt.subplot2grid((3,2), (2,1))

    axs = [ax1,ax2,ax3,ax4,ax5,ax6]

    #for i,x in enumerate(xvalues):
    for axnum, x in enumerate(mywaves):
        plt.sca(axs[axnum])
        i = np.argmax(xvalues == x)
        print i,x
        dx_wave = dx[:,i]
        dy_wave = dy[:,i]

        # estimate the 2D histogram
        H, xedges, yedges = np.histogram2d(dx_wave, dy_wave, bins=edges)

        # H needs to be rotated and flipped
        H = np.rot90(H)
        H = np.flipud(H)      

        # plot 2D histogram using pcolor
        # explore different color map options http://matplotlib.org/users/colormaps.html
        plt.pcolormesh(xedges, yedges, H,
            cmap='gist_ncar_r', norm=mpl.colors.LogNorm(), vmin=10)#, vmax=10**5)
        cbar = plt.colorbar()

        # set axis limits
        plt.xlim([edges[0], edges[-1]])
        plt.ylim([edges[0], edges[-1]])
        plt.gca().set_aspect('equal')

        # label plot
        plt.xlabel('dx (arcseconds)')
        plt.ylabel('dy (arcseconds)')
        plt.text(0.05, 0.95, (r'$\lambda = %s\AA$' % x), transform=plt.gca().transAxes, 
            va='top', ha='left', bbox=dict(boxstyle='round', facecolor='white'))
        plt.grid(True)

        # add legend with stat summary
        #add_stat_legend(np.sqrt(dx_wave**2 + dy_wave**2))

        # add circle indicating fiber size
        fiber_circle = plt.Circle((0, 0), 1, color='k', fill=False)
        plt.gca().add_artist(fiber_circle)

    # save and close figure
    fig.savefig(args.output+'.pdf', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()
