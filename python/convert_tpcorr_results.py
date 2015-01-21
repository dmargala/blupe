#!/usr/bin/env python
import argparse

import glob

import numpy as np

import h5py

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--output", type=str, default=None,
        help="output file")
    parser.add_argument("-i", "--input", type=str, default=None,
        help="required input file")
    args = parser.parse_args()

    filenames = glob.glob(args.input)

    nfiles = len(filenames)

    # the first is the fiberid and the next two columns are xfocal and yfocal positions of the target
    nidtokens = 3
    # the rest are the tabulated throughput correction values
    npoints = 71

    # the throughput correction vectors span the range 3500A to 10500A
    xvalues = np.linspace(3500, 10500, npoints, endpoint=True)

    outfile = h5py.File(args.output, 'w')

    outfile.create_dataset('wave', data=xvalues)

    for i,filename in enumerate(sorted(filenames)):
        plate, mjd = filename.split('.')[0].split('-')[-2:]

        print 'Reading %s...' % filename
        data = np.loadtxt(filename, ndmin=2)

        if plate not in outfile.keys():
            outfile.create_group(plate)

        if mjd not in outfile[plate].keys():
            outfile[plate].create_group(mjd)

        nentries, ntokens = data.shape

        assert ntokens == npoints + nidtokens

        for row in data:
            fiberid, xfocal, yfocal = row[0:nidtokens]

            # skip unplugged fibers
            if int(fiberid) == -1:
                continue

            tpcorr = row[nidtokens:]

            try:
                dset = outfile.create_dataset('/'.join([plate,mjd,str(int(fiberid))]), data=tpcorr, dtype='f4')
            except RuntimeError, e:
                print plate, mjd, fiberid, e
            # dset.attrs['xfocal'] = xfocal
            # dset.attrs['yfocal'] = yfocal


if __name__ == '__main__':
    main()
