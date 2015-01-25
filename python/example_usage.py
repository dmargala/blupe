#!/usr/bin/env python
"""
Demonstrates how to apply throughput corrections to individual BOSS spectrum.

There a few different ways to specify a spectrum. In order of precedence:
    1) using the spec filename
    2) a target id string (must also provide top-level directory containing spec files)
    3) individual plate, mjd, and fiberid (must also provide top-level directory containing spec files)

Usage: 

python example_usage.py --tpcorr tpcorr.hdf5 --spec-file ~/data/boss/v5_7_0/spectra/lite/3615/spec-3615-55445-0008.fits 

python example_usage.py --tpcorr tpcorr.hdf5 --spec-dir ~/data/boss/v5_7_0/spectra/lite --target 3615-55445-8

python example_usage.py --tpcorr tpcorr.hdf5 --spec-dir ~/data/boss/v5_7_0/spectra/lite --plate 3615 --mjd 55445 --fiberid 8

"""
import argparse
import os

import h5py
import numpy as np

from astropy.io import fits
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default=None,
        help='throughput correction filename, required')
    parser.add_argument('--spec-file', type=str, default=None,
        help='full path of a spec file')
    parser.add_argument('--spec-dir', type=str, default=None,
        help='path to directory containing individual spec files')
    parser.add_argument('--target', type=str, default=None,
        help='target id string (plate-mjd-fiberid), must specify --spec-dir')
    parser.add_argument('--plate', type=int, default=None,
        help='plate id, must specify --spec-dir')
    parser.add_argument('--mjd', type=int, default=None,
        help='mjd id, must specify --spec-dir')
    parser.add_argument('--fiberid', type=int, default=None,
        help='fiber id, must specify --spec-dir')
    args = parser.parse_args()

    # Open the throughput correction file
    tpcorr = h5py.File(args.tpcorr, 'r')
    tpcorr_wave = tpcorr['wave'].value

    # Parse arguments: extract plate, mjd, fiberid as integers and the spec filename
    if args.spec_file:
        spec_filename = args.spec_file
        spec_basename = os.path.basename(args.spec_file)
        plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
    elif args.spec_dir:
        if args.target:
            plate, mjd, fiberid = [int(field) for field in args.target.split('-')]
        elif args.plate and args.mjd and args.fiberid:
            plate, mjd, fiberid = (args.plate, args.mjd, args.fiberid)
        else:
            parser.error('Must specify a target using the --target option or the --plate, --mjd, and --fiberid options')
        spec_filename = os.path.join(args.spec_dir, str(plate), 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid))
    else:
        parser.error('Must specify either a spec file using the --spec-file option or path to directory containing spec files and a target.')

    # Load the target's combined spectrum
    spec = fits.open(spec_filename)
    flux = spec[1].data.field('flux')
    ivar = spec[1].data.field('ivar')
    wavelength = np.power(10, spec[1].data.field('loglam'))
    spec.close()

    # Read the target's throughput correction vector
    tpcorr_key = '%s/%s/%s' % (plate, mjd, fiberid)
    correction = tpcorr[tpcorr_key].value

    # Create an interpolated correction function
    correction_interp = interp1d(tpcorr_wave, correction, kind='linear')

    # Sample the interpolated correction using the observation's wavelength grid
    resampled_correction = correction_interp(wavelength)

    # Apply the correction to the observed flux and ivar
    corrected_flux = flux*resampled_correction
    corrected_ivar = ivar/resampled_correction**2

    # Plot the original and corrected spectrum
    fig = plt.figure(figsize=(12,4))

    plt.plot(wavelength, flux, c='red', label='BOSS')
    plt.plot(wavelength, corrected_flux, c='blue', label='Corrected BOSS')

    ymin = min(0, np.percentile(flux,  1))
    ymax = 1.25*np.percentile(flux, 99)

    plt.xlim(min(wavelength), max(wavelength))
    plt.ylim(ymin, ymax)

    plt.title('%04d-%5d-%04d' % (plate, mjd, fiberid))
    plt.xlabel('Observed Wavelength $(\AA)$')
    plt.ylabel('Flux $(10^{-17} erg/cm^2/s/\AA)$')

    plt.legend()
    plt.grid()
    plt.show()

    tpcorr.close()


if __name__ == '__main__':
    main()