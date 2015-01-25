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

from astropy import units as u
from astropy import constants as const

# Units
FNU = u.erg / (u.cm**2 * u.s * u.Hz)
FLAM = u.erg / (u.cm**2 * u.s * u.AA)


def blackbody_nu(in_x, temperature):
    """Calculate blackbody flux per steradian, :math:`B_{\\nu}(T)`.

    .. note::

        Use `numpy.errstate` to suppress Numpy warnings, if desired.

    .. warning::

        Output values might contain ``nan`` and ``inf``.

    Parameters
    ----------
    in_x : number, array-like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Hz.

    temperature : number or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in
        :math:`erg \\; cm^{-2} s^{-1} Hz^{-1} sr^{-1}`.

    Raises
    ------
    ValueError
        Invalid temperature.

    ZeroDivisionError
        Wavelength is zero (when converting to frequency).

    """
    # Convert to units for calculations, also force double precision
    with u.add_enabled_equivalencies(u.spectral() + u.temperature()):
        freq = u.Quantity(in_x, u.Hz, dtype=np.float64)
        temp = u.Quantity(temperature, u.K, dtype=np.float64)

    # Check if input values are physically possible
    if temp < 0:
        raise ValueError('Invalid temperature {0}'.format(temp))
    if np.any(freq <= 0):  # pragma: no cover
        warnings.warn('Input contains invalid wavelength/frequency value(s)',
                      AstropyUserWarning)

    # Calculate blackbody flux
    bb_nu = (2.0 * const.h * freq ** 3 /
             (const.c ** 2 * np.expm1(const.h * freq / (const.k_B * temp))))
    flux = bb_nu.to(FNU, u.spectral_density(freq))

    return flux / u.sr  # Add per steradian to output flux unit

def blackbody_lambda(in_x, temperature):
    """Like :func:`blackbody_nu` but for :math:`B_{\\lambda}(T)`.

    Parameters
    ----------
    in_x : number, array-like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Angstrom.

    temperature : number or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in
        :math:`erg \\; cm^{-2} s^{-1} \\AA^{-1} sr^{-1}`.

    """
    if getattr(in_x, 'unit', None) is None:
        in_x = u.Quantity(in_x, u.AA)

    bb_nu = blackbody_nu(in_x, temperature) * u.sr  # Remove sr for conversion
    flux = bb_nu.to(FLAM, u.spectral_density(in_x))

    return flux / u.sr  # Add per steradian to output flux unit

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

    import scipy.optimize as optimization

    def func(x, a, b):
        return (a*(blackbody_lambda(x*u.AA, b*u.K))).value

    fit_flux = optimization.curve_fit(func, wavelength, flux, p0=[10**-5, 7000], sigma=None)
    fit_corrected_flux = optimization.curve_fit(func, wavelength, corrected_flux, p0=[10**-5, 7000], sigma=None)

    print fit_flux[0]
    print fit_corrected_flux[0]

    # Plot the original and corrected spectrum
    fig = plt.figure(figsize=(12,4))

    plt.plot(wavelength, flux, c='red', label='BOSS')
    plt.plot(wavelength, corrected_flux, c='blue', label='Corrected BOSS')
    plt.plot(wavelength, func(wavelength, *fit_flux[0]), c='black')
    plt.plot(wavelength, func(wavelength, *fit_corrected_flux[0]), c='black')

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