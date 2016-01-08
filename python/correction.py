#!/usr/bin/env python

import argparse

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

# import warnings
# warnings.filterwarnings('error', category=FutureWarning)
# np.seterr(all='raise')

import matplotlib.pyplot as plt

import scipy.interpolate

import bossdata
import specsim

from astropy.table import Table
import astropy.time
import astropy.coordinates
import astropy.units as u
import astropy.units.imperial
import astropy.units.cds
astropy.units.imperial.enable()
# astropy.units.cds.enable()
from astropy.coordinates import Angle

import galsim

import h5py

spAll = bossdata.meta.Database(lite=False)

# ## Optical Distortion Model

# The chromatic distortion is described [here](https://trac.sdss3.org/tracmailman/browser/private/sdss3-infrastructure/all/505.html) with corrected tables appearing [here](https://trac.sdss3.org/tracmailman/browser/private/sdss3-infrastructure/all/946.html). The table of *principal ray heights* are tabulated in [platedesign svn repo](https://trac.sdss.org/browser/repo/sdss/platedesign/trunk/data/sdss/image-heights.txt) and used in [this IDL code](https://trac.sdss3.org/browser/repo/platedesign/trunk/pro/plate/apo_rdistort.pro). Note that the relative distortion number only have one significant digit, so are of limited use.

# https://trac.sdss.org/browser/repo/sdss/platedesign/trunk/data/sdss/plParam.par
optical_distortion_coefs = np.array([
        -0.000137627, -0.00125238, 1.5447e-09, 8.23673e-08, -2.74584e-13, -1.53239e-12, 6.04194e-18,
        1.38033e-17, -2.97064e-23, -3.58767e-23, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

def get_optical_distortion_model(r_max, platescale, npts=100):
    r = np.linspace(0, r_max.to(u.mm).value, npts)*u.mm
    asin_r = np.arcsin((r / platescale).to(u.rad).value) * u.rad
    dr = np.polynomial.polynomial.polyval((asin_r * platescale).to(u.mm).value, optical_distortion_coefs)
    interpolator = scipy.interpolate.interp1d(asin_r.value, dr, copy=True, kind='cubic')
    def model(r):
        # Input r should be a Quantity with length units.
        return interpolator(np.arcsin((r / platescale).to(u.rad).value)) * u.mm
    return model

def plot_optical_distortion(r_max = 325 * u.mm, platescale = 217.7358 * u.mm / u.deg):
    model = get_optical_distortion_model(r_max, platescale)
    r = np.linspace(0, r_max, 100)
    dr = (model(r) / platescale).to(u.arcsec)
    plt.plot(r, dr)
    plt.grid()
    plt.xlabel('Undistorted radius [mm]')
    plt.ylabel('Radial distortion at 5000A [arcsec]')


distortion_wlen = np.array([
    4000., 5300., 5500., 6000., 8000., 10000., 15350., 15950., 16550. ]) * u.Angstrom

distortions_at_5300 = np.array([
    -0.00,  36.26,  72.53, 108.84, 145.18, 181.53, 217.90, 254.29, 290.77, 327.44
]) * u.mm

distortions_relative_to_5300 = np.array([
    [  0.000, -0.002, -0.003, -0.004, -0.005, -0.005, -0.005, -0.004, -0.002,  0.003 ],
    # Added row of zeros for 5300A
    [  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000 ],
    [ -0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.000 ],
    [  0.000,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001, -0.001 ],
    [  0.000,  0.001,  0.003,  0.003,  0.004,  0.004,  0.004,  0.003,  0.002, -0.003 ],
    [  0.000,  0.002,  0.004,  0.005,  0.005,  0.005,  0.005,  0.005,  0.003, -0.004 ],
    [ -0.000,  0.003,  0.006,  0.007,  0.008,  0.008,  0.008,  0.008,  0.004, -0.006 ],
    [ -0.000,  0.003,  0.006,  0.008,  0.008,  0.009,  0.009,  0.008,  0.004, -0.006 ],
    [  0.000,  0.004,  0.006,  0.008,  0.009,  0.009,  0.009,  0.008,  0.004, -0.007 ]
]) * u.mm


def get_chromatic_distortion_model(platescale):
    """
    https://trac.sdss3.org/browser/repo/platedesign/trunk/pro/plate/apo_rdistort.pro
    """
    # Chromatic distortions are tabulated at 10 radii, corresponding to
    # sin(r) = 0', 10', ..., 90' (with arcmins converted to radians).
    sin_r_table = (np.arange(10) * 10 * u.arcmin).to(u.rad)
    # Calculate the corresponding radii in mm.
    r_table = (np.arcsin(sin_r_table.value) * u.rad * platescale).to(u.mm)
    
    # Calculate fractional distortions relative to 5300A
    d5300 = distortions_relative_to_5300 / distortions_at_5300 
    # We are copying the original IDL code here, but setting these to zero might
    # make more sense.
    d5300[:, 0] = d5300[:, 1]
    
    # Calculate additive distortions relative to 5500A.
    assert distortion_wlen[2] == 5500. * u.Angstrom
    d5500 = d5300 - d5300[2]
    
    # Build a linear interpolator in wavelength of additive distortions relative to 5500A.
    wlen_interpolator = scipy.interpolate.interp1d(distortion_wlen, d5500, axis=0,
                                                   kind='linear', copy=True, bounds_error=True)
    
    def model(r, wlen):
        # Clip wavelengths to the tabulated limits.
        wlen = np.clip(wlen, distortion_wlen[0], distortion_wlen[-1])
        # Calculate the additive relative distortion at each wavelength.
        radial_distortion = wlen_interpolator(wlen)
        r_interpolator = scipy.interpolate.interp1d(r_table, radial_distortion, kind='cubic',
                                                    copy=False, bounds_error=True)
        return r * r_interpolator(r)
    
    return model


def plot_chromatic_distortion(wlen=4000*u.Angstrom, r_max = 325 * u.mm, platescale = 217.7358 * u.mm / u.deg):
    model = get_chromatic_distortion_model(platescale)
    r = np.linspace(0, r_max, 100)
    dr = (model(r, wlen.to(u.Angstrom).value) / platescale).to(u.arcsec)
    plt.plot(r, dr)
    plt.grid(True)
    plt.xlabel('Undistorted radius [mm]')
    plt.ylabel('Radial distortion relative to 5500A [arcsec]')


## Atmospheric Refraction Model

class Pointing(object):
    """Represents an observer pointing at a fixed sky location.
    """
    def __init__(self, ra_center, dec_center):
        self.plate_center = astropy.coordinates.ICRS(ra=ra_center, dec=dec_center)
        self.plate_fid = astropy.coordinates.ICRS(ra=ra_center, dec=dec_center + 1.5 * u.deg)
        self.platescale = 217.7358 * u.mm / u.deg
        self.wlen0 = 5500 * u.Angstrom
        self.relhum = 0.2
        self.where = specsim.transform.observatories['APO']
        # self.where = astropy.coordinates.EarthLocation.from_geodetic(lat=32.7797556*u.deg, lon='-105d49m13s', height=2797.*u.m)
        # self.where = astropy.coordinates.EarthLocation.of_site('apo')
#         print self.where.to_geodetic()
        self.distortion_model = get_optical_distortion_model(330.0 * u.mm, self.platescale)
        self.chromatic_model = get_chromatic_distortion_model(self.platescale)

    def transform(self, targets, tai, wlen, temperature, pressure):
        """Transform from sky coordinates to focal plane coordinates.
        
        Args:
            targets: astropy.coordinates.SkyCoord
            tai: float or numpy.ndarray
            wlen: astropy.units.quantity.Quantity
            temperature: astropy.units.quantity.Quantity
            pressure: astropy.units.quantity.Quantity or None to calculate based
                on temperature and elevation
            
        Returns:
            tuple x, y of astropy.units.quantity.Quantity objects giving focal plane
            positions in length units, broadcast over the input args.
        """
        # apo = astropy.coordinates.EarthLocation.from_geodetic(lat=32.7797556*u.deg, lon='-105d50m00s', height=2797.*u.m)

        # Initialize time objects from the input TAI values in MJD seconds.
        when = astropy.time.Time(tai/86400., format='mjd', scale='tai', location=self.where)
        # when = astropy.time.Time(56383.2034742*np.ones_like(tai), format='mjd', scale='tai', location=apo)
                
#         temperature = 5*u.deg_C
#         pressure = 71.890241*u.kPa

        # Calculate the Alt-Az path of the telescope boresight at the plate design wavelength (5500A).
        obs_model0 = specsim.transform.create_observing_model(
                    where=self.where, when=when, wavelength=self.wlen0, pressure=pressure,
                    temperature=temperature, relative_humidity=self.relhum)
        altaz0 = specsim.transform.sky_to_altaz(self.plate_center, obs_model0)
        alt0, az0 = altaz0.alt, altaz0.az

#         print 'Plate Center (alt, az): %r, %r'  % (alt0.value, az0.value)
    
        # Calculate the Alt-Az paths of each target over the input wavelength grid.
        obs_model = specsim.transform.create_observing_model(
            where=self.where, when=when, wavelength=wlen, pressure=pressure,
            temperature=temperature, relative_humidity=self.relhum)
        altaz = specsim.transform.sky_to_altaz(targets, obs_model)

        # Convert each target's Alt-Az into local X, Y focal plane coordinates.
        x, y = specsim.transform.altaz_to_focalplane(
            altaz.alt, altaz.az, altaz0.alt, altaz0.az, self.platescale)
        # Flip y to match the handedness of the XFOCAL, YFOCAL coordinate system.
        y = -y

        # Rotate the focal plane so that +y points towards a point that is offset from
        # the plate center along DEC by +1.5 degrees.
        altaz_fid = specsim.transform.sky_to_altaz(self.plate_fid, obs_model0)
        x_fid, y_fid = specsim.transform.altaz_to_focalplane(
            altaz_fid.alt, altaz_fid.az, altaz0.alt, altaz0.az, self.platescale)
        angle = np.arctan2(x_fid.si, -y_fid.si)
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        x_rot = x * cos_angle - y * sin_angle
        y_rot = x * sin_angle + y * cos_angle

        # Apply radial optical distortions.
        r = np.sqrt(x_rot**2 + y_rot**2)
        distortion = ((r + self.distortion_model(r)) / r).si

        # chromatic_distortion = self.chromatic_model(r*distortion, wlen if wlen.isscalar else wlen[:, np.newaxis])
        # distortion = ((r + chromatic_distortion) / r).si

        x_dist = distortion * x_rot
        y_dist = distortion * y_rot
        
        return x_dist, y_dist, altaz.alt, altaz.az
    
    def hour_angle(self, tai):
        """Convert TAI to the hour angle of this plate's RA"""
        when = astropy.time.Time(tai/86400., format='mjd', scale='tai', location=self.where)
        return when.sidereal_time('apparent') - self.plate_center.ra

## Ideal Guiding Model

class Guider(object):
    """Calculates optimum guider corrections.
    """
    def __init__(self, x0, y0, x, y):
        """
        Find the scale, rotation and offsets that minimizes the residuals between x,y and x0,y0.
        """
        assert x0.shape == y0.shape, 'x0,y0 have different shapes.'
        assert x.shape == y.shape, 'x,y have different shapes.'
        assert len(x.shape) == 2, 'x,y have unexpected shape.'
        assert x.shape[0] == x0.shape[0], 'x,y and x0,y0 have different lengths.'

        nxy, nt = x.shape
        A = np.empty((2 * nxy, 4))
        self.scale = np.empty((nt,))
        self.rotation = np.empty((nt,)) * u.rad
        self.dx = np.empty((nt,)) * u.m
        self.dy = np.empty((nt,)) * u.m
        self.nt = nt
        self.platescale = 217.7358 * u.mm / u.deg
        
        xy0 = np.concatenate([x0.si.value.flat, y0.si.value.flat])
        for it in range(nt):
            xt = x[:, it].si.value.flatten()
            yt = y[:, it].si.value.flatten()
            xyt = np.concatenate([xt, yt])
            for i in range(nxy):
                A[i, :] = (xt[i], -yt[i], 1. ,0.)
                A[nxy + i, :] = (yt[i], xt[i], 0., 1.)
            params, xy_residuals, rank, sing = scipy.linalg.lstsq(A, xy0)
            
            self.scale[it] = np.sqrt(params[0]**2 + params[1]**2)
            self.rotation[it] = np.arctan2(params[1], params[0]) * u.rad
            self.dx[it] = params[2] * u.m
            self.dy[it] = params[3] * u.m
            
        #print np.mean(self.rotation), np.mean(self.scale), np.mean(self.dx), np.mean(self.dy)

        self.scos = self.scale * np.cos(self.rotation)
        self.ssin = self.scale * np.sin(self.rotation)
        
        self.x0 = x0
        self.y0 = y0
        self.x_before = x
        self.y_before = y
        self.x_after, self.y_after = self.correct(x, y)

    def correct(self, x, y):
        """Apply guiding corrections to the focal plane positions x, y.
        """
        assert x.shape == y.shape, 'x,y have different shapes.'
        assert x.shape[-1] == self.nt, 'x,y have unexpected shape.'
        
        return self.scos * x - self.ssin * y + self.dx, self.ssin * x + self.scos * y + self.dy
    
    def plot(self, tai, zoom=3000, field_radius=None, fiber_radius=None, save=None):
        """
        """
        plt.figure(figsize=(12, 12))
        assert len(tai.shape) == 1 and len(tai) == self.nt, 'tai has unexpected shape.'
        
        assert field_radius is not None
        rmax = field_radius.to(u.mm).value            
        plt.xlim(-1.01 * rmax, +1.01 * rmax)
        plt.ylim(-1.01 * rmax, +1.01 * rmax)
        plt.axis('off')
        outline = plt.Circle((0,0), rmax, edgecolor='black', facecolor='none')
        plt.gca().add_artist(outline)
        plt.gca().set_aspect(1.0)
        
        # Draw the nominal guide hole centers.
        plt.scatter(self.x0.to(u.mm).value, self.y0.to(u.mm).value, marker='+', s=100, color='k')

        for i, (x0, y0) in enumerate(zip(self.x0, self.y0)):
            if fiber_radius is not None:
                fiber = plt.Circle(
                    (x0.to(u.mm).value, y0.to(u.mm).value),
                    zoom * fiber_radius.to(u.mm).value, edgecolor='k', facecolor='none', ls='dotted')
                plt.gca().add_artist(fiber)
            plt.plot((x0 + zoom * (self.x_before[i] - x0)).to(u.mm).value,
                     (y0 + zoom * (self.y_before[i] - y0)).to(u.mm).value, 'b-',
                     label=(None if i else 'No Guiding'))
            plt.plot((x0 + zoom * (self.x_after[i] - x0)).to(u.mm).value,
                     (y0 + zoom * (self.y_after[i] - y0)).to(u.mm).value, 'r-',
                    label=(None if i else 'Ideal Guiding'))
        plt.legend(loc='lower left')
        
        plt.tight_layout()
        if save:
            plt.savefig(save)

## Fiber Acceptance Model

class Telescope(object):
    """
    Represents a telescope.
    """
    def __init__(self,diameter=3.80*u.m, obscuration_area_fraction=0.25, 
                 throughput=0.95*0.77, plate_scale=67.40*u.um/u.arcsec):
        self.diameter = diameter
        self.obscuration_area_fraction = obscuration_area_fraction
        self.throughput = throughput
        self.plate_scale = plate_scale
        self.effective_area = np.pi*diameter**2/4.*(1-obscuration_area_fraction)
    def get_optical_psf(self, wavelength):
        #Convert dimensionless lam/D to arcsec units.
        lam_over_diam_arcsec = ((wavelength/self.diameter)*u.rad).to(u.arcsec)
        # Airy requires floats as inputs, not numpy scalars.
        return galsim.Airy(lam_over_diam=float(lam_over_diam_arcsec.value),
            obscuration=float(np.sqrt(self.obscuration_area_fraction)))
    def get_atmospheric_psf(self, wavelength, fwhm5400, gauss=False):
        wlen_ratio = (wavelength/(5400*u.Angstrom)).si
        assert wlen_ratio == wlen_ratio.value,'wavelength has invalid units.'
        fwhm = fwhm5400.to(u.arcsec).value*wlen_ratio**(-0.2)
        # Kolmogorov requires floats as inputs, not numpy scalars.
        if gauss:
            return galsim.Gaussian(fwhm=float(fwhm))
        else:
            return galsim.Kolmogorov(fwhm=float(fwhm))
    def get_psf(self, wavelength, fwhm5400, rms_jitter=0.1*u.arcsec, gauss=False):
        components = [ self.get_atmospheric_psf(wavelength, fwhm5400, gauss=gauss),self.get_optical_psf(wavelength) ]
        # Include a Gaussian pointing jitter, if requested.
        if rms_jitter is not None:
            components.append(galsim.Gaussian(sigma=rms_jitter.to(u.arcsec).value))
        return galsim.Convolve(components)

def calculate_fiber_acceptance(fiber_diameter, psf, sampling=100, max_offset=2):
    """
    Calculate the fiber acceptance fraction versus offset for a specified PSF.
    
    Args:
        fiber_diameter: Diameter of the fiber to use with explicit angular units.
        psf: PSF model to use, assumed to be specified in arcseconds.
        sampling: Sampling to use for the calculation. Higher samplings take longer
            but give more accurate results.
        max_offset: Maximum centroid offset to calculate, as a ratio of the
            fiber diameter.
    
    Returns:
        tuple: Tuple (offset,acceptance) where offset is given as a ratio of fiber
            diameter and acceptance is a fraction from 0-1.
    """
    # Render the PSF to an image with size fiber_diameter by (max_offset+1)*fiber_diameter.
    diam_arcsec = (fiber_diameter.to(u.arcsec)).value
    width = 2*sampling+1
    height = int(np.ceil((max_offset+1)*width))
    image = galsim.Image(width,height,scale=diam_arcsec/width)
    psf.shift(dx=0.,dy=-0.5*diam_arcsec*max_offset).drawImage(image=image)
    # Prepare a boolean mask of pixels falling inside the fiber aperture.
    xy = np.arange(width) - 0.5*(width-1)
    x,y = np.meshgrid(xy,xy)
    mask = (x**2 + y**2 < (0.5*width)**2)
    # Loop over centroid offsets.
    offset = np.arange(height-width+1)/float(width)
    acceptance = np.empty_like(offset)
    for dy in range(height-width):
        acceptance[dy] = np.sum(image.array[dy:dy+width]*mask)
    #return offset,acceptance
    return scipy.interpolate.interp1d(offset, acceptance, kind='linear', copy=True, bounds_error=True)

class AcceptanceModel(object):
    def __init__(self, seeing_fwhm, fiber_diameter=2*u.arcsec, sampling=100, max_offset=0.75):
        """
        Calculate the fiber acceptance fraction versus offset for a specified PSF.

        Args:
            fiber_diameter: Diameter of the fiber to use with explicit angular units.
            psf: PSF model to use, assumed to be specified in arcseconds.
            sampling: Sampling to use for the calculation. Higher samplings take longer
                but give more accurate results.
            max_offset: Maximum centroid offset to calculate, as a ratio of the
                fiber diameter.
        """
        # Render the PSF to an image with size fiber_diameter by (max_offset+1)*fiber_diameter.
        diam_arcsec = (fiber_diameter.to(u.arcsec)).value
        width = 2 * sampling + 1
        pixel_scale = diam_arcsec / width
        height = int(np.ceil((max_offset + 1) * width))
        sigma_arcsec = seeing_fwhm.to(u.arcsec).value / (2. * np.sqrt(2. * np.log(2)))
        xx = pixel_scale * (np.arange(width) - 0.5 * (width - 1))
        yy = pixel_scale * (np.arange(height) - 0.5 * (width - 1))
        x, y = np.meshgrid(xx, yy, sparse=True, copy=False)
        image = np.exp(-(x**2 + y**2) / (2. * sigma_arcsec**2)) * pixel_scale**2 / (2 * np.pi * sigma_arcsec**2)
        # Prepare a boolean mask of pixels falling inside the fiber aperture.
        mask = xx**2 + xx[:, np.newaxis]**2 < (0.5 * diam_arcsec)**2
        # Loop over centroid offsets.
        offset = np.arange(height - width + 1) * pixel_scale
        acceptance = np.zeros_like(offset)
        for dy in range(height - width + 1):
            acceptance[dy] = np.sum(image[dy:dy + width] * mask)
        ##plt.plot(offset, acceptance)
        self.interpolation = scipy.interpolate.interp1d(
            offset, acceptance, kind='linear', copy=True, bounds_error=True)

    def __call__(self, centroid_offsets):
        return self.interpolation(centroid_offsets.to(u.arcsec).value)


# ## Throughput Corrections

def normalize_angle(angle):
    while angle <= -180:
        angle += 360
    while angle > 180:
        angle -= 360
    return angle

def calculate_target_offsets(plate, mjd, guide_wlen=5400*u.Angstrom, std_wlen=5400.*u.Angstrom,
                             offset_wlen=4000*u.Angstrom, steps_per_exposure=5, wlen_grid_steps=15):
    print 'Calculating corrections for {} observed on MJD {}'.format(plate, mjd)
    
    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()

    # Get the list of exposures used in this observation's coadd from a spec lite file.
    spec_name = finder.get_spec_path(plate, mjd, fiber=1, lite=True)
    spec_file = bossdata.spec.SpecFile(mirror.get(spec_name))

    # Read the first b1 raw science exposure to find this plate's plug map.
    raw = spec_file.get_raw_image(0, 'blue', finder=finder, mirror=mirror)
    plug_map = raw.read_plug_map()

    # Look up the plate design pointing from the raw header and convert to
    # an index A,B,C,... -> 0,1,2,...
    pointing_label = raw.header['POINTING'].strip()
    pointing_index = ord(pointing_label) - ord('A')
    
    # Initialize a pointing object for this plate's sky location.
    ra_center = float(plug_map['raCen']) * u.deg
    dec_center = float(plug_map['decCen']) * u.deg
    print 'Plate center is RA={:.3f}, DEC={:.3f} for {}-{}'.format(ra_center, dec_center, plate, pointing_label)
    pointing = Pointing(ra_center, dec_center)
    
    # Find the nominal observing temperature and time that this plate's holes are drilled for.
    apo = specsim.transform.observatories['APO']
    # apo = astropy.coordinates.EarthLocation.of_site('apo')

    design_temp = float(plug_map['temp'])*u.deg_C
    design_pressure = None # Calculate based on elevation and temperature
    design_ha = float(plug_map['ha'].split()[pointing_index]) * u.deg
    midnight = astropy.time.Time(mjd, format='mjd', scale='tai', location=apo)
    design_time = specsim.transform.adjust_time_to_hour_angle(midnight, ra_center, design_ha)
    design_tai = design_time.mjd * 86400.
    print 'Holes drilled for T={:.1f} and HA={:.1f} (TAI={:.1f})'.format(design_temp, design_ha, design_tai)
    
    # Find this plate's guide stars.
    plugging = plug_map['PLUGMAPOBJ']
    guide_fibers = plugging['holeType'] == 'GUIDE'
    guide_ra, guide_dec = plugging['ra'][guide_fibers], plugging['dec'][guide_fibers]
    guide_targets = astropy.coordinates.ICRS(guide_ra * u.deg, guide_dec * u.deg)
    
    # Calculate the nominal guide fiber positions.
    guide_x0, guide_y0, _, _ = pointing.transform(guide_targets, design_tai, guide_wlen,
                                                  design_temp, design_pressure)
    
    # Find this plate's offset fibers. We have to use spAll for this since the plug map does
    # not record the design wavelengths.
    offset_fibers = spAll.select_all(
        where='PLATE={} and MJD={} and LAMBDA_EFF={}'
        .format(plate, mjd, offset_wlen.to(u.Angstrom).value),
        what='FIBER,OBJTYPE,PLUG_RA,PLUG_DEC,XFOCAL,YFOCAL')
    offset_xfocal = offset_fibers['XFOCAL'] * u.mm
    offset_yfocal = offset_fibers['YFOCAL'] * u.mm
    offset_targets = astropy.coordinates.ICRS(
        ra=offset_fibers['PLUG_RA'] * u.deg,
        dec=offset_fibers['PLUG_DEC'] * u.deg)
    print 'Plate has {:d} guide fibers and {:d} offset targets.'.format(
        len(guide_targets), np.count_nonzero(offset_targets))

    # Calculate the nominal science fiber positions. These will not match XFOCAL, YFOCAL
    # exactly since we do not exactly replicate the IDL transforms, but they should be
    # close (within ~0.2 arcsec) and we only use offsets calculated consistently with
    # transform() in the following.
    offset_x0, offset_y0, offset_alt, offset_az = pointing.transform(
        offset_targets, design_tai, offset_wlen, design_temp, design_pressure)
    
    # Calculate where the offset target fibers would have been positioned if they were
    # designed for the same wavelength as the standard stars.
    offset_x0_std, offset_y0_std, _, _ = pointing.transform(
        offset_targets, design_tai, std_wlen, design_temp, design_pressure)
    
    # Initialize the wavelength grid to use for calculating corrections.
    wlen_grid = np.linspace(3500., 10500., wlen_grid_steps)[:, np.newaxis] * u.Angstrom
    
    # Initialize guided target centroid list
    guided_centroids = []

    # Initialize seeing array
    obs_seeing = []

    # Precompute the conversion from inches of Hg to kPa.
    pconv = (1 * u.cds.mmHg * u.imperial.inch / u.mm).to(u.kPa).value
    
    # Loop over exposures
    for exp_index in range(spec_file.num_exposures):

        # Open the b1 cframe for this exposure, to access its metadata.
        b1_cframe_name = finder.get_plate_path(
            plate, spec_file.get_exposure_name(exp_index, 'blue', 'spFrame'))
        b1_cframe = bossdata.plate.FrameFile(mirror.get(b1_cframe_name))
        exp_id = b1_cframe.exposure_id

        # Lookup this exposure's observing time, seeing, and temperature.
        tai_beg, tai_end = b1_cframe.header['TAI-BEG'], b1_cframe.header['TAI-END']
        tai_mid = 0.5 * (tai_beg + tai_end)
        
        obs_ha = normalize_angle(pointing.hour_angle(tai_mid).to(u.deg).value)*u.deg
        
        seeing = b1_cframe.header['SEEING50'] * u.arcsec
        try:
            temperature = b1_cframe.header['AIRTEMP'] * u.deg_C
        except ValueError,e:
            temperature = design_temp
        try:
            pressure = b1_cframe.header['PRESSURE'] * pconv * u.kPa
        except ValueError,e:
            pressure = design_pressure

        # print 'Exp[{:02d}] #{:08d} seeing {:.3f}, T={:+5.1f}, P={:.1f}, TAI {:.1f} ({:+7.3f} days, HA {:+.1f})'.format(
        #     exp_index, exp_id, seeing, temperature, pressure, tai_mid, (tai_mid - design_tai)/86400., obs_ha.to(u.deg))
        
        obs_seeing.append(seeing)
        
        # Create time steps covering this exposure.
        tai_steps = np.linspace(tai_beg, tai_end, steps_per_exposure)

        # Calculate the actual guide target positions on the focal plane without any guiding.
        guide_x, guide_y, _, _ = pointing.transform(
            guide_targets[:, np.newaxis], tai_steps, guide_wlen, temperature, pressure)
        
        # Solve for the optimal guider corrections.
        guider = Guider(guide_x0, guide_y0, guide_x, guide_y)
        guider.plot(tai_steps, field_radius=340 * u.mm, zoom=5000.,
                    fiber_radius=0.1 * u.arcsec * pointing.platescale,
                    save='plot/guide-{}-{}-{}.png'.format(plate, mjd, exp_index))
        plt.close()

        # Calculate the offset target paths on the focal plane without any guiding, for the
        # actual observing conditions.
        offset_x, offset_y, _, _ = pointing.transform(
            offset_targets[:, np.newaxis, np.newaxis], tai_steps, wlen_grid, temperature, pressure)
        
        # Apply guiding corrections to estimate the actual offset target paths during the exposure.
        guided_centroids.append(guider.correct(offset_x, offset_y))

    results = {
        'fibers':offset_fibers['FIBER'], 'num_exposures':spec_file.num_exposures,
        'obs_seeing':obs_seeing, 'wlen_grid':wlen_grid, 'steps_per_exposure':steps_per_exposure,
        'guided_centroids':guided_centroids,
        'offset_x0':offset_x0, 'offset_y0':offset_y0,
        'offset_x0_std':offset_x0_std, 'offset_y0_std':offset_y0_std
        }

    return results

def calculate_corrections(offsets, seeing_wlen=5400.*u.Angstrom, platescale=217.7358*u.mm/u.deg):

    # Initialize telescope model
    sdss_25m = Telescope(diameter=2.5*u.m, obscuration_area_fraction=0.27, plate_scale=platescale)

    # Extract precomputed centroid positions
    offset_x0 = offsets['offset_x0']
    offset_y0 = offsets['offset_y0']
    offset_x0_std = offsets['offset_x0_std']
    offset_y0_std = offsets['offset_y0_std']
    guided_centroids = offsets['guided_centroids']

    correction = np.empty(
        (offsets['num_exposures'], len(offsets['fibers']), len(offsets['wlen_grid']), offsets['steps_per_exposure']),
        dtype=float)

    obs_seeing = offsets['obs_seeing']

    for exp_index in range(offsets['num_exposures']):

        seeing = obs_seeing[exp_index]

        psf = sdss_25m.get_atmospheric_psf(seeing_wlen, seeing, gauss=False)
        acceptance_model = calculate_fiber_acceptance(2*u.arcsec, psf)

        guided_x, guided_y = guided_centroids[exp_index]

        # Calculate centroid offsets for each offset target, relative to its nominal fiber center.
        offset = np.sqrt(
            (guided_x - offset_x0[:, np.newaxis, np.newaxis])**2 +
            (guided_y - offset_y0[:, np.newaxis, np.newaxis])**2)
        
        # Calculate centroid offsets for each offset target, relative to where its fiber center would
        # be if it were designed for the same wavelength as the standard stars.
        offset_std = np.sqrt(
            (guided_x - offset_x0_std[:, np.newaxis, np.newaxis])**2 +
            (guided_y - offset_y0_std[:, np.newaxis, np.newaxis])**2)

        # Calculate the acceptance fractions for both sets of centroid offsets.
        acceptance = acceptance_model(0.5*(offset / platescale).to(u.arcsec).value)
        acceptance_std = acceptance_model(0.5*(offset_std / platescale).to(u.arcsec).value)
        
        # Calculate the acceptance fraction ratios, tabulated for each offset target, wavelength and time.
        # The ratio calculated this way gives the correction of eqn (13).
        correction[exp_index] = acceptance_std / acceptance

    # Average the correction over each exposure time slice.
    avg_correction = np.mean(np.mean(correction, axis=-1), axis=0)

    return avg_correction

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--plate', type=str, default='6641',
        help='plate id')
    parser.add_argument('-m', '--mjd', type=str, default='56383',
        help='observation mjd')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
    args = parser.parse_args()

    filename = 'corrections-%s-%s.hdf5' % (args.plate, args.mjd) if args.output is None else args.output
    outfile = h5py.File(filename, 'w')

    offsets = calculate_target_offsets(int(args.plate), int(args.mjd), wlen_grid_steps=15)
    corrections = calculate_corrections(offsets)

    outfile.create_dataset('wave', data=offsets['wlen_grid'])
    outfile.create_group(args.plate)
    outfile[args.plate].create_group(args.mjd)

    for fiber, correction in zip(offsets['fibers'], corrections):
        dset = outfile.create_dataset('/'.join([args.plate,args.mjd,str(fiber)]), data=correction, dtype='f4')

    outfile.close()


if __name__ == "__main__":
    main()

