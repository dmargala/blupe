#!/usr/bin/env python
#
# Created 19-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# ./python/exposure_query.py --verbose --input work/blue-plate-mjd-list.txt --output work/plate-mjd-ha-psf.txt

# Import libraries
import argparse
import math
from datetime import datetime

import numpy as np

import os
import yanny
from astropy.io import fits

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import Angle

# APO Geographical Location
apolat = Angle('32d46m49s')
apolon = Angle('-105d49m13s')
apo = EarthLocation.from_geodetic(apolon, apolat)

def equatorial_to_horizontal(ra, dec, lat, ha):
    sin_alt = np.sin(dec)*np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha)
    alt = np.arcsin(sin_alt)
    cos_az = (np.sin(dec) - sin_alt*np.sin(lat))/(np.cos(alt)*np.cos(lat))
    if np.sin(ha) < 0:
        az = np.arccos(cos_az)
    else:
        az = 2*np.pi*u.radian - np.arccos(cos_az)
    return alt, az

def examine_exposure(info, cframe, cframe_keys):
    # Process CFrame header keywords
    for keyword in cframe_keys:
        info[keyword] = cframe.header[keyword]

    obs_ra = cframe.header['RADEG']
    obs_dec = cframe.header['DECDEG']
    taibeg = cframe.header['TAI-BEG']
    taiend = cframe.header['TAI-END']

    taimid = 0.5*(taibeg+taiend)
    dec = Angle(obs_dec, u.degree)
    ra = Angle(obs_ra, u.degree)

    time = Time(taimid/86400.0, format='mjd', scale='tai', location=apo)
    lst = time.sidereal_time('apparent')
    ha = (lst - ra)

    if ha > np.pi*u.radian:
        ha -= 2*np.pi*u.radian
    elif ha < -np.pi*u.radian:
        ha += 2*np.pi*u.radian
    info['mean_ha'] = ha.to(u.degree).value

    alt, az = equatorial_to_horizontal(ra, dec, apolat, ha)
    info['mean_alt'] = alt.to(u.degree).value


class Fluxcalib(object):
    def __init__(self, name):
        hdulist = fits.open(name)
        self.header = hdulist[0].header
        self.data = hdulist[0].data
        hdulist.close()
    
class CFrame(object):
    def __init__(self, name):
        hdulist = fits.open(name)
        self.header = hdulist[0].header
        self.flux = hdulist[0].data
        self.ivar = hdulist[1].data
        self.mask = hdulist[2].data
        self.loglam = hdulist[3].data
        self.wdisp = hdulist[4].data
        self.sky = hdulist[6].data
        hdulist.close()

class Plate(object):
    def __init__(self, name):
        hdulist = fits.open(name)
        self.header = hdulist[0].header
        self.flux = hdulist[0].data
        self.ivar = hdulist[1].data
        self.andmask = hdulist[2].data
        self.ormask = hdulist[3].data
        self.wdisp = hdulist[4].data
        self.sky = hdulist[6].data
        hdulist.close() 

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v','--verbose', action='store_true',
        help='provide more verbose output')
    parser.add_argument('-p', '--plate', type=str, default=None,
        help='specify plate number for single plate-mjd query')
    parser.add_argument('-m', '--mjd', type=str, default=None,
        help='specify mjd for single plate-mjd query')
    parser.add_argument('-i', '--input', type=str, default=None,
        help='specify file with plate,mjd entries to query')
    parser.add_argument('--bossdir', type=str, default='/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/v5_7_0', 
        help='path to 2D reduction dir ($BOSS_SPECTRO_REDUX/$RUN2D)')
    parser.add_argument('--speclog', type=str, default='/home/boss/products/NULL/speclog/trunk',
        help='path to speclog base dir ($SPECLOG_DIR)')
    parser.add_argument('-d', '--delim', type=str, default=' ',
        help='output token deliminator (useful for making tables)')
    parser.add_argument('-o', '--output', type=str, default=None,
        help='output filename')
    args = parser.parse_args()

    if not (args.input or (args.plate and args.mjd)):
        print 'Must specify a plate list file to read or a specific plate and mjd pair!'
        return -1

    plate_keys = ['SEEING50', 'RMSOFF50', 'AIRMASS', 'ALT', 'BESTEXP', 'NSTD']
    plugmap_keys = ['haMin']
    new_plate_keys = ['plate', 'mjd', 'mapmjd', 'design_alt']

    cframe_keys = ['MJD', 'SEEING50', 'RMSOFF50', 'AIRMASS', 'ALT']
    new_exp_keys = ['id', 'mean_alt', 'mean_ha']

    plate_mjd_list = []
    if args.input:
        with open(args.input) as platelist:
            for line in platelist:
                plate_mjd_list.append(line.strip().split())
    else:
        plate_mjd_list.append((args.plate, args.mjd))

    plate_info_list = []
    
    for plate, mjd in plate_mjd_list: 
        plate_name = os.path.join(args.bossdir, str(plate), 'spPlate-%s-%s.fits' % (plate, mjd))
        spPlate = Plate(plate_name)

        plan_name = os.path.join(args.bossdir, str(plate), 'spPlancomb-%s-%s.par' % (plate, mjd))
        plan = yanny.yanny(plan_name)
        unique_mapnames = set(zip(plan['SPEXP']['mjd'], plan['SPEXP']['mapname']))
        (mapmjd, mapname) = unique_mapnames.pop()
        plugmap_name = os.path.join(args.speclog, str(mapmjd), 'plPlugMapM-%s.par' % mapname)
        plugmap = yanny.yanny(plugmap_name)

        # Construct list of exposure IDs
        exposures = list()
        nexp = spPlate.header['NEXP']
        if (nexp % 4) != 0:
            print 'Warning, unexpected number of exposures...'
        for i in range(nexp/4):
            exposures.append('-'.join(spPlate.header['EXPID%02d'%(i+1)].split('-')[0:2]))

        plate_info = {}
        plate_info['plate'] = plate
        plate_info['mjd'] = mjd
        plate_info['mapmjd'] = mapmjd
        plate_info['mapname'] = mapname
        plate_info['exposures'] = exposures
        # Process Plate header keywords
        for keyword in plate_keys:
            plate_info[keyword] = spPlate.header[keyword]
        # Process plugmap header keywords
        for keyword in plugmap_keys:
            plate_info[keyword] = plugmap[keyword]

        design_ra = Angle(float(plugmap['raCen']), unit=u.degree)
        design_dec = Angle(float(plugmap['decCen']), unit=u.degree)
        design_ha = Angle(float(plugmap['haMin'])%360, unit=u.degree)
        design_alt, design_az = equatorial_to_horizontal(design_ra, design_dec, apolat, design_ha)
        plate_info['design_alt'] = design_alt.to(u.degree).value

        # Iterate over exposures and gather information of interest
        psf_fwhm_list = []
        ha_list = []
        alt_list = []
        rmsoff_list = []

        if args.verbose:
            print args.delim.join([str(plate_info[key]) for key in (new_plate_keys + plate_keys + plugmap_keys)])

        for exposure in exposures:
            exposure_info = dict()
            exposure_info['id'] = exposure

            cframe = CFrame(os.path.join(args.bossdir, str(plate), 'spCFrame-%s.fits' % exposure))

            examine_exposure(exposure_info, cframe, cframe_keys)

            if args.verbose:
                print '\t', args.delim.join([str(exposure_info[key]) for key in (new_exp_keys + cframe_keys)])

            psf_fwhm_list.append(exposure_info['SEEING50'])
            ha_list.append(exposure_info['mean_ha'])
            alt_list.append(exposure_info['mean_alt'])
            rmsoff_list.append(exposure_info['RMSOFF50'])

        if args.verbose:
            print np.mean(psf_fwhm_list), np.mean(ha_list), np.mean(alt_list)

        plate_info['mean_psf_fwhm'] = np.mean(psf_fwhm_list)
        plate_info['mean_ha'] = np.mean(ha_list)
        plate_info['mean_alt'] = np.mean(alt_list)
        plate_info['mean_rmsoff'] = np.mean(rmsoff_list)

        plate_info_list.append(plate_info)

    if args.output:
        outkeys = ['plate', 'mjd', 'mapmjd', 'mapname', 'mean_ha', 'mean_psf_fwhm', 'mean_alt', 'mean_rmsoff', 'design_alt']
        outfmt = '%s %s %s %s %.4f %.4f %.4f %.8f %.4f\n'
        with open(args.output, 'w') as output:
            for plate_info in plate_info_list:
                values = [plate_info[key] for key in outkeys]
                output.write(outfmt % tuple(values))

if __name__ == '__main__':
    main()







