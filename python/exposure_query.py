#!/usr/bin/env python
#
# Created 19-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# ./python/exposure_query.py --verbose --input work/blue-plate-mjd-list.txt --output work/plate-mjd-ha-psf.txt

# Import libraries
import argparse
import math
import json

import numpy as np

import os
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
        if az == np.nan:
            print ra, dec, lat, ha, sin_alt, alt, cos_az
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

def add_exp_info(plate_info, exp_keys):
    plate = plate_info['plate']

    # Iterate over exposures and gather information of interest
    psf_fwhm_list = []
    rmsoff_list = []
    ha_list = []
    alt_list = []
    nexp = 0

    for exposure in plate_info['exposures']:
        exposure_info = dict()

        cframe = CFrame(os.path.join(plate_info['bossdir'], str(plate), 'spCFrame-%s.fits' % exposure))
        examine_exposure(exposure_info, cframe, exp_keys)

        plate_info[exposure] = exposure_info

        if exposure_info['SEEING50'] > 0:
            nexp += 1
            psf_fwhm_list.append(exposure_info['SEEING50'])
            rmsoff_list.append(exposure_info['RMSOFF50'])
        ha_list.append(exposure_info['mean_ha'])
        alt_list.append(exposure_info['mean_alt'])

    plate_info['mean_psf_fwhm'] = np.mean(psf_fwhm_list) if nexp > 0 else 0
    plate_info['mean_rmsoff'] = np.mean(rmsoff_list) if nexp > 0 else 0
    plate_info['mean_ha'] = np.mean(ha_list)
    plate_info['mean_alt'] = np.mean(alt_list)

def add_plate_info(plate_info, plate_keys):
    plate = plate_info['plate']
    mjd = plate_info['mjd']
    plate_name = os.path.join(plate_info['bossdir'], str(plate), 'spPlate-%s-%s.fits' % (plate, mjd))
    spPlate = Plate(plate_name)
    # Construct list of exposure IDs
    exposures = list()
    nexp = spPlate.header['NEXP']
    if (nexp % 4) != 0:
        print 'Warning, unexpected number of exposures...'
    for i in range(nexp/4):
        exposures.append('-'.join(spPlate.header['EXPID%02d'%(i+1)].split('-')[0:2]))
    plate_info['exposures'] = exposures

    # Process Plate header keywords
    for keyword in plate_keys:
        plate_info[keyword] = spPlate.header[keyword]

def add_plugmap_info(plate_info, plugmap_keys):
    import yanny
    plate = plate_info['plate']
    mjd = plate_info['mjd']

    plan_name = os.path.join(plate_info['bossdir'], str(plate), 'spPlancomb-%s-%s.par' % (plate, mjd))
    plan = yanny.yanny(plan_name)

    # look up plugmap filename
    unique_mapnames = set(zip(plan['SPEXP']['mjd'], plan['SPEXP']['mapname']))
    (mapmjd, mapname) = unique_mapnames.pop()
    plugmap_name = os.path.join(plate_info['speclog'], str(mapmjd), 'plPlugMapM-%s.par' % mapname)
    plugmap = yanny.yanny(plugmap_name)

    plate_info['mapmjd'] = mapmjd
    plate_info['mapname'] = mapname

    # Process plugmap header keywords
    for keyword in plugmap_keys:
        plate_info[keyword] = plugmap[keyword]

    design_ra = Angle(float(plugmap['raCen']), unit=u.degree)
    design_dec = Angle(float(plugmap['decCen']), unit=u.degree)
    design_ha = Angle(float(plugmap['haMin']), unit=u.degree)
    design_alt, design_az = equatorial_to_horizontal(design_ra, design_dec, apolat, design_ha)

    plate_info['design_ra'] = design_ra.to(u.degree).value
    plate_info['design_dec'] = design_dec.to(u.degree).value
    plate_info['design_ha'] = design_ha.to(u.degree).value
    plate_info['design_alt'] = design_alt.to(u.degree).value
    plate_info['design_az'] = design_az.to(u.degree).value

    
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
        help='specify plate id')
    parser.add_argument('-m', '--mjd', type=str, default=None,
        help='specify observation mjd')
    parser.add_argument('--bossdir', type=str, default='/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/v5_7_0', 
        help='path to 2D reduction dir ($BOSS_SPECTRO_REDUX/$RUN2D)')
    parser.add_argument('--speclog', type=str, default='/home/boss/products/NULL/speclog/trunk',
        help='path to speclog base dir ($SPECLOG_DIR)')
    parser.add_argument('-d', '--delim', type=str, default=' ',
        help='output token deliminator (useful for making tables)')
    parser.add_argument('-o', '--outdir', type=str, default=None,
        help='output directory')
    parser.add_argument('--skip-exp', action='store_true',
        help='skip spCFrame reading')
    parser.add_argument('--skip-plugmap', action='store_true',
        help='skip plugMapM reading')
    args = parser.parse_args()

    if not (args.plate and args.mjd):
        print 'Must specify a plate id and observation mjd!'
        return -1

    plate_info = {}
    plate_info['plate'] = args.plate
    plate_info['mjd'] = args.mjd
    plate_info['bossdir'] = args.bossdir
    plate_info['speclog'] = args.speclog

    summary_keys = ['plate', 'mjd']
    summary_fmts = ['%s', '%s']
    plate_keys = ['SEEING50', 'RMSOFF50', 'AIRMASS', 'ALT', 'BESTEXP', 'NSTD']
    add_plate_info(plate_info, plate_keys)

    plugmap_keys = ['haMin']
    if not args.skip_plugmap:
        add_plugmap_info(plate_info, plugmap_keys)
        summary_keys += ['mapmjd', 'mapname']
        summary_fmts += ['%s', '%s']

    exp_keys = ['MJD', 'SEEING50', 'RMSOFF50', 'AIRMASS', 'ALT']
    if not args.skip_exp:
        add_exp_info(plate_info, exp_keys)
        summary_keys += ['mean_ha', 'mean_psf_fwhm']
        summary_fmts += ['%.4f', '%.4f']

    header = args.delim.join(summary_keys)
    summary = args.delim.join([fmt % plate_info[key] for key,fmt in zip(summary_keys, summary_fmts)])
    print header
    print summary
    if args.outdir:
        summary_filename = os.path.join(args.outdir, 'expinfo-summary-%s-%s.txt' % (plate_info['plate'], plate_info['mjd']))
        with open(summary_filename, 'w') as outfile:
            outfile.write(header)
            outfile.write('\n')
            outfile.write(summary)
            outfile.write('\n')
        json_filename = os.path.join(args.outdir, 'expinfo-%s-%s.json' % (plate_info['plate'], plate_info['mjd']))
        with open(json_filename, 'w') as outfile:
            json.dump(plate_info, outfile, sort_keys=True, indent=2)
    if args.verbose:
        print json.dumps(plate_info, sort_keys=True, indent=2)

if __name__ == '__main__':
    main()







