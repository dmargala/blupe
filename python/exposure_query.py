#!/usr/bin/env python
#
# Created 19-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

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

    alt, az = equatorial_to_horizontal(ra, dec, apolat, ha)

    info['mean_alt'] = alt.to(u.degree).value
    if ha > np.pi*u.radian:
        ha -= 2*np.pi*u.radian
    info['mean_ha'] = ha.to(u.degree).value


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
    parser.add_argument('-v','--verbose', action='store_true', help='provide more verbose output')
    parser.add_argument('-p', '--plate', type=str, help='Plate')
    parser.add_argument('-m', '--mjd', type=str, help='MJD')
    parser.add_argument('--bossdir', type=str, help='Input default boss reduction $BOSS_SPECTRO_REDUX/$RUN2D')
    parser.add_argument('-i', '--input', type=str, help='Input plate list file to read')
    parser.add_argument('-d', '--delim', type=str, default=' ', help='Table deliminator')
    parser.add_argument('--speclog', type=str, help='Speclog base dir $SPECLOG_DIR')
    parser.add_argument('-o', '--output', type=str, default='', help='Output prefix')

    args = parser.parse_args()

    if not (args.input or (args.plate and args.mjd)):
        print 'Must specify a plate list file to read or a specific plate and mjd pair!'
        return -1

    plate_keys = []#,'SEEING50','RMSOFF50']
    plugmap_keys = ['haMin']#,'cartridgeId']#,'raCen','decCen']
    cframe_keys = ['MJD', 'SEEING50', 'RMSOFF50', 'AIRMASS', 'ALT']
    new_exp_keys = ['id', 'mean_alt', 'mean_ha']
    new_plate_keys = ['plate', 'mjd', 'design_alt']

    plate_mjd_list = []
    if args.input:
        with open(args.input) as platelist:
            for line in platelist:
                plate_mjd_list.append(line.strip().split())
    else:
        plate_mjd_list.append((args.plate, args.mjd))

    mean_psf_fwhm = []
    mean_ha = []
    mean_alt = []
    
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

        print args.delim.join([str(plate_info[key]) for key in (new_plate_keys + plate_keys + plugmap_keys)])

        for exposure in exposures:
            exposure_info = dict()
            exposure_info['id'] = exposure

            cframe = CFrame(os.path.join(args.bossdir, str(plate), 'spCFrame-%s.fits' % exposure))

            examine_exposure(exposure_info, cframe, cframe_keys)

            print '\t', args.delim.join([str(exposure_info[key]) for key in (new_exp_keys + cframe_keys)])

            psf_fwhm_list.append(exposure_info['SEEING50'])
            ha_list.append(exposure_info['mean_ha'])
            alt_list.append(exposure_info['mean_alt'])

        print np.mean(psf_fwhm_list), np.mean(ha_list), np.mean(alt_list)

        mean_psf_fwhm.append(np.mean(psf_fwhm_list))
        mean_ha.append(np.mean(ha_list))
        mean_alt.append(np.mean(alt_list))

    if args.output:
        with open(args.output, 'w') as output:
            for i, (plate, mjd) in enumerate(plate_mjd_list):
                output.write('%s %s %.4f %.4f %.4f\n' % (plate, mjd, mean_ha[i], mean_psf_fwhm[i], mean_alt[i]))

if __name__ == '__main__':
    main()







