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

from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u
#from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle

def equatorial_to_horizontal(ra, dec, lat, ha):
    sin_alt = np.sin(dec)*np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha)
    alt = np.arcsin(sin_alt)
    cos_az = (np.sin(dec) - sin_alt*np.sin(lat))/(np.cos(alt)*np.cos(lat))
    if np.sin(ha) < 0:
        az = np.arccos(cos_az)
    else:
        az = 2*np.pi*u.radian - np.arccos(cos_az)
    return alt, az

class Fluxcalib(object):
    def __init__(self, fluxcalibfile):
        hdulist = fits.open(fluxcalibfile)
        self.header = hdulist[0].header
        self.data = hdulist[0].data
        hdulist.close()
    
class CFrame(object):
    def __init__(self, cframefile):
        hdulist = fits.open(cframefile)
        self.header = hdulist[0].header
        self.flux = hdulist[0].data
        self.ivar = hdulist[1].data
        self.mask = hdulist[2].data
        self.loglam = hdulist[3].data
        self.wdisp = hdulist[4].data
        self.sky = hdulist[6].data
        hdulist.close()

class Plate(object):
    def __init__(self, platefile):
        hdulist = fits.open(platefile)
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
    parser.add_argument('-k', '--keyword', type=str, help='header keyword')
    parser.add_argument('-p', '--plate', type=str, help='Plate')
    parser.add_argument('-m', '--mjd', type=str, help='MJD')
    parser.add_argument('--bossdir', type=str, help='Input default boss reduction $BOSS_SPECTRO_REDUX/$RUN2D')
    parser.add_argument('-i', '--input', type=str, help='Input plate list file to read')
    parser.add_argument('-d', '--delim', type=str, default=' ', help='Table deliminator')
    parser.add_argument('--speclog', type=str, help='Speclog base dir $SPECLOG_DIR')
    parser.add_argument('-o', '--output', type=str, help='Output base directory')
    parser.add_argument('-a','--alt', action='store_true', help='Perform exposure alt calculation')

    args = parser.parse_args()

    if not (args.input or (args.plate and args.mjd)):
        print 'Must specify a plate list file to read or a specific plate and mjd pair!'
        return -1

    keywords = ['MJD', 'SEEING50', 'RMSOFF50', 'AIRMASS', 'ALT']
    if args.keyword:
        keywords.append(args.keyword)

    plate_keywords = []#,'SEEING50','RMSOFF50']
    plugmap_keywords = ['haMin']#,'cartridgeId']#,'raCen','decCen']

    # APO Geographical Location
    apolat = Angle('32d46m49s')
    apolon = Angle('-105d49m13s')
    apo = EarthLocation.from_geodetic(apolon, apolat)

    keys = ['plate', 'mjd', 'id']
    keys += plate_keywords
    keys += keywords
    if args.alt:
        keys += ['mean_alt', 'design_alt', 'mean_ha']
    keys += plugmap_keywords

    print args.delim.join([key for key in keys])

    def examine_exposures(datadir, speclog, plate, mjd, keywords, outputdir):
        plate_name = os.path.join(datadir, str(plate), 'spPlate-%s-%s.fits' % (plate, mjd))
        spPlate = Plate(plate_name)

        plan_name = os.path.join(datadir, str(plate), 'spPlancomb-%s-%s.par' % (plate, mjd))
        plan = yanny.yanny(plan_name)
        unique_mapnames = set(zip(plan['SPEXP']['mjd'], plan['SPEXP']['mapname']))
        (mapmjd, mapname) = unique_mapnames.pop()
        plugmap_name = os.path.join(speclog, str(mapmjd), 'plPlugMapM-%s.par' % mapname)
        plugmap = yanny.yanny(plugmap_name)

        # Construct list of exposure IDs
        exposures = list()
        nexp = spPlate.header['NEXP']
        if (nexp % 4) != 0:
            print 'Warning, unexpected number of exposures...'
        for i in range(nexp/4):
            exposures.append('-'.join(spPlate.header['EXPID%02d'%(i+1)].split('-')[0:2]))

        # Iterate over exposures and gather information of interest
        for exposure in exposures:
            info = dict()
            info['plate'] = plate
            info['mjd'] = mjd
            info['id'] = exposure
            # Process Plate header keywords
            for keyword in plate_keywords:
                info[keyword] = spPlate.header[keyword]
            # Process CFrame header keywords
            cframe = CFrame(os.path.join(datadir, plate, 'spCFrame-%s.fits'%exposure))
            for keyword in keywords:
                info[keyword] = str(cframe.header[keyword])
            # Process plugmap header keywords
            for keyword in plugmap_keywords:
                info[keyword] = str(plugmap[keyword])
            # Do fancy stuff
            if args.alt:
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

                info['mean_alt'] = '%.3f' % (alt.to(u.degree)).degree
                if ha > np.pi*u.radian:
                    ha -= 2*np.pi*u.radian
                info['mean_ha'] = '%.4f' % ha.degree

                design_ra = Angle(float(plugmap['raCen']), unit=u.degree)
                design_dec = Angle(float(plugmap['decCen']), unit=u.degree)
                design_ha = Angle(float(plugmap['haMin'])%360, unit=u.degree)
                design_alt, design_az = equatorial_to_horizontal(design_ra, design_dec, apolat, design_ha)
                info['design_alt'] = '%.3f' % (design_alt.to(u.degree)).degree

            print args.delim.join([str(info[key]) for key in keys])

    if args.input:
        with open(args.input) as platelist:
            for line in platelist:
                (plate, mjd) = line.strip().split()
                examine_exposures(datadir=args.bossdir, speclog=args.speclog, 
                    plate=plate, mjd=mjd, keywords=keywords, outputdir=args.output)
    else:
        examine_exposures(datadir=args.bossdir, speclog=args.speclog, 
            plate=args.plate, mjd=args.mjd, keywords=keywords, outputdir=args.output)

if __name__ == '__main__':
    main()







