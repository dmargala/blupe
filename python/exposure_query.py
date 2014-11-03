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
#import yanny
from astropy.io import fits

from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.coordinates import SkyCoord


import sidereal

def tai2utc(tai):
    return tai-(40587*86400)

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
    apolat = (32. + 46/60. + 49/3600.)*math.pi/180.
    apolon = -(105. + 49/60. + 13/3600.)*math.pi/180.
    #apo = sidereal.LatLon(apolat,apolon)
    apolat_str = '32d46m49s'
    apolon_str = '-105d49m13s'
    apo = EarthLocation.from_geodetic(apolon_str, apolat_str)

    keys = ['plate', 'mjd', 'id']
    keys += plate_keywords
    keys += keywords
    if args.alt:
        keys += ['meanAlt', 'designAlt', 'meanHa']
    keys += plugmap_keywords

    print args.delim.join([key for key in keys])

    def examine_plate(path, plate, mjd):
        plate_name = os.path.join(path, plate, 'spPlate-%s-%s.fits' % (plate, mjd))
        plate = Plate(plate_name)

        taibeg = plate.header['TAI-BEG']
        taiend = plate.header['TAI-END']
        taimid = 0.5*(taibeg+taiend)
        print taimid
        print taimid/86400

        obs_ra = plate.header['RADEG']
        obs_dec = plate.header['DECDEG']

        ut = datetime.utcfromtimestamp(tai2utc(taimid))

        pointing = sidereal.RADec(obs_ra*math.pi/180.,obs_dec*math.pi/180.)
        ha = pointing.hourAngle(ut, apolon) # radians

        print ha*180./math.pi

        time = Time(tai2utc(taimid), format='unix', scale='tai', location=apo)
        print time.tai
        print time.mjd
        lst = time.sidereal_time('mean')

        c = SkyCoord(ra=obs_ra*u.degree, dec=obs_dec*u.degree)
        lha = (lst - c.ra)
        print lha.degree

    examine_plate(args.bossdir, args.plate, args.mjd)

    def examine_exposures(datadir, speclog, plate, mjd, keywords, outputdir):
        plate_name = os.path.join(datadir, plate, 'spPlate-%s-%s.fits' % (plate, mjd))
        plate = Plate(plate_name)

        plan_name = os.path.join(datadir, plate, 'spPlancomb-%s-%s.par' % (plate, mjd))
        plan = yanny.yanny(plan_name)
        spexp = plan['SPEXP']
        unique_mapnames = set(zip(spexp['mjd'],spexp['mapname']))
        (mapmjd, mapname) = unique_mapnames.pop()
        plugmap_name = os.path.join(speclog, str(mapmjd), 'plPlugMapM-%s.par' % mapname)
        plugmap = yanny.yanny(plugmap_name)

        # Construct list of exposure IDs
        exposures = list()
        nexp = plate.header['NEXP']
        if (nexp % 4) != 0:
            print 'Warning, unexpected number of exposures...'
        for i in range(nexp/4):
            exposures.append('-'.join(plate.header['EXPID%02d'%(i+1)].split('-')[0:2]))

        # Iterate over exposures and gather information of interest
        for exposure in exposures:
            info = dict()
            info['plate'] = plate
            info['mjd'] = mjd
            info['id'] = exposure
            # Process Plate header keywords
            for keyword in plate_keywords:
                info[keyword] = plate.header[keyword]
            # Process CFrame header keywords
            cframe = CFrame(os.path.join(datadir, plate, 'spCFrame-%s.fits'%exposure))
            for keyword in keywords:
                info[keyword] = str(cframe.header[keyword])
            # Process plugmap header keywords
            for keyword in plugmap_keywords:
                info[keyword] = str(plugmap[keyword])
            # Do fancy stuff
            if args.alt:
                ra = cframe.header['RADEG']
                dec = cframe.header['DECDEG']
                taibeg = cframe.header['TAI-BEG']
                taiend = cframe.header['TAI-END']
                taimid = 0.5*(taibeg+taiend)
                ut = datetime.utcfromtimestamp(tai2utc(taimid))

                time = Time(taibeg, format='unix', scale='tai', location=apo)
                lst = time.sidereal_time('mean')

                c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
                lha = lst.hour - c.hour

                pointing = sidereal.RADec(ra*math.pi/180.,dec*math.pi/180.)
                ha = pointing.hourAngle(ut, apolon) # radians

                print lha, ha

                altaz = pointing.altAz(ha, apolat) # AltAz instance
                info['meanAlt'] = '%.3f' % (altaz.alt*180./math.pi) # degrees
                if ha > math.pi:
                    ha -= 2*math.pi
                info['meanHa'] = '%.4f' % (ha*180./math.pi) # degrees

                design_ra = float(plugmap['raCen'])*math.pi/180.
                design_dec = float(plugmap['decCen'])*math.pi/180.
                design_ha = (float(plugmap['haMin'])%360)*math.pi/180.
                desing_pointing = sidereal.RADec(design_ra, design_dec)
                # The sidereal altAz calculation bombs for designHa == 0, add a tiny smudge factor
                if design_ha == 0:
                    design_ha = 1e-7
                design_altaz = design_pointing.altAz(design_ha, apolat)
                info['designAlt'] = '%.3f' % (design_altaz.alt*180./math.pi)

            print args.delim.join([str(info[key]) for key in keys])

    # if args.input:
    #     with open(args.input) as platelist:
    #         for line in platelist:
    #             (plate, mjd) = line.strip().split()
    #             examine_exposures(datadir=args.datadir, speclog=args.speclog, 
    #                 plate=plate, mjd=mjd, keywords=keywords, outputdir=args.output)
    # else:
    #     examine_exposures(datadir=args.datadir, speclog=args.speclog, 
    #         plate=args.plate, mjd=args.mjd, keywords=keywords, outputdir=args.output)

if __name__ == '__main__':
    main()







