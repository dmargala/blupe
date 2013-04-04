#!/usr/bin/env python
#
# Created 19-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# Import libraries
import numpy as np
import os
import yanny
import pyfits as pf
import sys

import math
from datetime import datetime

import sidereal

def tai2utc(tai):
    return tai-(40587*86400)

class Fluxcalib(object):
    def __init__(self, fluxcalibfile):
        hdulist = pf.open(fluxcalibfile)
        self.header = hdulist[0].header
        self.data = hdulist[0].data
        hdulist.close()
    
class CFrame(object):
    def __init__(self, cframefile):
        hdulist = pf.open(cframefile)
        self.header = hdulist[0].header
        self.flux = hdulist[0].data
        self.ivar = hdulist[1].data
        self.mask = hdulist[2].data
        self.loglam = hdulist[3].data
        self.wdisp  = hdulist[4].data
        self.sky    = hdulist[6].data
        hdulist.close()

class Plate(object):
    def __init__(self, platefile):
        hdulist = pf.open(platefile)
        self.header = hdulist[0].header
        self.flux = hdulist[0].data
        self.ivar = hdulist[1].data
        self.andmask = hdulist[2].data
        self.ormask = hdulist[3].data
        self.wdisp  = hdulist[4].data
        self.sky    = hdulist[6].data
        hdulist.close() 

import optparse

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-k", "--keyword",  type="string", help="header keyword")
parser.add_option('-p', "--plate", type="string", help="Plate")
parser.add_option('-m', "--mjd", type="string", help="MJD")
parser.add_option("--datadir", type="string", help="Input blue reduction $BOSS_SPECTRO_REDUX/$RUN2D")
parser.add_option("--bossdir", type="string", help="Input default boss reduction $BOSS_SPECTRO_REDUX/$RUN2D")
parser.add_option('-i', "--input", type="string", help="Input plate list file to read")
parser.add_option('-d', "--delim", type="string", help="Table deliminator", default=" ")
parser.add_option('--speclog', type="string", help="Speclog base dir $SPECLOG_DIR")
parser.add_option('-o', "--output", type="string", help="Output base directory")
parser.add_option('-a',"--alt",action="store_true",help="Perform exposure alt calculation",default=False)
parser.add_option('-c',"--calib",action="store_true",help="Print Fluxcalib vector to text file",default=False)

opts, args = parser.parse_args()

if not (opts.input or (opts.plate and opts.mjd)):
    print "Must specify a plate list file to read or a specific plate and mjd pair!"
    sys.exit(1)

keywords = []
if opts.keyword:
    keywords.append(opts.keyword)
else:
    keywords = ["AIRMASS","ALT"]

plateKeywords = ["PLATEID","MJD","SEEING50","RMSOFF50"]
plugmapKeywords = ["haMin","cartridgeId"]#,"raCen","decCen"]

# APO Geographical Location
apolat = (32. + 46/60. + 49/3600.)*math.pi/180.
apolon = -(105. + 49/60. + 13/3600.)*math.pi/180.
apo = sidereal.LatLon(apolat,apolon)

keys = ['id']
keys += plateKeywords
keys += keywords
if opts.alt:
    keys += ['meanAlt','meanHa','designAlt']
keys += plugmapKeywords

print opts.delim.join([key for key in keys])

def examineExposures(datadir,speclog,plate,mjd,keywords,outputdir):
    spPlate = Plate(os.path.join(datadir,plate,"spPlate-%s-%s.fits"%(plate,mjd)))

    plan2dfile = os.path.join(datadir,plate,"spPlan2d-%s-%s.par"%(plate,mjd))
    plan2d = yanny.yanny(plan2dfile)
    spexp = plan2d['SPEXP']
    uniqueMapNames = set(zip(spexp['mjd'],spexp['mapname']))
    if len(uniqueMapNames) != 1:
        print "warning unique plugmap names != 1"
        sys.exit(1)
    (mapmjd,mapname) = uniqueMapNames.pop()
    plugmapfile = os.path.join(speclog,str(mapmjd),'plPlugMapM-%s.par'%mapname)
    plugmap = yanny.yanny(plugmapfile)

    exposureIDs = list()
    for i in range(len(spexp['flavor'])):
        if spexp['flavor'][i] == 'science':
                exposureIDs.append('-'.join(spexp['name'][i][0].split('-')[1:3]).split('.')[0])
    
    exposureData = list()
    for exposure in exposureIDs:

        info = dict()
        info['id'] = exposure
        # Process Plate header keywords
        for keyword in plateKeywords:
            info[keyword] = spPlate.header[keyword]
        # Process CFrame header keywords
        cframe = CFrame(os.path.join(datadir,plate,"spCFrame-%s.fits"%exposure))
        for keyword in keywords:
            info[keyword] = str(cframe.header[keyword])
        # Process plugmap header keywords
        for keyword in plugmapKeywords:
            info[keyword] = str(plugmap[keyword])
        # Do fancy stuff
        if opts.alt:
            ra = cframe.header['RADEG']
            dec = cframe.header['DECDEG']
            taibeg = cframe.header['TAI-BEG']
            taiend = cframe.header['TAI-END']
            ut = datetime.utcfromtimestamp(tai2utc((taibeg+taiend)/2.))

            pointing = sidereal.RADec(ra*math.pi/180.,dec*math.pi/180.)
            ha = pointing.hourAngle(ut, apolon) # radians
            altaz = pointing.altAz(ha, apolat) # AltAz instance
            info['meanAlt'] = "%.3f" % (altaz.alt*180./math.pi) # degrees
            if ha > math.pi:
                ha -= 2*math.pi
            info['meanHa'] = "%.4f" % (ha*180./math.pi) # degrees

            designRA = float(plugmap['raCen'])*math.pi/180.
            designDec = float(plugmap['decCen'])*math.pi/180.
            designHa = (float(plugmap['haMin'])%360)*math.pi/180.
            designPointing = sidereal.RADec(designRA,designDec)
            # The sidereal altAz calculation bombs for designHa == 0, add a tiny smudge factor
            if designHa == 0:
                designHa = 1e-7
            designAltAz = designPointing.altAz(designHa, apolat)
            info['designAlt'] = "%.3f" % (designAltAz.alt*180./math.pi)

        print opts.delim.join([str(info[key]) for key in keys])
        #exposureData.append(info)

        if opts.calib:
            BSR = '/home/dmargala/BSR'
            bossFluxcalibname = os.path.join(BSR,'v5_6_0',plate,'spFluxcalib-%s.fits.gz'%exposure)
            blueFluxcalibname = os.path.join(datadir,plate,'spFluxcalib-%s.fits.gz'%exposure)

            bossFluxcalib = Fluxcalib(bossFluxcalibname)
            blueFluxcalib = Fluxcalib(blueFluxcalibname)

            if not os.path.exists(os.path.join(outputdir,plate)):
                os.makedirs(os.path.join(outputdir,plate))
            textfile = open(os.path.join(outputdir,plate,'spFluxcalibRatio-%s.txt'%exposure),'w')
            for pixel in range(len(bossFluxcalib.data[0])):
                value = 0
                if blueFluxcalib.data[0][pixel] > 0:
                    value = bossFluxcalib.data[0][pixel]/blueFluxcalib.data[0][pixel]
                textfile.write("%d %f\n"%(10**cframe.loglam[0][pixel],value))
            textfile.close()
    return exposureData

exposureData = list()
if opts.input:
    platelist = open(opts.input)
    for line in platelist:
        (plate,mjd) = line.split()
        try:
            exposureData += examineExposures(datadir=opts.datadir,speclog=opts.speclog,\
                plate=plate,mjd=mjd,keywords=keywords,outputdir=opts.output)
        except ValueError:
            print 'error procressing exposures for %s-%s...' % (plate,mjd)
else:
    exposureData += examineExposures(opts.datadir,opts.speclog,opts.plate,opts.mjd,keywords)







