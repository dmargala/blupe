#!/usr/bin/env python
#
# Created 19-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# Import libraries
import numpy as np
import os
#import yanny
import pyfits as pf
import sys
import re

import math


class CAMERA: pass
class RED_CAMERA (CAMERA): pass
class BLUE_CAMERA (CAMERA): pass

class guideDerivs(object):
    def __init__(self, filename):
        prefix = "~/Cosmo/LyAlpha/BlueStandards/guiding/plate"

class spec(object):
    def __init__(self, filename):
        hdulist = pf.open(filename,memmap=True)
        self._header = hdulist[0].header
        self._combined = hdulist[1].data
        self._spAll = hdulist[2].data
        self._spZline = hdulist[3].data
        self._numExposures = self._header['NEXP']
        self._spCFrames = []
        for i in range(self._numExposures):
            self._spCFrames.append(hdulist[4+i].data)
        hdulist.close()
    def header(self, keyword):
        return self._header[keyword]
    def combined(self, keyword):
        return self._combined[keyword]
    def spAll(self, keyword):
        return self._spAll[keyword][0]
    def spZline(self, keyword):
        return self._spZline[keyword]
    def getNumExposures(self):
        return self._numExposures
    def getFrames(self):
        return self._spCFrames


class spAll(object):
    def __init__(self, filename):
        hdulist = pf.open(filename, memmap=True)
        self.header = hdulist[0].header
        hdulist.close()

class spPlate(object):
    def __init__(self, platefile):
        hdulist = pf.open(platefile)
        self.header = hdulist[0].header
        self.flux = hdulist[0].data
        self.ivar = hdulist[1].data
        self.andmask = hdulist[2].data
        self.ormask = hdulist[3].data
        self.wdisp = hdulist[4].data
        self.sky = hdulist[6].data
        hdulist.close()

    def getFrameIDs(self, fiber, cameraMask):
        frameIDs = []
        specNum = fiber.spectrographNumber()
        numExposures = self.header["NEXP"]
        # Target regex pattern
        expIDpattern = re.compile(r'([rb])([12])-([0-9]{8})-[0-9]{8}-[0-9]{8}')
        for expNum in range(1,numExposures+1):
            # Read next exposure ID
            expID = self.header["EXPID%02d" % expNum]
            # Parse the exposure ID
            m = expIDpattern.match(expID)
            if m:
                if int(expID[1]) != specNum: continue
                if (not (cameraMask is RED_CAMERA) and expID[1] == 'r'): continue
                if (not (cameraMask is BLUE_CAMERA) and expOD[1] == 'b'): continue
                frameIDs.append("%s%s-%s" % (m.groups()))
            else:
                raise RuntimeError('getFrameIDs: invalid exposure ID: %s' % expID)
        return frameIDs

class Fiber(object):
    Min = 1
    Max = 1000
    def __init__(self, fiber):
        self.fiber = int(fiber)
        assert self.fiber >= Fiber.Min
        assert self.fiber <= Fiber.Max
    def spectrographNumber(self):
        if self.fiber <= 500:
            return 1
        else:
            return 2
    def __str__(self):
        return '%d' % self.fiber
    def __call__(self):
        return self.fiber
    def __int__(self):
        return self.fiber

class MJD(object):
    Min = 50000
    Max = 60000
    def __init__(self, mjd):
        "Initialize MJD object."
        self.mjd = int(mjd)
        assert self.mjd >= MJD.Min
        assert self.mjd <= MJD.Max
    def __str__(self):
        return '%d' % self.mjd
    def __call__(self):
        return self.mjd
    def __int__(self):
        return self.mjd

class Plate(object):
    Min = 1
    Max = 10000
    def __init__(self, plate):
        "Initialize Plate object."
        self.plate = int(plate)
        assert self.plate >= Plate.Min
        assert self.plate <= Plate.Max
    def __str__(self):
        return '%d' % self.plate
    def __call__(self):
        return self.plate
    def __int__(self):
        return self.plate

class Target(object):
    def __init__(self):
        "Initialize Target object."
        self.plate = 0
        self.mjd = 0
        self.fiber = 0

    @classmethod
    def fromTuple(cls, target):
        "Initialize Target object."
        myTarget = cls()
        myTarget.plate = target[0]
        myTarget.mjd = target[1]
        myTarget.fiber = target[2]

        return myTarget

    @classmethod
    def fromPlateMJDFiber(cls, plate, mjd, fiber):
        "Initialize Target object."
        return cls.fromTuple((plate,mjd,fiber))

    @classmethod
    def fromTargetString(cls, targetString):
        "Initialize Target object from target string."
        target = targetString.strip().split("-")
        return cls.fromTuple(target)

    def isInitialized(self):
        return (self.plate != 0) and (self.mjd != 0) and (self.fiber != 0)

    def getPlate(self):
        return Plate(self.plate)

    def getMJD(self):
        return MJD(self.mjd)

    def getFiber(self):
        return Fiber(self.fiber)

    def __str__(self):
        "Human readable target string."
        return '%s-%s-%s' % (self.plate,self.mjd,self.fiber)

    def __call__(self):
        return (self.getPlate(),self.getMJD(),self.getFiber())


class TargetSet(object):
    def __init__(self):
        self.targets = []
        self.index = -1

    def __iter__(self):
        return self

    def next(self):
        if self.index < len(self.targets)-1:
            self.index += 1
            return self.targets[self.index]
        else:
            self.index = -1
            raise StopIteration

    def readTargetSet(self,filename):
        """
        Reads targets as "plate-mjd-fiber" lines from the specified file name and appends
        them to the specified target set. Each target should be on its own line, input
        following the target on each line is ignored.
        """
        count = 0
        targetfile = open(filename)
        # Target regex pattern
        p = re.compile(r'([0-9]+)-(5[0-9]{4})-([0-9]{1,4})')
        for line in targetfile:
            # Look for a match
            m = p.match(line)
            if m:
                # Add this new target to our set
                self.targets.append(Target.fromTargetString(m.group()))
                count += 1
            else:
                raise RuntimeError('readTargetSet: invalid target on line: %d' % (count+1))
        if len(self.targets) != len(set(self.targets)):
            raise RuntimeError('readTargetSet: Duplicate targets found!')
        targetfile.close()
        return count

    def saveTargetSetSpecList(self,filename):
        outfile = open(filename,'w')
        for target in self.targets:
            outfile.write('%4d/spec-%4d-%5d-%04d.fits\n' % (target.getPlate(),target.getPlate(),target.getMJD(),target.getFiber()))
        outfile.close()

class PathFinder(object):
    def __init__(self, boss_root, boss_version):
        self.boss_root = boss_root
        self.boss_version = boss_version

    def spPlate(self, target):
        return os.path.join(self.boss_root, self.boss_version, target.plate, 'spPlate-%d-%d.fits' % (target.getPlate(), target.getMJD()))

    def spCFrame(self, target, frameID):
        return os.path.join(self.boss_root, self.boss_version, target.plate, 'spCFrame-%s.fits' % frameID)

    def spAll(self, extension='fits'):
        return os.path.join(self.boss_root, self.boss_version, 'spAll-%s.%s' % (self.boss_version, extension))

    def spec(self, target):
        return os.path.join(self.boss_root, self.boss_version, 'spectra', '%04d' % target.getPlate(), 'spec-%04d-%5d-%04d.fits' % target())


