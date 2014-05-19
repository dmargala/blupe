#!/usr/bin/env python
#
# Created 11-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# Import libraries
import numpy as np
import os
import shutil
import yanny
import argparse

def swapStandards(bluePlateMJDPairDict, speclog_src, speclog_new, bsr_run2d, verbose=False, clobber=False):
    for plateMJDPair in sorted(bluePlateMJDPairDict.keys()):
        (plate, mjd) = plateMJDPair
        # build path to spPlan for current plate-mjd
        spPlanFilename = os.path.join(bsr_run2d, plate, 'spPlancomb-%s-%s.par' % plateMJDPair)
        # read spPlan file
        spPlan = yanny.yanny(spPlanFilename)
        # get mjd list of observations for this plate-mjd
        uniqueMapNames = set(zip(spPlan['SPEXP']['mjd'],spPlan['SPEXP']['mapname']))
        # The plugmaps are saved in files corresponding to the actual 
        # observation nights, not necessarily the last night
        for night, mapname in uniqueMapNames:
            night = str(night)
            if verbose:
                print plate, mjd, night, mapname
            # build path to plPlugMagM for unique mapname
            plPlugMapFilename = os.path.join(speclog_src, night, 'plPlugMapM-%s.par' % mapname)
            # read the plugmap file
            try:
                plPlugMap = yanny.yanny(plPlugMapFilename)
            except:
                print "Could not open plPlugMapFilename: %s" % plPlugMapFilename
                return -1
            plugMapOBJs = plPlugMap['PLUGMAPOBJ']
            # hide the 'standard' boss spectrophoto standards
            bossStandards = np.where(np.array(plugMapOBJs['objType']) == 'SPECTROPHOTO_STD')[0]
            if verbose and len(bossStandards) != 20:
                print 'Found %d boss standards on plate %s (expected 20)' % (len(bossStandards), plate)
            for i in bossStandards:
                plugMapOBJs['objType'][i] = 'NA'
            # set the ancillary targets' objType as "SPECTROPHOTO_STD" by iterating 
            # over plugmap objects and matching fiberids to ancillary target fiberids
            # also make sure to check holetype is correct (guide stars and light traps have 
            # overlapping fiberids with the boss 'science' targets)
            nBlueStandards = 0
            for i in range(len(plugMapOBJs['objType'])):
                if str(plugMapOBJs['fiberId'][i]) in set(bluePlateMJDPairDict[plateMJDPair]) and plugMapOBJs['holeType'][i] == 'OBJECT':
                    # We found a blue standard!! 
                    nBlueStandards += 1
                    plugMapOBJs['objType'][i] = 'SPECTROPHOTO_STD'
            if verbose and  nBlueStandards != len(bluePlateMJDPairDict[plateMJDPair]):
                print 'Found %d blue standards on plate %s (expected %d)' % (nBlueStandards, plate, len(bluePlateMJDPairDict[plateMJDPair]))
            # save the plug map for this night's observation
            if not os.path.exists(os.path.join(speclog_new, night)):
                os.makedirs(os.path.join(speclog_new, night))
            newPlugMapFilename = os.path.join(speclog_new, night, 'plPlugMapM-%s.par' % mapname)
            if verbose:
                print newPlugMapFilename
            # yanny will not overwrite files so remove this by hand, in case we are running
            if os.path.isfile(newPlugMapFilename):
                if clobber:
                    os.remove(newPlugMapFilename)
                else: 
                    print 'Modified plugmap file already exists: %s'% newPlugMapFilename
                    print ' specify --clobber option to overwrite...' 
                    continue
            plPlugMap.set_filename(newPlugMapFilename)
            plPlugMap.write()
            # copy sdHdrFix file over
            sdHdrFixFilename = os.path.join(night, 'sdHdrFix-%s.par' % night)
            if os.path.isfile(os.path.join(speclog_src, sdHdrFixFilename)):
                if verbose:
                    print os.path.join(speclog_new, sdHdrFixFilename)
                shutil.copy(os.path.join(speclog_src, sdHdrFixFilename), os.path.join(speclog_new, sdHdrFixFilename))
            # copy guiderMon file over
            guiderMonFilename = os.path.join(night, 'guiderMon-%s.par' % night)
            if os.path.isfile(os.path.join(speclog_src, guiderMonFilename)):
                if verbose:
                    print os.path.join(speclog_new, guiderMonFilename)
                shutil.copy(os.path.join(speclog_src, guiderMonFilename), os.path.join(speclog_new, guiderMonFilename))
    # copy speclog/opfiles
    if not os.path.exists(os.path.join(speclog_new, "opfiles")):
        os.makedirs(os.path.join(speclog_new, "opfiles"))
    shutil.copy(os.path.join(speclog_src, "opfiles", "spPlateMinSN2.par"), os.path.join(speclog_new, "opfiles", "spPlateMinSN2.par"))
    shutil.copy(os.path.join(speclog_src, "opfiles", "spPlateZrange.par"), os.path.join(speclog_new, "opfiles", "spPlateZrange.par"))
    shutil.copy(os.path.join(speclog_src, "opfiles", "spPlateList.par"), os.path.join(speclog_new, "opfiles", "spPlateList.par"))
    shutil.copy(os.path.join(speclog_src, "spPlateList.par"), os.path.join(speclog_new, "spPlateList.par"))

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v","--verbose", action = "store_true",
        help = "provide more verbose output")
    parser.add_argument("-i","--input", type=str, default ="/home/dmargala/blue-standards.txt",
        help = "input list of ancillary program targets (plate-mjd-fiber)")
    parser.add_argument("-o","--output", type=str, default="bluePlateMJDList.txt",
        help = "save plate mjd list to specified file")
    parser.add_argument("--min-stds", type=int, default=10,
        help = "Minimum number of ancillary targets observed to consider plate")
    parser.add_argument("--clobber", action = "store_true",
        help = "overwrite existing files")
    parser.add_argument("--speclog-src", type=str, default="/home/boss/products/NULL/speclog/trunk",
        help = "source speclog directory")
    parser.add_argument("--speclog-new", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/test/dmargala/speclog",
        help = "new speclog directory to write modified plugmap files to")
    parser.add_argument("--bsr-run2d", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/v5_7_0",
        help = "path to reduction to mirror")
    parser.add_argument("--no-swap", action="store_true",
        help = "skip modification of plugmap and file copies, only create filtered plate-mjd list")
    args = parser.parse_args()

    # open the observed ancillary target list
    try:
        targetListFile = open(args.input,'r')
    except IOError:
        print 'Could not open target list file: %s' % args.input
        return -1
    # create a dictionary of observed ancillary targets, keyed by plate-mjd pairs
    nTargets = 0
    bluePlateMJDPairDict = {}
    for line in targetListFile:
        plateMJDFiber = line.strip().split('-')
        if len(plateMJDFiber) != 3:
            print 'Invalid target: %s' % line
            return -1
        (plate,mjd,fiberid) = plateMJDFiber
        plateMJDPair = (plate,mjd)
        # If we haven't seen this plate-mjd yet, create a new list entry
        if not (plateMJDPair in bluePlateMJDPairDict.keys()):
            bluePlateMJDPairDict[plateMJDPair] = []
        # Add this target to the plate-mjd entry
        bluePlateMJDPairDict[plateMJDPair].append(fiberid)
        nTargets += 1
    targetListFile.close()
    if args.verbose:
        "Successfully read a total of %d targets on %d plate-mjd pairs" % (nTargets, len(bluePlateMJDPairDict.keys()))  

    # Iterate over plate-mjd pairs and create new plugmap with 
    # blue standards swapped with spectrophoto standards
    try:
        plateMJDListFile = open(args.output,'w')
    except IOError:
        print "Could not open output file for plate-mjd list: %s" % args.output
        return -1
    # create a trimmed dictionary containing plate-mjd pairs that have enough ancillary targets
    trimmedPlateMJDPairDict = {}
    for plateMJDPair in sorted(bluePlateMJDPairDict.keys()):
        numFound = len(bluePlateMJDPairDict[plateMJDPair])
        if numFound < args.min_stds:
            continue
        # write plate,mjd pair to output file
        plateMJDListFile.write('%s %s\n' % plateMJDPair)
        if args.verbose:
            print '%s %s %d' % (plateMJDPair[0], plateMJDPair[1], numFound)
        trimmedPlateMJDPairDict[plateMJDPair] = bluePlateMJDPairDict[plateMJDPair]
    plateMJDListFile.close()

    if not args.no_swap:
        # Modify speclog and copy relevant files to new dir
        if args.verbose:
            print "Mirroring reduction at: %s" % args.bsr_run2d
            print "Creating modified speclog with with blue standards and spectrophoto standards swapped..."
            print " source speclog: %s" % args.speclog_src
            print " new speclog: %s" % args.speclog_new
        swapStandards(trimmedPlateMJDPairDict, speclog_src=args.speclog_src, 
            speclog_new=args.speclog_new, bsr_run2d=args.bsr_run2d, verbose=args.verbose, clobber=args.clobber)

if __name__ == "__main__":
    main()
    
