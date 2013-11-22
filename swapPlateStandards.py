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

# Set up appropriate paths
# os.environ['BOSS_SPECTRO_REDUX']
BSR = '/clusterfs/riemann/raid006/bosswork/boss/spectro/redux'

# os.environ['SPECLOG_DIR']
origSpecLogDir = '/home/boss/products/NULL/speclog/trunk'
newSpecLogDir = os.path.join('/data/dmargala/speclog')

bossVersion = 'v5_6_0'
blueVersion = 'test/dmargala/redux/v5_6_5'

targetListFilename = '/home/dmargala/blue-standards.txt'

verbose = True
force = True

def parseTargetList(targetListFilename):
    # Read in list of plates we are interested in
    bluePlateMJDPairDict = {}
    try:
        targetList = open(targetListFilename,'r')
    except:
        print 'Could not open target list file: %s' % targetListFilename
        return -1
    nTargets = 0
    for line in targetList:
        words = line.strip().split('-')
        if len(words) != 3:
            print 'Invalid target: %s' % line
            return -1
        (plate,mjd,fiberid) = (words[0],words[1],words[2])
        plateMJDPair = (plate,mjd)
        if not (plateMJDPair in bluePlateMJDPairDict.keys()):
            bluePlateMJDPairDict[plateMJDPair] = []
        bluePlateMJDPairDict[plateMJDPair].append(fiberid)
        nTargets += 1
    targetList.close()
    if verbose:
        "Successfully read a total of %d targets on %d plate-mjd pairs" % (nTargets, len(bluePlateMJDPairDict.keys()))  
    return bluePlateMJDPairDict

def swapStandards(bluePlateMJDPairDict):
    for plateMJDPair in sorted(bluePlateMJDPairDict.keys()):
        (plate, mjd) = plateMJDPair
        # Set path to spPlan for current plate-mjd
        spPlanFilename = os.path.join(BSR,bossVersion,plate,'spPlancomb-%s-%s.par'%(plate,mjd))
        # Read spPlan file
        spPlan = yanny.yanny(spPlanFilename)
        # Get mjd list of observations for this plate-mjd
        print plate, mjd
        print spPlanFilename
        uniqueMapNames = set(zip(spPlan['SPEXP']['mjd'],spPlan['SPEXP']['mapname']))
        #print len(bluePlateMJDPairDict[plateMJDPair]) # debug
        #print plateMJDPair # debug
        # The plugmaps are saved in files corresponding to the actual 
        # observation nights, not necessarily the last night
        for night, mapname in uniqueMapNames:
            #print night # debug
            night = str(night)
            print plate, mjd, night, mapname
            # do stuff for unique mjd path
            plPlugMapFilename = os.path.join(origSpecLogDir,night,'plPlugMapM-%s.par'%mapname)
            #print plPlugMapFilename # debug
            try:
                plPlugMap = yanny.yanny(plPlugMapFilename)
            except:
                print "Could not open plPlugMapFilename: %s" % plPlugMapFilename
                return -1
            plugMapOBJs = plPlugMap['PLUGMAPOBJ']
            # Hide the 'standard' boss spectrophoto standards
            bossStandards = np.where(np.array(plugMapOBJs['objType']) == 'SPECTROPHOTO_STD')[0]
            if verbose and len(bossStandards) != 20:
                print 'Hey! There are %d boss standards on plate %s' % (len(bossStandards), plate)
            for i in bossStandards:
                plugMapOBJs['objType'][i] = 'NA'
            ### Set blue standards' objType as "SPECTROPHOTO_STD"
            nPlugMapObjs = len(plugMapOBJs['objType'])
            # Iterate over plugmap objects and match fiberids to blue standard fiberids
            # also make sure check holetype is correct (guide stars and light traps have 
            # overlapping fiberids with the boss 'science' targets)
            nBlueStandards = 0
            for i in range(nPlugMapObjs):
                if str(plugMapOBJs['fiberId'][i]) in set(bluePlateMJDPairDict[plateMJDPair]) \
                and plugMapOBJs['holeType'][i] == 'OBJECT':
                    # We found a blue standard!! 
                    nBlueStandards += 1
                    plugMapOBJs['objType'][i] = 'SPECTROPHOTO_STD'
            if verbose:
                print 'Hey! We found %d blue standards on plate %s' % (nBlueStandards, plate)
            # Save the plug map for this night's observation
            if not os.path.exists(os.path.join(newSpecLogDir,night)):
                os.makedirs(os.path.join(newSpecLogDir,night))
            newFilename = os.path.join(newSpecLogDir,night,'plPlugMapM-%s.par'%mapname)
            # yanny will not overwrite files so let's remove this by hand, in case we are running
            if force and os.path.isfile(newFilename):
                os.remove(newFilename)
            elif not force and os.path.isfile(newFilename):
                print 'New plugmap file (%s) already exists run with force option set to overwrite...' % newFilename
                return -1
            plPlugMap.set_filename(newFilename)
            plPlugMap.write()
            # Copy sdHdrFix file and any other files that need to be copied over
            sdHdrFixFilename = os.path.join(origSpecLogDir,night,'sdHdrFix-%s.par'%night)
            if os.path.isfile(sdHdrFixFilename):
                shutil.copy(sdHdrFixFilename,os.path.join(newSpecLogDir,night,'sdHdrFix-%s.par'%night))
            guiderMonFilename = os.path.join(origSpecLogDir,night,'guiderMon-%s.par'%night)
            if os.path.isfile(guiderMonFilename):
                shutil.copy(guiderMonFilename,os.path.join(newSpecLogDir,night,'guiderMon-%s.par'%night))

if __name__ == "__main__":
    bluePlateMJDPairDict = parseTargetList(targetListFilename=targetListFilename)
    # Iterate over plate-mjd pairs and create new plugmap with 
    # blue standards swapped with spectrophoto standards
    try:
        plateMJDListFile = open('bluePlateMJDList.txt','w')
    except IOError:
        print "Could not open output file for plate-MJD list."
    finalPlateMJDPairDict = {}
    minBlueStandards = 10
    for plateMJDPair in sorted(bluePlateMJDPairDict.keys()):
        (plate, mjd) = plateMJDPair
        # We are only interested in plates with at least 10 blue standards that survived target selection
        if len(bluePlateMJDPairDict[plateMJDPair]) < minBlueStandards:
            continue
        # Write plate,mjd pair to output file
        plateMJDListFile.write('%s %s\n' % plateMJDPair)
        print ('%s %s' % plateMJDPair) + (' %d' % len(bluePlateMJDPairDict[plateMJDPair]))
        finalPlateMJDPairDict[plateMJDPair] = bluePlateMJDPairDict[plateMJDPair]
    plateMJDListFile.close()

    swapStandards(bluePlateMJDPairDict=finalPlateMJDPairDict)
    
