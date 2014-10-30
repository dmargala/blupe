#!/usr/bin/env python
#
# Created 19-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage:

# Import Libraries
import bosslya as lya
import optparse
import sys
import numpy as np

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("--boss-root", type="string", help="Specify boss root")
parser.add_option("--boss-version", type="string", help="Specify boss version")
parser.add_option('-i', "--input", type="string", help="Input target list.")
parser.add_option('-o', "--output", type="string", help="Output filename.")

opts, args = parser.parse_args()

if not opts.input:
    print "Must specify an input target list to read!"
    sys.exit(1)

if not (opts.boss_version and opts.boss_root):
    print "Must specify boss_root and boss_version!"
    sys.exit(1)

outfile = open(opts.output,'w')

targetFile = open(opts.input)

finder = lya.PathFinder(opts.boss_root,opts.boss_version)

numPlatesOpened = 0
numTargetsRead = 0

plateFile = ''

targetSet = lya.TargetSet()
numTargetsRead = lya.TargetSet.readTargetSet(targetSet, opts.input)

camera = lya.BLUE_CAMERA

# loop over targets
for target in targetSet:
    # Open this target's combined spectrum file.
    # if finder.spPlate(target) != plateFile:
    #     plateFile = finder.spPlate(target)
    #     myPlate = spPlate(plateFile)
    #     numPlatesOpened += 1
    # frameIDs = myPlate.getFrameIDs(target.getFiber(),camera)
    # for exposure in range(0,len(frameIDs)):
    #     print frameIDs[exposure],
    # print

    myspec = lya.spec(finder.spec(target))
    numExposures = myspec.getNumExposures()
    psfmag = myspec.spAll('PSFMAG')
    spectroflux = myspec.spAll('SPECTROFLUX')

    # spectromag = 22.5 - 2.5*np.log10(spectroflux)

    # myobjid = int(myspec.spAll('OBJID'))

    # print int(bin(myobjid)[-48:-32],2), int(bin(myobjid)[-59:-48],2), int(bin(myobjid)[-32:-29],2), int(bin(myobjid)[-28:-16],2), int(bin(myobjid)[-16:],2)

    myrun = myspec.spAll('RUN')
    myrerun = myspec.spAll('RERUN')
    mycamcol = myspec.spAll('CAMCOL')
    myfield = myspec.spAll('FIELD')
    myid = myspec.spAll('ID')

    #print myrun, myrerun, mycamcol, myfield, myid
    outfile.write('%s %s %s %s %s\n' % (myrun, myrerun, mycamcol, myfield, myid))
    # for mag in psfmag: 
    #     if not np.isnan(mag):
    #         outfile.write('%.4f ' % mag)
    #     else:
    #         outfile.write('%.4f ' % 0)
    #         print "invalid value encounted on target %s: %s" % (target,spectroflux)
    # for mag in spectromag:
    #     if not np.isnan(mag):
    #         outfile.write('%.4f ' % mag)
    #     else:
    #         outfile.write('%.4f ' % 0)
    #         print "invalid value encounted on target %s: %s" % (target,spectroflux)
    # outfile.write('\n')

outfile.close()

print 'Number of targets read: %d' % numTargetsRead
print 'Number of plates files opened: %d' % numPlatesOpened


print 'Done!'
