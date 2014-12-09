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
import cPickle as pickle
import os

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("--boss-root", type="string", help="Specify boss root")
parser.add_option("--boss-version", type="string", help="Specify boss version")
parser.add_option('-i', "--input", type="string", help="Input target list.")
parser.add_option('-o', "--output", type="string", help="Output filename.")
parser.add_option('-j', "--jnput", type="string", help="DR8 Targetlist")

opts, args = parser.parse_args()

if not opts.input:
    print "Must specify an input target list to read!"
    sys.exit(1)

if not (opts.boss_version and opts.boss_root):
    print "Must specify boss_root and boss_version!"
    sys.exit(1)

#outfile = open(opts.output,'w')

targetFile = open(opts.input)

finder = lya.PathFinder(opts.boss_root,opts.boss_version)
dr8finder = lya.PathFinder(opts.boss_root,'dr8')
bluefinder = lya.PathFinder(opts.boss_root,'test/dmargala/redux/v5_6_5')


numPlatesOpened = 0
numTargetsRead = 0

plateFile = ''

targetSet = lya.TargetSet()
numTargetsRead = lya.TargetSet.readTargetSet(targetSet, opts.input)

otherTargetSet = lya.TargetSet()
otherTargetsRead = lya.TargetSet.readTargetSet(otherTargetSet, opts.jnput)
camera = lya.BLUE_CAMERA

targetmap = {}
force = False

if (not os.path.isfile("targetmap.p")) or force:
    for target in targetSet:
        if int(target()[0]) >= 6307:
            continue    

        myspec = lya.spec(finder.spec(target))

        numExposures = myspec.getNumExposures()
        spectroflux = myspec.spAll('SPECTROFLUX')
        ra = myspec.spAll('PLUG_RA')
        dec = myspec.spAll('PLUG_DEC')

        # myrun = myspec.spAll('RUN')
        # myrerun = myspec.spAll('RERUN')
        # mycamcol = myspec.spAll('CAMCOL')
        # myfield = myspec.spAll('FIELD')
        # myid = myspec.spAll('ID')

        print len(otherTargetSet.targets)
        for other in otherTargetSet:
            otherspec = lya.spec(dr8finder.spec(other))
            otherRA = otherspec.spAll('PLUG_RA')
            otherDec = otherspec.spAll('PLUG_DEC')

            # otherobjid = otherspec.spAll('OBJID')
            # otherrun = int(otherobjid[0])
            # otherrerun = int(otherobjid[1])
            # othercamcol = int(otherobjid[2])
            # otherfield = int(otherobjid[3])
            # otherid = int(otherobjid[4])
            #if myrun == otherrun and myrerun == otherrerun and mycamcol == othercamcol \
            #and myfield == otherfield and myid == otherid:
            if np.abs(otherRA-ra) < .05 and np.abs(otherDec-dec) < .05:
                print "\t", other, otherRA, otherDec
                otherTargetSet.targets.remove(other)

                if str(target) in targetmap.keys():
                    targetmap[str(target)].append(other)
                else:
                    targetmap[str(target)] = [other]

    pickle.dump(targetmap, open("targetmap.p","wb"))
else:
    print "opening pickle!"
    targetmap = pickle.load(open("targetmap.p","rb"))


writeOut = True

if writeOut:
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

        if int(target()[0]) >= 6307:
            continue

        # Standard BOSS Spectrum
        myspec = lya.spec(finder.spec(target))
        numExposures = myspec.getNumExposures()
        #psfmag = myspec.spAll('PSFMAG')
        spectroflux = myspec.spAll('SPECTROFLUX')
        ra = myspec.spAll('PLUG_RA')
        dec = myspec.spAll('PLUG_DEC')
        xfocal = myspec.spAll('XFOCAL')
        yfocal = myspec.spAll('YFOCAL')
        lameff = myspec.spAll('LAMBDA_EFF')

        flux = myspec.combined("flux")
        loglam = myspec.combined("loglam")

        print target, ra, dec, xfocal, yfocal, lameff, numExposures #myrun, myrerun, mycamcol, myfield, myid

        # Open output file for this target
        outfile = open("%s.dat" % target,"w")

        # Save standard BOSS combined spectrum data
        outfile.write("combined %d " % lameff)
        for index in range(len(loglam)):
            outfile.write("%f %f " % (loglam[index],flux[index]))
        outfile.write("\n")
        # Save standard BOSS exposure spectrum data
        frames = myspec.getFrames()
        for frame in frames:
            frameflux = frame['flux']
            frameloglam = frame['loglam']
            outfile.write("exp %d " % lameff)
            for index in range(len(frameloglam)):
                outfile.write("%f %f " % (frameloglam[index],frameflux[index]))
            outfile.write("\n")

        # Blue Reduction Spectrum
        bluespec = lya.spec(bluefinder.spec(target))
        blueflux = bluespec.combined("flux")
        blueloglam = bluespec.combined("loglam")
        # Save blue reduction combined spectrum data
        outfile.write("combined %d " % lameff)
        for index in range(len(blueloglam)):
            outfile.write("%f %f " % (blueloglam[index],blueflux[index]))
        outfile.write("\n")
        # Save blue reduction exposure spectrum data
        for frame in bluespec.getFrames():
            frameflux = frame["flux"]
            frameloglam = frame["loglam"]
            outfile.write("exp %d " % lameff)
            for index in range(len(frameloglam)):
                outfile.write("%f %f " % (frameloglam[index],frameflux[index]))
            outfile.write("\n")

        # DR8 Spectrum
        others = []
        if str(target) in targetmap.keys():
            others = targetmap[str(target)]
        for other in others:
            otherspec = lya.spec(dr8finder.spec(other))
            otherRA = otherspec.spAll('PLUG_RA')
            otherDec = otherspec.spAll('PLUG_DEC')

            print "\t", other, otherRA, otherDec

            otherflux = otherspec.combined("flux")
            otherloglam = otherspec.combined("loglam")

            # Save DR8 combined spectrum
            outfile.write("combined %d " % lameff)
            for index in range(len(otherloglam)):
                outfile.write("%f %f " % (otherloglam[index],otherflux[index]))
            outfile.write("\n")
            # Save DR8 exposure spectra
            for frame in otherspec.getFrames():
                frameflux = frame["flux"]
                frameloglam = frame["loglam"]
                outfile.write("exp %d " % lameff)
                for index in range(len(frameloglam)):
                    outfile.write("%f %f " % (frameloglam[index],frameflux[index]))
                outfile.write("\n")

        outfile.close()

