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
parser.add_option('-i', "--input", type="string", help="Input target list.")
parser.add_option('-o', "--output", type="string", help="Output filename.")

opts, args = parser.parse_args()

targetFile = open(opts.input)

outfile = open(opts.output,'w')

for line in targetFile:
	(plate,mjd,fiber) = line.strip().split('-')

	plate = int(plate)
	mjd = int(mjd)
	fiber = int(fiber)

	outfile.write('%04d/spec-%04d-%5d-%04d.fits\n' % (plate,plate,mjd,fiber))



#'%d' % target.getPlate(), 'spec-%4d-%5d-%04d.fits' % target())