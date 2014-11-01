#!/usr/bin/env python
#
# Created 11-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# Import libraries
import argparse
import os

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v","--verbose", action="store_true",
        help="provide more verbose output")
    parser.add_argument("-i","--input", type=str, default="",
        help= "input list of ancillary program targets (plate-mjd-fiber)")
    parser.add_argument("--min-stds", type=int, default=10,
        help="Minimum number of ancillary targets observed to consider plate")
    parser.add_argument("--work-dir", type=str, default="",
        help="Working directory")

    args = parser.parse_args()

    with open(args.input,'r') as targetlist:
        # create a dictionary of observed ancillary targets, keyed by plate-mjd pairs
        nTargets = 0
        bluePlateMJDPairDict = {}
        for line in targetlist:
            plateMJDFiber = line.strip().split('-')
            if len(plateMJDFiber) != 3:
                print 'Invalid target: %s' % line
                return -1
            (plate, mjd, fiberid) = plateMJDFiber
            plateMJDPair = (plate, mjd)
            # If we haven't seen this plate-mjd yet, create a new list entry
            if not (plateMJDPair in bluePlateMJDPairDict.keys()):
                bluePlateMJDPairDict[plateMJDPair] = []
            # Add this target to the plate-mjd entry
            bluePlateMJDPairDict[plateMJDPair].append(int(fiberid))
            nTargets += 1

    if args.verbose:
        print "Read a total of %d targets on %d plate-mjd pairs from: %s" % (nTargets, len(bluePlateMJDPairDict.keys()), args.input)  

    # create a trimmed dictionary containing plate-mjd pairs that have enough ancillary targets
    trimmedPlateMJDPairDict = {}
    for plateMJDPair, fiberids in sorted(bluePlateMJDPairDict.iteritems()):
        count2 = sum(i > 500 for i in fiberids)
        count1 = len(fiberids) - count2
        if (count1 >= args.min_stds) and (count2 >= args.min_stds):
            trimmedPlateMJDPairDict[plateMJDPair] = fiberids
            if args.verbose:
                print plateMJDPair, count1, count2

    def savePerPlateFibers(plate, mjd, fiberids):
        plate_dir = os.path.join(args.work_dir, plate)
        if not os.path.exists(plate_dir):
            os.makedirs(plate_dir)
        with open(os.path.join(plate_dir, 'blue-stds-%s-%s.txt' % (plate, mjd)), 'w') as output:
            for fiberid in fiberids:
                output.write('%s-%s-%s\n' % (plate, mjd, fiberid))

    # save plate, mjd pair list and per plate target lists
    if not os.path.exists(args.work_dir):
        os.makedirs(args.work_dir)
    plate_mjd_list_name = os.path.join(args.work_dir, 'blue-plate-mjd-list.txt')
    if args.verbose:
        print 'Saving plate, mjd list to: %s' % plate_mjd_list_name
    with open(plate_mjd_list_name, 'w') as output:
        for plateMJDPair, fiberids in sorted(trimmedPlateMJDPairDict.iteritems()):
            (plate, mjd) = plateMJDPair
            output.write('%s %s\n' % plateMJDPair)
            savePerPlateFibers(plate, mjd, fiberids)


if __name__ == "__main__":
    main()
    
