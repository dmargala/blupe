#!/usr/bin/env python
#
# Created 12-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# Import libraries
import os
import shutil
import fnmatch
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v","--verbose", action="store_true",
        help="provide more verbose output")
    parser.add_argument("-i","--input", type=str, default="bluePlateMJDList.txt",
        help="input list plate-mjd pairs to copy")
    parser.add_argument("--speclog", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/test/dmargala/speclog",
        help="speclog dir containing modified plugmap files")
    parser.add_argument("--bsr-from", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux",
        help="path to original reduction")
    parser.add_argument("--bsr-to", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/test/dmargala/redux",
        help="path to new reduction")
    parser.add_argument("--run2d", type=str, default="v5_7_0",
        help="run2d version")
    parser.add_argument("--qsub", action="store_true",
        help="submit jobs to queue")
    parser.add_argument("--dry-run", action="store_true",
        help="dry run")
    parser.add_argument("--clobber", action="store_true",
        help="overwrite existing files")
    args = parser.parse_args()

    path_from = os.path.join(args.bsr_from, args.run2d)
    path_to = os.path.join(args.bsr_to, args.run2d)
    if args.verbose:
        print 'Creating mirror of %s at %s' % (path_from, path_to)

    # read in list of plates to copy
    for line in open(args.input,'r'):
        (plate, mjd) = line.strip().split()
        
        plate_path_from = os.path.join(path_from, plate)
        plate_path_to = os.path.join(path_to, plate)
        if not os.path.exists(plate_path_to):
            os.makedirs(plate_path_to)

        # file patterns to copy over to new working directory
        patterns = ['spFlat*', 'spArc*', 'spFrame*', 'photo*', '*.par', 'redux-%s-%s' % (plate, mjd)]
        for filename in os.listdir(plate_path_from):
            if any(fnmatch.fnmatch(filename, p) for p in patterns):
                fullname_from = os.path.join(plate_path_from, filename)
                fullname_to = os.path.join(plate_path_to, filename)

                if not args.clobber and os.path.isfile(fullname_to):
                    raise IOError('File already exists, remove it or run with --clobber: %s' % fullname_to)
                if args.verbose:
                    print fullname_to
                if not args.dry_run:
                    shutil.copy(fullname_from, fullname_to)

        # change into working directory
        os.chdir(plate_path_to)
        # modify the submission script so it doesn't overwrite custom env variables
        reduxfilename = 'redux-%s-%s' % (plate, mjd)
        if args.verbose:
            print 'Modifying submission script: %s' % reduxfilename
        with open(reduxfilename,'r') as reduxfile:
            output = []
            for reduxline in reduxfile:
                output.append(reduxline)
                if 'setup' in reduxline:
                    output.append('export SPECLOG_DIR=%s\n' % args.speclog)
        # save modified submission script
        if not args.dry_run:
            with open(reduxfilename, 'w') as reduxfile:
                reduxfile.writelines(output)

        # submit job if specified
        if args.verbose:
            subprocess.call(['pwd'])
            subprocess.call(['echo','qsub','-q','batch','redux-%s-%s' % (plate, mjd)])
        if args.qsub and not args.dry_run:
            subprocess.call(['qsub','-q','batch','redux-%s-%s' % (plate, mjd)])

if __name__ == '__main__':
    main()
