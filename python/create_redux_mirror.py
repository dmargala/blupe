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
    parser.add_argument("--skip-existing", action="store_true",
        help="skip if file already exits")
    args = parser.parse_args()

    path_from = os.path.join(args.bsr_from, args.run2d)
    path_to = os.path.join(args.bsr_to, args.run2d)
    if args.verbose:
        print 'Creating mirror of %s at %s' % (path_from, path_to)

    # file patterns to copy over to new working directory
    patterns = ['spFlat*', 'spArc*', 'spFrame*', 'photo*', '*.par']

    # read in list of plate-mjd observations to mirror
    for line in open(args.input,'r'):
        (plate, mjd) = line.strip().split()

        plate_path_from = os.path.join(path_from, plate)
        plate_path_to = os.path.join(path_to, plate)
        if not os.path.exists(plate_path_to):
            os.makedirs(plate_path_to)

        # copy pipeline input files
        for filename in os.listdir(plate_path_from):
            if any(fnmatch.fnmatch(filename, p) for p in patterns):
                fullname_from = os.path.join(plate_path_from, filename)
                fullname_to = os.path.join(plate_path_to, filename)

                if args.verbose:
                    print fullname_to

                if os.path.isfile(fullname_to):
                    if args.skip_existing:
                        continue
                    elif not args.clobber:
                        raise IOError('File already exists, remove it or run with --clobber')

                if not args.dry_run:
                    shutil.copy(fullname_from, fullname_to)

        # modify the submission script so it doesn't overwrite custom env variables
        redux_from = os.path.join(plate_path_from, 'redux-%s-%s' % (plate, mjd))
        redux_to = os.path.join(plate_path_to, 'redux-%s-%s' % (plate, mjd))
        with open(redux_from, 'r') as reduxfile:
            output = []
            for reduxline in reduxfile:
                output.append(reduxline)
                if 'setup' in reduxline:
                    output.append('export SPECLOG_DIR=%s\n' % args.speclog)
        # save modified submission script
        if args.verbose:
            print 'Saving modified submission script: %s' % redux_to
        if not args.dry_run:
            with open(redux_to, 'w') as reduxfile:
                reduxfile.writelines(output)

        # submit job if specified
        os.chdir(plate_path_to)
        if args.verbose:
            subprocess.call(['pwd'])
            subprocess.call(['echo','qsub','-q','batch','redux-%s-%s' % (plate, mjd)])
        if args.qsub and not args.dry_run:
            subprocess.call(['qsub','-q','batch','redux-%s-%s' % (plate, mjd)])

if __name__ == '__main__':
    main()
