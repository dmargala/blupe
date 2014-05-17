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
    parser.add_argument("-v","--verbose", action = "store_true",
        help = "provide more verbose output")
    parser.add_argument("-i","--input", type=str, default="bluePlateMJDList.txt",
        help = "input list plate-mjd pairs to copy")
    parser.add_argument("--clobber", action = "store_true",
        help = "overwrite existing files")
    parser.add_argument("--speclog", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/test/dmargala/speclog",
        help = "speclog dir containing modified plugmap files")
    parser.add_argument("--bsr-orig", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux",
        help = "path to original reduction")
    parser.add_argument("--bsr-new", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/test/dmargala/redux",
        help = "path to new reduction")
    parser.add_argument("--run2d", type=str, default="v5_7_0",
        help = "run2d version")
    parser.add_argument("--qsub", action="store_true",
        help = "submit jobs to queue")
    args = parser.parse_args()

    origDir = os.path.join(args.bsr_orig,args.run2d)
    newDir = os.path.join(args.bsr_new,args.run2d)
    if args.verbose:
        print 'Creating mirror of %s at %s' % (origDir, newDir)

    # read in list of plates to copy
    for line in open(args.input,'r'):
        (plate, mjd) = line.strip().split()
        origDirPlate = os.path.join(origDir,plate)
        newDirPlate = os.path.join(newDir,plate)
        # file patterns to copy over to new working directory
        patterns = ['spFlat*','spArc*','spFrame*','photo*','*.par','redux-%s-%s'%(plate,mjd)]
        if not os.path.exists(newDirPlate):
            os.makedirs(newDirPlate)
        for filename in os.listdir(origDirPlate):
            if any(fnmatch.fnmatch(filename, p) for p in patterns):
                if not args.clobber and os.path.isfile(os.path.join(newDirPlate,filename)):
                    print 'The following file already exists, please remove it or run with --clobber to continue:'
                    print os.path.join(newDirPlate,filename)
                    continue
                print os.path.join(newDirPlate,filename)
                shutil.copy(os.path.join(origDirPlate,filename),os.path.join(newDirPlate,filename)) 
        # change into working directory
        os.chdir(newDirPlate)
        # Modify the submission script so it doesn't overwrite custom env variables
        reduxfilename = 'redux-%s-%s'%(plate,mjd)
        reduxfile = open(reduxfilename,'r')
        output = []
        for reduxline in reduxfile:
            output.append(reduxline)
            if "setup" in reduxline:
                output.append('export SPECLOG_DIR=%s\n'%speclog_dir)
        reduxfile.close()
        reduxfile = open(reduxfilename, 'w')
        reduxfile.writelines(output)
        reduxfile.close()
        # print working directory and qsub command
        if args.verbose:
            subprocess.call(['pwd'])
            subprocess.call(['echo','qsub','-q','batch','redux-%s-%s'%(plate,mjd)])
        # submit job if specified
        if args.qsub:
            subprocess.call(['qsub','-q','batch','redux-%s-%s'%(plate,mjd)])

if __name__ == '__main__':
    main()
