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

# Set up appropriate paths
origBSRPath = '/clusterfs/riemann/raid006/bosswork/boss/spectro/redux/v5_6_0'
newBSRPath = '/data/dmargala/redux/test/dmargala/redux/v5_6_5'
speclog_dir = '/data/dmargala/speclog'

force = False

def copyPipelineFiles():
    # Read in list of plates we are interested in
    for line in open('bluePlateMJDList.txt','r'):
        (plate, mjd) = line.strip().split()
        origDir = os.path.join(origBSRPath,plate)
        newDir = os.path.join(newBSRPath,plate)
        patterns = ['spFlat*','spArc*','spFrame*','photo*','*.par','redux-%s-%s'%(plate,mjd)]
        #patterns = ['*.par','redux-%s-%s'%(plate,mjd)]
        if not os.path.exists(newDir):
            os.makedirs(newDir)
        for filename in os.listdir(origDir):
            if any(fnmatch.fnmatch(filename, p) for p in patterns):
                if not force and os.path.isfile(os.path.join(newDir,filename)):
                    #print 'The following file already exists, please remove it or run with force to continue:'
                    #print os.path.join(newDir,filename)
                    continue
                shutil.copy(os.path.join(origDir,filename),os.path.join(newDir,filename)) 
        os.chdir(newDir)
        reduxfilename = 'redux-%s-%s'%(plate,mjd)
        reduxfile = open(reduxfilename,'r')
        output = []

        for reduxline in reduxfile:
            if "setup" in reduxline:
                output.append('export SPECLOG_DIR=%s\n'%speclog_dir)
            elif 'v5_6_0' in reduxline:
                parts = reduxline.split("\"")
                newparts = []
                for part in parts:
                    if part == 'v5_6_0':
                        newparts.append("v5_6_5")
                    else:
                        newparts.append(part)
                output.append("\"".join(newparts))
            else:
                output.append(reduxline)
        reduxfile.close()
        reduxfile = open(reduxfilename, 'w')
        reduxfile.writelines(output)
        reduxfile.close()

        subprocess.call(['pwd'])
        subprocess.call(['echo','qsub','redux-%s-%s'%(plate,mjd)])
        #subprocess.call(['qsub','redux-%s-%s'%(plate,mjd)])

if __name__ == '__main__':
    copyPipelineFiles()
