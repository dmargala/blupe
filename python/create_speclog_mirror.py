#!/usr/bin/env python
#
# Created 11-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# Modify speclog and copy relevant files to new dir


# Import libraries
import numpy as np
import os
import shutil
import yanny
import argparse

def swap_plugmap_stds(speclog_from, speclog_to, plate, night, mapname, fiberids, dry_run, clobber, verbose):
    # build path to plPlugMagM for unique mapname
    plugmap_name = os.path.join(speclog_from, night, 'plPlugMapM-%s.par' % mapname)
    new_plugmap_name = os.path.join(speclog_to, night, 'plPlugMapM-%s.par' % mapname)

    # Check to see if file exists already
    # yanny will not overwrite files so remove this by hand, in case we are running
    if os.path.isfile(new_plugmap_name):
        if clobber:
            os.remove(new_plugmap_name)
        else: 
            print 'Modified plugmap file already exists: %s'% new_plugmap_name
            print ' use --clobber option to overwrite...' 
            return -1

    # read the plugmap file
    try:
        plugmap = yanny.yanny(plugmap_name)
    except:
        print 'Could not open plugmap file: %s' % plugmap_name
        return -1
    plugmap_objs = plugmap['PLUGMAPOBJ']
    # hide the 'standard' boss spectrophoto standards
    boss_stds = np.where(np.array(plugmap_objs['objType']) == 'SPECTROPHOTO_STD')[0]
    if verbose and len(boss_stds) != 20:
        print 'Found %d boss standards on plate %s (expected 20)' % (len(boss_stds), plate)
    for i in boss_stds:
        plugmap_objs['objType'][i] = 'NA'
    # set the ancillary targets' objType as "SPECTROPHOTO_STD" by iterating 
    # over plugmap objects and matching fiberids to ancillary target fiberids
    # also make sure to check holetype is correct (guide stars and light traps have 
    # overlapping fiberids with the boss 'science' targets)
    nblue_stds = 0
    for i in range(len(plugmap_objs['objType'])):
        if str(plugmap_objs['fiberId'][i]) in set(fiberids) and plugmap_objs['holeType'][i] == 'OBJECT':
            # We found a blue standard, change its OBJTYPE
            plugmap_objs['objType'][i] = 'SPECTROPHOTO_STD'
            nblue_stds += 1
    if verbose and nblue_stds != len(fiberids):
        print 'Found %d blue standards on plate %s (expected %d)' % (nblue_stds, plate, len(fiberids))
    # save the plug map for this night's observation
    if verbose:
        print 'Saving modified plugmap file: %s' % new_plugmap_name
    plugmap.set_filename(new_plugmap_name)
    if not dry_run:
        plugmap.write()

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v","--verbose", action = "store_true",
        help="provide more verbose output")
    parser.add_argument("-i","--input", type=str, default="",
        help="input plate, mjd list")
    parser.add_argument("--clobber", action="store_true",
        help="overwrite existing files")
    parser.add_argument("--dry-run", action="store_true",
        help="dry run, no copies")
    parser.add_argument("--speclog-from", type=str, default="/home/boss/products/NULL/speclog/trunk",
        help="source speclog directory")
    parser.add_argument("--bsr-run2d", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/v5_7_0",
        help="path to reduction to mirror")
    parser.add_argument("--speclog-to", type=str, default="/clusterfs/riemann/raid008/bosswork/boss/spectro/redux/test/dmargala/speclog",
        help="destination speclog directory to write modified plugmap files to")
    parser.add_argument("--work-dir", type=str, default="",
        help="Working directory")
    args = parser.parse_args()

    # Read list of plate, mjd pairs to process
    plate_mjd_file = open(args.input, 'r')
    plate_mjd_list = []
    for line in plate_mjd_file:
        plate, mjd = line.split()
        plate_mjd_list.append((plate, mjd))

    if args.verbose:
        print 'Read %d from file: %s' % (len(plate_mjd_list), args.input)
        print 'Mirroring reduction at: %s' % args.bsr_run2d
        print 'Creating modified speclog with with blue standards and spectrophoto standards swapped...'
        print ' source speclog: %s' % args.speclog_from
        print ' new speclog: %s' % args.speclog_to

    def copy_speclog_file(name):
        from_name = os.path.join(args.speclog_from, name)
        to_name = os.path.join(args.speclog_to, name)
        if os.path.isfile(from_name):
            if not args.clobber and os.path.isfile(to_name):
                raise ValueError('Destination file already exists: %s' % to_name)
            if args.verbose:
                print os.path.join(to_name)
            if not args.dry_run:
                shutil.copy(from_name, to_name)
        else:
            raise ValueError('Source file does not exist: %s' % from_name)

    copy_speclog_file('spPlateList.par')

    # copy opfiles
    if not os.path.exists(os.path.join(args.speclog_to, 'opfiles')):
        os.makedirs(os.path.join(args.speclog_to, 'opfiles'))
    opfiles = ['spPlateMinSN2.par', 'spPlateZrange.par', 'spPlateList.par']
    for opfile in opfiles:
        copy_speclog_file(os.path.join('opfiles', opfile))

    for plate_mjd_pair in sorted(plate_mjd_list):
        (plate, mjd) = plate_mjd_pair
        # read blue standards list
        target_file = os.path.join(args.work_dir, plate, 'blue-stds-%s-%s.txt' % plate_mjd_pair)
        fiberids = []
        with open(target_file,'r') as targetlist:
            for line in targetlist:
                plate_mjd_fiber = line.strip().split('-')
                assert (plate == plate_mjd_fiber[0]) and (mjd == plate_mjd_list[1])
                fiberids.append(plate_mjd_fiber[2])
        # build path to spPlan for current plate-mjd
        plan_name = os.path.join(args.bsr_run2d, plate, 'spPlancomb-%s-%s.par' % plate_mjd_pair)
        # read spPlan file
        plan = yanny.yanny(plan_name)
        # get mjd list of observations for this plate-mjd
        unique_mapnames = set(zip(plan['SPEXP']['mjd'], plan['SPEXP']['mapname']))

        # The plugmaps are saved in files corresponding to the actual 
        # observation nights, not necessarily the last night (which is part of the target identifier)
        for night, mapname in unique_mapnames:
            night = str(night)
            if verbose:
                print plate, mjd, night, mapname
            if not os.path.exists(os.path.join(args.speclog_to, night)):
                os.makedirs(os.path.join(args.speclog_to, night))
            # copy sdHdrFix file over
            copy_speclog_file(os.path.join(night, 'sdHdrFix-%s.par' % night))
            # copy guiderMon file over
            copy_speclog_file(os.path.join(night, 'guiderMon-%s.par' % night))
            # copy plPlugMapM file and swap standards' OBJTYPE
            swap_plugmap_stds(args.speclog_from, args.speclog_to, plate, night, mapname, fiberids, args.dry_run, args.clobber, args.verbose)


if __name__ == '__main__':
    main()
    
