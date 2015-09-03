Modify and rerun pipeline
#########################

Modify the pipeline input files so that the ancillary targets are labeled as spectrophotometric standards. Run the pipeline on the ancillary plates. Results are available at blahblahblah.

* Copy speclog and opfiles files and to new speclog dir
* Swap spectrophotometric standards in plugmap files
* Copy pipeline input files `['spFlat*', 'spArc*', 'spFrame*', 'photo*', '*.par']` into new redux dir
* Set modified speclog and redux environment vars `SPECLOG_DIR` and `BOSS_SPECTRO_REDUX`
* Copy and modify the job submission script so it doesn't overwrite custom env variables

Compare magnitudes
##################

# pre selection
./examples/filter.py -i /dm/all/data/boss/spAll-v5_7_0.fits -s "['OBJTYPE'] == 'SPECTROPHOTO_STD'" --verbose --output compare_tpcorr/spAll-tpcorr-stds.fits --plates compare_tpcorr/tpcorr_plates.txt

# save annotated target list
./examples/filter.py -i compare_tpcorr/spAll-tpcorr-stds.fits --verbose --annotate 'psfmag' --save compare_tpcorr/tpcorr-stds.txt

# extract columns of interest
cat compare_tpcorr/tpcorr-stds.txt | cut -d' ' -f 3,4,5 > compare_tpcorr/std-mags.txt

# calculate synthetic magnitudes
./examples/calc_mags.py --targets compare_tpcorr/tpcorr-stds.txt --output compare_tpcorr/std-syn-mags.txt

# copy to local system
scp hpc:~/source/qusp/compare_tpcorr/std*mags.txt compare_tpcorr/

# compare mags
./examples/compare_mags.py -y compare_tpcorr/std-mags.txt -x compare_tpcorr/std-syn-mags.txt --output compare_tpcorr/std-syn-vs-psf

########## offset stds
./examples/filter.py -i /dm/all/data/boss/spAll-v5_7_0.fits -s "['ANCILLARY_TARGET2'] == 1<<20" --verbose --output compare_tpcorr/spAll-tpcorr-offset-stds.fits --plates compare_tpcorr/tpcorr_plates.txt

./examples/filter.py -i compare_tpcorr/spAll-tpcorr-offset-stds.fits --verbose --annotate 'psfmag' --save compare_tpcorr/tpcorr-offset-stds.txt

./examples/calc_mags.py --targets compare_tpcorr/tpcorr-offset-stds.txt --output compare_tpcorr/offset-std-syn-mags.txt

./examples/calc_mags.py --targets compare_tpcorr/tpcorr-offset-stds.txt --output compare_tpcorr/offset-std-syn-mags-anc.txt --boss-version 'test/dmargala/redux/v5_7_0'

./examples/calc_mags.py --targets compare_tpcorr/tpcorr-offset-stds.txt --output compare_tpcorr/offset-std-syn-mags-tpcorr.txt --tpcorr 'data/tpcorr/tpcorr*'

scp hpc:~/source/qusp/compare_tpcorr/offset-std*mags*.txt compare_tpcorr/

##########

./examples/compare_mags.py -y compare_tpcorr/std-mags.txt -x compare_tpcorr/std-syn-mags.txt --output compare_tpcorr/std-syn-vs-psf
./examples/compare_mags.py -y compare_tpcorr/offset-std-mags.txt -x compare_tpcorr/offset-std-syn-mags.txt --output compare_tpcorr/offset-std-syn-vs-psf
./examples/compare_mags.py -y compare_tpcorr/offset-std-mags.txt -x compare_tpcorr/offset-std-syn-mags-anc.txt --output compare_tpcorr/offset-std-syn-anc-vs-psf
./examples/compare_mags.py -y compare_tpcorr/offset-std-mags.txt -x compare_tpcorr/offset-std-syn-mags-tpcorr.txt --output compare_tpcorr/offset-std-syn-tpcorr-vs-psf

###########################
###########################

./examples/filter.py --verbose -i /dm/all/data/boss/spAll-v5_7_0.fits -s "(['LAMBDA_EFF'] == 4000) & ((['ZWARNING'] & 1<<7) == 0)" --save data/target_lists/offset-targets.txt --output data/target_lists/offset-targets.fits
./examples/filter.py --verbose -i /dm/all/data/boss/spAll-v5_7_0.fits -s "(['OBJTYPE'] == 'SPECTROPHOTO_STD') & (['CLASS'] == 'STAR') & (['LAMBDA_EFF'] == 5400) & ((['ZWARNING'] & 1<<7) == 0)" --save data/target_lists/specphoto-stds.txt --output data/target_lists/specphoto-stds.fits
./examples/filter.py --verbose -i data/target_lists/offset-targets.fits -s "(['OBJTYPE'] == 'QSO') & (['CLASS'] == 'QSO') & (['LAMBDA_EFF'] == 4000) & ((['ZWARNING'] & 1<<7) == 0)" --save data/target_lists/quasars.txt --output data/target_lists/quasars.fits
./examples/filter.py --verbose -i data/target_lists/offset-targets.fits -s "(['OBJTYPE'] == 'QSO') & (['CLASS'] == 'STAR') & (['LAMBDA_EFF'] == 4000) & ((['ZWARNING'] & 1<<7) == 0)" --save data/target_lists/failed-quasars.txt --output data/target_lists/failed-quasars.fits
./examples/filter.py --verbose -i data/target_lists/offset-targets.fits -s "(['ANCILLARY_TARGET2'] == 1<<20) & (['CLASS'] == 'STAR') & (['LAMBDA_EFF'] == 4000) & ((['ZWARNING'] & 1<<7) == 0)" --save data/target_lists/offset-stds.txt --output data/target_lists/offset-stds.fits
