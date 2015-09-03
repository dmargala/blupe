Generate tpcorr
###############

* Create plate, mjd list for observations with offset targets: `work/plate-mjd-lambdaeff-4000.txt`
* Run `expsure_query.py` all plate, mjd pairs:

.. code-block::

    cat work/plate-mjd-lambdaeff-4000.txt | \
    awk 'print "python python/exposure_query.py --plate "$1" --mjd "$2" --outdir /data/dmargala/tpcorr/"$1""' > bash/run_info_new.sh

* add header to `bash/run_info_new.sh`:

.. code-block::

    #PBS -l nodes=1
    #PBS -l walltime=48:00:00
    #PBS -W umask=0022
    #PBS -V
    #PBS -j oe
    cd $PBS_O_WORKDIR

* Submit to queue: `qsub bash/run_info_new.sh`
* Check invalid azimuth calculation messages: fixed, floating point problem, \texttt{arccos(1.00001) -> NaN}
* Check sidereal time error: fixed, problem with astropy and conversion of recent times
* Collate results:

.. code-block::

    cat /data/dmargala/tpcorr/????/expinfo-summary-*.txt | \
    awk 'NR % 2 == 0' > expinfo-summary-all.txt

* make plot of `psf_fwhm` distribution: medisn psf value is 1.49

* Fix 79 plates with `psf_fwhm == 0`:

.. code-block::

    awk '{if($6 == 0) {$6=1.49}; print}' work/expinfo-summary-all.txt > work/expinfo-summary-all-fixed.txt

* Run `plate_tpcorr.pro` on DR12:

.. code-block::

    cat ../work/expinfo-summary-all-fixed.txt | \
    awk '{print "idl -e \"plate_tpcorr, "$1", "$2", "$3", '\''"$4"'\'', "$5", "$6", \
    '\'~/blupe/work/tpcorr/tpcorr-'"$1"-"$2"'.txt\''\""}' > run_tpcorr_new.sh
    split -l 100 -d run_tpcorr_new.sh; \
    for x in $(ls -1 x??); \
        do cat qsub_header.txt $x > $x.sh; chmod +x $x.sh; rm $x; qsub $x.sh; \
    done

* Run `fit_tpcorr.py` on `plate_tpcorr.pro` results:

.. code-block::

    for f in $(ls -1 work/tpcorr/tpcorr-*.txt);
        do echo ./python/fit_tpcorr.py --skip-plots -i $f -o work/tpcorr/fit-tpcorr --scipy;
    done > run_fit_tpcorr.sh
    split -l 100 -d run_fit_tpcorr.sh;
    for x in $(ls -1 x??);
        do cat pro/qsub_header.txt $x > $x.sh; chmod +x $x.sh; rm $x; qsub $x.sh;
    done

* Summarize model fit results:

.. code-block::

    ./python/plot_fit_tpcorr_results.py -i 'work/tpcorr/fit-tpcorr-*txt' -o work/fit-tpcorr-summary.png

* Convert `plate_tpcorr.pro` results to release format: 539 targets with `fiberid == -1`

.. code-block::

    ./python/convert_tpcorr_results.py --input 'work/tpcorr/tpcorr-*.txt' --output release.hdf5
