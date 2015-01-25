#!/bin/bash

TPCORRFILE=$1
TPCORRDIR=$2

if [ -d $2 ]
then
    echo "Error: specified directory "$TPCORRDIR" already exists."
    exit 1
fi

mkdir $TPCORRDIR

echo "= Improved Calibration of the SDSS-III Quasar Sample =\n\n" > $TPCORRDIR/README.txt
echo "For description on how to use the correction, see http://darkmatter.ps.uci.edu/tpcorr/\n" >> $TPCORRDIR/README.txt

cp $TPCORRFILE $TPCORRDIR/
cp python/example_usage.py $TPCORRDIR/example_usage.py

tar czfv $TPCORRDIR.tar.gz $TPCORRDIR

rm -r $TPCORRDIR