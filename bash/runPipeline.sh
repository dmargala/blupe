export MJD=55445
export PLATESTART=3615
export PLATEEND=3615

idl -e "spplan2d, mjd='$MJD', /clobber"
idl -e "spplan1d, mjd='$MJD', /clobber"
idl -e "batchpbs, mjd='$MJD', /clobber, /zcode"
idl -e "platelist, /create, run2d='$RUN2D', run1d='$RUN1D'"
idl -e "platemerge, run2d='$RUN2D'"
