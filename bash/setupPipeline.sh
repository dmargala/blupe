#- Setup pipeline test environment
export VERS=v5_6_5

setup -r $HOME/svn/idlspec2d/$VERS idlspec2d
export BOSS_SPECTRO_REDUX=/data/dmargala/spectro/redux/test
export RUN1D=$VERS
export RUN2D=$VERS

