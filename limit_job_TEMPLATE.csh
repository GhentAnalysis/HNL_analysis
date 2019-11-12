#!/bin/csh

set label='TEMPLABEL'
set cmsdir='/storage_mnt/storage/user/trocino/Analysis/HNL/Displaced/2019-03-11_LegacyRelease/CMSSW_10_2_9/src/'
set wrkdir='TEMPDIR'
set curdir=`pwd`
cd $cmsdir
source /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scramv1 runtime -csh`
cd $wrkdir
make -f jobs/makefiles/make_${label}
jobs/executables/analisi_${label}
touch jobs/completed/done_${label}
