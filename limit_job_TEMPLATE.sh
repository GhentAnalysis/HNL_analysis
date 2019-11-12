#!/bin/bash                                                                                                                                           

label='TEMPLABEL'
cmsdir='/storage_mnt/storage/user/trocino/Analysis/HNL/Displaced/2019-03-11_LegacyRelease/CMSSW_10_2_9/src/'
wrkdir='TEMPDIR'
curdir=`pwd`
cd $cmsdir
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $wrkdir
make -f jobs/makefiles/make_${label}
jobs/executables/analisi_${label}
touch jobs/completed/done_${label}
