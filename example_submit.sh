cd /storage_mnt/storage/user/mvit/CMSSW_9_4_4
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
cd /user/mvit/CMSSW_9_4_4/src/HNL_analysis
make -f makeFiles/makeAnalisi
./analisi_hnl
