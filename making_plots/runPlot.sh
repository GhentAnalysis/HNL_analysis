#!/bin/bash

sr="/user/kskovpen/analysis/HNL/CMSSW_10_2_13/src/Limits/input/3L/majorana"
cr="/user/kskovpen/analysis/HNL/CMSSW_10_2_13/src/Limits/input/3L/majorana/CR"

python \
readDatacard.py \
${sr}/M-2.0_V-0.0089106678_mu_muo_datacard \
${sr}/M-5.0_V-0.0044721360_mu_muo_datacard \
${sr}/M-9.0_V-0.0011224972_mu_muo_datacard \
${sr}/M-9.0_V-0.0008910668_mu_muo_datacard \
muon SR 2GeV 5GeV 9GeV 8_0=5 2_0=5 1_2=6 Run2 \
${cr}/M-2.0_V-0.0089106678_mu_muo_datacard

#python \
#readDatacard.py \
#${sr}/M-1_V-0.022383_mu_muo_datacard \
#${sr}/M-2_V-0.00447214_mu_muo_datacard \
#${sr}/M-3_V-0.0022383_mu_muo_datacard \
#${sr}/M-3_V-0.0022383_mu_muo_datacard \
#muon SR 1GeV 2GeV 3GeV 1_6=4 1_8=5 5_2=6 16 \
#${cr}/M-1_V-0.212367605816_mu_muo_datacard
