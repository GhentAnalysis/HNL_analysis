#!/bin/bash

hnl="majorana"
#hnl="dirac"

# use for SR
sr="/user/kskovpen/analysis/HNLAnalysis/CMSSW_9_4_4/src/HNL_analysis/SR/${hnl}"
cr="/user/kskovpen/analysis/HNLAnalysis/CMSSW_9_4_4/src/HNL_analysis/SR/${hnl}/CR"

# use for Fig. 2
#sr="/user/kskovpen/analysis/HNLAnalysis/CMSSW_9_4_4/src/HNL_analysis/test"
#cr="/user/kskovpen/analysis/HNLAnalysis/CMSSW_9_4_4/src/HNL_analysis/test/CR"

chan=("ele" "muon")
#year=("16" "17" "18" "Run2")
year=("Run2")
mode="M" # SR, M, D, M3

for ch in ${chan[@]}
do

  if [[ ${ch} == "ele" ]]; then
    chname="e_ele"
    if [[ ${hnl} == "dirac" ]]; then
      chname="e_Dirac_ele"
    fi
  else
    chname="mu_muo"
    if [[ ${hnl} == "dirac" ]]; then
      chname="mu_Dirac_muo"
    fi
  fi
  
  for y in ${year[@]}
  do
    # use for SR figures
    python \
    readDatacard.py \
    ${sr}/M-2.0_V-0.0089106678_${chname}_datacard \
    ${sr}/M-6.0_V-0.0011224972_${chname}_datacard \
    ${sr}/M-12.0_V-0.0010000000_${chname}_datacard \
    ${sr}/M-12.0_V-0.0010000000_${chname}_datacard \
    ${ch} ${mode} 2 6 12 0_8=4 1_3=6 1_0=6 ${y} \
    ${cr}

    # use for Fig. 2
#    coup="0.00707813534767"
#    python \
#    readDatacard.py \
#    ${sr}/M-3_V-${coup}_${chname}_datacard \
#    ${sr}/M-3_V-${coup}_${chname}_datacard \
#    ${sr}/M-3_V-${coup}_${chname}_datacard \
#    ${sr}/M-3_V-${coup}_${chname}_datacard \
#    ${ch} ${mode} 2 6 12 0_8=4 1_3=6 1_0=6 ${y} \
#    ${cr}
  done
  
done
