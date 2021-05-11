#!/bin/bash

sr="martinaCards"
cr="../dataCards_shapeRoot_CR"

python \
readDatacard.py \
${sr}/M-1_V-0.022383_mu_muo_datacard \
${sr}/M-2_V-0.00447214_mu_muo_datacard \
${sr}/M-3_V-0.0022383_mu_muo_datacard \
${sr}/M-3_V-0.0022383_mu_muo_datacard \
muon SR 1GeV 2GeV 3GeV 1_6=4 1_8=5 5_2=6 16 \
${cr}/M-1_V-0.212367605816_mu_muo_datacard
