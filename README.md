# HNL_analysis

This repository has the code to make SR plots-tables and produce dataCards and shapeFile. 

To compile:
           
           make -f makeFiles/makeAnalisi

to run:
           
           ./analisi_hnl

# analisi.C
To set: year and flags regarding the decision about: which samples to run on, which txt use, which samples' folder, produce limits, produce plots. 

    unsigned year = 0;  // 2016: 0; 2017: 1; 2018: 2;
    //                              skipData, skipSignal, skipBackground, skipPlotting, skipLimits
    all.analisi(basename.c_str(), false    , false     , false          , false        , false    );

--> Analysis_mc.cc:     this is the "main" where all the other classes are called

# Re-weighting procedure, ONLY on m2, using condor. 

/submit_re_weighting_condor/sum_1*_*.sh

it produces ~100 data/shape output files for each mass in /dataCards_shapeRoot. They run over only signal samples so the obs/nonpromptDF/SF/Xgamma are ==0. 
They have to be merged --> mergeDatacardsRootfiles.py with the complete output files. 



# making_plots

To make plots and tables 

Input files should be in martinaCards/

example:             
                              
                 python readDatacard_martina.py martinaCards/M-1_V-0.022383_mu_muo_datacard martinaCards/M-2_V-0.00447214_mu_muo_datacard martinaCards/M-3_V-0.0022383_mu_muo_datacard martinaCards/M-3_V-0.0022383_mu_muo_datacard muon SR 1GeV 2GeV 3GeV 1_6=4 1_8=5 5_2=6 16
                 (4 signal samples
                 muon/ele coupling
                 SR,M,D type of distribution
                 1_6=4 coupling for legend --> 1.6 *10^-4)
