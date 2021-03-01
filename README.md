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

# Analysis_mc.cc
this is the "main" where all the other classes are called



