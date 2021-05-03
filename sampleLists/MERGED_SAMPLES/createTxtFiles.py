import sys, os, math
#import ROOT as r
#from array import array

file_xs_ctau = [
    ['1.0', -5, +0, [
        ['HeavyNeutrino_trilepton_M-1.0_V-0.022360679775_e_massiveAndCKM_LO_2016_MERGED.root'     , '2.130e+00', '1276.0800'],
        ['HeavyNeutrino_trilepton_M-1.0_V-0.022360679775_mu_massiveAndCKM_LO_2016_MERGED.root'    , '2.130e+00', '1276.0800'],
        ['HeavyNeutrino_trilepton_M-1.0_V-0.022360679775_e_massiveAndCKM_LO_2017_MERGED.root'     , '2.130e+00', '1276.0800'],
        ['HeavyNeutrino_trilepton_M-1.0_V-0.022360679775_mu_massiveAndCKM_LO_2017_MERGED.root'    , '2.130e+00', '1276.0800'],
        ['HeavyNeutrino_trilepton_M-1.0_V-0.022360679775_e_massiveAndCKM_LO_2018_MERGED.root'     , '2.130e+00', '1276.0800'],
        ['HeavyNeutrino_trilepton_M-1.0_V-0.022360679775_mu_massiveAndCKM_LO_2018_MERGED.root'    , '2.130e+00', '1276.0800']
    ]],
    ['2.0', -6, +0, [ 
        ['HeavyNeutrino_trilepton_M-2.0_V-0.0110905365064_e_massiveAndCKM_LO_2016_MERGED.root'    , '5.270e-01',  '156.2400'],
        ['HeavyNeutrino_trilepton_M-2.0_V-0.0110905365064_mu_massiveAndCKM_LO_2016_MERGED.root'   , '5.270e-01',  '156.2400'],
        ['HeavyNeutrino_trilepton_M-2.0_V-0.0110905365064_e_massiveAndCKM_LO_2017_MERGED.root'    , '5.270e-01',  '156.2400'],
        ['HeavyNeutrino_trilepton_M-2.0_V-0.0110905365064_mu_massiveAndCKM_LO_2017_MERGED.root'   , '5.270e-01',  '156.2400'],
        ['HeavyNeutrino_trilepton_M-2.0_V-0.0110905365064_e_massiveAndCKM_LO_2018_MERGED.root'    , '5.270e-01',  '156.2400'],
        ['HeavyNeutrino_trilepton_M-2.0_V-0.0110905365064_mu_massiveAndCKM_LO_2018_MERGED.root'   , '5.270e-01',  '156.2400']
    ]],
    ['3.0', -6, +0, [
        ['HeavyNeutrino_trilepton_M-3.0_V-0.00244948974278_e_massiveAndCKM_LO_2016_MERGED.root'   , '2.412e-02',  '380.3636'],
        ['HeavyNeutrino_trilepton_M-3.0_V-0.00244948974278_mu_massiveAndCKM_LO_2016_MERGED.root'  , '2.412e-02',  '380.3636'],
        ['HeavyNeutrino_trilepton_M-3.0_V-0.00707813534767_e_massiveAndCKM_LO_2017_MERGED.root'   , '2.014e-01',   '45.5525'],
        ['HeavyNeutrino_trilepton_M-3.0_V-0.00707813534767_mu_massiveAndCKM_LO_2017_MERGED.root'  , '2.014e-01',   '45.5525'],
        ['HeavyNeutrino_trilepton_M-3.0_V-0.00707813534767_e_massiveAndCKM_LO_2018_MERGED.root'   , '2.014e-01',   '45.5525'],
        ['HeavyNeutrino_trilepton_M-3.0_V-0.00707813534767_mu_massiveAndCKM_LO_2018_MERGED.root'  , '2.014e-01',   '45.5525']
    ]],
    ['4.0', -7, -1, [
        ['HeavyNeutrino_trilepton_M-4.0_V-0.0013_e_massiveAndCKM_LO_2016_MERGED.root'             , '6.707e-03',  '287.0136'],
        ['HeavyNeutrino_trilepton_M-4.0_V-0.0013_mu_massiveAndCKM_LO_2016_MERGED.root'            , '6.707e-03',  '287.0136'],
        ['HeavyNeutrino_trilepton_M-4.0_V-0.00290516780927_e_massiveAndCKM_LO_2017_MERGED.root'   , '3.350e-02',   '57.4706'],
        ['HeavyNeutrino_trilepton_M-4.0_V-0.00290516780927_mu_massiveAndCKM_LO_2017_MERGED.root'  , '3.350e-02',   '57.4706'],
        ['HeavyNeutrino_trilepton_M-4.0_V-0.00290516780927_e_massiveAndCKM_LO_2018_MERGED.root'   , '3.350e-02',   '57.4706'],
        ['HeavyNeutrino_trilepton_M-4.0_V-0.00290516780927_mu_massiveAndCKM_LO_2018_MERGED.root'  , '3.350e-02',   '57.4706']
    ]],
    ['5.0', -7, -1, [ 
        ['HeavyNeutrino_trilepton_M-5.0_V-0.000316227766017_e_massiveAndCKM_LO_2016_MERGED.root'  , '3.981e-04', '1476.7352'],
        ['HeavyNeutrino_trilepton_M-5.0_V-0.000316227766017_mu_massiveAndCKM_LO_2016_MERGED.root' , '3.981e-04', '1476.7352'],
        ['HeavyNeutrino_trilepton_M-5.0_V-0.000316227766017_e_massiveAndCKM_LO_2017_MERGED.root'  , '3.981e-04', '1476.7352'],
        ['HeavyNeutrino_trilepton_M-5.0_V-0.000316227766017_mu_massiveAndCKM_LO_2017_MERGED.root' , '3.981e-04', '1476.7352'],
        ['HeavyNeutrino_trilepton_M-5.0_V-0.000316227766017_e_massiveAndCKM_LO_2018_MERGED.root'  , '3.981e-04', '1476.7352'],
        ['HeavyNeutrino_trilepton_M-5.0_V-0.000316227766017_mu_massiveAndCKM_LO_2018_MERGED.root' , '3.981e-04', '1476.7352']
    ]],
    ['6.0', -7, -1, [
        ['HeavyNeutrino_trilepton_M-6.0_V-0.000370135110466_e_massiveAndCKM_LO_2016_MERGED.root'  , '5.560e-04',  '409.0738'],
        ['HeavyNeutrino_trilepton_M-6.0_V-0.000370135110466_mu_massiveAndCKM_LO_2016_MERGED.root' , '5.560e-04',  '409.0738'],
        ['HeavyNeutrino_trilepton_M-6.0_V-0.000370135110466_e_massiveAndCKM_LO_2017_MERGED.root'  , '5.560e-04',  '409.0738'],
        ['HeavyNeutrino_trilepton_M-6.0_V-0.000370135110466_mu_massiveAndCKM_LO_2017_MERGED.root' , '5.560e-04',  '409.0738'],
        ['HeavyNeutrino_trilepton_M-6.0_V-0.000370135110466_e_massiveAndCKM_LO_2018_MERGED.root'  , '5.560e-04',  '409.0738'],
        ['HeavyNeutrino_trilepton_M-6.0_V-0.000370135110466_mu_massiveAndCKM_LO_2018_MERGED.root' , '5.560e-04',  '409.0738']
    ]],
    ['7.0', -7, -2, [
        ['HeavyNeutrino_trilepton_M-7.0_V-0.000316227766_e_massiveAndCKM_LO_2016_MERGED.root'     , '4.093e-04',  '249.9849'],
        ['HeavyNeutrino_trilepton_M-7.0_V-0.000316227766_mu_massiveAndCKM_LO_2016_MERGED.root'    , '4.093e-04',  '249.9849'],
        ['HeavyNeutrino_trilepton_M-7.0_V-0.000316227766_e_massiveAndCKM_LO_2017_MERGED.root'     , '4.093e-04',  '249.9849'],
        ['HeavyNeutrino_trilepton_M-7.0_V-0.000316227766_mu_massiveAndCKM_LO_2017_MERGED.root'    , '4.093e-04',  '249.9849'],
        ['HeavyNeutrino_trilepton_M-7.0_V-0.000316227766_e_massiveAndCKM_LO_2018_MERGED.root'     , '4.093e-04',  '249.9849'],
        ['HeavyNeutrino_trilepton_M-7.0_V-0.000316227766_mu_massiveAndCKM_LO_2018_MERGED.root'    , '4.093e-04',  '249.9849']
    ]],
    ['8.0', -7, -2, [
        ['HeavyNeutrino_trilepton_M-8.0_V-0.000151327459504_e_massiveAndCKM_LO_2016_MERGED.root'  , '9.396e-05',  '545.1831'],
        ['HeavyNeutrino_trilepton_M-8.0_V-0.000151327459504_mu_massiveAndCKM_LO_2016_MERGED.root' , '9.396e-05',  '545.1831'],
        ['HeavyNeutrino_trilepton_M-8.0_V-0.000316227766017_e_massiveAndCKM_LO_2017_MERGED.root'  , '4.096e-04',  '124.8462'],
        ['HeavyNeutrino_trilepton_M-8.0_V-0.000547722557505_mu_massiveAndCKM_LO_2017_MERGED.root' , '1.229e-03',   '41.7314'],
        ['HeavyNeutrino_trilepton_M-8.0_V-0.000316227766017_e_massiveAndCKM_LO_2018_MERGED.root'  , '4.096e-04',  '124.8462'],
        ['HeavyNeutrino_trilepton_M-8.0_V-0.000547722557505_mu_massiveAndCKM_LO_2018_MERGED.root' , '1.229e-03',   '41.7314']
    ]],
    ['9.0', -7, -2, [ 
        ['HeavyNeutrino_trilepton_M-9.0_V-0.000316227766_e_massiveAndCKM_LO_2016_MERGED.root'     , '4.133e-04',   '68.0172'],
        ['HeavyNeutrino_trilepton_M-9.0_V-0.000316227766_mu_massiveAndCKM_LO_2016_MERGED.root'    , '4.133e-04',   '68.0172'],
        ['HeavyNeutrino_trilepton_M-9.0_V-0.000316227766_e_massiveAndCKM_LO_2017_MERGED.root'     , '4.133e-04',   '68.0172'],
        ['HeavyNeutrino_trilepton_M-9.0_V-0.000316227766_mu_massiveAndCKM_LO_2017_MERGED.root'    , '4.133e-04',   '68.0172'],
        ['HeavyNeutrino_trilepton_M-9.0_V-0.000316227766_e_massiveAndCKM_LO_2018_MERGED.root'     , '4.133e-04',   '68.0172'],
        ['HeavyNeutrino_trilepton_M-9.0_V-0.000316227766_mu_massiveAndCKM_LO_2018_MERGED.root'    , '4.133e-04',   '68.0172']
    ]],
    ['10.0', -7, -3, [
        ['HeavyNeutrino_trilepton_M-10.0_V-7.56967634711e-05_e_massiveAndCKM_LO_2016_MERGED.root' , '2.363e-05',  '690.2373'],
        ['HeavyNeutrino_trilepton_M-10.0_V-7.56967634711e-05_mu_massiveAndCKM_LO_2016_MERGED.root', '2.363e-05',  '690.2373'],
        ['HeavyNeutrino_trilepton_M-10.0_V-0.000316227766017_e_massiveAndCKM_LO_2017_MERGED.root' , '4.118e-04',   '39.5504'],
        ['HeavyNeutrino_trilepton_M-10.0_V-0.000316227766017_mu_massiveAndCKM_LO_2017_MERGED.root', '4.118e-04',   '39.5504'],
        ['HeavyNeutrino_trilepton_M-10.0_V-0.000316227766017_e_massiveAndCKM_LO_2018_MERGED.root' , '4.118e-04',   '39.5504'],
        ['HeavyNeutrino_trilepton_M-10.0_V-0.000316227766017_mu_massiveAndCKM_LO_2018_MERGED.root', '4.118e-04',   '39.5504']
    ]],
    ['11.0', -7, -3, [
        ['HeavyNeutrino_trilepton_M-11.0_V-0.000316227766_e_massiveAndCKM_LO_2016_MERGED.root'    , '4.165e-04',   '24.2603'],
        ['HeavyNeutrino_trilepton_M-11.0_V-0.000316227766_mu_massiveAndCKM_LO_2016_MERGED.root'   , '4.165e-04',   '24.2603'],
        ['HeavyNeutrino_trilepton_M-11.0_V-0.000316227766_e_massiveAndCKM_LO_2017_MERGED.root'    , '4.165e-04',   '24.2603'],
        ['HeavyNeutrino_trilepton_M-11.0_V-0.000316227766_mu_massiveAndCKM_LO_2017_MERGED.root'   , '4.165e-04',   '24.2603'],
        ['HeavyNeutrino_trilepton_M-11.0_V-0.000316227766_e_massiveAndCKM_LO_2018_MERGED.root'    , '4.165e-04',   '24.2603'],
        ['HeavyNeutrino_trilepton_M-11.0_V-0.000316227766_mu_massiveAndCKM_LO_2018_MERGED.root'   , '4.165e-04',   '24.2603']
    ]],
    ['12.0', -7, -3, [
        ['HeavyNeutrino_trilepton_M-12.0_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.305e-03',    '0.7755'],
        ['HeavyNeutrino_trilepton_M-12.0_V-0.000836660026_mu_massiveAndCKM_LO_2016_MERGED.root'   , '2.913e-03',    '2.2204'],
        ['HeavyNeutrino_trilepton_M-12.0_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.305e-03',    '0.7755'],
        ['HeavyNeutrino_trilepton_M-12.0_V-0.000836660026_mu_massiveAndCKM_LO_2017_MERGED.root'   , '2.913e-03',    '2.2204'],
        ['HeavyNeutrino_trilepton_M-12.0_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.305e-03',    '0.7755'],
        ['HeavyNeutrino_trilepton_M-12.0_V-0.000836660026_mu_massiveAndCKM_LO_2018_MERGED.root'   , '2.913e-03',    '2.2204']
    ]],
    ['12.2', -7, -3, [
        ['HeavyNeutrino_trilepton_M-12.2_V-0.000836660026_mu_massiveAndCKM_LO_2016_MERGED.root'   , '2.907e-03',    '2.0371'],
        ['HeavyNeutrino_trilepton_M-12.2_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.300e-03',    '0.7130'],
        ['HeavyNeutrino_trilepton_M-12.2_V-0.000836660026_mu_massiveAndCKM_LO_2017_MERGED.root'   , '2.907e-03',    '2.0371'],
        ['HeavyNeutrino_trilepton_M-12.2_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.300e-03',    '0.7130'],
        ['HeavyNeutrino_trilepton_M-12.2_V-0.000836660026_mu_massiveAndCKM_LO_2018_MERGED.root'   , '2.907e-03',    '2.0371']
    ]],
    ['12.4', -7, -3, [
        ['HeavyNeutrino_trilepton_M-12.4_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.294e-03',    '0.6559'],
        ['HeavyNeutrino_trilepton_M-12.4_V-0.000836660026_mu_massiveAndCKM_LO_2016_MERGED.root'   , '2.902e-03',    '1.8747'],
        ['HeavyNeutrino_trilepton_M-12.4_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.294e-03',    '0.6559'],
        ['HeavyNeutrino_trilepton_M-12.4_V-0.000836660026_mu_massiveAndCKM_LO_2017_MERGED.root'   , '2.902e-03',    '1.8747'],
        ['HeavyNeutrino_trilepton_M-12.4_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.294e-03',    '0.6559'],
        ['HeavyNeutrino_trilepton_M-12.4_V-0.000836660026_mu_massiveAndCKM_LO_2018_MERGED.root'   , '2.902e-03',    '1.8747']
    ]],
    ['12.6', -7, -3, [
        ['HeavyNeutrino_trilepton_M-12.6_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.291e-03',    '0.6043'],
        ['HeavyNeutrino_trilepton_M-12.6_V-0.000836660026_mu_massiveAndCKM_LO_2016_MERGED.root'   , '2.902e-03',    '1.7269'],
        ['HeavyNeutrino_trilepton_M-12.6_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.291e-03',    '0.6043'],
        ['HeavyNeutrino_trilepton_M-12.6_V-0.000836660026_mu_massiveAndCKM_LO_2017_MERGED.root'   , '2.902e-03',    '1.7269'],
        ['HeavyNeutrino_trilepton_M-12.6_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.291e-03',    '0.6043'],
        ['HeavyNeutrino_trilepton_M-12.6_V-0.000836660026_mu_massiveAndCKM_LO_2018_MERGED.root'   , '2.902e-03',    '1.7269']
    ]],
    ['12.8', -7, -3, [
        ['HeavyNeutrino_trilepton_M-12.8_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.277e-03',    '0.5573'],
        ['HeavyNeutrino_trilepton_M-12.8_V-0.000836660026_mu_massiveAndCKM_LO_2016_MERGED.root'   , '2.897e-03',    '1.5926'],
        ['HeavyNeutrino_trilepton_M-12.8_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.277e-03',    '0.5573'],
        ['HeavyNeutrino_trilepton_M-12.8_V-0.000836660026_mu_massiveAndCKM_LO_2017_MERGED.root'   , '2.897e-03',    '1.5926'],
        ['HeavyNeutrino_trilepton_M-12.8_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.277e-03',    '0.5573'],
        ['HeavyNeutrino_trilepton_M-12.8_V-0.000836660026_mu_massiveAndCKM_LO_2018_MERGED.root'   , '2.897e-03',    '1.5926']
    ]],
    ['13.0', -7, -4, [
        ['HeavyNeutrino_trilepton_M-13.0_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.264e-03',    '0.5148'],
        ['HeavyNeutrino_trilepton_M-13.0_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.264e-03',    '0.5148'],
        ['HeavyNeutrino_trilepton_M-13.0_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.264e-03',    '0.5148'],
        ['HeavyNeutrino_trilepton_M-13.0_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.264e-03',    '0.5148'],
        ['HeavyNeutrino_trilepton_M-13.0_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.264e-03',    '0.5148'],
        ['HeavyNeutrino_trilepton_M-13.0_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.264e-03',    '0.5148']
    ]],
    ['13.2', -7, -4, [
        ['HeavyNeutrino_trilepton_M-13.2_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.256e-03',    '0.4760'],
        ['HeavyNeutrino_trilepton_M-13.2_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.256e-03',    '0.4760'],
        ['HeavyNeutrino_trilepton_M-13.2_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.256e-03',    '0.4760'],
        ['HeavyNeutrino_trilepton_M-13.2_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.256e-03',    '0.4760'],
        ['HeavyNeutrino_trilepton_M-13.2_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.256e-03',    '0.4760']
    ]],
    ['13.4', -7, -4, [
        ['HeavyNeutrino_trilepton_M-13.4_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.255e-03',    '0.4406'],
        ['HeavyNeutrino_trilepton_M-13.4_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.255e-03',    '0.4406'],
        ['HeavyNeutrino_trilepton_M-13.4_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.255e-03',    '0.4406'],
        ['HeavyNeutrino_trilepton_M-13.4_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.255e-03',    '0.4406'],
        ['HeavyNeutrino_trilepton_M-13.4_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.255e-03',    '0.4406'],
        ['HeavyNeutrino_trilepton_M-13.4_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.255e-03',    '0.4406']
    ]],
    ['13.6', -7, -4, [
        ['HeavyNeutrino_trilepton_M-13.6_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.247e-03',    '0.4083'],
        ['HeavyNeutrino_trilepton_M-13.6_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.247e-03',    '0.4083'],
        ['HeavyNeutrino_trilepton_M-13.6_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.247e-03',    '0.4083'],
        ['HeavyNeutrino_trilepton_M-13.6_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.247e-03',    '0.4083'],
        ['HeavyNeutrino_trilepton_M-13.6_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.247e-03',    '0.4083'],
        ['HeavyNeutrino_trilepton_M-13.6_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.247e-03',    '0.4083']
    ]],
    ['13.8', -7, -4, [
        ['HeavyNeutrino_trilepton_M-13.8_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.231e-03',    '0.3789'],
        ['HeavyNeutrino_trilepton_M-13.8_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.231e-03',    '0.3789'],
        ['HeavyNeutrino_trilepton_M-13.8_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.231e-03',    '0.3789'],
        ['HeavyNeutrino_trilepton_M-13.8_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.231e-03',    '0.3789'],
        ['HeavyNeutrino_trilepton_M-13.8_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.231e-03',    '0.3789'],
        ['HeavyNeutrino_trilepton_M-13.8_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.231e-03',    '0.3789']
    ]],
    ['14.0', -7, -4, [
        ['HeavyNeutrino_trilepton_M-14.0_V-0.00141421356237_e_massiveAndCKM_LO_2016_MERGED.root'  , '8.219e-03',    '0.3519'],
        ['HeavyNeutrino_trilepton_M-14.0_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.219e-03',    '0.3519'],
        ['HeavyNeutrino_trilepton_M-14.0_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.219e-03',    '0.3519'],
        ['HeavyNeutrino_trilepton_M-14.0_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.219e-03',    '0.3519'],
        ['HeavyNeutrino_trilepton_M-14.0_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.219e-03',    '0.3519'],
        ['HeavyNeutrino_trilepton_M-14.0_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.219e-03',    '0.3519']
    ]],
    ['14.2', -7, -4, [
        ['HeavyNeutrino_trilepton_M-14.2_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.192e-03',    '0.3272'],
        ['HeavyNeutrino_trilepton_M-14.2_V-0.00141421356237_e_massiveAndCKM_LO_2017_MERGED.root'  , '8.192e-03',    '0.3272'],
        ['HeavyNeutrino_trilepton_M-14.2_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.192e-03',    '0.3272'],
        ['HeavyNeutrino_trilepton_M-14.2_V-0.00141421356237_e_massiveAndCKM_LO_2018_MERGED.root'  , '8.192e-03',    '0.3272'],
        ['HeavyNeutrino_trilepton_M-14.2_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.192e-03',    '0.3272']
    ]],
    ['14.4', -7, -4, [
        ['HeavyNeutrino_trilepton_M-14.4_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.194e-03',    '0.3045'],
        ['HeavyNeutrino_trilepton_M-14.4_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.194e-03',    '0.3045'],
        ['HeavyNeutrino_trilepton_M-14.4_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.194e-03',    '0.3045']
    ]],
    ['14.6', -7, -4, [
        ['HeavyNeutrino_trilepton_M-14.6_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.196e-03',    '0.2837'],
        ['HeavyNeutrino_trilepton_M-14.6_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.196e-03',    '0.2837'],
        ['HeavyNeutrino_trilepton_M-14.6_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.196e-03',    '0.2837']
    ]],
    ['14.8', -7, -4, [
        ['HeavyNeutrino_trilepton_M-14.8_V-0.00141421356237_mu_massiveAndCKM_LO_2016_MERGED.root' , '8.171e-03',    '0.2646'],
        ['HeavyNeutrino_trilepton_M-14.8_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.171e-03',    '0.2646'],
        ['HeavyNeutrino_trilepton_M-14.8_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.171e-03',    '0.2646']
    ]],
    ['15.0', -7, -4, [
        ['HeavyNeutrino_trilepton_M-15.0_V-2.13775583264e-05_e_massiveAndCKM_LO_2016_MERGED.root' , '1.852e-06', '1080.7379'],
        ['HeavyNeutrino_trilepton_M-15.0_V-2.13775583264e-05_mu_massiveAndCKM_LO_2016_MERGED.root', '1.852e-06', '1080.7379'],
        ['HeavyNeutrino_trilepton_M-15.0_V-0.00141421356237_mu_massiveAndCKM_LO_2017_MERGED.root' , '8.154e-03',    '0.2470'],
        ['HeavyNeutrino_trilepton_M-15.0_V-0.00141421356237_mu_massiveAndCKM_LO_2018_MERGED.root' , '8.154e-03',    '0.2470']
    ]]
]

line_m = '%-28s %-90s %-12s\n'
yrs = ['2016', '2017', '2018']

bkgs = {'2016' :
        [
            line_m % ('Obs'   , 'data_2016_2.root'                                                         , ''       ),
            line_m % ('Xgamma', 'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root', '18610'  ),
            line_m % ('Xgamma', 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root'    , '6077.22'),
            line_m % ('Xgamma', 'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root'                    , '0.2086' ),
            line_m % ('Xgamma', 'WZTo3LNu_mllmin01_13TeV-powheg-pythia8_ext1_Summer16.root'                , '4.42965'),
            line_m % ('Xgamma', 'WZG_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root'                    , '0.04123'),
            line_m % ('Xgamma', 'WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root'                    , '0.2147' ),
            line_m % ('Xgamma', 'WWToLNuQQ_13TeV-powheg_Summer16.root'                                     , '49.997' ),
            line_m % ('Xgamma', 'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root'                    , '0.05565'),
            line_m % ('Xgamma', 'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root'                    , '0.01398'),
            line_m % ('Xgamma', 'ZZTo4L_13TeV-amcatnloFXFX-pythia8_Summer16.root'                          , '1.256'  ),
            line_m % ('Xgamma', 'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root'            , '489'    ),
            line_m % ('Xgamma', 'ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root'     , '50.2'   ),
            line_m % ('Xgamma', 'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root'                    , '0.00'   )
        ],
        '2017' :
        [
            line_m % ('Obs'   , 'data_2017.root'                                                                 , ''       ),
            line_m % ('Xgamma', 'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v10_Fall17.root', '18610'  ),
            line_m % ('Xgamma', 'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_realistic_v14_Fall17.root'   , '6077.22'),
            line_m % ('Xgamma', 'WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_realistic_v11_Fall17.root'                , '0.1651' ),
            line_m % ('Xgamma', 'WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_realistic_v11_Fall17.root'                , '0.2086' ),
            line_m % ('Xgamma', 'WZG_TuneCP5_13TeV-amcatnlo-pythia8_realistic_v10_Fall17.root'                   , '0.04123'),
            line_m % ('Xgamma', 'WGGJets_TuneCP5_13TeV_madgraphMLM_pythia8_realistic_v11_Fall17.root'            , '0.2147' ),
            line_m % ('Xgamma', 'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_realistic_v10_Fall17.root'       , '12.178' ),
            line_m % ('Xgamma', 'WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_realistic_v10_Fall17.root'       , '49.997' ),
            line_m % ('Xgamma', 'WZZ_TuneCP5_13TeV-amcatnlo-pythia8_realistic_v11_Fall17.root'                   , '0.05565'),
            line_m % ('Xgamma', 'ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_realistic_v11_Fall17.root'                   , '0.01398'),
            line_m % ('Xgamma', 'ZZTo4L_13TeV_powheg_pythia8_realistic_v14_Fall17.root'                          , '1.256'  ),
            line_m % ('Xgamma', 'WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v14_Fall17.root'           , '489'    ),
            line_m % ('Xgamma', 'ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_realistic_v14_Fall17.root'    , '50.2'   ),
            line_m % ('Xgamma', 'ZZTo4L_13TeV_powheg_pythia8_realistic_v14_Fall17.root'                          , '0.00'   )
        ],
        '2018' :
        [
            line_m % ('Obs'   ,  'data2018.root'                                                        , ''       ),
            line_m % ('Xgamma',  'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_Autumn18.root'  , '18610'  ),
            line_m % ('Xgamma',  'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_Autumn18.root'      , '6077.22'),
            line_m % ('Xgamma',  'WWZ_TuneCP5_13TeV-amcatnlo-pythia8_Autumn18.root'                     , '0.1651' ),
            line_m % ('Xgamma',  'WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_Autumn18.root'                  , '0.2086' ),
            line_m % ('Xgamma',  'WZTo3LNu_TuneCP5_13TeV-powheg-pythia8_Autumn18.root'                  , '4.42965'),
            line_m % ('Xgamma',  'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Autumn18.root'            , '5.595'  ),
            line_m % ('Xgamma',  'WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_Autumn18.root'         , '49.997' ),
            line_m % ('Xgamma',  'WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Autumn18.root'         , '49.997' ),
            line_m % ('Xgamma',  'WZZ_TuneCP5_13TeV-amcatnlo-pythia8_Autumn18.root'                     , '0.05565'),
            line_m % ('Xgamma',  'ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_Autumn18.root'                     , '0.01398'),
            line_m % ('Xgamma',  'ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_Autumn18.root'   , '4.04'   ),
            line_m % ('Xgamma',  'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Autumn18.root'            , '3.22'   ),
            line_m % ('Xgamma',  'ZZTo4L_TuneCP5_13TeV_powheg_pythia8_Autumn18.root'                    , '1.256'  ),
            line_m % ('Xgamma',  'TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_Autumn18.root'     , '3.697'  ),
            line_m % ('Xgamma',  'TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_Autumn18.root', '0.2043' ),
            line_m % ('Xgamma',  'TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_Autumn18.root'        , '0.2728' ),
            line_m % ('Xgamma',  'WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_Autumn18.root'     , '489'    ),
            line_m % ('Xgamma',  'ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_Autumn18.root'      , '50.2'   ),
            line_m % ('Xgamma',  'ZZTo4L_TuneCP5_13TeV_powheg_pythia8_Autumn18.root'                    , '0.00'   )
        ]
}

nstr = int(line_m.split('-')[1].split('s')[0])
addspc = 6

line_m = line_m.split('\n')[0]+' %-12s %-10s\n'
line_d = line_m.replace(str(nstr), str(nstr+addspc), 1)

pts = ['1.00', '1.26', '1.58', '2.00', '2.51', '3.16', '3.98', '5.01', '6.31', '7.94']

for mass, emin, emax, smps in file_xs_ctau:
    for yr in yrs:
        #outf_m = open(yr+'_signal_merged_'      +mass+'.txt' , 'w')
        outf_d = open(yr+'_signal_Dirac_merged_'+mass+'.txt' , 'w')
        datasmpl = bkgs[yr][0].split()[0]
        #outf_m.write(bkgs[yr][0])
        outf_d.write(bkgs[yr][0].replace(datasmpl, datasmpl+(' '*addspc)))
        srcs = []
        tmplist = [smp for smp in smps if yr+'_MERGED' in smp[0] and '_e_' in smp[0]]
        if len(tmplist)>0:
            srcs.append(['e' , tmplist[0][0], tmplist[0][1], tmplist[0][2]])
        tmplist = [smp for smp in smps if yr+'_MERGED' in smp[0] and '_mu_' in smp[0]]
        if len(tmplist)>0:
            srcs.append(['mu', tmplist[0][0], tmplist[0][1], tmplist[0][2]])
        for src in srcs:
            for e in range(emin, emax):
                for pt in pts:
                    v2 = pt+'e'+str(e)
                    pntname = 'M-%s_V-%12.10f_%s'       % (mass, math.sqrt(float(v2)), src[0])
                    #outf_m.write(line_m % (pntname         , src[1], src[2], src[3], v2))
                    outf_d.write(line_d % (pntname+'_Dirac', src[1], src[2], src[3], v2))
        for i in range(1, len(bkgs[yr])):
            bkg = bkgs[yr][i].split()[0]
            #outf_m.write(bkgs[yr][i])
            outf_d.write(bkgs[yr][i].replace(bkg, bkg+(' '*addspc)))
        #outf_m.close()
        outf_d.close()


