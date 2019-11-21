import ROOT as r
import sys, os
from array import array

#indir    = '/user/trocino/Analysis/HNL/Displaced/2019-10-21_ETH_Synchronization/HNL_analysis/dataCards_shapeRoot/'
indir    = '/user/mvit/CMSSW_9_4_4/src/HNL_analysis/dataCards_shapeRoot/'
outdir   = 'merged_datacards_rootfiles_2016/'
inbkgdir = '/user/mvit/CMSSW_9_4_4/src/HNL_analysis/dataCards_shapeRoot/'

if not os.path.exists(outdir): os.mkdir(outdir)

inbkgdata_e = open(inbkgdir+'M-4_V-0.00290516780927_e_ele_datacard.txt' , 'r')
inbkgdata_m = open(inbkgdir+'M-4_V-0.00290516780927_mu_muo_datacard.txt', 'r')

obs_e = ''
bkg_e = ''
obs_m = ''
bkg_m = ''

for line in inbkgdata_e.readlines():
    if 'observation' in line:
        obs_e = line.split()[1]
    elif 'rate' in line:
        div = line.split()[2]
        pos = line.find(div)
        bkg_e = line[pos:]
inbkgdata_e.close()

for line in inbkgdata_m.readlines():
    if 'observation' in line:
        obs_m = line.split()[1]
    elif 'rate' in line:
        div = line.split()[2]
        pos = line.find(div)
        bkg_m = line[pos:]
inbkgdata_m.close()

indatacards = [f for f in os.listdir(indir) if '_datacard'   in f]

for f in indatacards:
    inf  = open( indir+f, 'r' )
    outf = open(outdir+f, 'w+')
    obs = ''
    bkg = ''
    if '_e_ele' in f:
        obs = obs_e
        bkg = bkg_e
    else:
        obs = obs_m
        bkg = bkg_m
    for line in inf.readlines():
        if 'observation' in line:
            outf.write('observation '+obs+'\n')
        elif 'rate' in line:
            div = line.split()[2]
            pos = line.find(div)
            outf.write(line[:pos]+bkg)
        else:
            outf.write(line)
    outf.close()
    inf.close()


inrootfiles = [f for f in os.listdir(indir) if 'shape_file_' in f]

for f in inrootfiles:
    inf  = r.TFile( indir+f, 'read'    )
    outf = r.TFile(outdir+f, 'recreate')
    for key in inf.GetListOfKeys():
        if 'signal' in key.GetName():
            hist = inf.Get(key.GetName())
            outf.cd()
            hist.Write(key.GetName())
    inf.Close()
    inbkgroot = r.TFile()
    if '_e_ele' in f:
        inbkgroot = r.TFile(inbkgdir+'shape_file_M-4_V-0.00290516780927_e_ele.root' , 'read')
    else:
        inbkgroot = r.TFile(inbkgdir+'shape_file_M-4_V-0.00290516780927_mu_muo.root', 'read')
    for key in inbkgroot.GetListOfKeys():
        if 'signal' in key.GetName(): continue
        hist = inbkgroot.Get(key.GetName())
        outf.cd()
        hist.Write(key.GetName())
    outf.Close()
    inbkgroot.Close()
