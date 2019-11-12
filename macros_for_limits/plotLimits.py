import ROOT as r
import sys, os
from array import array

r.gROOT.LoadMacro('/afs/cern.ch/user/t/trocino/Utilities/danystyle.C') 
r.setTDRStyle()

# Input files
indir  = '2019-09-27_prova/'
infile = open(indir+'list_datacards.txt', 'r')
inlist = infile.readlines()

# Output files
#fout = r.TFile('plots_3l_4l.root', 'recreate')

# Plot name
couplings = ['e_ele', 'mu_muo']
#couplings = ['mu_muo']

ncoupl = len(couplings)

xbins = array('f', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.0, 9.0, 11.0])
##
## Y axis with 60 bins
# ybins = array('f', [ 1.00e-07 , 1.26e-07 , 1.58e-07 , 2.00e-07 , 2.51e-07 , 3.16e-07 , 3.98e-07 , 5.01e-07 , 6.31e-07 , 7.94e-07 ,\
#                      1.00e-06 , 1.26e-06 , 1.58e-06 , 2.00e-06 , 2.51e-06 , 3.16e-06 , 3.98e-06 , 5.01e-06 , 6.31e-06 , 7.94e-06 ,\
#                      1.00e-05 , 1.26e-05 , 1.58e-05 , 2.00e-05 , 2.51e-05 , 3.16e-05 , 3.98e-05 , 5.01e-05 , 6.31e-05 , 7.94e-05 ,\
#                      1.00e-04 , 1.26e-04 , 1.58e-04 , 2.00e-04 , 2.51e-04 , 3.16e-04 , 3.98e-04 , 5.01e-04 , 6.31e-04 , 7.94e-04 ,\
#                      1.00e-03 , 1.26e-03 , 1.58e-03 , 2.00e-03 , 2.51e-03 , 3.16e-03 , 3.98e-03 , 5.01e-03 , 6.31e-03 , 7.94e-03 ,\
#                      1.00e-02 , 1.26e-02 , 1.58e-02 , 2.00e-02 , 2.51e-02 , 3.16e-02 , 3.98e-02 , 5.01e-02 , 6.31e-02 , 7.94e-02 ,\
#                      1.00e-01 ])
##
## Y axis with 30 bins
##   --> mid-points:      1.26e-07     2e-07    3.16e-07   5.01e-07   7.94e-07
##                          /   \      /   \      /   \      /   \      /   
ybins = array('f', [ 1.00e-07 , 1.58e-07 , 2.51e-07 , 3.98e-07 , 6.31e-07 ,\
                     1.00e-06 , 1.58e-06 , 2.51e-06 , 3.98e-06 , 6.31e-06 ,\
                     1.00e-05 , 1.58e-05 , 2.51e-05 , 3.98e-05 , 6.31e-05 ,\
                     1.00e-04 , 1.58e-04 , 2.51e-04 , 3.98e-04 , 6.31e-04 ,\
                     1.00e-03 , 1.58e-03 , 2.51e-03 , 3.98e-03 , 6.31e-03 ,\
                     1.00e-02 , 1.58e-02 , 2.51e-02 , 3.98e-02 , 6.31e-02 ,\
                     1.00e-01 ])
hs_m2s = {}
hs_m1s = {}
hs_ctr = {}
hs_p1s = {}
hs_p2s = {}
hs_obs = {}
for ic in couplings:
    cstr = 'e'
    if 'mu_muo' in ic: cstr = '#mu'
    hs_m2s[ic] = r.TH2F('hlimits_m2s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #minus 2s.d.)/#sigma_{SM}', 9, xbins, 30, ybins)
    hs_m1s[ic] = r.TH2F('hlimits_m1s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #minus s.d.)/#sigma_{SM}' , 9, xbins, 30, ybins)
    hs_ctr[ic] = r.TH2F('hlimits_ctr_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected #sigma_{95% CL}/#sigma_{SM}'               , 9, xbins, 30, ybins)
    hs_p1s[ic] = r.TH2F('hlimits_p1s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #plus s.d.)/#sigma_{SM}'  , 9, xbins, 30, ybins)
    hs_p2s[ic] = r.TH2F('hlimits_p2s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #plus 2 s.d.)/#sigma_{SM}', 9, xbins, 30, ybins)
    hs_obs[ic] = r.TH2F('hlimits_obs_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Observed #sigma_{95% CL}/#sigma_{SM}'               , 9, xbins, 30, ybins)

# Get plots
# leg = r.TLegend(0.75, 0.65, 0.95, 0.90)
# leg.SetBorderSize(0)
# leg.SetLineWidth(0)
# leg.SetLineColor(0)
# leg.SetFillStyle(0)
# leg.SetFillColor(0)

c = r.TCanvas('cc', 'cc', 700, 600)
c.SetLeftMargin(0.14)
c.SetRightMargin(0.20)
c.Draw()
for cpl in couplings:
    for inf in inlist:
        if cpl not in inf: continue
        inn = inf.strip()
        name = indir+'higgsCombine'+inn+'.AsymptoticLimits.mH120.root'
        mass = int(inf.split('_')[0][2:])
        vsqr = float(inf.split('_')[1][2:])
        vsqr = vsqr*vsqr
        f = r.TFile(name)
        # t = f.Get('limit')
        # lim = array('f', [0.0])
        # t.SetBranchAddress('limit', lim)
        t = f.limit
        #
        # expected -2sd
        t.GetEntry(0)
        hs_m2s[cpl].Fill(mass, vsqr, t.limit)
        print str(t.limit)
        #
        # expected -1sd
        t.GetEntry(1)
        hs_m1s[cpl].Fill(mass, vsqr, t.limit)
        print str(t.limit)
        #
        # expected central
        t.GetEntry(2)
        hs_ctr[cpl].Fill(mass, vsqr, t.limit)
        print str(t.limit)
        #
        # expected +1sd
        t.GetEntry(3)
        hs_p1s[cpl].Fill(mass, vsqr, t.limit)
        print str(t.limit)
        #
        # expected +2sd
        t.GetEntry(4)
        hs_p2s[cpl].Fill(mass, vsqr, t.limit)
        print str(t.limit)
        #
        # observed central
        t.GetEntry(5)
        hs_obs[cpl].Fill(mass, vsqr, t.limit)
        print str(t.limit)
    c.cd()
    c.SetLogy()
    hs_m2s[cpl].Draw('colz')
    hs_m2s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_m2s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_m2s[cpl].GetZaxis().SetLabelSize(0.03)
    c.SaveAs('signalStrength_'+cpl+'_m2s.png')
    hs_m1s[cpl].Draw('colz')
    hs_m1s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_m1s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_m1s[cpl].GetZaxis().SetLabelSize(0.03)
    c.SaveAs('signalStrength_'+cpl+'_m1s.png')
    hs_ctr[cpl].Draw('colz')
    hs_ctr[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_ctr[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_ctr[cpl].GetZaxis().SetLabelSize(0.03)
    c.SaveAs('signalStrength_'+cpl+'_ctr.png')
    hs_p1s[cpl].Draw('colz')
    hs_p1s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_p1s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_p1s[cpl].GetZaxis().SetLabelSize(0.03)
    c.SaveAs('signalStrength_'+cpl+'_p1s.png')
    hs_p2s[cpl].Draw('colz')
    hs_p2s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_p2s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_p2s[cpl].GetZaxis().SetLabelSize(0.03)
    c.SaveAs('signalStrength_'+cpl+'_p2s.png')
    hs_obs[cpl].Draw('colz')
    hs_obs[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_obs[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_obs[cpl].GetZaxis().SetLabelSize(0.03)
    c.SaveAs('signalStrength_'+cpl+'_obs.png')
