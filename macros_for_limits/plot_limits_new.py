import ROOT as r
import sys, os
from array import array
import math

r.gROOT.LoadMacro('/afs/cern.ch/user/t/trocino/Utilities/danystyle.C') 
r.setTDRStyle()

def interpolateLimit(hist):
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    masses = []
    v2s = []
    for x in range(1, hist.GetNbinsX()+1):
        mass = round(hist.GetXaxis().GetBinCenter(x))
        masses.append(mass)
        prevv2 = -1.
        prevmu = -1.
        foundv2 = False
        for y in range(1, hist.GetNbinsY()+1):
            v2 = hist.GetYaxis().GetBinCenter(y)
            mu = hist.GetBinContent(x, y)
            if mu<1e-9: continue
            if mu<1. and prevv2<0.:
                print ' >>> WARNING: no lower limit for histogram '+hist.GetName()+' and mass = '+str(mass)+' <<< '
                v2s.append(-1.)
                break
            if prevmu>1. and mu<=1:
                limv2 = prevv2 + ((1-prevmu)/(mu-prevmu))*(v2-prevv2)
                #print '   -- Found crossing point: mass = '+str(mass)+', v2 = '+str(limv2)+' ('+str(prevv2)+' -- '+str(v2)+'), ('+str(prevmu)+' -- '+str(mu)+')'
                v2s.append(limv2)
                foundv2 = True
                break
            prevv2 = v2
            prevmu = mu
        if not foundv2: v2s.append(-1.)
    return masses,v2s

# List of original param-space points
orig_mass = array('f', [1.     , 1.    , 2.      , 2.      , 3.      , 4.      , 5.      , 6.     , 8.      , 10.     ])
orig_vsqr = array('f', [0.00902, 0.0451, 0.000123, 0.000617, 5.01e-05, 8.44e-06, 2.12e-06, 4.1e-06, 2.29e-06, 5.73e-07])
#orig_v = [0.0949736805647, 0.212367605816 , 0.0110905365064, 0.0248394846967, 0.00707813534767, 0.00290516780927, 0.00145602197786, 0.00202484567313, 0.00151327459504, 0.000756967634711]

new_mass = array('f',
                 [1.   , 1.   , 2.   , 2.   , 3.   , 3.   , 4.   , 4.   , 5.   , 5.   , 6.   , 6.   , 7.   , 7.   , 7.   , 8.   , 8.   , 9.   , 9.   , 10.   , 10.   , 11.   , 11.   , 12.   , 12.   ]
)
new_vsqr = array('f',
                 [5.e-1, 5.e-4, 5.e-1, 5.e-2, 5.e-1, 5.e-2, 5.e-2, 5.e-3, 5.e-2, 5.e-3, 1.e-2, 1.e-3, 1.e-2, 1.e-3, 5.e-6, 1.e-2, 1.e-3, 1.e-3, 1.e-5,  1.e-4,  1.e-5,  1.e-4,  1.e-5,  1.e-4,  1.e-5]
)

# Print out new samples
print '--------------------------------------------------------------'
print ' Old samples'
print ('%8s %12s %12s' % ('Mass', 'V^2', 'V'))
for i in range(0, len(orig_mass)):
    print ('%8d %12.2e %12.4e' % (orig_mass[i], orig_vsqr[i], math.sqrt(orig_vsqr[i])))
print '--------------------------------------------------------------'
print ' New samples'
print ('%8s %12s %12s' % ('Mass', 'V^2', 'V'))
for i in range(0, len(new_mass)):
    print ('%8d %12.2e %12.4e' % (new_mass[i], new_vsqr[i], math.sqrt(new_vsqr[i])))
print '--------------------------------------------------------------'

# Input files
indir  = '/afs/cern.ch/work/m/mvit/public/limiti_per_daniele/'
#infile = open(indir+'list_datacards.txt', 'r')
infile = open('list_datacards.txt', 'r')
inlist = infile.readlines()

# Plot name
couplings = ['e_ele', 'mu_muo']
#couplings = ['mu_muo']

ncoupl = len(couplings)

xbins = array('f', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.0, 11.0])
nxbins = len(xbins) - 1
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
nybins = len(ybins) - 1


# Extended axes
ext_xbins = xbins[:]
ext_ybins = ybins[:]
ext_xbins.insert(0, 0.)
ext_xbins.append(12.5)
ext_ybins.append(1.)

## Plotting
hs_m2s = {}
hs_m1s = {}
hs_ctr = {}
hs_p1s = {}
hs_p2s = {}
hs_obs = {}
fakeh2 = {}
for ic in couplings:
    cstr = 'e'
    if 'mu_muo' in ic: cstr = '#mu'
    hs_m2s[ic] = r.TH2F('hlimits_m2s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #minus 2s.d.)/#sigma_{SM}', nxbins, xbins, nybins, ybins)
    hs_m1s[ic] = r.TH2F('hlimits_m1s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #minus s.d.)/#sigma_{SM}' , nxbins, xbins, nybins, ybins)
    hs_ctr[ic] = r.TH2F('hlimits_ctr_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected #sigma_{95% CL}/#sigma_{SM}'               , nxbins, xbins, nybins, ybins)
    hs_p1s[ic] = r.TH2F('hlimits_p1s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #plus s.d.)/#sigma_{SM}'  , nxbins, xbins, nybins, ybins)
    hs_p2s[ic] = r.TH2F('hlimits_p2s_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected (#sigma_{95% CL} #plus 2 s.d.)/#sigma_{SM}', nxbins, xbins, nybins, ybins)
    hs_obs[ic] = r.TH2F('hlimits_obs_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Observed #sigma_{95% CL}/#sigma_{SM}'               , nxbins, xbins, nybins, ybins)
    #
    fakeh2[ic] = r.TH2F('fakeh2_'+ic, ';M_{N} [GeV];|V_{'+cstr+'}|^{2};Expected #sigma_{95% CL}/#sigma_{SM}', len(ext_xbins)-1, ext_xbins, len(ext_ybins)-1, ext_ybins)


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
        if not os.path.isfile(name):
            print ' >>> WARNING: there is no file '+name
            continue
        mass = int(inf.split('_')[0][2:])
        vsqr = float(inf.split('_')[1][2:])
        vsqr = vsqr*vsqr
        f = r.TFile(name)
        # t = f.Get('limit')
        # lim = array('f', [0.0])
        # t.SetBranchAddress('limit', lim)
        if not f.GetListOfKeys().Contains('limit'):
            print ' >>> WARNING: file '+name+' has no TTree named "limit"'
            continue
        t = f.limit
        #
        # expected -2sd
        t.GetEntry(0)
        hs_m2s[cpl].Fill(mass, vsqr, t.limit)
        #print str(t.limit)
        #
        # expected -1sd
        t.GetEntry(1)
        hs_m1s[cpl].Fill(mass, vsqr, t.limit)
        #print str(t.limit)
        #
        # expected central
        t.GetEntry(2)
        hs_ctr[cpl].Fill(mass, vsqr, t.limit)
        #print str(t.limit)
        #
        # expected +1sd
        t.GetEntry(3)
        hs_p1s[cpl].Fill(mass, vsqr, t.limit)
        #print str(t.limit)
        #
        # expected +2sd
        t.GetEntry(4)
        hs_p2s[cpl].Fill(mass, vsqr, t.limit)
        #print str(t.limit)
        #
        # observed central
        t.GetEntry(5)
        hs_obs[cpl].Fill(mass, vsqr, t.limit)
        #print str(t.limit)
    c.cd()
    c.SetLogy()
    hs_m2s[cpl].Draw('colz')
    hs_m2s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_m2s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_m2s[cpl].GetZaxis().SetLabelSize(0.03)
    hs_m2s[cpl].GetZaxis().SetRangeUser(0., 5.)
    c.SaveAs('signalStrength_'+cpl+'_m2s.png')
    hs_m1s[cpl].Draw('colz')
    hs_m1s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_m1s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_m1s[cpl].GetZaxis().SetLabelSize(0.03)
    hs_m1s[cpl].GetZaxis().SetRangeUser(0., 5.)
    c.SaveAs('signalStrength_'+cpl+'_m1s.png')
    hs_ctr[cpl].Draw('colz')
    hs_ctr[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_ctr[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_ctr[cpl].GetZaxis().SetLabelSize(0.03)
    hs_ctr[cpl].GetZaxis().SetRangeUser(0., 5.)
    c.SaveAs('signalStrength_'+cpl+'_ctr.png')
    hs_p1s[cpl].Draw('colz')
    hs_p1s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_p1s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_p1s[cpl].GetZaxis().SetLabelSize(0.03)
    hs_p1s[cpl].GetZaxis().SetRangeUser(0., 5.)
    c.SaveAs('signalStrength_'+cpl+'_p1s.png')
    hs_p2s[cpl].Draw('colz')
    hs_p2s[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_p2s[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_p2s[cpl].GetZaxis().SetLabelSize(0.03)
    hs_p2s[cpl].GetZaxis().SetRangeUser(0., 5.)
    c.SaveAs('signalStrength_'+cpl+'_p2s.png')
    hs_obs[cpl].Draw('colz')
    hs_obs[cpl].GetYaxis().SetTitleOffset(1.10)
    hs_obs[cpl].GetZaxis().SetTitleOffset(1.15)
    hs_obs[cpl].GetZaxis().SetLabelSize(0.03)
    hs_obs[cpl].GetZaxis().SetRangeUser(0., 5.)
    c.SaveAs('signalStrength_'+cpl+'_obs.png')

##
## Interpolate limits
##
for cpl in couplings:
    #if 'mu' in cpl: continue
    ms, lims = interpolateLimit(hs_ctr[cpl])
    #print '   ==== Limits for '+cpl+': '
    #for ii in range(0, len(ms)):
    #    print '('+str(ms[ii])+', '+str(lims[ii])+')'
    masspoints = array('f', ms)
    vsqrpoints = array('f', lims)
    gr = r.TGraph(len(ms), masspoints, vsqrpoints);
    gr.SetMarkerSize(1.0)
    gr.SetMarkerColor(r.kRed)
    gr.SetLineColor(r.kRed)
    gr.SetLineWidth(3)
    c.cd()
    hs_ctr[cpl].Draw('colz')
    gr.Draw('L')
    c.SaveAs('signalStrength_'+cpl+'_ctr_limits.png')
    ## Add original grid points and new ones
    orig_gr = r.TGraph(len(orig_mass), orig_mass, orig_vsqr);
    orig_gr.SetMarkerStyle(20)
    orig_gr.SetMarkerSize(1.5)
    orig_gr.SetMarkerColor(r.kOrange)
    new_gr = r.TGraph(len(new_mass), new_mass, new_vsqr);
    new_gr.SetMarkerStyle(21)
    new_gr.SetMarkerSize(1.5)
    new_gr.SetMarkerColor(r.kGreen)
    c.cd()
    fakeh2[cpl].Draw()
    fakeh2[cpl].GetYaxis().SetTitleOffset(1.10)
    fakeh2[cpl].GetZaxis().SetTitleOffset(1.15)
    fakeh2[cpl].GetZaxis().SetLabelSize(0.03)
    hs_ctr[cpl].Draw('colzsame')
    gr.Draw('L')
    orig_gr.Draw('P')
    new_gr.Draw('P')
    c.SaveAs('signalStrength_'+cpl+'_ctr_limits_orig_new.png')
    ##
    ## 1-sigma band
    ms, lims1 = interpolateLimit(hs_m1s[cpl])
    ms, lims2 = interpolateLimit(hs_p1s[cpl])
    ms_ctr = []
    ms_err = []
    lims_ctr = []
    lims_err = []
    #print ms
    #print lims1
    #print lims2
    for ii in range(0, len(ms)):
        if lims1[ii]<0. or lims2[ii]<0.: continue
        ms_ctr.append(ms[ii])
        ms_err.append(0.)
        lims_ctr.append((lims1[ii]+lims2[ii])/2.)
        lims_err.append(abs(lims2[ii]-lims1[ii])/2.)
    masspoints = array('f', ms_ctr)
    masserrors = array('f', ms_err)
    vsqrpoints = array('f', lims_ctr)
    vsqrerrors = array('f', lims_err)
    gr1s = r.TGraphErrors(len(masspoints), masspoints, vsqrpoints, masserrors, vsqrerrors)
    gr1s.SetFillColor(r.kGreen)
    ##
    ## 2-sigma band
    del ms_ctr[:]
    del ms_err[:]
    del lims_ctr[:]
    del lims_err[:]
    ms, lims1 = interpolateLimit(hs_m2s[cpl])
    ms, lims2 = interpolateLimit(hs_p2s[cpl])
    #print ms
    #print lims1
    #print lims2
    for ii in range(0, len(ms)):
        if lims1[ii]<0. or lims2[ii]<0.: continue
        ms_ctr.append(ms[ii])
        ms_err.append(0.)
        lims_ctr.append((lims1[ii]+lims2[ii])/2.)
        lims_err.append(abs(lims2[ii]-lims1[ii])/2.)
    masspoints = array('f', ms_ctr)
    masserrors = array('f', ms_err)
    vsqrpoints = array('f', lims_ctr)
    vsqrerrors = array('f', lims_err)
    gr2s = r.TGraphErrors(len(masspoints), masspoints, vsqrpoints, masserrors, vsqrerrors)
    gr2s.SetFillColor(r.kYellow)
    c.cd()
    gr2s.Draw('a3')
    gr1s.Draw('3same')
    gr.Draw('L')
    gr2s.GetXaxis().SetLimits(xbins[0], xbins[nxbins])
    gr2s.GetYaxis().SetRangeUser(ybins[0], ybins[nybins])
    c.SaveAs('signalStrength_'+cpl+'_limits.png')
