import ROOT as r
import sys, os

r.gROOT.LoadMacro('/Applications/root/test/danystyle.C') 
r.setTDRStyle()

# Output file
#fout = r.TFile('plots_3l_4l.root', 'recreate')

# Plot name
#plots = ['e_ele', 'mu_muo']
plots = ['mu_muo']

nplot = len(plots)

# Labels
labels = [
    'M_{HNL} = 2 GeV, |V_{N#mu}|^{2} = 6\times10^{#minus4}'
]

mcnames = {}
mcnames['M_{HNL} = 2 GeV, |V_{N#mu}|^{2} = 6\times10^{#minus4}'] = 'M-2_V-0.0248394846967'

systs = [
    'pu',
    'qcd',
    'pdf',
    'pEle',
    'pMuo',
    'npEle',
    'npMuo',
    'jec',
    'jer',
    'btag'
]


# Get plots
iplot=0
leg = r.TLegend(0.75, 0.65, 0.95, 0.90)
leg.SetBorderSize(0)
leg.SetLineWidth(0)
leg.SetLineColor(0)
leg.SetFillStyle(0)
leg.SetFillColor(0)
for plot in plots:
    for label in labels:
        name = '/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/NOTES/AN-18-014/tmp/PLOT_AN/shape_signal/shape_file_'+mcnames[label]+'_'+plot+'.root'
        f = r.TFile(name)
        h = f.Get('signal')
        valct = h.Integral(0, -1)
        h.SetLineWidth(2)
        c = r.TCanvas('cc', 'cc', 600, 700)
        p1 = r.TPad('p1', 'p1', 0., 0., 1., 0.25)
        p2 = r.TPad('p2', 'p2', 0., 0.25, 1., 1.)
        p1.SetBottomMargin(0.34)
        p1.SetTopMargin(0.02)
        p2.SetBottomMargin(0.01)
        c.Draw()
        p1.Draw()
        p2.Draw()
        p2.cd()
        for syst in systs:
            h.SetTitle(label+' - Syst: '+syst)
            hdn = f.Get('signal_'+syst+'Down')
            hup = f.Get('signal_'+syst+'Up')
            hdn.SetName('signal_'+mcnames[label]+'_'+syst+'Down')
            hup.SetName('signal_'+mcnames[label]+'_'+syst+'Up')
            valdn = hdn.Integral(0, -1)
            valup = hup.Integral(0, -1)
            errstr = str((valdn+valup)/(2*valct)-1)
            #print ' >> ' + syst + ': ' + str(valdn/valct-1) + ' ' + str(valup/valct-1)
            hdn.SetLineWidth(2)
            hup.SetLineWidth(2)
            hdn.SetLineColor(r.kBlue)
            hup.SetLineColor(r.kRed)
            hdn.SetLineStyle(2)
            hup.SetLineStyle(7)
            if leg.GetNRows()==0:
                leg.AddEntry(h  , 'central', 'l')
                leg.AddEntry(hdn, 'central #minus #sigma', 'l')
                leg.AddEntry(hup, 'central #plus #sigma', 'l')
            p2.cd()
            h.Draw('hist')
            hdn.Draw('histsame')
            hup.Draw('histsame')
            leg.Draw()
            hratiodn = hdn.Clone('ratio_'+mcnames[label]+'_'+syst+'Down')
            hratioup = hup.Clone('ratio_'+mcnames[label]+'_'+syst+'Up')
            hratiodn.Divide(h)
            hratioup.Divide(h)
            hratiodn.SetLineWidth(2)
            hratioup.SetLineWidth(2)
            hratiodn.SetLineColor(r.kBlue)
            hratioup.SetLineColor(r.kRed)
            hratiodn.SetLineStyle(2)
            hratioup.SetLineStyle(7)
            for x in range(1, h.GetNbinsX()+1):
                if h.GetBinContent(x)==0:
                    if hdn.GetBinContent(x)==0:
                        hratiodn.SetBinContent(x, 1.0)
                    if hup.GetBinContent(x)==0:
                        hratioup.SetBinContent(x, 1.0)
            p1.cd()
            hratiodn.Draw('hist')
            hratioup.Draw('histsame')
            hratiodn.GetYaxis().SetRangeUser(0.46, 1.54)
            p1.SetGridx()
            p1.SetGridy()
            hratiodn.GetXaxis().SetTitleSize(0.15)
            hratiodn.GetXaxis().SetLabelSize(0.10)
            hratiodn.GetYaxis().SetTitle('rel. unc.')
            hratiodn.GetYaxis().SetTitleSize(0.14)
            hratiodn.GetYaxis().SetTitleOffset(0.50)
            hratiodn.GetYaxis().SetLabelSize(0.10)
            c.SaveAs('systVar_'+mcnames[label]+'_'+syst+'.png')
