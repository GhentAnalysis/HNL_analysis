import sys, os
import math as m
import ROOT as r
from array import array
from urllib import urlopen

r.gROOT.LoadMacro('/storage_mnt/storage/user/trocino/Utilities/danystyle.C')
r.setTDRStyle()
r.gStyle.SetPalette(r.kBird)
r.gROOT.SetBatch(True)

def sortByMassV2(e):
  return [float(e[2]), float(e[3])]

## Final state
fs_e = 'e'
fs_m = 'mu'
fs = fs_e

## Year
year_16 = '2016'
year_17 = '2017'
year_18 = '2018'
year = year_16

## Masses
masses_all = ['1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0', '11.0', '12.0', '12.2', '12.4', '12.6', '12.8', '13.0', '13.2', '13.4', '13.6', '13.8', '14.0', '14.2', '14.4', '14.6', '14.8', '15.0']
#
masses_lo = [ '1.0',  '2.0',  '3.0',  '4.0',  '5.0',  '6.0']
masses_md = [ '7.0',  '8.0',  '9.0', '10.0', '11.0', '12.0']
masses_hi = ['12.2', '12.4', '12.6', '12.8', '13.0', '13.2', '13.4', '13.6', '13.8', '14.0', '14.2', '14.4', '14.6', '14.8', '15.0']
#
masses_nonew = [ '1.0', '2.0', '3.0', '4.0', '5.0', '8.0', '10.0']
masses_new = ['6.0', '7.0', '9.0', '11.0', '12.0', '12.2', '12.4', '12.6', '12.8', '13.0', '13.2', '13.4', '13.6', '13.8', '14.0', '14.2', '14.4', '14.6', '14.8', '15.0']
#
masses = masses_nonew

## Output dir
outputdir = '/storage_mnt/storage/user/trocino/Analysis/HNL/Displaced/2020-12-15_LeptonNumberViolation/CMSSW_10_2_22/src/heavyNeutrino/multilep/test/MERGED_SAMPLES/'
outputdir += ('' if outputdir[-1]=='/' else '/')
if not os.path.exists(outputdir): os.mkdir(outputdir)

coupls = ['e', 'mu']

####
## Private samples
####

## Ntuple location
indir_pvt = '/pnfs/iihe/cms/store/user/mvit/samples/FINAL/'+year+'/new_signal/'
indir_pvt += ('' if indir_pvt[-1]=='/' else '/')

campaign = ''
if   year=='2016': campaign = 'Moriond17_aug2018_miniAODv3'
elif year=='2017': campaign = 'Fall17'
elif year=='2018': campaign = 'Autumn18'

sampletempl = 'majorana*'+campaign+'%2Fdisplaced*_trilepton_*_'+fs+'_'
tomspage = urlopen('https://tomc.web.cern.ch/tomc/heavyNeutrino/?match='+sampletempl)
#tomslist = [line.split() for line in tomspage.readlines() if 'pnfs' in line and line.split()[0]=='*' and float(line.split()[2])<19.5 and int(line.split()[6])>50000]
tomslist = [line.split() for line in tomspage.readlines() if len(line.split())==12 and 'pnfs' in line.split()[11] and line.split()[0]=='*' and line.split()[2] in masses and line.split()[4]!='-' and line.split()[5]!='-' and line.split()[6]!='?' and int(line.split()[6])>50000 and line.split()[8]=='+-']
tomslength = len(tomslist)

# print ' ---- Tom\'s list size: ' + str(tomslength)
# print tomslist

#tomslist.sort(key=sortByMassV2)

####
## Centrally produced samples
####

## Ntuple location
indir_ctr = '/storage_mnt/storage/user/trocino/Analysis/HNL/Displaced/2020-12-15_LeptonNumberViolation/CMSSW_10_2_22/src/heavyNeutrino/multilep/test/CENTRAL_SAMPLES/'
indir_ctr += ('' if indir_ctr[-1]=='/' else '/')

if   year=='2016': campaign = 'Summer16'
#elif year=='2017': campaign = 'Fall17'
#elif year=='2018': campaign = 'Autumn18'

basilespage = open('/storage_mnt/storage/user/trocino/Analysis/HNL/Displaced/2020-12-15_LeptonNumberViolation/CMSSW_10_2_22/src/heavyNeutrino/multilep/test/CENTRAL_SAMPLES/DAS_parsing_scripts/availableHeavyNeutrinoSamples_official.txt', 'r')
basileslist = [line.split() for line in basilespage.readlines() if len(line.split())==12 and 'HeavyNeutrino_' in line.split()[11] and '_'+fs+'_' in line.split()[11] and campaign in line.split()[11] and line.split()[2] in masses and line.split()[4]!='-' and line.split()[5]!='-' and line.split()[6]!='?' and int(line.split()[6])>50000 and line.split()[8]=='+-']
basileslength = len(basileslist)
basilespage.close()

# print ' ---- Basile\'s list size: ' + str(basileslength)
# print basileslist

#basileslist.sort(key=sortByMassV2)

####
## New centrally produced samples 
####

## Ntuple location
indir_new = '/storage_mnt/storage/user/trocino/Analysis/HNL/Displaced/2020-12-15_LeptonNumberViolation/CMSSW_10_2_22/src/heavyNeutrino/multilep/test/NEW_SAMPLES/'
indir_new += ('' if indir_new[-1]=='/' else '/')

if   year=='2016': campaign = 'Summer16'
#elif year=='2017': campaign = 'Fall17'
#elif year=='2018': campaign = 'Autumn18'

kirillspage = open('/storage_mnt/storage/user/trocino/Analysis/HNL/Displaced/2020-12-15_LeptonNumberViolation/CMSSW_10_2_22/src/heavyNeutrino/multilep/test/NEW_SAMPLES/DAS_parsing_scripts/availableHeavyNeutrinoSamples_official.txt', 'r')
kirillslist = [line.split() for line in kirillspage.readlines() if len(line.split())==12 and 'HeavyNeutrino_' in line.split()[11] and '_'+fs+'_' in line.split()[11] and campaign in line.split()[11] and line.split()[2] in masses and line.split()[4]!='-' and line.split()[5]!='-' and line.split()[6]!='?' and int(line.split()[6])>50000 and line.split()[8]=='+-']
kirillslength = len(kirillslist)
kirillspage.close()

# print ' ---- Kirill\'s list size: ' + str(kirillslength)
# print kirillslist

#kirillslist.sort(key=sortByMassV2)

####
## Merge lists
####

#if tomslength==0 and basileslength==0: sys.exit(0)
if tomslength==0 and basileslength==0 and kirillslength==0: sys.exit(0)

listlength = tomslength + basileslength + kirillslength
#if listlength==0: sys.exit(0)

if      tomslength>0: linelength = len(   tomslist[0])
elif basileslength>0: linelength = len(basileslist[0])
elif kirillslength>0: linelength = len(kirillslist[0])

## Initialize merged list 
totlist = [ [ '' for J in range(linelength+1) ] for I in range(listlength) ]

for i in range(0, tomslength):
  for j in range(0, linelength-1):
    #print ' (%2d, %2d) %s' % (i, j, tomslist[i][j])
    totlist[i][j] = tomslist[i][j]
  #print ' (%2d, %2d) %s' % (i, linelength-1, tomslist[i][linelength-1].split('/')[-1]+'_'+year+'.root')
  totlist[i][linelength-1] = tomslist[i][linelength-1].split('/')[-1]+'_'+year+'.root'
  # mstr = totlist[i][linelength-1].split('_')[2]
  # if not '.' in mstr:
  #   totlist[i][linelength-1] = totlist[i][linelength-1].replace(mstr, mstr+'.0')
  #print ' (%2d, %2d) %s' % (i, linelength, indir_pvt)
  totlist[i][linelength]   = indir_pvt

for i in range(0, basileslength):
  for j in range(0, linelength-1):
    #print ' (%2d, %2d) %s' % (tomslength+i, j, basileslist[i][j])
    totlist[tomslength+i][j] = basileslist[i][j]
  #print ' (%2d, %2d) %s' % (tomslength+i, linelength-1, basileslist[i][linelength-1].split('/')[1].split('_LO')[0]+'_LO_'+year+'.root')
  totlist[tomslength+i][linelength-1] = basileslist[i][linelength-1].split('/')[1].split('_LO')[0]+'_LO_'+year+'.root'
  totlist[tomslength+i][linelength-1] = totlist[tomslength+i][linelength-1].replace('_V-0_', '_V-0.')
  # mstr = totlist[tomslength+i][linelength-1].split('_')[2]
  # if not '.' in mstr:
  #   totlist[tomslength+i][linelength-1] = totlist[tomslength+i][linelength-1].replace(mstr, mstr+'.0')
  #print ' (%2d, %2d) %s' % (tomslength+i, linelength, indir_ctr)
  totlist[tomslength+i][linelength]   = indir_ctr

for i in range(0, kirillslength):
  for j in range(0, linelength-1):
    #print ' (%2d, %2d) %s' % (tomslength+basileslength+i, j, kirillslist[i][j])
    totlist[tomslength+basileslength+i][j] = kirillslist[i][j]
  #print ' (%2d, %2d) %s' % (tomslength+basileslength+i, linelength-1, kirillslist[i][linelength-1].split('/')[1].split('_'+fs+'_')[0]+'_'+fs+'_massiveAndCKM_LO_'+year+'.root')
  totlist[tomslength+basileslength+i][linelength-1] = kirillslist[i][linelength-1].split('/')[1].split('_'+fs+'_')[0]+'_'+fs+'_massiveAndCKM_LO_'+year+'.root'
  totlist[tomslength+basileslength+i][linelength-1] = totlist[tomslength+basileslength+i][linelength-1].replace('_V-0p', '_V-0.')
  mstr = totlist[tomslength+basileslength+i][linelength-1].split('_')[2]
  totlist[tomslength+basileslength+i][linelength-1] = totlist[tomslength+basileslength+i][linelength-1].replace(mstr, mstr.replace('p', '.'))
  #print ' (%2d, %2d) %s' % (tomslength+basileslength+i, linelength, indir_new)
  totlist[tomslength+basileslength+i][linelength]   = indir_new

## Sort final list by (mass, V2)
totlist.sort(key=sortByMassV2)

# print ' ---- Final list size: ' + str(len(totlist)) + ' (' + str(listlength) + ')'
# for i in range(0, len(totlist)):
#   print totlist[i]

# n = 0
# for smpl in tomslist:
#     n += 1
#     #print '%3d. %-140s %7s' % (n, smpl[11], smpl[6])
#     filename = smpl[11].split('/')[-1]+'_'+year+'.root'
#     if os.path.exists(indir+filename):
#       tmpfile = r.TFile.Open(indir+filename)
#       tmptree = tmpfile.Get('blackJackAndHookers/blackJackAndHookersTree')
#       #print '%3d. %-140s %7s' % (n, indir+filename, tmptree.GetEntries())
#     # h2.Fill(float(smpl[2]), float(smpl[3]), float(smpl[6])/1000.)
#     tmpfile.Close()
#     #print '------------------------------------------------------------'


treesAndExps = {}
nGenTot = 0.
nRealGenTot = 0.
nEntriesTot = 0.

## Collect exponential parameters 
for cpl in coupls:
  if cpl != fs: continue
  for mass in masses:
    label = cpl+':'+mass
    treesAndExps[label] = [[], [], [], [], []]
    #for smpl in tomslist:
    for smpl in totlist:
      #indir = smpl[12] + ('' if smpl[12][-1]=='/' else '/')
      #filename = smpl[11].split('/')[-1]+'_'+year+'.root'
      if mass==smpl[2] and '_'+cpl+'_massive' in smpl[11]:
        if not os.path.exists(smpl[12]+smpl[11]):
          print ' >>> WARNING: file '+smpl[12]+smpl[11]+' does not exist <<<'
          continue
        treesAndExps[label][0].append(float(smpl[3])) # V^2
        treesAndExps[label][1].append(float(smpl[4])) # ctau, given in mm!!!
        treesAndExps[label][2].append(float(smpl[6])) # number of generated events
        #treesAndExps[label][3].append(smpl[11].split('/')[-1]+'_'+year+'.root')
        treesAndExps[label][3].append(smpl[11])
        treesAndExps[label][4].append(smpl[12])

    for ii in range(len(treesAndExps[label][3])): print treesAndExps[label][4][ii]+treesAndExps[label][3][ii]
    print '-----------------------'

print '------------------------------------------------------------'

# ## Build exponentials 
# for cpl in coupls:
#   for mass in masses:
#     label = cpl+':'+mass
#     nfiles = len(treesAndExps[label][0])
#     if nfiles==0: continue
#     fittformula = ''
#     for i in range(0, nfiles):
#       if i>0: fittformula += '+'
#       fittformula += ('['+str(2*i)+']*exp(-x/['+str(2*i+1)+'])')
#     fitt = r.TF1('fitt', fittformula      , 0., 1000.)
#     fit0 = r.TF1('fit0', '[0]*exp(-x/[1])', 0., 1000.)
#     ngentot = 0
#     for i in range(0, nfiles):
#       tmpngen = treesAndExps[label][2][i]
#       tmpctau = treesAndExps[label][1][i]
#       fitt.FixParameter(2*i  , tmpngen/tmpctau)
#       fitt.FixParameter(2*i+1, tmpctau)
#       ngentot += tmpngen
#     fit0.FixParameter(0, ngentot/treesAndExps[label][1][0])
#     fit0.FixParameter(1, treesAndExps[label][1][0])
#     #fit0.GetParameters()
#     #fitt.GetParameters()

# print '------------------------------------------------------------'

## Build new tree
histdir = 'blackJackAndHookers'
histnames = ['nVertices', 'hCounter', 'lheCounter', 'psCounter', 'tauCounter', 'nTrueInteractions']
nhists = len(histnames)
for cpl in coupls:
  if cpl != fs: continue
  for mass in masses:
    label = cpl+':'+mass
    nfiles = len(treesAndExps[label][0])
    if nfiles==0: continue
    nGenTot = 0.
    nRealGenTot = 0.
    nEntriesTot = 0.
    outputfile = treesAndExps[label][3][0].replace('.root', '_MERGED.root')
    mstr = outputfile.split('_')[2]
    if not '.' in mstr:
      outputfile = outputfile.replace(mstr, mstr+'.0')
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'Creating '+outputdir+outputfile
    fnew = r.TFile(outputdir+outputfile, 'recreate')
    #print 'New file '+fnew.GetName()+':  '+str(fnew)+' --- '+hex(id(fnew))
    dnew = fnew.mkdir(histdir)
    #print 'New dir '+dnew.GetName()+':  '+str(dnew)+' --- '+hex(id(dnew))
    dnew.cd()
    tchain = r.TChain(histdir+'/blackJackAndHookersTree')
    #tlist = r.TList()
    #print 'New chain '+tchain.GetName()+':  '+str(tchain)+' --- '+hex(id(tchain))
    hists = []
    for i in range(0, nfiles):
      #indir = treesAndExps[label][4][i]
      #indir += ('' if indir[-1]=='/' else '/')
      ftmp = r.TFile.Open(treesAndExps[label][4][i]+treesAndExps[label][3][i])
      print '  ~~ file '+str(i)+': '+ftmp.GetName()
      for j in range(0, nhists):
        htmp = r.TH1D(ftmp.Get(histdir+'/'+histnames[j]))
        #if j==1: print '  + %8d' % (htmp.GetBinContent(1))
        if i==0:
          dnew.cd()
          #print hists[j].InheritsFrom('TH1'),hists[j].GetName()
          hists.append(htmp.Clone(htmp.GetName()))
        else:
          #print hists[j].InheritsFrom('TH1'),hists[j].GetName()
          hists[j].Add(htmp)
        if j==1:
          genevts      = treesAndExps[label][2][i]
          realgenevts  = htmp.GetBinContent(1)
          nGenTot     += genevts
          nRealGenTot += realgenevts
          treesAndExps[label][2][i] = realgenevts
          #print '          + %8d   (= %8d)' % (realgenevts, nRealGenTot)
          print '         --> %8d   vs %8d' % (realgenevts, treesAndExps[label][2][i])
          if not (realgenevts>0. and abs(genevts-realgenevts)/realgenevts<0.01):
            print '             WATCH OUT: N. gen evts = %f vs %f' % (genevts, realgenevts)
      ttmp = ftmp.Get(histdir+'/blackJackAndHookersTree')
      nEntriesTot += ttmp.GetEntries()
      #tlist.Add(ttmp)
      ftmp.Close()
      tchain.Add(treesAndExps[label][4][i]+treesAndExps[label][3][i])
    print ' TOTAL   ==> %8d   vs %8d' % (nRealGenTot, nGenTot)
    dnew.cd()
    for h in hists: h.Write()
    #
    # Histogram with exopnential parameters
    dnew.cd()
    hexp = r.TH1D('mergingWeightParams', ';parameter index;parameter value', 80, 0., 80.)
    hexp.SetBinContent(1, 2*nfiles + 2)
    hexp.GetXaxis().SetBinLabel(1, 'N pars')
    hexp.SetBinContent(2, nRealGenTot/treesAndExps[label][1][0])
    #hexp.SetBinContent(2, nGenTot/treesAndExps[label][1][0])
    hexp.GetXaxis().SetBinLabel(2, 'par0')
    hexp.SetBinContent(3, treesAndExps[label][1][0])
    hexp.GetXaxis().SetBinLabel(3, 'par1')
    #for x in range(0, len(treesAndExps[label][1])):
    for x in range(0, nfiles):
      hexp.SetBinContent(2*x + 4, treesAndExps[label][2][x]/treesAndExps[label][1][x])
      hexp.GetXaxis().SetBinLabel(2*x + 4, 'par'+str(2*x + 2))
      hexp.SetBinContent(2*x + 5, treesAndExps[label][1][x])
      hexp.GetXaxis().SetBinLabel(2*x + 5, 'par'+str(2*x + 3))
      hexp.GetXaxis().SetBinLabel(1, 'N pars')
    dnew.cd()
    hexp.Write()
    print 'Cloning tree into '+fnew.GetName()+'...'
    dnew.cd()
    #tnew = tchain.CloneTree()
    tchain.Merge(fnew, 0, 'keep')
    #tnew = r.TTree.MergeTrees(tlist)
    #tnew.SetName('blackJackAndHookersTree')
    print '... cloned!'
    tnew = dnew.Get('blackJackAndHookersTree')
    print '  ===> Cloned tree entries = %d vs Sum(original trees) = %d' % (tnew.GetEntries(), nEntriesTot)
    #tnew.Write()
    #print '... and written!'
    #print 'New output tree '+tnew.GetName()+':  '+str(tnew)+' --- '+hex(id(tnew))
    # expweight = array('d', [0.])
    # print 'Exp weight array expweight:  '+str(expweight)+' --- '+hex(id(expweight))
    # mergeWeight = tnew.Branch('_mergeWeight', expweight, '_mergeWeight/D')
    # print 'New branch '+mergeWeight.GetName()+':  '+str(mergeWeight)+' --- '+hex(id(mergeWeight))
    # ntreeentries = tnew.GetEntries()
    # cnt = 0
    # nsteps = 50
    # print '[' + ' '*nsteps + ']'
    # for i in range(0, ntreeentries):
    #   if int(nsteps*i/ntreeentries)>cnt:
    #     cnt = int(nsteps*i/ntreeentries)
    #     print '[' + '-'*cnt + ' '*(nsteps-cnt) + ']'
    #   tnew.GetEntry(i)
    #   ctau = tnew.GetLeaf('_ctauHN').GetValue()
    #   expweight[0] = fit0.Eval(ctau)/fitt.Eval(ctau)
    #   mergeWeight.Fill()
    # print '[' + '-'*nsteps + ']'
    # print 'Writing tree to new file '+fnew.GetName()+'...'
    #tnew.Write('', r.TObject.kOverwrite)
    #print '... written!'
    print 'Closing file '+fnew.GetName()+'...'
    fnew.Close()
    print '... closed!'

## Central samples (from DAS)
##  command:
##    dasgoclient --query="dataset dataset=/HeavyNeutrino_trilepton*_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" | gawk '{print "dasgoclient --query=\"file dataset="$1" | sum(file.nevents)\""}' | csh 
# centrSamplesE = [
#   ['/HeavyNeutrino_trilepton_M-1_V-0_0949736805647_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'    , 250000],
#   ['/HeavyNeutrino_trilepton_M-1_V-0_212367605816_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'     , 242000],
#   ['/HeavyNeutrino_trilepton_M-2_V-0_0110905365064_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'    , 250000],
#   ['/HeavyNeutrino_trilepton_M-2_V-0_0248394846967_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'    , 250000],
#   ['/HeavyNeutrino_trilepton_M-3_V-0_00707813534767_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 241000],
#   ['/HeavyNeutrino_trilepton_M-4_V-0_00290516780927_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-5_V-0_00145602197786_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-6_V-0_00202484567313_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-8_V-0_00151327459504_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-10_V-0_000756967634711_e_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM' , 250000]
# ]

# centrSamplesMu = [
#   ['/HeavyNeutrino_trilepton_M-1_V-0_0949736805647_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'    , 250000],
#   ['/HeavyNeutrino_trilepton_M-1_V-0_212367605816_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'     , 250000],
#   ['/HeavyNeutrino_trilepton_M-2_V-0_0110905365064_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'    , 241000],
#   ['/HeavyNeutrino_trilepton_M-2_V-0_0248394846967_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'    , 250000],
#   ['/HeavyNeutrino_trilepton_M-3_V-0_00707813534767_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-4_V-0_00290516780927_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-5_V-0_00145602197786_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-6_V-0_00202484567313_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-8_V-0_00151327459504_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'   , 250000],
#   ['/HeavyNeutrino_trilepton_M-10_V-0_000756967634711_mu_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM' , 250000]
# ]

# centralSamples = {'e' : centrSamplesE, 'mu' : centrSamplesMu}

# for sample,nevents in centralSamples[fs]:
#   mass = float(sample[sample.index('_M-')+3:sample.index('_V-')])
#   vsqr = float(sample[sample.index('_V-')+3:sample.index('_'+fs+'_massive')].replace('_', '.'))**2
#   # h2.Fill(mass, vsqr, nevents/1000.)

# cc = r.TCanvas('cc', 'cc', 720, 600)
# cc.SetTopMargin(0.03)
# cc.SetLeftMargin(0.15)
# cc.SetRightMargin(0.18)
# cc.SetLogy()
# cc.SetGridx()
# cc.SetGridy()
# cc.cd()
# h2.Draw('COLZTEXT45')
# h2.GetYaxis().SetTitleOffset(1.10)
# h2.SetMarkerSize(1.5)
# cc.SaveAs('ll_hnl_samples_'+fs+'.png')
# cc.SaveAs('ll_hnl_samples_'+fs+'.pdf')
#cc.SaveAs('dummy.root')

