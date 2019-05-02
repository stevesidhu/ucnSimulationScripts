import ROOT
import math
import scipy.optimize
from ROOT import TCanvas, TGraph, TGraphErrors
from array import array
from ROOT import TF1
from ROOT import TH1

#productionRateReal = 1.8e7
#productionRateReal = 2.62e7
degaussingTime = 1800.
stableHours = 16.14
fillsPerCycle = 8
polarityTime = 150.
hbar = 6.71e-16
polarization = 0.5
EField = 12000.
tau_Xe = 6566.92
detEff = 0.9
alpha0 = 0.95
T2 = 500.
tWait = 2.
tPulse = 2.
T1 = 1000.
depolCollection = 0.95
pAnalyzer = 0.9
simTime = 200.
#simTime = 300.

# valveClosedTime: duration of irradiation while IV1 is closed
# fillTime: durattion of irradiation while IV1 is open (filling of cell)
# storageTime: duration of ramsey cycle
# emptyingTime: duration of collecting

def FillingSpec(histogram, valveClosedTime, fillTime):
  startBin = histogram.GetXaxis().FindBin(100. - valveClosedTime)
  fillBin = histogram.GetXaxis().FindBin(100. + fillTime)
  endBin = histogram.GetYaxis().FindBin(100. + fillTime)
  histogram.GetXaxis().SetRange(startBin, fillBin)
  histogram.GetYaxis().SetRange(endBin, endBin)
  h = histogram.Project3D('z')
  h.Sumw2()
  return h


def StorageProbSpec(histogram, storageTime):
  prob = histogram.Clone()
  for b in range(prob.GetNbinsX() + 2):
    tau = histogram.GetBinContent(b)
    if storageTime <= 0. or tau <= 0.:
      prob.SetBinContent(b, 0.)
      prob.SetBinError(b, 0.)
    else:
      p = math.exp(-storageTime/tau)
      prob.SetBinContent(b, p)
      prob.SetBinError(b, histogram.GetBinError(b)*p*storageTime/tau**2)
  prob.Sumw2()
  return prob


def EmptyingProbSpec(histogram, emptyingTime):
  emptyBin = histogram.GetXaxis().FindBin(emptyingTime)
  h = histogram.ProjectionY('collProb', 0, emptyBin)
  h.Sumw2()
  return h


#def daysToReach(histograms, valveClosedTime, fillTime, storageTime, emptyingTime):
def daysToReach(histograms, valveClosedTime, fillTime, storageTime, emptyingTime,productionRateReal = 1.8E7):
  Espec = {}
  for cell in ['top', 'bottom']:
    Espec[cell] = FillingSpec(histograms['fill_' + cell], valveClosedTime, fillTime)    
    storageHist = StorageProbSpec(histograms['storage_' + cell], storageTime)
    collProb = EmptyingProbSpec(histograms['emptying_' + cell], emptyingTime)
    Espec[cell].Multiply(storageHist)
    Espec[cell].Multiply(collProb)

  Espec['both'] = Espec['top']
  Espec['both'].Add(Espec['bottom'])
  dNcoll = ROOT.Double(0.)
  Ncoll = Espec['both'].IntegralAndError(0, -1, dNcoll)
  #dNcoll = math.sqrt(Ncoll)
  #print("Number collected: ", Ncoll, dNcoll)

  Ndet = Ncoll * detEff * polarization * productionRateReal /histograms['productionRate'] * math.exp(-storageTime/tau_Xe)
  if Ndet == 0.:
    return float('inf'), float('inf')
  dNdet = dNcoll/Ncoll*Ndet
  #print("Number detected: ", Ndet, dNdet)

  tEDM = storageTime - tWait - 2.*tPulse
  alpha = alpha0 * math.exp(-tEDM/T2 - (tWait + 2.*tPulse)/T1) * depolCollection * pAnalyzer
#  print(alpha)

  sensitivityPerFill = hbar/(2.*tEDM*EField*math.sqrt(Ndet)*alpha)
  dsensitivityPerFill = 0.5*dNdet/Ndet*sensitivityPerFill
#  print(sensitivityPerFill, dsensitivityPerFill)
  cycleTime = valveClosedTime + fillTime + storageTime + emptyingTime
  cyclesPerDay = stableHours*3600./(fillsPerCycle*cycleTime + 2.*polarityTime + degaussingTime/10.)
#  print(cyclesPerDay)
  cyclesToReach = (sensitivityPerFill/math.sqrt(fillsPerCycle)/1.e-27)**2
  dcyclesToReach = 2.*dsensitivityPerFill/sensitivityPerFill*cyclesToReach
  #print("Cycles to reach: ", cyclesToReach, dcyclesToReach)
  days = cyclesToReach/cyclesPerDay
  ddays = dcyclesToReach/cyclesPerDay
#  print(valveClosedTime, fillTime, storageTime, emptyingTime, days)
  return days, ddays, int(Ncoll*polarization * productionRateReal /histograms['productionRate'])


def readHistograms(fillingFile, topStorageFile, bottomStorageFile, topEmptyingFile, bottomEmptyingFile, cellCenter, cellCut):
  histograms = {}

  c = ROOT.TCanvas('c', 'c')
  fillFile = ROOT.TFile(fillingFile)
  fillTree = fillFile.Get('neutronsnapshot')

  if cellCut ==0:
      fillTree.Draw('Hend : tend : tstart >> topCell(200, 0, 200, 200, 0, 200, 50, 0, 250e-9)', 'xend>5.3 && zend<0.2 && zend>0')
     # fillTree.Draw('Hend : tend : tstart >> topCell(300, 0, 300, 300, 0, 300, 50, 0, 250e-9)', 'xend>5.3 && zend<0.2 && zend>0')
  if cellCut == 1:
      fillTree.Draw("Hend- {0}*1.025E-7: tend: tstart>>topCell(200,0,200, 200,0,200, 50,0,250e-9)".format(cellCenter), 'solidend>274 && solidend<276 && zend>{0}'.format(cellCenter))
     # fillTree.Draw("Hend- {0}*1.025E-7: tend: tstart>>topCell(300,0,300, 300,0,300, 50,0,250e-9)".format(cellCenter), 'solidend>274 && solidend<276 && zend>{0}'.format(cellCenter))
  histograms['fill_top'] = ROOT.gDirectory.Get('topCell')

  if cellCut ==0:
      fillTree.Draw('Hend : tend : tstart >> bottomCell(200, 0, 200, 200, 0, 200, 50, 0, 250e-9)', 'xend>5.3 && zend> -0.2 && zend<0')
     # fillTree.Draw('Hend : tend : tstart >> bottomCell(300, 0, 300, 300, 0, 300, 50, 0, 250e-9)', 'xend>5.3 && zend> -0.2 && zend<0')
  if cellCut == 1:
      fillTree.Draw("Hend-{0}*1.025E-7: tend: tstart>>bottomCell(200,0,200, 200,0,200, 50,0,250e-9)".format(cellCenter), 'solidend>274 && solidend<276 && zend<{0}'.format(cellCenter))
     # fillTree.Draw("Hend-{0}*1.025E-7: tend: tstart>>bottomCell(300,0,300, 300,0,300, 50,0,250e-9)".format(cellCenter), 'solidend>274 && solidend<276 && zend<{0}'.format(cellCenter))
  histograms['fill_bottom'] = ROOT.gDirectory.Get('bottomCell')

  tStorageFile = ROOT.TFile(topStorageFile)
  histograms['storage_top'] = tStorageFile.Get('lifetime_1')
  bStorageFile = ROOT.TFile(bottomStorageFile)
  histograms['storage_bottom'] = bStorageFile.Get('lifetime_1')

  tEmptyingFile = ROOT.TFile(topEmptyingFile)
  histograms['emptying_top'] = tEmptyingFile.Get('emptyEff')
  bEmptyingFile = ROOT.TFile(bottomEmptyingFile)
  histograms['emptying_bottom'] = bEmptyingFile.Get('emptyEff')

  for hist in histograms:
    histograms[hist].SetDirectory(0)
    histograms[hist].Sumw2()
  histograms['productionRate'] = fillFile.Get('neutronend').GetEntries()/simTime

  return histograms

def importantValues(fillingFile,valveClosedTime,fillingTime,cellStorageTime,emptyingTime):
 
  c = ROOT.TCanvas('c', 'c')
  fillFile = ROOT.TFile(fillingFile)
  fillTree = fillFile.Get('neutronsnapshot')
  
  fillTree.Draw("tend-tstart>>lifetime(90,5,95)","tstart<5")
  lifetime =ROOT.gDirectory.Get("lifetime")
  g2= ROOT.TF1("fit","expo",5,95)  #exponential fit to calculate tau
  lifetime.Fit(g2,"R")
  sourceTau = -1./g2.GetParameter(1)
  sourceTauError = g2.GetParError(1)/(math.pow(g2.GetParameter(1),2))

  
  fillTree.Draw("tend-tstart>>sysLifetime","tstart>100 && tstart<105")
  sysLifetime =ROOT.gDirectory.Get("sysLifetime")
  g3= ROOT.TF1("fit","expo",0,95)  #exponential fit to calculate tau
  sysLifetime.Fit(g3,"R")
  sysTau = -1./g3.GetParameter(1)
  sysTauError = g3.GetParError(1)/(math.pow(g3.GetParameter(1),2))
  return sourceTau, sourceTauError, sysTau, sysTauError

#'''
#Files stored on angerona_data, all simulations had a smaller cell for filling, but using a slightly larger one here for storage and emptying
#'''

#V3 cryostat; standard (0.05 mm thick foil; solidend ==3)
extractionHeight=[
           {'label':'Center of guide from bottle of He bottle (cm)','parameter': 0.5,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/low_extraction.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 7.5,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/middle_extraction.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 7.5,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/double_kink.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':0, 'filenames':  ('/ucn/angerona_data/newCryostatSims/testAnalysis/standard.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#V3 cryostat; standard (0.05 mm thick foil; solidend == 3)
sourceTube=[
           {'label':'Source tube diamater (cm)','parameter': 10,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/10cmTube.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 12.5,'cellCut':0,'filenames':('/ucn/angerona_data/newCryostatSims/testAnalysis/12.5cmTube.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/standard.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 18,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/18cmTube.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#V3 cryostat; standard (0.05 mm thick foil; solidend == 3)
HEXdiameter=[
           {'label':'Diameter of HEX (cm)','parameter': 12.,'cellCut':0, 'filenames':('/ucn/angerona_data/newCryostatSims/testAnalysis/HEX12cm_75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/standard.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 18,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/HEX18cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 20,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/HEX20.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 12.5,'cellCut':0,'filenames':('/ucn/angerona_data/newCryostatSims/testAnalysis/HEX_12cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#kink height study: 10cm overfill and 50cm long HEX
KinkHeight=[
           {'label':'Kink height (cm)','parameter': 35,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/newKinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/newKinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/newKinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/newKinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/newKinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('newKinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('newKinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#noBfield but foil vs kink height: 10cm overfill and 50cm long HEX
noBKinkHeight=[
           {'label':'Kink height (cm)','parameter': 35,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 80,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight80.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBkinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#no bfield or foil vs kink height:10cm overfill and 50cm long HEX
noBnoFoilKinkHeight=[
           {'label':'Kink height (cm)','parameter': 35,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 80,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight80.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/noBFkinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#double bfield coild diamater vs kink height: 10cm overf/ucn/angerona_data/newCryostatSims/testAnalysis/ill and 50cm long HEX
doubleBKinkHeight=[
           {'label':'Kink height (cm)','parameter': 35,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/2BkinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/2BkinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/2BkinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/2BkinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/2BkinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/2BkinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/2BkinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#Hard bore material V5 Cryostat
hardBore=[
           {'label':'Diameter of Ni58 bore guidie (cm)','parameter': 7,'cellCut':0, 'filenames':  ('/ucn/angerona_data/newCryostatSims/testAnalysis/70HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 8.5,'cellCut':0, 'filenames':('/ucn/angerona_data/newCryostatSims/testAnalysis/85HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 10,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/100HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 12,'cellCut':0, 'filenames': ('/ucn/angerona_data/newCryostatSims/testAnalysis/125HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]

#'''
#Files stored on orithiya_data
#'''
boreDiameter=[
           {'label':'Diameter of Ni58 bore guide (mm)','parameter': 60,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/boreDiameter60.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/boreDiameter70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/boreDiameter75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/boreDiameter85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/boreDiameter100.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/boreDiameter125.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 150,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/boreDiameter150.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
doubleBboreDiameter=[
           {'label':'Diameter of Ni58 bore guide (mm)','parameter': 60,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/2BboreDiameter60.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/2BboreDiameter70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/2BboreDiameter75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames':  ('/ucn/orithyia_data/ssidhu2/cryoSimulations/2BboreDiameter85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/2BboreDiameter100.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/2BboreDiameter125.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 150,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/2BboreDiameter150.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
angTransferGuide=[
           {'label':'Height of cell-center wrt bottom of He bottle (cm)','parameter': 45,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain45cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.18},  
           {'parameter': 45,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain45cmKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain65cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.38},  
           {'parameter': 65,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain65cmKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain75cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.48},  
           {'parameter': 75,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain75cmKinkHt75cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain95cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.68},  
           {'parameter': 95,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/angTransGain95cmKinkHt95cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
funnelShapes=[
           {'label':'Funnel Shape','parameter': 1,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFunKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 2,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFunKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 3,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/eccFunKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 4,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/eccFunKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 5,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/newKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 6,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/newKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
kinkRadius1=[
           {'label':'Radius of kink (cm)','parameter':  8,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad8cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad15cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 22,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad22cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 29,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad29cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 36,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad36cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad43cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad50cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 57,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad57cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
kinkRadius2=[
           {'label':'Radius of kink (cm)','parameter':  8,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V8cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V15cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 22,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V22cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 29,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V29cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 36,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V36cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V43cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V50cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 57,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkRad2V57cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
kinkAngle1=[
           {'label':'Kink angle (deg)','parameter': 45,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 55,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 60,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng60deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 65,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 70,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 75,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 80,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 85,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 95,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           ]
kinkAngle2=[
           {'label':'Kink Angle (deg)','parameter': 45,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 55,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 60,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V60deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 65,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 70,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 75,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 80,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 85,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 95,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           ]
minSepFeeders1=[
           {'label':'Kink height (cm)','parameter': 35,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepKinkHt2V35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptTop_DetEff.root','/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptBottom_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepKinkHt2V45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptTop_DetEff.root','/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 55,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepKinkHt2V55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptTop_DetEff.root','/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 65,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepKinkHt2V65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptTop_DetEff.root','/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 75,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepKinkHt2V75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptTop_DetEff.root','/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 85,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepKinkHt2V85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptTop_DetEff.root','/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 95,'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepKinkHt2V95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptTop_DetEff.root','/ucn/orithyia_data/ssidhu2/cryoSimulations/minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           ]
asymmetricFeeders1=[
           {'label':'Offset of top cell from cell center (cm)','parameter':  0,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset0cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root',  '/ucn/orithyia_data/ssidhu2/cryoSimulations/kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter':  5,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset5cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root',  '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptTop5cm_DetEff.root',   '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptBottom5cm_DetEff.root'), 'cellCenter': 0.05},
           {'parameter': 10,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset10cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptTop10cm_DetEff.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptBottom10cm_DetEff.root'),'cellCenter': 0.1},
           {'parameter': 15,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset15cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptTop15cm_DetEff.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptBottom15cm_DetEff.root'),'cellCenter': 0.15},
           {'parameter': 20,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset20cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptTop20cm_DetEff.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptBottom20cm_DetEff.root'),'cellCenter': 0.2},
           {'parameter': 25,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset25cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptTop25cm_DetEff.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptBottom25cm_DetEff.root'),'cellCenter': 0.25},
           {'parameter': 30,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset30cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptTop30cm_DetEff.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptBottom30cm_DetEff.root'),'cellCenter': 0.3},
           {'parameter': 35,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedOffset35cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptTop35cm_DetEff.root', '/ucn/orithyia_data/ssidhu2/cryoSimulations/asymFeedEmptBottom35cm_DetEff.root'),'cellCenter': 0.35},
           ]
camSpider =[
           {'label':'Spider set up','parameter': 1, 'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/noSpider.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 2, 'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/camSpider.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 3, 'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/camSpider90deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 4, 'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/camSpider90degFurther.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 5, 'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/camSpider0cmCap.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 6, 'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/camSpider5cmCap.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 7, 'cellCut':0, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/camSpider10cmCap.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0}
           ]

guideDiameter =[
           {'label':'Diamter of UCN guides (mm)','parameter': 85, 'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/guideDiameterBigger85mm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95, 'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/guideDiameterBigger95mm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100, 'cellCut':1, 'filenames':('/ucn/orithyia_data/ssidhu2/cryoSimulations/guideDiameterBigger100mm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125, 'cellCut':1, 'filenames':('/ucn/orithyia_data/ssidhu2/cryoSimulations/guideDiameterBigger125mm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
               ] 

#'''
#Files stored on \ucn\data\ssidhu
#'''
cellRadius =[
           {'label':'Radius of cells (cm)','parameter': 18, 'cellCut':1, 'filenames': ('cellRad18cm.root', 'cellRadStoreTop18cm_hist.root', 'cellRadStoreBottom18cm_hist.root', 'cellRadEmptTop18cm_DetEff.root', 'cellRadEmptBottom18cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 20, 'cellCut':1, 'filenames': ('cellRad20cm.root', 'cellRadStoreTop20cm_hist.root', 'cellRadStoreBottom20cm_hist.root', 'cellRadEmptTop20cm_DetEff.root', 'cellRadEmptBottom20cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 22, 'cellCut':1, 'filenames': ('cellRad22cm.root', 'cellRadStoreTop22cm_hist.root', 'cellRadStoreBottom22cm_hist.root', 'cellRadEmptTop22cm_DetEff.root', 'cellRadEmptBottom22cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 25, 'cellCut':1, 'filenames': ('cellRad25cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 28, 'cellCut':1, 'filenames': ('cellRad28cm.root', 'cellRadStoreTop28cm_hist.root', 'cellRadStoreBottom28cm_hist.root', 'cellRadEmptTop28cm_DetEff.root', 'cellRadEmptBottom28cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 30, 'cellCut':1, 'filenames': ('cellRad30cm.root', 'cellRadStoreTop30cm_hist.root', 'cellRadStoreBottom30cm_hist.root', 'cellRadEmptTop30cm_DetEff.root', 'cellRadEmptBottom30cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 40, 'cellCut':1, 'filenames': ('cellRad40cm.root', 'cellRadStoreTop40cm_hist.root', 'cellRadStoreBottom40cm_hist.root', 'cellRadEmptTop40cm_DetEff.root', 'cellRadEmptBottom40cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50, 'cellCut':1, 'filenames': ('cellRad50cm.root', 'cellRadStoreTop50cm_hist.root', 'cellRadStoreBottom50cm_hist.root', 'cellRadEmptTop50cm_DetEff.root', 'cellRadEmptBottom50cm_DetEff.root'), 'cellCenter': 0},  
           ]

moderator =[
           {'label':'Moderator simulation','parameter': 11, 'cellCut':1, 'filenames': ('moderatorSimsAl2219cylAlBeMet.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':17.9E6},  
           {'parameter': 10, 'cellCut':1, 'filenames': ('moderatorSimsAl2219cyl.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':16.0E6},
           {'parameter': 4, 'cellCut':1, 'filenames': ('moderatorSimsAl6061cylAlBeMetCloser.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':19.1E6},
           {'parameter': 2, 'cellCut':1, 'filenames': ('moderatorSimsAl6061cylAlBeMet.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':17.8E6},
           {'parameter': 1, 'cellCut':1, 'filenames': ('moderatorSimsAl6061cyl.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':15.7E6},
           {'parameter': 5, 'cellCut':1, 'filenames': ('moderatorSimsAl6061sph.root','kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':15.4E6},
           {'parameter': 6, 'cellCut':1, 'filenames': ('moderatorSimsAl6061sphAlBeMet.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':17.3E6},
           {'parameter': 7, 'cellCut':1, 'filenames': ('moderatorSimsAl6061sphAlBeMetCloser.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':18.6E6},
           {'parameter': 3, 'cellCut':1, 'filenames': ('moderatorSimsAl2219cylCloser.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':17E6},
           {'parameter': 9, 'cellCut':1, 'filenames': ('moderatorSimsAl6061sphAlBeMetFurther.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':15.6E6},
           {'parameter': 8, 'cellCut':1, 'filenames': ('moderatorSimsAl6061sphMoreLead.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0, 'ucnProduction':14E6}
             ]

guideDiameterBiggerLD3 =[
           {'label':'Diameter of UCN guides (mm)','parameter': 85, 'cellCut':1, 'filenames': ('guideDiamBigLamb385mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95, 'cellCut':1, 'filenames': ('guideDiamBigLamb395mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100, 'cellCut':1, 'filenames': ('guideDiamBigLamb3100mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125, 'cellCut':1, 'filenames': ('guideDiamBigLamb3125mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
               ] 
guideDiameterBiggerLD1 =[
           {'label':'Diameter of UCN guides (mm)','parameter': 85, 'cellCut':1, 'filenames': ('guideDiamBigLamb185mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95, 'cellCut':1, 'filenames': ('guideDiamBigLamb195mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100, 'cellCut':1, 'filenames': ('guideDiamBigLamb1100mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125, 'cellCut':1, 'filenames': ('guideDiamBigLamb1125mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
               ] 
guideDiameterBiggerLD5 =[
           {'label':'Diameter of UCN guides (mm)','parameter': 85, 'cellCut':1, 'filenames': ('guideDiamBigLamb585mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95, 'cellCut':1, 'filenames': ('guideDiamBigLamb595mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100, 'cellCut':1, 'filenames': ('guideDiamBigLamb5100mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125, 'cellCut':1, 'filenames': ('guideDiamBigLamb5125mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root','cellRadEmptTop25cm_DetEff.root ', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
               ] 

kinkRiseBiggerLD3=[
           {'label':'Kink rise above LHe level (cm)','parameter': 18,'cellCut':1, 'filenames': ('kinkHtBigger18cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 23,'cellCut':1, 'filenames': ('kinkHtBigger23cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 28,'cellCut':1, 'filenames': ('kinkHtBigger28cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 33,'cellCut':1, 'filenames': ('kinkHtBigger33cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 38,'cellCut':1, 'filenames': ('kinkHtBigger38cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('kinkHtBigger43cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
          # {'parameter':48,'cellCut':1, 'filenames': ('kinkHtBigger48cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           ]
kinkRiseBiggerLD1=[
           {'label':'Kink rise above LHe level (cm)','parameter': 18,'cellCut':1, 'filenames': ('kinkHtBiggerLamb118cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 23,'cellCut':1, 'filenames': ('kinkHtBiggerLamb123cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 28,'cellCut':1, 'filenames': ('kinkHtBiggerLamb128cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 33,'cellCut':1, 'filenames': ('kinkHtBiggerLamb133cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 38,'cellCut':1, 'filenames': ('kinkHtBiggerLamb138cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('kinkHtBiggerLamb143cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           ]
kinkRiseBiggerLD5=[
           {'label':'Kink rise above LHe level (cm)','parameter': 18,'cellCut':1, 'filenames': ('kinkHtBiggerLamb518cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 23,'cellCut':1, 'filenames': ('kinkHtBiggerLamb523cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 28,'cellCut':1, 'filenames': ('kinkHtBiggerLamb528cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 33,'cellCut':1, 'filenames': ('kinkHtBiggerLamb533cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 38,'cellCut':1, 'filenames': ('kinkHtBiggerLamb538cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('kinkHtBiggerLamb543cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           ]

kinkRadiusBiggerLD3=[
           {'label':'Kink radius (cm)','parameter': 18,'cellCut':1, 'filenames': ('kinkHtBigger18cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 23,'cellCut':1, 'filenames': ('kinkHtBigger23cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 28,'cellCut':1, 'filenames': ('kinkHtBigger28cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 33,'cellCut':1, 'filenames': ('kinkHtBigger33cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 38,'cellCut':1, 'filenames': ('kinkHtBigger38cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('kinkHtBigger43cm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           ]

KG35HEXDiameter=[
           {'label':'35KG HEX diameter (mm)','parameter': 125,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/HEXDiameter125mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 148,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/HEXDiameter148mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 180,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/HEXDiameter180mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 200,'cellCut':1, 'filenames': ('/ucn/orithyia_data/ssidhu2/cryoSimulations/HEXDiameter200mm.root', 'cellRadStoreTop25cm_hist.root', 'cellRadStoreBottom25cm_hist.root', 'cellRadEmptTop25cm_DetEff.root', 'cellRadEmptBottom25cm_DetEff.root'), 'cellCenter': 0}  
           ]

ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kError

#Make a dictionary here called studyList: key = studayName:

studyList1 = {'ExtractionHeight':extractionHeight, 'SourceTubeDiameter':sourceTube, 'HEXDiameter':HEXdiameter, 'BoreDiameter' :boreDiameter, 'DoubleBBoreDiameter':doubleBboreDiameter, 'AngledTransferGuide':angTransferGuide, 'FunnelShapes':funnelShapes, 'KinkRadiusV1':kinkRadius1, 'KinkRadiusV2':kinkRadius2, 'KinkAngleV1':kinkAngle1, 'KinkAngleV2':kinkAngle2, 'MinimalSeperationFeeders':minSepFeeders1, 'AsymmetricFeeders':asymmetricFeeders1, 'CamSpider':camSpider, 'KinkRiseBiggerLD3':kinkRiseBiggerLD3,'KinkRiseBiggerLD1':kinkRiseBiggerLD1,'KinkRiseBiggerLD5':kinkRiseBiggerLD5,'KinkRadiusBiggerLD3':kinkRadiusBiggerLD3,'GuideDiameterBiggerLD3':guideDiameterBiggerLD3,'GuideDiameterBiggerLD1':guideDiameterBiggerLD1,'GuideDiameterBiggerLD5':guideDiameterBiggerLD5}
#, 'Kink Height':KinkHeight, 'NoBKinkHeight':noBKinkHeight, 'NoBnoFKinkHeight':noBnoFoilKinkHeight, 'DoubleBKinkHeight':doubleBKinkHeight, 'HardBore':hardBore}
studyList2 = {'CellRadius':cellRadius}  #change simNumber and limits on optimization as well as binning on Hend histograms
studyList3 = {'Moderator':moderator}   #change days function parameters to include productio number
studyList4 = {'KG35HEXDiameter':KG35HEXDiameter}

for studyKey,studyValue in studyList4.items():
    studyDays, studyDaysErrors, studyParameters, zeroArray  = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
    studies = studyList4[studyKey]
    f = open(studyKey+'.txt', mode='w')
    f.write("FileName,TimeBeforeValveOpen,FillingTime,TotalSourcePumpingTime,CellStorageTime,EmptyingTime,CycleTime,DaysToReach,ErrorDaysToReach,NumberCollected,SourceTau,SourceTauError,SystemTau,SystemTauError\n")
    
    for study in studies:
      histograms = readHistograms(study['filenames'][0], study['filenames'][1], study['filenames'][2], study['filenames'][3], study['filenames'][4], study['cellCenter'], study['cellCut'])
     #for cell radius studies change to 300 s 
      optimized = scipy.optimize.differential_evolution(lambda x: daysToReach(histograms, x[0], x[1], x[2], x[3])[0], [[0., 100], [0., 100.], [0., 300.], [0., 250.]], polish=False)
     # optimized = scipy.optimize.differential_evolution(lambda x: daysToReach(histograms, x[0], x[1], x[2], x[3])[0], [[0., 100], [0., 200.], [0., 300.], [0., 250.]], polish=False)
      xmin = optimized.x
      nit = optimized.nit
      print(study['parameter'], xmin, daysToReach(histograms, xmin[0], xmin[1], xmin[2], xmin[3]), nit)
     # daysResult = daysToReach(histograms, xmin[0], xmin[1], xmin[2], xmin[3],study['ucnProduction'])
      daysResult = daysToReach(histograms, xmin[0], xmin[1], xmin[2], xmin[3])
      studyDays.append(daysResult[0])
      print(daysResult[0])
      studyDaysErrors.append(daysResult[1])
      zeroArray.append(0.)
      studyParameters.append(study['parameter'])
      studyProperties = importantValues(study['filenames'][0], xmin[0],xmin[1],xmin[2],xmin[3])  
      print("Source Tau: ",studyProperties[0], " +/- ", studyProperties[1], " System Tau: " ,studyProperties[2], " +/- ", studyProperties[3])
      f.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}\n".format(study['filenames'][0], xmin[0], xmin[1], xmin[0]+xmin[1] , xmin[2] , xmin[3] , +xmin[0] + xmin[1] + xmin[2] + xmin[3] ,daysResult[0],daysResult[1],daysResult[2],studyProperties[0],studyProperties[1],studyProperties[2],studyProperties[3]))

    f.close()
    c = ROOT.TCanvas('c','c')
    gr = TGraphErrors(len(studyParameters), studyParameters, studyDays, zeroArray, studyDaysErrors);
    gr.SetMarkerStyle(20);
    gr.SetTitle("Days to reach vs " + studyKey) ; 
    gr.GetXaxis().SetTitle(studyValue[0]['label']);
    gr.GetYaxis().SetTitle("Days to reach");
    gr.GetXaxis().SetTitleSize(0.05);
    gr.GetXaxis().SetTitleOffset(0.75);
    gr.GetYaxis().SetTitleSize(0.05);
    gr.GetYaxis().SetTitleOffset(0.9);
    gr.Draw("Ap");
    c.Print(studyKey+"DaysToReachVsP.png");
    

'''
c = ROOT.TCanvas('c','c')
f1 = ROOT.TF2('closed-fill', lambda x: daysToReach(histograms, x[0], x[1], xmin[2], xmin[3])[0], 0, 100, 10, 100)
f1.SetMaximum(1000)
f1.Draw('COLZ')
c.Print('days.pdf(')
f2 = ROOT.TF2('fill-storage', lambda x: daysToReach(histograms, xmin[0], x[0], x[1], xmin[3])[0], 10, 100, 10, 300)
f2.SetMaximum(1000)
f2.Draw('COLZ')
c.Print('days.pdf')
f3 = ROOT.TF2('storage-empty', lambda x: daysToReach(histograms, xmin[0], xmin[1], x[0], x[1])[0], 10, 300, 10, 250)
f3.SetMaximum(1000)
f3.Draw('COLZ')
c.Print('days.pdf')
f4 = ROOT.TF2('fill-empty', lambda x: daysToReach(histograms, xmin[0], x[0], xmin[2], x[1])[0], 10, 100, 10, 250)
f4.SetMaximum(1000)
f4.Draw('COLZ')
c.Print('days.pdf)')

for cell in ['top','bottom']:
  for valveClosedTime in range(100):
    FillingSpec(histograms['fill_' + cell], valveClosedTime, xmin[1]).Draw()
    if valveClosedTime == 0:
      c.Print(cell + '_valve.pdf(')
    elif valveClosedTime == 99:
      c.Print(cell + '_valve.pdf)')
    else:
      c.Print(cell + '_valve.pdf')
  for fillTime in range(100):
    FillingSpec(histograms['fill_' + cell], xmin[0], fillTime).Draw()
    if fillTime == 0:
      c.Print(cell + '_fill.pdf(')
    elif fillTime == 99:
      c.Print(cell + '_fill.pdf)')
    else:
      c.Print(cell + '_fill.pdf')
  for storageTime in range(300):
    f = FillingSpec(histograms['fill_' + cell], xmin[0], xmin[1])
    f.Multiply(StorageProbSpec(histograms['storage_' + cell], storageTime))
    f.Draw()
    if storageTime == 0:
      c.Print(cell + '_storage.pdf(')
    elif storageTime == 299:
      c.Print(cell + '_storage.pdf)')
    else:
      c.Print(cell + '_storage.pdf')
  for emptyingTime in range(250):
    f = FillingSpec(histograms['fill_' + cell], xmin[0], xmin[1])
    f.Multiply(StorageProbSpec(histograms['storage_' + cell], storageTime))
    f.Multiply(EmptyingProbSpec(histograms['emptying_' + cell], emptyingTime))
    if emptyingTime == 0:
      c.Print(cell + '_emptying.pdf(')
    elif emptyingTime == 249:
      c.Print(cell + '_emptying.pdf)')
    else:
      c.Print(cell + '_emptying.pdf')
'''
