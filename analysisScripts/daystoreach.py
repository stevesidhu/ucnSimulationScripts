import ROOT
import math
import scipy.optimize
from ROOT import TCanvas, TGraph, TGraphErrors
from array import array

productionRateReal = 1.8e7
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


def daysToReach(histograms, valveClosedTime, fillTime, storageTime, emptyingTime):
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
#  print(Ncoll, dNcoll)

  Ndet = Ncoll * detEff * polarization * productionRateReal/histograms['productionRate'] * math.exp(-storageTime/tau_Xe)
  if Ndet == 0.:
    return float('inf'), float('inf')
  dNdet = dNcoll/Ncoll*Ndet
#  print(Ndet, dNdet)

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
#  print(cyclesToReach, dcyclesToReach)
  days = cyclesToReach/cyclesPerDay
  ddays = dcyclesToReach/cyclesPerDay
#  print(valveClosedTime, fillTime, storageTime, emptyingTime, days)
  return days, ddays


def readHistograms(fillingFile, topStorageFile, bottomStorageFile, topEmptyingFile, bottomEmptyingFile, cellCenter, cellCut):
  histograms = {}

  c = ROOT.TCanvas('c', 'c')


  fillFile = ROOT.TFile(fillingFile)
  fillTree = fillFile.Get('neutronsnapshot')

  if cellCut ==0:
      fillTree.Draw('Hend : tend : tstart >> topCell(200, 0, 200, 200, 0, 200, 50, 0, 250e-9)', 'xend>5.3 && zend<0.2 && zend>0')
  if cellCut == 1:
      fillTree.Draw("Hend- {0}*1.025E-7: tend: tstart>>topCell(200,0,200, 200,0,200, 50,0,250e-9)".format(cellCenter), 'solidend>274 && solidend<276 && zend>{0}'.format(cellCenter))
  histograms['fill_top'] = ROOT.gDirectory.Get('topCell')

  if cellCut ==0:
      fillTree.Draw('Hend : tend : tstart >> bottomCell(200, 0, 200, 200, 0, 200, 50, 0, 250e-9)', 'xend>5.3 && zend> -0.2 && zend<0')
  if cellCut == 1:
      fillTree.Draw("Hend-{0}*1.025E-7: tend: tstart>>bottomCell(200,0,200, 200,0,200, 50,0,250e-9)".format(cellCenter), 'solidend>274 && solidend<276 && zend<{0}'.format(cellCenter))
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
  histograms['productionRate'] = fillFile.Get('neutronend').GetEntries()/200.

  return histograms


'''
Files stored on angerona_data, all simulations had a smaller cell for filling, but using a slightly larger one here for storage and emptying
'''
#V3 cryostat; standard (0.05 mm thick foil; solidend ==3)
extractionHeight=[
           {'parameter': 0.5,'cellCut':0, 'filenames': ('low_extraction.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 7.5,'cellCut':0, 'filenames': ('middle_extraction.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 7.5,'cellCut':0, 'filenames': ('double_kink.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':0, 'filenames': ('standard.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#V3 cryostat; standard (0.05 mm thick foil; solidend == 3)
sourceTube=[
           {'parameter': 10,'cellCut':0, 'filenames': ('10cmTube.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 12.5,'cellCut':0, 'filenames': ('12.5cmTube.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':0, 'filenames': ('standard.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 18,'cellCut':0, 'filenames': ('18cmTube.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#V3 cryostat; standard (0.05 mm thick foil; solidend == 3)
HEXdiameter=[
           {'parameter': 12.,'cellCut':0, 'filenames': ('HEX12cm_75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':0, 'filenames': ('standard.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 18,'cellCut':0, 'filenames': ('HEX18cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 20,'cellCut':0, 'filenames': ('HEX20.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 12.5,'cellCut':0, 'filenames': ('HEX_12cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#kink height study: 10cm overfill and 50cm long HEX
KinkHeight=[
           {'parameter': 35,'cellCut':0, 'filenames': ('newKinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('newKinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('newKinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('newKinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('newKinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('newKinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('newKinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#noBfield but foil vs kink height: 10cm overfill and 50cm long HEX
noBKinkHeight=[
           {'parameter': 35,'cellCut':0, 'filenames': ('noBkinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('noBkinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('noBkinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('noBkinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames': ('noBkinkHeight70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('noBkinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 80,'cellCut':0, 'filenames': ('noBkinkHeight80.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('noBkinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('noBkinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#no bfield or foil vs kink height:10cm overfill and 50cm long HEX
noBnoFoilKinkHeight=[
           {'parameter': 35,'cellCut':0, 'filenames': ('noBFkinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('noBFkinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('noBFkinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('noBFkinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames': ('noBFkinkHeight70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('noBFkinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 80,'cellCut':0, 'filenames': ('noBFkinkHeight80.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('noBFkinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('noBFkinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#double bfield coild diamater vs kink height: 10cm overfill and 50cm long HEX
doubleBKinkHeight=[
           {'parameter': 35,'cellCut':0, 'filenames': ('2BkinkHeight35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('2BkinkHeight45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 55,'cellCut':0, 'filenames': ('2BkinkHeight55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':0, 'filenames': ('2BkinkHeight65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('2BkinkHeight75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('2BkinkHeight85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':0, 'filenames': ('2BkinkHeight95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
#Hard bore material V5 Cryostat
hardBore=[
           {'parameter': 7,'cellCut':0, 'filenames': ('70HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 8.5,'cellCut':0, 'filenames': ('85HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 10,'cellCut':0, 'filenames': ('100HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 12,'cellCut':0, 'filenames': ('125HardBore.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]

'''
Files stored on orithiya
'''
boreDiameter=[
           {'parameter': 60,'cellCut':0, 'filenames': ('boreDiameter60.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames': ('boreDiameter70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('boreDiameter75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('boreDiameter85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100,'cellCut':0, 'filenames': ('boreDiameter100.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125,'cellCut':0, 'filenames': ('boreDiameter125.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 150,'cellCut':0, 'filenames': ('boreDiameter150.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
doubleBboreDiameter=[
           {'parameter': 60,'cellCut':0, 'filenames': ('2BboreDiameter60.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 70,'cellCut':0, 'filenames': ('2BboreDiameter70.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':0, 'filenames': ('2BboreDiameter75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 85,'cellCut':0, 'filenames': ('2BboreDiameter85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 100,'cellCut':0, 'filenames': ('2BboreDiameter100.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 125,'cellCut':0, 'filenames': ('2BboreDiameter125.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 150,'cellCut':0, 'filenames': ('2BboreDiameter150.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
angTransferGuide=[
           {'parameter': 45,'cellCut':1, 'filenames': ('angTransGain45cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.18},  
           {'parameter': 45,'cellCut':1, 'filenames': ('angTransGain45cmKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 65,'cellCut':1, 'filenames': ('angTransGain65cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.38},  
           {'parameter': 65,'cellCut':1, 'filenames': ('angTransGain65cmKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 75,'cellCut':1, 'filenames': ('angTransGain75cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.48},  
           {'parameter': 75,'cellCut':1, 'filenames': ('angTransGain75cmKinkHt75cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 95,'cellCut':1, 'filenames': ('angTransGain95cmKinkHt27cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0.68},  
           {'parameter': 95,'cellCut':1, 'filenames': ('angTransGain95cmKinkHt95cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
funnelShapes=[
           {'parameter': 1,'cellCut':0, 'filenames': ('asymFunKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 2,'cellCut':0, 'filenames': ('asymFunKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 3,'cellCut':0, 'filenames': ('eccFunKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 4,'cellCut':0, 'filenames': ('eccFunKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 5,'cellCut':0, 'filenames': ('newKinkHt45cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 6,'cellCut':0, 'filenames': ('newKinkHt65cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
kinkRadius1=[
           {'parameter':  8,'cellCut':1, 'filenames': ('kinkRad8cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':1, 'filenames': ('kinkRad15cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 22,'cellCut':1, 'filenames': ('kinkRad22cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 29,'cellCut':1, 'filenames': ('kinkRad29cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 36,'cellCut':1, 'filenames': ('kinkRad36cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('kinkRad43cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('kinkRad50cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 57,'cellCut':1, 'filenames': ('kinkRad57cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
kinkRadius2=[
           {'parameter':  8,'cellCut':1, 'filenames': ('kinkRad2V8cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 15,'cellCut':1, 'filenames': ('kinkRad2V15cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 22,'cellCut':1, 'filenames': ('kinkRad2V22cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 29,'cellCut':1, 'filenames': ('kinkRad2V29cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 36,'cellCut':1, 'filenames': ('kinkRad2V36cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 43,'cellCut':1, 'filenames': ('kinkRad2V43cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('kinkRad2V50cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 57,'cellCut':1, 'filenames': ('kinkRad2V57cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           ]
kinkAngle1=[
           {'parameter': 45,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 55,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 60,'cellCut':1, 'filenames': ('kinkAng60deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 65,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 70,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 75,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 80,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 85,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 95,'cellCut':1, 'filenames': ('kinkAng45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           ]
kinkAngle2=[
           {'parameter': 45,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 50,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 55,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 60,'cellCut':1, 'filenames': ('kinkAng2V60deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 65,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 70,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 75,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 80,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 85,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 95,'cellCut':1, 'filenames': ('kinkAng2V45deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           ]
minSepFeeders1=[
           {'parameter': 35,'cellCut':0, 'filenames': ('minSepKinkHt2V35.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'minSepEmptTop_DetEff.root','minSepEmptBottom_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 45,'cellCut':0, 'filenames': ('minSepKinkHt2V45.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'minSepEmptTop_DetEff.root','minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 55,'cellCut':0, 'filenames': ('minSepKinkHt2V55.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'minSepEmptTop_DetEff.root','minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 65,'cellCut':0, 'filenames': ('minSepKinkHt2V65.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'minSepEmptTop_DetEff.root','minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 75,'cellCut':0, 'filenames': ('minSepKinkHt2V75.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'minSepEmptTop_DetEff.root','minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 85,'cellCut':0, 'filenames': ('minSepKinkHt2V85.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'minSepEmptTop_DetEff.root','minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           {'parameter': 95,'cellCut':0, 'filenames': ('minSepKinkHt2V95.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'minSepEmptTop_DetEff.root','minSepEmptBottom_DetEff.root'), 'cellCenter': 0},
           ]
asymmetricFeeders1=[
           {'parameter':  0,'cellCut':1, 'filenames': ('asymFeedOffset0cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter':  5,'cellCut':1, 'filenames': ('asymFeedOffset5cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'asymFeedEmptTop5cm_DetEff.root', 'asymFeedEmptBottom5cm_DetEff.root'), 'cellCenter': 0.05},
           {'parameter': 10,'cellCut':1, 'filenames': ('asymFeedOffset10cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'asymFeedEmptTop10cm_DetEff.root', 'asymFeedEmptBottom10cm_DetEff.root'),'cellCenter': 0.1},
           {'parameter': 15,'cellCut':1, 'filenames': ('asymFeedOffset15cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'asymFeedEmptTop15cm_DetEff.root', 'asymFeedEmptBottom15cm_DetEff.root'),'cellCenter': 0.15},
           {'parameter': 20,'cellCut':1, 'filenames': ('asymFeedOffset20cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'asymFeedEmptTop20cm_DetEff.root', 'asymFeedEmptBottom20cm_DetEff.root'),'cellCenter': 0.2},
           {'parameter': 25,'cellCut':1, 'filenames': ('asymFeedOffset25cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'asymFeedEmptTop25cm_DetEff.root', 'asymFeedEmptBottom25cm_DetEff.root'),'cellCenter': 0.25},
           {'parameter': 30,'cellCut':1, 'filenames': ('asymFeedOffset30cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'asymFeedEmptTop30cm_DetEff.root', 'asymFeedEmptBottom30cm_DetEff.root'),'cellCenter': 0.3},
           {'parameter': 35,'cellCut':1, 'filenames': ('asymFeedOffset35cm.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'asymFeedEmptTop35cm_DetEff.root', 'asymFeedEmptBottom35cm_DetEff.root'),'cellCenter': 0.35},
           ]
camSpider =[
           {'parameter': 1, 'cellCut':0, 'filenames': ('noSpider.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},  
           {'parameter': 2, 'cellCut':0, 'filenames': ('camSpider.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 3, 'cellCut':0, 'filenames': ('camSpider90deg.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 4, 'cellCut':0, 'filenames': ('camSpider90degFurther.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 5, 'cellCut':0, 'filenames': ('camSpider0cmCap.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 6, 'cellCut':0, 'filenames': ('camSpider5cmCap.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0},
           {'parameter': 7, 'cellCut':0, 'filenames': ('camSpider10cmCap.root', 'kinkHtTopTau2V_hist.root', 'kinkHtBottomTau2V_hist.root', 'kinkHtTopEmptying2V_DetEff.root', 'kinkHtBottomEmptying2V_DetEff.root'), 'cellCenter': 0}
           ]



ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kError

#Make a dictionary here called studyList: key = studayName: value = 
studyList1 = {'ExtractionHeight':extractionHeight, 'SourceTubeDiameter':sourceTube, 'HEXdiameter':HEXdiameter, 'KinkHeight':KinkHeight, 'NoBKinkHeight':noBKinkHeight, 'NoBnoFKinkHeight':noBnoFoilKinkHeight, 'DoubleBKinkHeight':doubleBKinkHeight, 'HardBore':hardBore}
studyList2 = {'BoreDiameter':boreDiameter, 'DoubleBBoreDiameter':doubleBboreDiameter, 'AngledTransferGuide':angTransferGuide, 'FunnelShapes':funnelShapes, 'KinkRadiusV1':kinkRadius1, 'KinkRadiusV2':kinkRadius2, 'KinkAngleV1':kinkAngle1, 'KinkAngleV2':kinkAngle2, 'MinimalSeperationFeeders':minSepFeeders1, 'AsymmetricFeeders':asymmetricFeeders1, 'CamSpider':camSpider} 

studyList = {'KinkHeight':KinkHeight}

for studyKey in studyList.keys():
    
    studyDays, studyDaysErrors, studyParameters, zeroArray  = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
    studies = studyList[studyKey]
    f = open(studyKey+'.txt', mode='w')
    f.write("FileName,TimeBeforeValveOpen,FillingTime,TotalSourcePumpingTime,CellStorageTime,EmptyingTime,CycleTime,DaysToReach,ErrorDaysToReach\n")
    
    for study in studies:
      histograms = readHistograms(study['filenames'][0], study['filenames'][1], study['filenames'][2], study['filenames'][3], study['filenames'][4], study['cellCenter'], study['cellCut'])
      optimized = scipy.optimize.differential_evolution(lambda x: daysToReach(histograms, x[0], x[1], x[2], x[3])[0], [[0., 100], [0., 100.], [0., 300.], [0., 250.]], polish=False)
      xmin = optimized.x
      nit = optimized.nit
      print(study['parameter'], xmin, daysToReach(histograms, xmin[0], xmin[1], xmin[2], xmin[3]), nit)
      tempDays = daysToReach(histograms, xmin[0], xmin[1], xmin[2], xmin[3])[0]
      studyDays.append(tempDays)
      tempErrorDays = daysToReach(histograms, xmin[0], xmin[1], xmin[2], xmin[3])[1]
      studyDaysErrors.append(tempErrorDays)
      zeroArray.append(0.)
      studyParameters.append(study['parameter'])
      f.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(study['filenames'][0], xmin[0], xmin[1], xmin[0]+xmin[1] , xmin[2] , xmin[3] , +xmin[0] + xmin[1] + xmin[2] + xmin[3] , tempDays ,tempErrorDays))

    f.close()
    c = ROOT.TCanvas('c','c')
    gr = TGraphErrors(len(studyParameters), studyParameters, studyDays, zeroArray, studyDaysErrors);
    gr.SetMarkerStyle(20);
    gr.SetTitle("Days to reach vs " + studyKey) ; 
    gr.GetXaxis().SetTitle("Parameter");
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
