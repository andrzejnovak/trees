#

import os
import sys
import math
from array import array
import ROOT
from ROOT import *
import Plotting_Header
from Plotting_Header import *


from optparse import OptionParser
parser = OptionParser()
parser.add_option('-n', '--name', metavar='NAME', type='string', dest='name', help="The name of the output file, minus the .root.")
parser.add_option('-f', '--files', metavar='FILES', type='string', dest='files', help="Location of the ntuples to run over.")
(options, args) = parser.parse_args()

OutPut = ROOT.TFile( "PU_weights_" + options.name + ".root", "recreate" )
OutPut.cd()

path1 = options.files+"0000/"
path2 = options.files+"0001/"

files = []
for file in os.listdir(path1):
	if file.endswith(".root"):
		files.append(path1+file)
try:
	for file in os.listdir(path2):
		if file.endswith(".root"):
			files.append(path2+file)
except: print "Less than 1000 files"


PU_hist = TH1F("DataPileupHist", "DataPileupHist", 100, 0, 100)
W_hist = TH1F("WeightHist", "WeightHist", 1, 0, 1)
for num, i in enumerate(files):
	print "Reading file:", num,"/", len(files)
	File = TFile(i)
	Tree = File.Get("B2GTTreeMaker/B2GTree")
	WeightHist = File.Get("EventCounter/totweight")
	W_hist.Add(WeightHist, 1.)
	n = Tree.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		if j % 50000 == 0:
  			percentDone = float(j) / float(n) * 100.0
     			print 'Processing {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(j, n, percentDone )
		Tree.GetEntry(j)
		PU_hist.Fill(Tree.pu_NtrueInt)
print W_hist.GetEntries()

PUF = TFile("MyDataPileupHistogram.root")
PUn = PUF.Get("pileup")
PUU = TFile("MyDataPileupHistogram_UP.root")
PUu = PUU.Get("pileup")
PUD = TFile("MyDataPileupHistogram_DOWN.root")
PUd = PUD.Get("pileup")


OutPut.cd()
pPUWN = PUn.Clone("nom_weight")
pPUWU = PUu.Clone("up_weight")
pPUWD = PUd.Clone("down_weight")
pPU_hist = PU_hist.Clone("DataPileupHist")
PUWN = PUn.Clone("nom_weight")
PUWU = PUu.Clone("up_weight")
PUWD = PUd.Clone("down_weight")

for i in [PUWN, PUWU, PUWD, PU_hist]:
	i.Scale(1/i.Integral())

nPU_hist = PU_hist.Clone("DataPileupHist")
nPUWN = PUWN.Clone("nom_weight")
nPUWU = PUWU.Clone("up_weight")
nPUWD = PUWD.Clone("down_weight")

for j in [PUWN, PUWU, PUWD]:
	j.GetXaxis().SetTitle("True Number of Primary Vertices")
	j.GetYaxis().SetTitle("Weight")
	j.SetStats(0)
	j.SetLineWidth(2)
	j.SetLineColor(kBlue)
	j.Divide(PU_hist)
PUWN.SetLineStyle(1)
PUWU.SetLineStyle(2)
PUWD.SetLineStyle(3)
PUWD.GetXaxis().SetRangeUser(0,50)
PUWD.SetTitle("")


OutPut.Write()
OutPut.Save()

leg = TLegend(0.6,0.6,0.89,0.89)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(PUWN, "Pileup weight (nominal)", "LP")
leg.AddEntry(PUWU, "Pileup weight (up)", "LP")
leg.AddEntry(PUWD, "Pileup weight (down)", "LP")

CMSLABL = TLatex()
CMSLABL.SetNDC()
CMSLABL.SetTextSize(0.035)

THILABL = TLatex()
THILABL.SetNDC()
THILABL.SetTextSize(0.04)

FindAndSetMax([PUWN,PUWD,PUWU])

C = TCanvas("PU_vizcheck_"+options.name, "", 800,600)
C.cd()
PUWD.Draw("hist")
PUWN.Draw("same hist")
PUWU.Draw("same hist")
leg.Draw()
CMSLABL.DrawLatex(0.135,0.85,"CMS Preliminary")
THILABL.DrawLatex(0.7,0.91,"#bf{12.9 fb^{-1}}, 13 TeV")
C.SaveAs("PU_diag_"+options.name+"final.png")

leg = TLegend(0.6,0.6,0.89,0.89)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(pPU_hist, "Pure npv", "LP")
leg.AddEntry(pPUWN, "Pure data (nominal)", "LP")
leg.AddEntry(pPUWU, "Pure data (up)", "LP")
leg.AddEntry(pPUWD, "Pure data (down)", "LP")
FindAndSetMax([pPUWN,pPUWD,pPUWU, pPU_hist])
D = TCanvas("PU_vizcheck_"+options.name, "", 800,600)
D.SetLineColor(kBlack)
D.cd()
for i, e in enumerate([pPUWD, pPUWN, pPUWU]):
	e.SetLineColor(kRed)
	e.SetLineStyle(i+1)
pPU_hist.SetLineColor(kBlack)
pPU_hist.Draw("hist")
pPUWD.Draw("same hist")
pPUWN.Draw("same hist")
pPUWU.Draw("same hist")
leg.Draw()
CMSLABL.DrawLatex(0.135,0.85,"CMS Preliminary")
THILABL.DrawLatex(0.7,0.91,"#bf{12.9 fb^{-1}}, 13 TeV")
D.SaveAs("PU_diag_"+options.name+"raw.png")

leg = TLegend(0.6,0.6,0.89,0.89)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(nPU_hist, "Normed npv", "LP")
leg.AddEntry(nPUWN, "Normed data (nominal)", "LP")
leg.AddEntry(nPUWU, "Normed data (up)", "LP")
leg.AddEntry(nPUWD, "Normed data (down)", "LP")
FindAndSetMax([nPUWN,nPUWD,nPUWU, nPU_hist])
E = TCanvas("PU_vizcheck_"+options.name, "", 800,600)
E.cd()
for i, e in enumerate([nPUWD, nPUWN, nPUWU]):
	e.SetLineColor(kRed)
	e.SetLineStyle(i+1)
nPU_hist.SetLineColor(kBlack)
nPU_hist.Draw("hist")
nPUWD.Draw("same hist")
nPUWN.Draw("same hist")
nPUWU.Draw("same hist")
leg.Draw()
CMSLABL.DrawLatex(0.135,0.85,"CMS Preliminary")
THILABL.DrawLatex(0.7,0.91,"#bf{12.9 fb^{-1}}, 13 TeV")
E.SaveAs("PU_diag_"+options.name+"normed.png")



