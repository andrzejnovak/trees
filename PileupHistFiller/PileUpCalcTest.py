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
parser.add_option('-c', '--choose', metavar='NAME', type='string', dest='choose', help="choochootrain")
(options, args) = parser.parse_args()

OutPut = ROOT.TFile( "PU_npv" + options.name + ".root", "recreate" )
OutPut.cd()

files = []
for file in os.listdir(options.files ):
	if file.endswith(options.choose+".root"): 
		files.append(file)
		print(file)

Nevents = 0
PU_hist = TH1F("DataPileupHist", "DataPileupHist", 100, 0, 100)
for i in files:
			print "Reading from " + i
			File = TFile(options.files  + "/" +i)
			Tree = File.Get("B2GTTreeMaker/B2GTree")
			n = Tree.GetEntries()
			Nevents += n
			for j in range(0, n): # Here is where we loop over all events.
				if j % 50000 == 0:
	      				percentDone = float(j) / float(n) * 100.0
	       				print 'Processing {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(j, n, percentDone )
				Tree.GetEntry(j)
				PU_hist.Fill(Tree.pu_NtrueInt)
print Nevents
OutPut.cd()

OutPut.Write()
OutPut.Save()

leg = TLegend(0.6,0.6,0.89,0.89)
leg.SetLineColor(0)
leg.SetFillColor(0)


CMSLABL = TLatex()
CMSLABL.SetNDC()
CMSLABL.SetTextSize(0.035)

THILABL = TLatex()
THILABL.SetNDC()
THILABL.SetTextSize(0.04)


C = TCanvas("PU_vizcheck_"+options.name, "", 800,600)
C.cd()
leg.Draw()
CMSLABL.DrawLatex(0.135,0.85,"CMS Preliminary")
THILABL.DrawLatex(0.7,0.91,"#bf{12.9 fb^{-1}}, 13 TeV")


