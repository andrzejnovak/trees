import os
import ROOT
from ROOT import *
from array import array
import math
from math import *
import sys


pufiles = [] 
for file in os.listdir("/uscms/home/anovak/FinalZPrime/CMSSW_8_0_20/src/Trees/PileupHistFiller/"):
	if file.startswith("PU_weights_W"):
		if file.endswith(".root"):
			pufiles.append(file)
print pufiles
for pu in pufiles:
	puf = ROOT.TFile(pu)
	totweight = puf.WeightHist.GetEntries()

	print pu, "Totweight:", totweight
	
	
	
