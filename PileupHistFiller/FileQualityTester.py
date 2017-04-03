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
parser.add_option('-f', '--files', metavar='FILES', type='string', dest='files', help="Location of the ntuples to run over.")
(options, args) = parser.parse_args()

files = []
for file in os.listdir(options.files ):
	if file.endswith(".root"): 
		files.append(file)
		print(file)

for i in files:
			print "Reading from " + i
			File = TFile(options.files  + "/" +i)
			Tree = File.Get("B2GTTreeMaker/B2GTree")
			n = Tree.GetEntries()
			#for j in range(0, n): # Here is where we loop over all events.
				#if j % 50000 == 0:
	      				#percentDone = float(j) / float(n) * 100.0
	       				#print 'Processing {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(j, n, percentDone )
				#Tree.GetEntry(j)


