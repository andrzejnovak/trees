import os
import sys
import math
from array import array
import ROOT
from ROOT import *
#import Zprime_to_tTprime_Treemaker
#from Zprime_to_tTprime_Treemaker import *
from SimpleTreemaker import *

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-n', '--name', metavar='NAME', type='string', dest='name', help="name")
parser.add_option('-c', '--choose', metavar='NAME', type='string', dest='choose', help="choochootrain")
parser.add_option('-f', '--files', metavar='FILES', type='string', dest='F', help="Location of the ntuples to run over.")
parser.add_option('-p', '--PU', metavar='FILES', type='string', dest='PU', help="PileUpFile")
parser.add_option("-w", '--weight', metavar="FILES", type='string', dest='weight', help="event weight")
parser.add_option("-s", '--saveto', metavar='NAME', type='string', dest='saveto', help="target dir")
parser.add_option("-t", '--trigg', metavar='NAME', type='string', dest='trigg', help="trigger")

#parser.add_option('-m', '--isMC', metavar='FILES', action='store_true', default=False, dest='mc', help="set MC vars?")
#parser.add_option('-t', '--isTT', metavar='FILES', action='store_true', default=True, dest='tt', help="set MC vars?")

(options, args) = parser.parse_args()

test = Zp_tTp_Treemaker(options.name, options.F, False, False , float(options.weight), options.choose, options.PU, options.saveto, options.trigg)
test.Fill("B2GTTreeMaker/B2GTree")
#def __init__(self, name, TreeFolder, mc, tt, weight, choose, PU): 
