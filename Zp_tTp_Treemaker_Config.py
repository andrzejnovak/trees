import os
import sys
import math
from array import array
import ROOT
from ROOT import *
#import Zprime_to_tTprime_Treemaker
from Zprime_to_tTprime_Treemaker import *
#from SimpleTreemaker import *

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-n', '--name', metavar='NAME', type='string', dest='name', help="name")
parser.add_option('-c', '--choose', metavar='NAME', type='string', dest='choose', help="choochootrain")
parser.add_option('-f', '--files', metavar='FILES', type='string', dest='F', help="Location of the ntuples to run over.")
parser.add_option('-p', '--PU', metavar='FILES', type='string', dest='PU', help="PileUpFile")
parser.add_option("-x", '--xs', metavar="FILES", type='string', dest='xs', help="cross section")
parser.add_option("-s", '--saveto', metavar='NAME', type='string', dest='saveto', help="target dir")
parser.add_option("-t", '--trigg', metavar='NAME', type='string', dest='trigg', help="trigger")
parser.add_option("-m", '--treemaker', metavar='NAME', type='string', dest='treemaker', help="treemaker")

(options, args) = parser.parse_args()

#exec("from Zprime_to_tTprime_Treemaker_N"+options.treemaker+" import *")

test = Zp_tTp_Treemaker(options.name, options.F, True, eval(options.xs), options.choose, options.PU, options.saveto)
test.Fill("B2GTTreeMaker/B2GTree")

#def __init__(self, name, TreeFolder, mc, tt, weight, choose, PU): 
