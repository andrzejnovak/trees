import os
import ROOT
from ROOT import *
from array import array
import math
from math import *
import sys

TF = "/eos/uscms/store/user/anovak/QCD/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_234704/0000/"

QCD = ["/eos/uscms/store/user/anovak/QCD/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_234704/0000/", "/eos/uscms/store/user/anovak/QCD/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_234825/0000/", "/eos/uscms/store/user/anovak/QCD/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_234946/0000/", "/eos/uscms/store/user/anovak/QCD/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/160806_235108/0000/","/eos/uscms/store/user/anovak/QCD/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v3/160806_235229/0000/", "/eos/uscms/store/user/anovak/QCD/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_235350/0000/"]

W = ["/eos/uscms/store/user/anovak/Wjets/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/160806_235511/0000/", "/eos/uscms/store/user/anovak/Wjets/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_235632/0000/", "/eos/uscms/store/user/anovak/Wjets/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_235755/0000/", "/eos/uscms/store/user/anovak/Wjets/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160806_235917/0000/" ,"/eos/uscms/store/user/anovak/Wjets/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/160807_000038/0000/", "/eos/uscms/store/user/anovak/Wjets/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/160807_000159/0000/" , "/eos/uscms/store/user/anovak/Wjets/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160807_000322/0000/",  "/eos/uscms/store/user/anovak/Wjets/WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160807_000443/0000/"]

sig = ["/eos/uscms/store/user/anovak/Signal/ZprimeToTprimeT_TprimeToWB_MZp-1500Nar_MTp-900Nar_LH_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/160807_111901/0000/", "/eos/uscms/store/user/anovak/Signal/ZprimeToTprimeT_TprimeToWB_MZp-2500Nar_MTp-1500Nar_LH_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/160807_113102/0000/" , "/eos/uscms/store/user/anovak/Signal/ZprimeToTprimeT_TprimeToWB_MZp-2000Nar_MTp-1200Nar_LH_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_80X_V2p0_PR53_Aug06_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/160807_112259/0000/"]


def Fill(TreeName): # Loop through events and fill them. Actual Fill step is done at the end, allowing us to make a few quality control cuts.
	total = 0
	for i in files:
		print i
		File = TFile("root://cmsxrootd.fnal.gov/"+TF  + "/" +i)
		Tree = File.Get(TreeName)
		n = Tree.GetEntries()
		total += n
	return total

def count(TF, end=".root"):
	files =[]
	for file in os.listdir(TF):
		if file.endswith(end): # select files, allows you to select a subset of files for running on very large directories
			files.append(file)

	TreeName = "B2GTTreeMaker/B2GTree"
	#TreeName = "tree_T1"
	total = 0
	for i in files:
		print i
		File = TFile("root://cmsxrootd.fnal.gov/"+TF  + "/" +i)
		Tree = File.Get(TreeName)
		n = Tree.GetEntries()
		total += n
	return total

direc = "/home/storage/andrzejnovak/T"

for i in sig:
	n = count(i)
	print n
	
