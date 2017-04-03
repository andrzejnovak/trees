#
import os
import ROOT
from ROOT import *
from array import array
import math
from math import *
import sys

def MakeW(met, lep): #both should be TLor vectors.
	newmet = ROOT.TLorentzVector()
	newmet_m = ROOT.TLorentzVector()
	newmet_p = ROOT.TLorentzVector()
	newmet.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
	newmet_m.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
	newmet_p.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
	phivec = [math.cos(met.Phi()), math.sin(met.Phi())]
	P_lep = math.sqrt((lep.Px()*lep.Px())+(lep.Py()*lep.Py())+(lep.Pz()*lep.Pz()))
	P_phi = (lep.Px()*phivec[0])+(lep.Py()*phivec[1])
	b = (80.4*80.4) + (P_lep*P_lep) - (lep.E()*lep.E()) + (2*met.Pt()*P_phi)
	arg = (lep.E()*lep.E()) * ((4*met.Pt()*met.Pt()*((lep.Pz()*lep.Pz())-(lep.E()*lep.E())))+(b*b))
	if arg <= 0:
		Pz_met = lep.Pz()*b/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
		newmet.SetPz(Pz_met)
		newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))
		return newmet+lep
	else:
		Pz_met_p = ((lep.Pz()*b)+math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
		Pz_met_m = ((lep.Pz()*b)-math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
		newmet_p.SetPz(Pz_met_p)
		newmet_p.SetE(math.sqrt(newmet_p.Px()*newmet_p.Px()+newmet_p.Py()*newmet_p.Py()+newmet_p.Pz()*newmet_p.Pz()))
		newmet_m.SetPz(Pz_met_m)
		newmet_m.SetE(math.sqrt(newmet_m.Px()*newmet_m.Px()+newmet_m.Py()*newmet_m.Py()+newmet_m.Pz()*newmet_m.Pz()))
		return newmet_p+lep

class Zp_tTp_Treemaker:
	def __init__(self, name, TreeFolder, mc, tt, weight, choose, PU, saveto): # Initialize: mc is an important option and data will crash if run with it. Weight should be xs/NevtGen.
		self.mc = mc
		self.tt = tt
		self.name = name
		self.w = weight
		self.saveto = saveto
		self.files = []
		self.TF = TreeFolder
		self.puf = ROOT.TFile(PU)
		self.puf_n = self.puf.Get("nom_weight")
		self.puf_u = self.puf.Get("up_weight")
		self.puf_d = self.puf.Get("down_weight")
		self.sadnumber = 0
		if self.mc or self.tt:
			# Set up b-tag Calibrations:
			self.BTagCalib = BTagCalibration("csvv2", "CSVv2_ichep.csv")
			# Central
			self.BTagCalib_mediumB = BTagCalibrationReader(1, "central")
			self.BTagCalib_mediumB.load(self.BTagCalib, 0, "comb")
			self.BTagCalib_mediumC = BTagCalibrationReader(1, "central")
			self.BTagCalib_mediumC.load(self.BTagCalib, 1, "comb")
			self.BTagCalib_mediumL = BTagCalibrationReader(1, "central")
			self.BTagCalib_mediumL.load(self.BTagCalib, 2, "incl")
			# Up
			self.BTagCalib_mediumBu = BTagCalibrationReader(1, "up")
			self.BTagCalib_mediumBu.load(self.BTagCalib, 0, "comb")
			self.BTagCalib_mediumCu = BTagCalibrationReader(1, "up")
			self.BTagCalib_mediumCu.load(self.BTagCalib, 1, "comb")
			self.BTagCalib_mediumLu = BTagCalibrationReader(1, "up")
			self.BTagCalib_mediumLu.load(self.BTagCalib, 2, "incl")
			# Down
			self.BTagCalib_mediumBd = BTagCalibrationReader(1, "down")
			self.BTagCalib_mediumBd.load(self.BTagCalib, 0, "comb")
			self.BTagCalib_mediumCd = BTagCalibrationReader(1, "down")
			self.BTagCalib_mediumCd.load(self.BTagCalib, 1, "comb")
			self.BTagCalib_mediumLd = BTagCalibrationReader(1, "down")
			self.BTagCalib_mediumLd.load(self.BTagCalib, 2, "incl")
		if self.mc:
			self.MuTrigFile = TFile("SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root")
			self.MuTrigHist1 = self.MuTrigFile.Get("Mu45_eta2p1_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt_DATA")
			self.MuTrigHist2 = self.MuTrigFile.Get("Mu45_eta2p1_PtEtaBins_Run274094_to_276097/efficienciesDATA/abseta_pt_DATA")
			self.MuIDFile = TFile("MuonID_Z_RunBCD_prompt80X_7p65.root")
			self.MuIDHist = self.MuIDFile.Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")
			self.ElIDFile = TFile("egammaEffi.txt_SF2D.root")
			self.ElIDHist = self.ElIDFile.Get("EGamma_SF2D")
		#print "Files in " + self.TF + " :"
		for file in os.listdir(self.TF):
			if file.endswith(choose): # select files, allows you to select a subset of files for running on very large directories
				self.files.append(file)
				##print(file)
		#print self.files
		self.__book__()
	def addBranch(self, name, var, T): # Simple script to keep things straight when creating the branches. Only advantage is that the names are aligned (reduces chance of bugs from copy pasting dozens of times for all the up/dn variations of variables).
		T.Branch(name, var, name+'/F')
	def __book__(self): # Create the trees/branches which we will fill

		self.f = ROOT.TFile(self.saveto + self.name + ".root", "recreate" )
        	self.f.cd()
		# Two main trees: one for each channel"\
		self.CutHist = TH1F("CutHist", "", 10, 0, 10)
		self.tree1 = ROOT.TTree("tree_T1", "tree_T1")
		self.tree2 = ROOT.TTree("tree_T2", "tree_T2")

		self.tree1jecU = ROOT.TTree("tree_T1jecU", "tree_T1jesU")
		self.tree1jecD = ROOT.TTree("tree_T1jecU", "tree_T1jesU")
		self.tree2jecU = ROOT.TTree("tree_T2jecU", "tree_T2jesU")
		self.tree2jecD = ROOT.TTree("tree_T2jecU", "tree_T2jesU")

		self.tree1jesU = ROOT.TTree("tree_T1jesU", "tree_T1jesU")
		self.tree1jesD = ROOT.TTree("tree_T1jesU", "tree_T1jesU")
		self.tree2jesU = ROOT.TTree("tree_T2jesU", "tree_T2jesU")
		self.tree2jesD = ROOT.TTree("tree_T2jesU", "tree_T2jesU")

		# Vars shared by both trees:		
		self.weight = array('f', [0.0])
		self.addBranch('weight', self.weight, self.tree1)
		self.addBranch('weight', self.weight, self.tree2)
		#
		self.LepIso = array('f', [-1.0]) 
		self.addBranch('LepMiniIso', self.LepIso, self.tree1)
		self.addBranch('LepMiniIso', self.LepIso, self.tree2)
		#
		self.LepTight = array('f', [-1.0]) # This will be 0, 1, 2 or 3 depending on the highest ID passed (Tight = 3, Medium = 2, etc)
		self.addBranch('LepTightness', self.LepTight, self.tree1)
		self.addBranch('LepTightness', self.LepTight, self.tree2)
		#
		self.Lep2D_dr = array('f', [0.0])
		self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree1)
		self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree2)
		#
		self.Lep2D_rel = array('f', [0.0])
		self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree1)
		self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree2)
		#
		self.LepType = array('f', [0.0]) # -1 = Electron, 1 = Muon
		self.addBranch('LepType', self.LepType, self.tree1)
		self.addBranch('LepType', self.LepType, self.tree2)
		#
		self.LepPt = array('f', [0.0]) # Lepton pT: highest el or mu
		self.addBranch('LepPt', self.LepPt, self.tree1)
		self.addBranch('LepPt', self.LepPt, self.tree2)
		#
		self.METPt = array('f', [0.0])
		self.addBranch('METPt', self.METPt, self.tree1)
		self.addBranch('METPt', self.METPt, self.tree2)
		#
		self.WPt = array('f', [0.0]) # Recreated from MET and LEP kinematics
		self.addBranch('WPt', self.WPt, self.tree1)
		self.addBranch('WPt', self.WPt, self.tree2)
		#
		self.TAGPt = array('f', [0.0])
		self.addBranch('TAGPt', self.TAGPt, self.tree1)
		self.addBranch('TAGPt', self.TAGPt, self.tree2)
		#
		self.TAGM = array('f', [0.0])
		self.addBranch('TAGM', self.TAGM, self.tree1)
		self.addBranch('TAGM', self.TAGM, self.tree2)
		#
		self.TAGTau32 = array('f', [-1.0])
		self.addBranch('TAGTau32', self.TAGTau32, self.tree1)
		self.addBranch('TAGTau32', self.TAGTau32, self.tree2)
		#
		self.TAGTau21 = array('f', [-1.0])
		self.addBranch('TAGTau21', self.TAGTau21, self.tree1)
		self.addBranch('TAGTau21', self.TAGTau21, self.tree2)
		#
		self.JLep1Pt = array('f', [0.0])
		self.JLep2Pt = array('f', [0.0])
		self.addBranch('lepJetPt', self.JLep1Pt, self.tree1)
		self.addBranch('lepJetPt', self.JLep2Pt, self.tree2)
		#
		self.JLep1CSV = array('f', [-1.0])
		self.JLep2CSV = array('f', [-1.0])
		self.addBranch('lepJetCSV', self.JLep1CSV, self.tree1)
		self.addBranch('lepJetCSV', self.JLep2CSV, self.tree2)
		#
		self.JHadPt = array('f', [-1.0])
		self.addBranch('hadJetPt', self.JHadPt, self.tree2)
		#
		self.JHadCSV = array('f', [-1.0])
		self.addBranch('hadJetCSV', self.JHadCSV, self.tree2)
		#
		self.Tp1M = array('f', [0.0])
		self.addBranch('TPRIMEM', self.Tp1M, self.tree1)
		#
		self.Zp1M = array('f', [0.0])
		self.addBranch('ZPRIMEM', self.Zp1M, self.tree1)
		#
		self.Tp2M = array('f', [0.0])
		self.addBranch('TPRIMEM', self.Tp2M, self.tree2)
		#
		self.TopM = array('f', [0.0])
		self.addBranch('lepTopM', self.TopM, self.tree2)
		#
		self.Zp2M = array('f', [0.0])
		self.addBranch('ZPRIMEM', self.Zp2M, self.tree2)
		if self.mc or self.tt:
			# Inert SFs:
			self.topSF = array('f', [0.0])
			self.addBranch('topSF', self.topSF, self.tree1)
			self.topSFu = array('f', [0.0])
			self.addBranch('topSF_up', self.topSFu, self.tree1)
			self.topSFd = array('f', [0.0])
			self.addBranch('topSF_down', self.topSFd, self.tree1)

			self.bSF = array('f', [0.0])
			self.addBranch('bSF', self.bSF, self.tree1)
			self.bSFu = array('f', [0.0])
			self.addBranch('bSF_up', self.bSFu, self.tree1)
			self.bSFd = array('f', [0.0])
			self.addBranch('bSF_down', self.bSFd, self.tree1)


			self.blSF = array('f', [0.0])
			self.addBranch('blSF', self.blSF, self.tree2)
			self.blSFu = array('f', [0.0])
			self.addBranch('blSF_up', self.blSFu, self.tree2)
			self.blSFd = array('f', [0.0])
			self.addBranch('blSF_down', self.blSFd, self.tree2)


			self.bhSF = array('f', [0.0])
			self.addBranch('bhSF', self.bhSF, self.tree2)
			self.bhSFu = array('f', [0.0])
			self.addBranch('bhSF_up', self.bhSFu, self.tree2)
			self.bhSFd = array('f', [0.0])
			self.addBranch('bhSF_down', self.bhSFd, self.tree2)
			if self.mc:	

				self.ttHT = array('f', [0.0])
				self.addBranch('ttHT', self.ttHT, self.tree1)
				self.addBranch('ttHT', self.ttHT, self.tree2)
				#
				self.LepTrigSF = array('f', [-1.0]) # Lepton pT: highest el or mu
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree1)
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree2)	
				#
				self.LepIDSF = array('f', [-1.0]) # Lepton pT: highest el or mu
				self.addBranch('LepIDSF', self.LepIDSF, self.tree1)
				self.addBranch('LepIDSF', self.LepIDSF, self.tree2)	


				#
				self.LepTrigSFU = array('f', [-1.0]) # Lepton pT: highest el or mu
				self.addBranch('LepTrigSFU', self.LepTrigSFU, self.tree1)
				self.addBranch('LepTrigSFU', self.LepTrigSFU, self.tree2)	
				#
				self.LepIDSFU = array('f', [-1.0]) # Lepton pT: highest el or mu
				self.addBranch('LepIDSFU', self.LepIDSFU, self.tree1)
				self.addBranch('LepIDSFU', self.LepIDSFU, self.tree2)	


				#
				self.LepTrigSFD = array('f', [-1.0]) # Lepton pT: highest el or mu
				self.addBranch('LepTrigSFD', self.LepTrigSFD, self.tree1)
				self.addBranch('LepTrigSFD', self.LepTrigSFD, self.tree2)	
				#
				self.LepIDSFD = array('f', [-1.0]) # Lepton pT: highest el or mu
				self.addBranch('LepIDSFD', self.LepIDSFD, self.tree1)
				self.addBranch('LepIDSFD', self.LepIDSFD, self.tree2)	
	
				self.PUweightN = array('f', [0.0])
				self.addBranch('PU_weight_N', self.PUweightN, self.tree1)
				self.addBranch('PU_weight_N', self.PUweightN, self.tree2)	
				self.PUweightU = array('f', [0.0])
				self.addBranch('PU_weight_U', self.PUweightU, self.tree1)
				self.addBranch('PU_weight_U', self.PUweightU, self.tree2)		
				self.PUweightD = array('f', [0.0])
				self.addBranch('PU_weight_D', self.PUweightD, self.tree1)
				self.addBranch('PU_weight_D', self.PUweightD, self.tree2)

				self.PDFweightU = array('f', [0.0])
				self.addBranch('PDF_weight_U', self.PDFweightU, self.tree1)
				self.addBranch('PDF_weight_U', self.PDFweightU, self.tree2)	
				self.PDFweightD = array('f', [0.0])
				self.addBranch('PDF_weight_D', self.PDFweightD, self.tree1)
				self.addBranch('PDF_weight_D', self.PDFweightD, self.tree2)	

	
				self.SCLweightU = array('f', [0.0])
				self.addBranch('SCL_weight_U', self.SCLweightU, self.tree1)
				self.addBranch('SCL_weight_U', self.SCLweightU, self.tree2)	
				self.SCLweightD = array('f', [0.0])
				self.addBranch('SCL_weight_D', self.SCLweightD, self.tree1)
				self.addBranch('SCL_weight_D', self.SCLweightD, self.tree2)	
				
				self.isLEP = array('f', [-1.0])
				self.addBranch('isLEP', self.isLEP, self.tree1)
				self.addBranch('isLEP', self.isLEP, self.tree2)
				self.addBranch('isLEP', self.isLEP, self.tree1jecU)
				self.addBranch('isLEP', self.isLEP, self.tree2jecU)
				self.addBranch('isLEP', self.isLEP, self.tree1jecD)
				self.addBranch('isLEP', self.isLEP, self.tree2jecD)
				self.addBranch('isLEP', self.isLEP, self.tree1jesU)
				self.addBranch('isLEP', self.isLEP, self.tree2jesU)
				self.addBranch('isLEP', self.isLEP, self.tree1jesD)
				self.addBranch('isLEP', self.isLEP, self.tree2jesD)


				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree1jecU)
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree2jecU)	
				self.addBranch('LepIDSF', self.LepIDSF, self.tree1jecU)
				self.addBranch('LepIDSF', self.LepIDSF, self.tree2jecU)	
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree1jecD)
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree2jecD)	
				self.addBranch('LepIDSF', self.LepIDSF, self.tree1jecD)
				self.addBranch('LepIDSF', self.LepIDSF, self.tree2jecD)	
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree1jesU)
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree2jesU)	
				self.addBranch('LepIDSF', self.LepIDSF, self.tree1jesU)
				self.addBranch('LepIDSF', self.LepIDSF, self.tree2jesU)
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree1jesD)
				self.addBranch('LepTrigSF', self.LepTrigSF, self.tree2jesD)	
				self.addBranch('LepIDSF', self.LepIDSF, self.tree1jesD)
				self.addBranch('LepIDSF', self.LepIDSF, self.tree2jesD)		
				# JEC UP
				self.addBranch('PU_weight_N', self.PUweightN, self.tree1jecU)
				self.addBranch('PU_weight_N', self.PUweightN, self.tree2jecU)
				self.addBranch('topSF', self.topSF, self.tree1jecU)
				self.addBranch('bSF', self.bSF, self.tree1jecU)
				self.addBranch('blSF', self.blSF, self.tree2jecU)
				self.addBranch('bhSF', self.bhSF, self.tree2jecU)
				self.addBranch('weight', self.weight, self.tree1jecU)
				self.addBranch('weight', self.weight, self.tree2jecU)
				self.addBranch('LepMiniIso', self.LepIso, self.tree1jecU)
				self.addBranch('LepMiniIso', self.LepIso, self.tree2jecU)
				self.addBranch('LepTightness', self.LepTight, self.tree1jecU)
				self.addBranch('LepTightness', self.LepTight, self.tree2jecU)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree1jecU)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree2jecU)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree1jecU)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree2jecU)
				self.addBranch('LepType', self.LepType, self.tree1jecU)
				self.addBranch('LepType', self.LepType, self.tree2jecU)
				self.addBranch('LepPt', self.LepPt, self.tree1jecU)
				self.addBranch('LepPt', self.LepPt, self.tree2jecU)
				self.addBranch('METPt', self.METPt, self.tree1jecU)
				self.addBranch('METPt', self.METPt, self.tree2jecU)
				self.addBranch('WPt', self.WPt, self.tree1jecU)
				self.addBranch('WPt', self.WPt, self.tree2jecU)
				self.addBranch('lepJetCSV', self.JLep1CSV, self.tree1jecU)
				self.addBranch('lepJetCSV', self.JLep2CSV, self.tree2jecU)
				self.addBranch('hadJetCSV', self.JHadCSV, self.tree2jecU)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree1jecU)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree2jecU)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree1jecU)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree2jecU)
				self.addBranch('TAGM', self.TAGM, self.tree1jecU)
				self.addBranch('TAGM', self.TAGM, self.tree2jecU)
				# NEW VARS:
				self.TAGPtjecU = array('f', [0.0])
				self.addBranch('TAGPt', self.TAGPtjecU, self.tree1jecU)
				self.addBranch('TAGPt', self.TAGPtjecU, self.tree2jecU)
				self.JLep1PtjecU = array('f', [0.0])
				self.JLep2PtjecU = array('f', [0.0])
				self.addBranch('lepJetPt', self.JLep1PtjecU, self.tree1jecU)
				self.addBranch('lepJetPt', self.JLep2PtjecU, self.tree2jecU)
				self.JHadPtjecU = array('f', [-1.0])
				self.addBranch('hadJetPt', self.JHadPtjecU, self.tree2jecU)
				self.Tp1MjecU = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp1MjecU, self.tree1jecU)
				self.Zp1MjecU = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp1MjecU, self.tree1jecU)
				self.Tp2MjecU = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp2MjecU, self.tree2jecU)
				self.TopMjecU = array('f', [0.0])
				self.addBranch('lepTopM', self.TopMjecU, self.tree2jecU)
				self.Zp2MjecU = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp2MjecU, self.tree2jecU)

				# JEC DN
				self.addBranch('PU_weight_N', self.PUweightN, self.tree1jecD)
				self.addBranch('PU_weight_N', self.PUweightN, self.tree2jecD)
				self.addBranch('topSF', self.topSF, self.tree1jecD)
				self.addBranch('bSF', self.bSF, self.tree1jecD)
				self.addBranch('blSF', self.blSF, self.tree2jecD)
				self.addBranch('bhSF', self.bhSF, self.tree2jecD)
				self.addBranch('weight', self.weight, self.tree1jecD)
				self.addBranch('weight', self.weight, self.tree2jecD)
				self.addBranch('LepMiniIso', self.LepIso, self.tree1jecD)
				self.addBranch('LepMiniIso', self.LepIso, self.tree2jecD)
				self.addBranch('LepTightness', self.LepTight, self.tree1jecD)
				self.addBranch('LepTightness', self.LepTight, self.tree2jecD)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree1jecD)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree2jecD)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree1jecD)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree2jecD)
				self.addBranch('LepType', self.LepType, self.tree1jecD)
				self.addBranch('LepType', self.LepType, self.tree2jecD)
				self.addBranch('LepPt', self.LepPt, self.tree1jecD)
				self.addBranch('LepPt', self.LepPt, self.tree2jecD)
				self.addBranch('METPt', self.METPt, self.tree1jecD)
				self.addBranch('METPt', self.METPt, self.tree2jecD)
				self.addBranch('WPt', self.WPt, self.tree1jecD)
				self.addBranch('WPt', self.WPt, self.tree2jecD)
				self.addBranch('lepJetCSV', self.JLep1CSV, self.tree1jecD)
				self.addBranch('lepJetCSV', self.JLep2CSV, self.tree2jecD)
				self.addBranch('hadJetCSV', self.JHadCSV, self.tree2jecD)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree1jecD)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree2jecD)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree1jecD)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree2jecD)
				self.addBranch('TAGM', self.TAGM, self.tree1jecD)
				self.addBranch('TAGM', self.TAGM, self.tree2jecD)
				# NEW VARS:
				self.TAGPtjecD = array('f', [0.0])
				self.addBranch('TAGPt', self.TAGPtjecD, self.tree1jecD)
				self.addBranch('TAGPt', self.TAGPtjecD, self.tree2jecD)
				self.JLep1PtjecD = array('f', [0.0])
				self.JLep2PtjecD = array('f', [0.0])
				self.addBranch('lepJetPt', self.JLep1PtjecD, self.tree1jecD)
				self.addBranch('lepJetPt', self.JLep2PtjecD, self.tree2jecD)
				self.JHadPtjecD = array('f', [-1.0])
				self.addBranch('hadJetPt', self.JHadPtjecD, self.tree2jecD)
				self.Tp1MjecD = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp1MjecD, self.tree1jecD)
				self.Zp1MjecD = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp1MjecD, self.tree1jecD)
				self.Tp2MjecD = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp2MjecD, self.tree2jecD)
				self.TopMjecD = array('f', [0.0])
				self.addBranch('lepTopM', self.TopMjecD, self.tree2jecD)
				self.Zp2MjecD = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp2MjecD, self.tree2jecD)

				# JES UP
				self.addBranch('PU_weight_N', self.PUweightN, self.tree1jesU)
				self.addBranch('PU_weight_N', self.PUweightN, self.tree2jesU)
				self.addBranch('topSF', self.topSF, self.tree1jesU)
				self.addBranch('bSF', self.bSF, self.tree1jesU)
				self.addBranch('blSF', self.blSF, self.tree2jesU)
				self.addBranch('bhSF', self.bhSF, self.tree2jesU)
				self.addBranch('weight', self.weight, self.tree1jesU)
				self.addBranch('weight', self.weight, self.tree2jesU)
				self.addBranch('LepMiniIso', self.LepIso, self.tree1jesU)
				self.addBranch('LepMiniIso', self.LepIso, self.tree2jesU)
				self.addBranch('LepTightness', self.LepTight, self.tree1jesU)
				self.addBranch('LepTightness', self.LepTight, self.tree2jesU)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree1jesU)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree2jesU)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree1jesU)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree2jesU)
				self.addBranch('LepType', self.LepType, self.tree1jesU)
				self.addBranch('LepType', self.LepType, self.tree2jesU)
				self.addBranch('LepPt', self.LepPt, self.tree1jesU)
				self.addBranch('LepPt', self.LepPt, self.tree2jesU)
				self.addBranch('METPt', self.METPt, self.tree1jesU)
				self.addBranch('METPt', self.METPt, self.tree2jesU)
				self.addBranch('WPt', self.WPt, self.tree1jesU)
				self.addBranch('WPt', self.WPt, self.tree2jesU)
				self.addBranch('lepJetCSV', self.JLep1CSV, self.tree1jesU)
				self.addBranch('lepJetCSV', self.JLep2CSV, self.tree2jesU)
				self.addBranch('hadJetCSV', self.JHadCSV, self.tree2jesU)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree1jesU)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree2jesU)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree1jesU)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree2jesU)
				self.addBranch('TAGM', self.TAGM, self.tree1jesU)
				self.addBranch('TAGM', self.TAGM, self.tree2jesU)
				# NEW VARS:
				self.TAGPtjesU = array('f', [0.0])
				self.addBranch('TAGPt', self.TAGPtjesU, self.tree1jesU)
				self.addBranch('TAGPt', self.TAGPtjesU, self.tree2jesU)
				self.JLep1PtjesU = array('f', [0.0])
				self.JLep2PtjesU = array('f', [0.0])
				self.addBranch('lepJetPt', self.JLep1PtjesU, self.tree1jesU)
				self.addBranch('lepJetPt', self.JLep2PtjesU, self.tree2jesU)
				self.JHadPtjesU = array('f', [-1.0])
				self.addBranch('hadJetPt', self.JHadPtjesU, self.tree2jesU)
				self.Tp1MjesU = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp1MjesU, self.tree1jesU)
				self.Zp1MjesU = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp1MjesU, self.tree1jesU)
				self.Tp2MjesU = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp2MjesU, self.tree2jesU)
				self.TopMjesU = array('f', [0.0])
				self.addBranch('lepTopM', self.TopMjesU, self.tree2jesU)
				self.Zp2MjesU = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp2MjesU, self.tree2jesU)

				# JES DN
				self.addBranch('PU_weight_N', self.PUweightN, self.tree1jesD)
				self.addBranch('PU_weight_N', self.PUweightN, self.tree2jesD)
				self.addBranch('topSF', self.topSF, self.tree1jesD)
				self.addBranch('bSF', self.bSF, self.tree1jesD)
				self.addBranch('blSF', self.blSF, self.tree2jesD)
				self.addBranch('bhSF', self.bhSF, self.tree2jesD)
				self.addBranch('weight', self.weight, self.tree1jesD)
				self.addBranch('weight', self.weight, self.tree2jesD)
				self.addBranch('LepMiniIso', self.LepIso, self.tree1jesD)
				self.addBranch('LepMiniIso', self.LepIso, self.tree2jesD)
				self.addBranch('LepTightness', self.LepTight, self.tree1jesD)
				self.addBranch('LepTightness', self.LepTight, self.tree2jesD)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree1jesD)
				self.addBranch('Lep2D_dR', self.Lep2D_dr, self.tree2jesD)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree1jesD)
				self.addBranch('Lep2D_rel', self.Lep2D_rel, self.tree2jesD)
				self.addBranch('LepType', self.LepType, self.tree1jesD)
				self.addBranch('LepType', self.LepType, self.tree2jesD)
				self.addBranch('LepPt', self.LepPt, self.tree1jesD)
				self.addBranch('LepPt', self.LepPt, self.tree2jesD)
				self.addBranch('METPt', self.METPt, self.tree1jesD)
				self.addBranch('METPt', self.METPt, self.tree2jesD)
				self.addBranch('WPt', self.WPt, self.tree1jesD)
				self.addBranch('WPt', self.WPt, self.tree2jesD)
				self.addBranch('lepJetCSV', self.JLep1CSV, self.tree1jesD)
				self.addBranch('lepJetCSV', self.JLep2CSV, self.tree2jesD)
				self.addBranch('hadJetCSV', self.JHadCSV, self.tree2jesD)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree1jesD)
				self.addBranch('TAGTau32', self.TAGTau32, self.tree2jesD)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree1jesD)
				self.addBranch('TAGTau21', self.TAGTau21, self.tree2jesD)
				self.addBranch('TAGM', self.TAGM, self.tree1jesD)
				self.addBranch('TAGM', self.TAGM, self.tree2jesD)
				# NEW VARS:
				self.TAGPtjesD = array('f', [0.0])
				self.addBranch('TAGPt', self.TAGPtjesD, self.tree1jesD)
				self.addBranch('TAGPt', self.TAGPtjesD, self.tree2jesD)
				self.JLep1PtjesD = array('f', [0.0])
				self.JLep2PtjesD = array('f', [0.0])
				self.addBranch('lepJetPt', self.JLep1PtjesD, self.tree1jesD)
				self.addBranch('lepJetPt', self.JLep2PtjesD, self.tree2jesD)
				self.JHadPtjesD = array('f', [-1.0])
				self.addBranch('hadJetPt', self.JHadPtjesD, self.tree2jesD)
				self.Tp1MjesD = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp1MjesD, self.tree1jesD)
				self.Zp1MjesD = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp1MjesD, self.tree1jesD)
				self.Tp2MjesD = array('f', [0.0])
				self.addBranch('TPRIMEM', self.Tp2MjesD, self.tree2jesD)
				self.TopMjesD = array('f', [0.0])
				self.addBranch('lepTopM', self.TopMjesD, self.tree2jesD)
				self.Zp2MjesD = array('f', [0.0])
				self.addBranch('ZPRIMEM', self.Zp2MjesD, self.tree2jesD)
			if self.tt:
				self.ttHT = array('f', [0.0])
				self.addBranch('ttHT', self.ttHT, self.tree1)
				self.addBranch('ttHT', self.ttHT, self.tree2)
	def Fill(self, TreeName): # Loop through events and fill them. Actual Fill step is done at the end, allowing us to make a few quality control cuts.
		total = 0
		print "filling..."
		for i in self.files:
			print "Reading from " + i
			File = TFile("root://cmsxrootd.fnal.gov/"+self.TF  + "/" +i)
			self.Tree = File.Get(TreeName)
			n = self.Tree.GetEntries()
			total += n-1
			for j in range(0, n): # Here is where we loop over all events.
				if j % 50000 == 0 or j == 1:
	      				percentDone = float(j) / float(n) * 100.0
	       				print 'Processing {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(j, n, percentDone )
				self.Tree.GetEntry(j)
				#print "EVENT STARTED"
				# Start finding variables and doing analysis things:
				self.weight[0] = self.w
				if self.tt:
					print "DO TTBAR as MC"
					
				if self.mc:
					#print "DOING MC"
					if self.doLepPart(): # Fills in basic leptonic variables (same for MC and data)
						# Do AK8 Jet: needs on jet with pT > 150 and mass > 50.
						if self.doAK8PartMC(): # MC version of the AK8 jet check. Uses Smeared quantities
							# Do MET/W: Needs 25 GeV of MET and a "good" W solution.
							if self.doWPart():
								# Do AK4 Jet:
								if self.doAK4JetsMC():
									# Split into two branches now (if two good jets are found). Drop events without a TPRIME candidate above the t-mass scale:
									if self.Type1:
										if self.doType1Vars(): # We can now start to fill variables:
											self.FillVarLeptons()
											self.FillVarAK8()
											self.FillVar1AK4()
											self.FillVar1Phys()
											self.FillVar1JEC()
											self.FillVar1JES()
											self.getBSF1()
											self.getMCVar()
											self.doMCUnc()
											self.doLepUnc()
											self.tree1.Fill()
											self.tree1jecU.Fill()
											self.tree1jecD.Fill()
											self.tree1jesU.Fill()
											self.tree1jesD.Fill()
									if len(self.jetList) > 1:
										if self.doType2Vars(): # We can now start to fill variables:
											self.FillVarLeptons()
											self.FillVarAK8()
											self.FillVar2AK4()
											self.FillVar2Phys()
											self.FillVar2JEC()
											self.FillVar2JES()
											self.getBSF2()
											self.getMCVar()
											self.doMCUnc()
											self.doLepUnc()
											self.tree2.Fill()
											self.tree2jecU.Fill()
											self.tree2jecD.Fill()
											self.tree2jesU.Fill()
											self.tree2jesD.Fill()
				else:
				    #print "DOING DATA"
				    if self.doTrigger():
					if self.doLepPart(): # Fills in basic leptonic variables (same for MC and data)
						# Do AK8 Jet: needs on jet with pT > 150 and mass > 50.
						if self.doAK8Part(): # unsmeared version of AK8 part
							# Do MET/W: Needs 25 GeV of MET and a "good" W solution.
							if self.doWPart():
								# Do AK4 Jet:
								if self.doAK4Jets():
									# Split into two branches now (if two good jets are found). Drop events without a TPRIME candidate above the t-mass scale:
									if self.Type1:
										if self.doType1Vars(): # We can now start to fill variables:
											self.FillVarLeptons()
											self.FillVarAK8()
											self.FillVar1AK4()
											self.FillVar1Phys()
											self.tree1.Fill()
									if len(self.jetList) > 1:
										if self.doType2Vars(): # We can now start to fill variables:
											self.FillVarLeptons()
											self.FillVarAK8()
											self.FillVar2AK4()
											self.FillVar2Phys()
											self.tree2.Fill()
			File.Close()
		self.f.cd()
	        self.f.Write()
	        self.f.Close()
		print str(total) + " events processed"

	def doLepUnc(self):
		if self.leptype == "mu":
			muIDbin = self.MuIDHist.FindBin(min(2.4,math.fabs(self.LEP.Eta())),min(200.,self.LEP.Pt()))
			muTRIGbin = self.MuTrigHist2.FindBin(min(2.,math.fabs(self.LEP.Eta())),min(500.,self.LEP.Pt()))
			self.LepTrigSF[0] = (self.MuTrigHist1.GetBinContent(muTRIGbin)*(0.05)) + (self.MuTrigHist2.GetBinContent(muTRIGbin)*(0.95))
			self.LepTrigSFU[0] = self.LepTrigSF[0] + 0.005 + (self.MuTrigHist1.GetBinError(muTRIGbin)*(0.05)) + (self.MuTrigHist2.GetBinError(muTRIGbin)*(0.95))# -0.5% recc by POG
			self.LepTrigSFD[0] = self.LepTrigSF[0] - 0.005 - (self.MuTrigHist1.GetBinError(muTRIGbin)*(0.05)) - (self.MuTrigHist2.GetBinError(muTRIGbin)*(0.95))
			self.LepIDSF[0] = self.MuIDHist.GetBinContent(muIDbin)
			self.LepIDSFU[0] = self.MuIDHist.GetBinContent(muIDbin)+0.01+ self.MuIDHist.GetBinError(muIDbin)
			self.LepIDSFD[0] = self.MuIDHist.GetBinContent(muIDbin)-0.01- self.MuIDHist.GetBinError(muIDbin) # -1% recc by POG
		if self.leptype == "el":
			elIDbin = self.ElIDHist.FindBin(min(2.4,math.fabs(self.LEP.Eta())),min(200.,self.LEP.Pt()))
			self.LepTrigSF[0] = 0.9598
			self.LepTrigSFU[0] = 0.9598 +0.0067
			self.LepTrigSFD[0] =0.9598 - 0.0067
			self.LepIDSF[0] = self.ElIDHist.GetBinContent(elIDbin)
			self.LepIDSFU[0] = self.ElIDHist.GetBinContent(elIDbin) + 0.01+ self.ElIDHist.GetBinError(elIDbin)
			self.LepIDSFD[0] = self.ElIDHist.GetBinContent(elIDbin) - 0.01 - self.ElIDHist.GetBinError(elIDbin)
	def doMCUnc(self):
		# Fill in MC variables:
		# PDF/ALPHA/SCALE:
		W = 0.
		for w in self.Tree.pdf_Weights:
			W = W + math.fabs((1-w))*math.fabs((1-w))
		size = max(self.Tree.pdf_size, 1.)
		W =math.sqrt(W/size)
		self.PDFweightU[0] = 1 + W
		self.PDFweightD[0] = 1 - W
		W = 0.
		for w in self.Tree.scale_Weights:
			W = W + math.fabs((1-w))*math.fabs((1-w))
		size = max(self.Tree.scale_size, 1.)
		W = math.sqrt(W/size)
		self.SCLweightU[0] = 1 + W
		self.SCLweightD[0] = 1 - W
		# PU:
		bin = self.puf_n.FindBin(self.Tree.evt_NGoodVtx)
		self.PUweightN[0] = self.puf_n.GetBinContent(bin)
		self.PUweightU[0] = self.puf_u.GetBinContent(bin)
		self.PUweightD[0] = self.puf_d.GetBinContent(bin)
	def getMCVar(self):
		lepiness = 0.
		for g in range(len(self.Tree.gen_ID)):
			if math.fabs(self.Tree.gen_ID[g]) == 24:
				if 10 < math.fabs(self.Tree.gen_Dau0ID[g]) < 19: 
					lepiness += 1.
		self.isLEP[0] = lepiness
		ttHT = 0.
		for g in range(len(self.Tree.gen_ID)):
			if math.fabs(self.Tree.gen_ID[g]) == 6:
				if (math.fabs(self.Tree.gen_Dau0ID[g]) == 5 and math.fabs(self.Tree.gen_Dau1ID[g]) == 24) or (math.fabs(self.Tree.gen_Dau1ID[g]) == 5 and math.fabs(self.Tree.gen_Dau0ID[g]) == 24):
					ttHT += self.Tree.gen_Pt[g]
		self.ttHT[0] = ttHT
	def getTTVar(self):
		if self.TAG.Pt() < 500.:
			self.topSF[0] = 1.15
			self.topSFu[0] = 1.15+0.19
			self.topSFd[0] = 1.15-0.19
			
		else:
			self.topSF[0] = 1.15
			self.topSFu[0] = 1.15+0.19
			self.topSFd[0] = 1.15-0.19
		ttHT = 0.
		for g in range(len(self.Tree.gen_ID)):
			if math.fabs(self.Tree.gen_ID[g]) == 6:
				if (math.fabs(self.Tree.gen_Dau0ID[g]) == 5 and math.fabs(self.Tree.gen_Dau1ID[g]) == 24) or (math.fabs(self.Tree.gen_Dau1ID[g]) == 5 and math.fabs(self.Tree.gen_Dau0ID[g]) == 24):
					ttHT += self.Tree.gen_Pt[g]
		self.ttHT[0] = ttHT
	def getBSF1(self):
		Flav = math.fabs(self.Tree.jetAK4Puppi_PartonFlavour[self.jetList[0][1]])
		if Flav > 4.:
			self.bSF[0] = self.BTagCalib_mediumB.eval_auto_bounds("central", 0, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
			self.bSFu[0] = self.BTagCalib_mediumBu.eval_auto_bounds("up", 0, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
			self.bSFd[0] = self.BTagCalib_mediumBd.eval_auto_bounds("down", 0, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
		elif 5. > Flav > 0.:
			self.bSF[0] = self.BTagCalib_mediumC.eval_auto_bounds("central", 1, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
			self.bSFu[0] = self.BTagCalib_mediumCu.eval_auto_bounds("up", 1, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
			self.bSFd[0] = self.BTagCalib_mediumCd.eval_auto_bounds("down", 1, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
		else:
			self.bSF[0] = self.BTagCalib_mediumL.eval_auto_bounds("central", 2, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
			self.bSFu[0] = self.BTagCalib_mediumLu.eval_auto_bounds("up", 2, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
			self.bSFd[0] = self.BTagCalib_mediumLd.eval_auto_bounds("down", 2, self.jetList[0][0].Eta(), self.jetList[0][0].Pt())
	def getBSF2(self):
		FlavL = math.fabs(self.Tree.jetAK4Puppi_PartonFlavour[self.LEPJETIND])
		FlavH = math.fabs(self.Tree.jetAK4Puppi_PartonFlavour[self.HADJETIND])
		if FlavL > 4.:
			self.blSF[0] = self.BTagCalib_mediumB.eval_auto_bounds("central", 0, self.LEPJET.Eta(), self.LEPJET.Pt())
			self.blSFu[0] = self.BTagCalib_mediumBu.eval_auto_bounds("up", 0, self.LEPJET.Eta(), self.LEPJET.Pt())
			self.blSFd[0] = self.BTagCalib_mediumBd.eval_auto_bounds("down", 0, self.LEPJET.Eta(), self.LEPJET.Pt())
		elif 5. > FlavL > 0.:
			self.blSF[0] = self.BTagCalib_mediumC.eval_auto_bounds("central", 1, self.LEPJET.Eta(), self.LEPJET.Pt())
			self.blSFu[0] = self.BTagCalib_mediumCu.eval_auto_bounds("up", 1, self.LEPJET.Eta(), self.LEPJET.Pt())
			self.blSFd[0] = self.BTagCalib_mediumCd.eval_auto_bounds("down", 1, self.LEPJET.Eta(), self.LEPJET.Pt())
		else:
			self.blSF[0] = self.BTagCalib_mediumL.eval_auto_bounds("central", 2, self.LEPJET.Eta(), self.LEPJET.Pt())
			self.blSFu[0] = self.BTagCalib_mediumLu.eval_auto_bounds("up", 2, self.LEPJET.Eta(), self.LEPJET.Pt())
			self.blSFd[0] = self.BTagCalib_mediumLd.eval_auto_bounds("down", 2, self.LEPJET.Eta(), self.LEPJET.Pt())
		if FlavH > 4.:
			self.bhSF[0] = self.BTagCalib_mediumB.eval_auto_bounds("central", 0, self.HADJET.Eta(), self.HADJET.Pt())
			self.bhSFu[0] = self.BTagCalib_mediumBu.eval_auto_bounds("up", 0, self.HADJET.Eta(), self.HADJET.Pt())
			self.bhSFd[0] = self.BTagCalib_mediumBd.eval_auto_bounds("down", 0, self.HADJET.Eta(), self.HADJET.Pt())
		elif 5. > FlavH > 0.:
			self.bhSF[0] = self.BTagCalib_mediumC.eval_auto_bounds("central", 1, self.HADJET.Eta(), self.HADJET.Pt())
			self.bhSFu[0] = self.BTagCalib_mediumCu.eval_auto_bounds("up", 1, self.HADJET.Eta(), self.HADJET.Pt())
			self.bhSFd[0] = self.BTagCalib_mediumCd.eval_auto_bounds("down", 1, self.HADJET.Eta(), self.HADJET.Pt())
		else:
			self.bhSF[0] = self.BTagCalib_mediumL.eval_auto_bounds("central", 2, self.HADJET.Eta(), self.HADJET.Pt())
			self.bhSFu[0] = self.BTagCalib_mediumLu.eval_auto_bounds("up", 2, self.HADJET.Eta(), self.HADJET.Pt())
			self.bhSFd[0] = self.BTagCalib_mediumLd.eval_auto_bounds("down", 2, self.HADJET.Eta(), self.HADJET.Pt())
	def FillVar1JES(self): # move the JES up and down
		TAGJESU = self.TAG*(self.Tree.jetAK8Puppi_JERSFUp[0]/self.Tree.jetAK8Puppi_JERSF[0])
		TAGJESD = self.TAG*(self.Tree.jetAK8Puppi_JERSFDown[0]/self.Tree.jetAK8Puppi_JERSF[0])
		JETJESU = self.jetList[0][0]*(self.Tree.jetAK4Puppi_JERSFUp[self.jetList[0][1]]/self.Tree.jetAK4Puppi_JERSF[self.jetList[0][1]])
		JETJESD = self.jetList[0][0]*(self.Tree.jetAK4Puppi_JERSFDown[self.jetList[0][1]]/self.Tree.jetAK4Puppi_JERSF[self.jetList[0][1]])

		TPU = JETJESU + self.W
		ZPU = TPU + TAGJESU

		TPD = JETJESD + self.W
		ZPD = TPD + TAGJESD

		self.TAGPtjesU[0] = TAGJESU.Pt() 
		self.JLep1PtjesU[0] = JETJESU.Pt()
		self.Tp1MjesU[0] = TPU.M()
		self.Zp1MjesU[0] = ZPU.M()

		self.TAGPtjesD[0] = TAGJESD.Pt()
		self.JLep1PtjesD[0] = JETJESD.Pt()
		self.Tp1MjesD[0] = TPD.M()
		self.Zp1MjesD[0] = ZPD.M()
	def FillVar2JES(self): # move the JES up and down

		TAGJESU = self.TAG*(self.Tree.jetAK8Puppi_JERSFUp[0]/self.Tree.jetAK8Puppi_JERSF[0])
		TAGJESD = self.TAG*(self.Tree.jetAK8Puppi_JERSFDown[0]/self.Tree.jetAK8Puppi_JERSF[0])
		LEPJESU = self.LEPJET*(self.Tree.jetAK4Puppi_JERSFUp[self.LEPJETIND]/self.Tree.jetAK4Puppi_JERSF[self.LEPJETIND])
		HADJESU = self.HADJET*(self.Tree.jetAK4Puppi_JERSFUp[self.LEPJETIND]/self.Tree.jetAK4Puppi_JERSF[self.LEPJETIND])
		LEPJESD = self.LEPJET*(self.Tree.jetAK4Puppi_JERSFDown[self.LEPJETIND]/self.Tree.jetAK4Puppi_JERSF[self.LEPJETIND])
		HADJESD = self.HADJET*(self.Tree.jetAK4Puppi_JERSFDown[self.LEPJETIND]/self.Tree.jetAK4Puppi_JERSF[self.LEPJETIND])

		TU = LEPJESU + self.W
		TPU = TAGJESU + HADJESU
		ZPU = TU + TPU

		TD = LEPJESD + self.W
		TPD = TAGJESD + HADJESD
		ZPD = TD + TPD

		self.TAGPtjesU[0] = TAGJESU.Pt() 
		self.JLep2PtjesU[0] = LEPJESU.Pt() 
		self.JHadPtjesU[0] = HADJESU.Pt()
		self.Tp2MjesU[0] = TPU.M()
		self.Zp2MjesU[0] = ZPU.M()
		self.TopMjesU[0] = TU.M()

		self.TAGPtjesD[0] = TAGJESD.Pt() 
		self.JLep2PtjesD[0] = LEPJESD.Pt() 
		self.JHadPtjesD[0] = HADJESD.Pt()
		self.Tp2MjesD[0] = TPD.M()
		self.Zp2MjesD[0] = ZPD.M()
		self.TopMjesD[0] = TD.M()
	def FillVar1JEC(self): # move the JEC up and down by its uncertainty and save new trees
		# might as well to Top Tag SF here: From AN-16-245
		if self.TAG.Pt() < 500.:
			self.topSF[0] = 1.15
			self.topSFu[0] = 1.15+0.19
			self.topSFd[0] = 1.15-0.19
			
		else:
			self.topSF[0] = 1.15
			self.topSFu[0] = 1.15+0.19
			self.topSFd[0] = 1.15-0.19
		TAGJECU = self.TAG*(1+self.Tree.jetAK8Puppi_jecUncertainty[0])
		TAGJECD = self.TAG*(1-self.Tree.jetAK8Puppi_jecUncertainty[0])
		JETJECU = self.jetList[0][0]*(1+self.Tree.jetAK4Puppi_jecUncertainty[self.jetList[0][1]])
		JETJECD = self.jetList[0][0]*(1-self.Tree.jetAK4Puppi_jecUncertainty[self.jetList[0][1]])

		TPU = JETJECU + self.W
		ZPU = TPU + TAGJECU

		TPD = JETJECD + self.W
		ZPD = TPD + TAGJECD

		self.TAGPtjecU[0] = TAGJECU.Pt() 
		self.JLep1PtjecU[0] = JETJECU.Pt()
		self.Tp1MjecU[0] = TPU.M()
		self.Zp1MjecU[0] = ZPU.M()

		self.TAGPtjecD[0] = TAGJECD.Pt()
		self.JLep1PtjecD[0] = JETJECD.Pt()
		self.Tp1MjecD[0] = TPD.M()
		self.Zp1MjecD[0] = ZPD.M()
	def FillVar2JEC(self): # move the JEC up and down by its uncertainty and save new trees
		TAGJECU = self.TAG*(1+self.Tree.jetAK8Puppi_jecUncertainty[0])
		TAGJECD = self.TAG*(1-self.Tree.jetAK8Puppi_jecUncertainty[0])
		LEPJECU = self.LEPJET*(1+self.Tree.jetAK4Puppi_jecUncertainty[self.LEPJETIND])
		HADJECU = self.HADJET*(1+self.Tree.jetAK4Puppi_jecUncertainty[self.HADJETIND])
		LEPJECD = self.LEPJET*(1-self.Tree.jetAK4Puppi_jecUncertainty[self.LEPJETIND])
		HADJECD = self.HADJET*(1-self.Tree.jetAK4Puppi_jecUncertainty[self.HADJETIND])

		TU = LEPJECU + self.W
		TPU = TAGJECU + HADJECU
		ZPU = TU + TPU

		TD = LEPJECD + self.W
		TPD = TAGJECD + HADJECD
		ZPD = TD + TPD

		self.TAGPtjecU[0] = TAGJECU.Pt() 
		self.JLep2PtjecU[0] = LEPJECU.Pt() 
		self.JHadPtjecU[0] = HADJECU.Pt()
		self.Tp2MjecU[0] = TPU.M()
		self.Zp2MjecU[0] = ZPU.M()
		self.TopMjecU[0] = TU.M()

		self.TAGPtjecD[0] = TAGJECD.Pt() 
		self.JLep2PtjecD[0] = LEPJECD.Pt() 
		self.JHadPtjecD[0] = HADJECD.Pt()
		self.Tp2MjecD[0] = TPD.M()
		self.Zp2MjecD[0] = ZPD.M()
		self.TopMjecD[0] = TD.M()
	def FillVar1Phys(self):
		self.Tp1M[0] = self.TPRIME1.M()
		self.Zp1M[0] = self.ZPRIME1.M()
	def FillVar2Phys(self):
		self.TopM[0] = self.LEPTOP.M()
		self.Tp2M[0] = self.TPRIME2.M()
		self.Zp2M[0] = self.ZPRIME2.M()	
	def FillVar1AK4(self):
		self.JLep1Pt[0] = self.jetList[0][0].Pt()
		self.JLep1CSV[0] = self.Tree.jetAK4Puppi_CSVv2[self.jetList[0][1]]
	def FillVar2AK4(self):
		self.JLep2Pt[0] = self.LEPJET.Pt()
		self.JLep2CSV[0] = self.Tree.jetAK4Puppi_CSVv2[self.LEPJETIND]
		self.JHadPt[0] = self.HADJET.Pt()
		self.JHadCSV[0] = self.Tree.jetAK4Puppi_CSVv2[self.HADJETIND]
	def FillVarAK8(self):
		self.TAGM[0] = self.Tree.jetAK8Puppi_softDropMass[0]
		self.TAGPt[0] = self.TAG.Pt()
		self.TAGTau32[0] = self.Tree.jetAK8Puppi_tau3[0]/self.Tree.jetAK8Puppi_tau2[0]
		self.TAGTau21[0] = self.Tree.jetAK8Puppi_tau2[0]/self.Tree.jetAK8Puppi_tau1[0]
	def FillVarLeptons(self):
		self.LepPt[0] = self.LEP.Pt()
		self.METPt[0] = self.MET.Pt()
		self.WPt[0]   = self.W.Pt()
		#print "--=--=--=--=--=--"
		#print "lep pT = " + str(self.LEP.Pt())
		if self.leptype == "mu":
			self.LepType[0] = 1.
			if self.Tree.mu_IsTightMuon[0] > 0. :
				self.LepTight[0] = 3.
			elif self.Tree.mu_IsMediumMuon[0] > 0.:
				self.LepTight[0] = 2.
			elif self.Tree.mu_IsLooseMuon[0] > 0.:
				self.LepTight[0] = 1.
			else :
				self.LepTight[0] = 0.
			self.LepIso[0] = self.Tree.mu_MiniIso[0]
			self.Lep2D_dr[0] = self.Tree.mu_AK4JetV2DR[0]
			self.Lep2D_rel[0] = self.Tree.mu_AK4JetV2PtRel[0]
			#print "2D cut DRs are " + str(self.Tree.mu_AK4JetV1DR[0]) + ", " + str(self.Tree.mu_AK4JetV2DR[0]) + ", " + str(self.Tree.mu_AK4JetV3DR[0])
			#print "2D cut RelPtss are " + str(self.Tree.mu_AK4JetV1PtRel[0]) + ", " + str(self.Tree.mu_AK4JetV2PtRel[0]) + ", " + str(self.Tree.mu_AK4JetV3PtRel[0])
		if self.leptype == "el":
			self.LepType[0] = -1.	
			if self.Tree.el_vidTight[0] > 0. :
				self.LepTight[0] = 3.
			elif self.Tree.el_vidMedium[0] > 0.:
				self.LepTight[0] = 2.
			elif self.Tree.el_vidLoose[0] > 0.:
				self.LepTight[0] = 1.
			else :
				self.LepTight[0] = 0.
			self.LepIso[0] = self.Tree.el_MiniIso[0]
			self.Lep2D_dr[0] = self.Tree.el_AK4JetV2DR[0]
			self.Lep2D_rel[0] = self.Tree.el_AK4JetV2PtRel[0]
			#print "2D cut DRs are " + str(self.Tree.el_AK4JetV1DR[0]) + ", " + str(self.Tree.el_AK4JetV2DR[0]) + ", " + str(self.Tree.el_AK4JetV3DR[0])
			#print "2D cut RelPtss are " + str(self.Tree.el_AK4JetV1PtRel[0]) + ", " + str(self.Tree.el_AK4JetV2PtRel[0]) + ", " + str(self.Tree.el_AK4JetV3PtRel[0])
			if self.Lep2D_dr[0] > 5 or self.Lep2D_rel[0] > 1000:
				print "============================="
				print str(self.Lep2D_dr[0]) + " (DR)"
				print str(self.LEP.Pt()) + " (pT)"
				print str(self.LEP.Eta()) + " (eta)"
				print str(self.Lep2D_rel[0]) + " (rel)"
				self.sadnumber += 1
	def doType1Vars(self):
		self.TPRIME1 = self.W + self.jetList[0][0] # make the T' candidate: no ambiguity for type 1 channel, the light jet must be the hardest jet (not including the AK8 jets)
		if self.TPRIME1.M() > 105: # quality cut, below 120 the TPRIME isn't even a top
			self.ZPRIME1 = self.TPRIME1 + self.TAG # Here is the ZPRIME candidate!
			return True
		return False
	def doType2Vars(self):
		# first pick a lepjet and a hadjet (based on distance to the lepton)
		if self.jetList[0][0].DeltaR(self.LEP) < self.jetList[1][0].DeltaR(self.LEP):
			self.LEPJET = self.jetList[0][0]
			self.LEPJETIND = self.jetList[0][1]
			self.HADJET = self.jetList[1][0]
			self.HADJETIND = self.jetList[1][1]
		if self.jetList[1][0].DeltaR(self.LEP) < self.jetList[0][0].DeltaR(self.LEP):
			self.LEPJET = self.jetList[1][0]
			self.LEPJETIND = self.jetList[1][1]
			self.HADJET = self.jetList[0][0]
			self.HADJETIND = self.jetList[0][1]
		# Now make the variables:
		self.TPRIME2 = self.HADJET + self.TAG
		self.LEPTOP = self.LEPJET + self.W
		if self.TPRIME2.M() > 105 and self.LEPTOP.M() > 105:
			self.ZPRIME2 = self.TPRIME2 + self.LEPTOP
			return True
		return False
	def doAK4JetsMC(self):
		self.jetList = []
		self.jetInd = []
		for j in range(len(self.Tree.jetAK4Puppi_Pt)):
			iJet = 	ROOT.TLorentzVector()
			iJet.SetPtEtaPhiE(self.Tree.jetAK4Puppi_SmearedPt[j],self.Tree.jetAK4Puppi_Eta[j],self.Tree.jetAK4Puppi_Phi[j],self.Tree.jetAK4Puppi_E[j])
			if iJet.Pt() > 25. and iJet.DeltaR(self.TAG) > 0.8:
				self.jetList.append([iJet,j])
		if len(self.jetList) > 0:
			return True
		return False
	def doAK4Jets(self):
		self.jetList = []
		for j in range(len(self.Tree.jetAK4Puppi_Pt)):
			iJet = 	ROOT.TLorentzVector()
			iJet.SetPtEtaPhiE(self.Tree.jetAK4Puppi_Pt[j],self.Tree.jetAK4Puppi_Eta[j],self.Tree.jetAK4Puppi_Phi[j],self.Tree.jetAK4Puppi_E[j])
			if iJet.Pt() > 15. and iJet.DeltaR(self.TAG) > 0.8:
				self.jetList.append([iJet,j])
		if len(self.jetList) > 0:
			return True
		return False
	def doWPart(self):
		self.MET = ROOT.TLorentzVector()
		self.MET.SetPtEtaPhiM(self.Tree.met_Pt[0],0.,self.Tree.met_Phi[0],0.)
		self.W = MakeW(self.MET,self.LEP)
		if self.W.M() > 80.3 and self.W.M() < 80.5:
			return True
		return False
	def doAK8PartMC(self):
		if self.Tree.jetAK8Puppi_softDropMass[0] > 50 and self.Tree.jetAK8Puppi_Pt[0] > 35 and self.Tree.jetAK8Puppi_tau2[0]>0 and self.Tree.jetAK8Puppi_tau1[0]>0:
			self.Type1 = False
			self.TAG = ROOT.TLorentzVector() # we will be using this to store the kinematic information of the AK8 jet
			self.TAG.SetPtEtaPhiE(self.Tree.jetAK8Puppi_SmearedPt[0], self.Tree.jetAK8Puppi_Eta[0], self.Tree.jetAK8Puppi_Phi[0], self.Tree.jetAK8Puppi_E[0]) # MC uses smeared values
			if self.Tree.jetAK8Puppi_softDropMass[0] > 70 and self.Tree.jetAK8Puppi_Pt[0] > 75:
				self.Type1 = True
			return True
		return False
	def doAK8Part(self):
		if self.Tree.jetAK8Puppi_softDropMass[0] > 50 and self.Tree.jetAK8Puppi_Pt[0] > 35 and self.Tree.jetAK8Puppi_tau2[0]>0 and self.Tree.jetAK8Puppi_tau1[0]>0:
			self.Type1 = False
			self.TAG = ROOT.TLorentzVector() # we will be using this to store the kinematic information of the AK8 jet
			self.TAG.SetPtEtaPhiE(self.Tree.jetAK8Puppi_Pt[0], self.Tree.jetAK8Puppi_Eta[0], self.Tree.jetAK8Puppi_Phi[0], self.Tree.jetAK8Puppi_E[0]) # data uses raw values
			if self.Tree.jetAK8Puppi_softDropMass[0] > 70 and self.Tree.jetAK8Puppi_Pt[0] > 75:
				self.Type1 = True
			return True
		return False
	def doLepPart(self): # Fill leptoni information and retain self.LEP for further analysis.
		self.leptype = "none"
		if len(self.Tree.mu_Pt) > 0 and not len(self.Tree.el_Pt) > 0: # Just one muon
			self.leptype = "mu"
		if len(self.Tree.el_Pt) > 0 and not len(self.Tree.mu_Pt) > 0: # Just one elec
			self.leptype = "el"
		if len(self.Tree.mu_Pt) > 0 and len(self.Tree.el_Pt) > 0: # Need to pick a lepton
			if self.Tree.mu_Pt[0] > self.Tree.el_Pt[0] and self.Tree.el_Pt[0] < 50: # It's a muon
				self.leptype = "mu"
			elif self.Tree.mu_Pt[0] < self.Tree.el_Pt[0] and self.Tree.mu_Pt[0] < 50: # It's an electron
				self.leptype = "el"
		return self.FillLepton(self.leptype)
	def FillLepton(self, leptype): # Cut leptons that don't have good pT or eta. Fill variables if they're ok.
		"""		
		if self.leptype == "mu" and self.Tree.mu_AK4JetV2DR[0] < 50. and self.Tree.mu_Pt[0] > 50 and self.Tree.mu_Pt[0] < 4000 and math.fabs(self.Tree.mu_Eta[0]) < 2.1:
			self.LEP = ROOT.TLorentzVector() # we will be using this to store the kinematic information of the lepton
			self.LEP.SetPtEtaPhiE(self.Tree.mu_Pt[0], self.Tree.mu_Eta[0], self.Tree.mu_Phi[0], self.Tree.mu_E[0])
			return True
		if self.leptype == "el" and self.Tree.el_AK4JetV2DR[0] < 50. and self.Tree.el_Pt[0] > 50 and self.Tree.el_Pt[0] < 4000 and math.fabs(self.Tree.el_Eta[0]) < 2.1:
			self.LEP = ROOT.TLorentzVector() # we will be using this to store the kinematic information of the lepton
			self.LEP.SetPtEtaPhiE(self.Tree.el_Pt[0], self.Tree.el_Eta[0], self.Tree.el_Phi[0], self.Tree.el_E[0])
			return True
		"""
		if self.leptype == "el" and self.Tree.el_Pt[0] > 32 and self.Tree.el_Pt[0] < 4000 and math.fabs(self.Tree.el_Eta[0]) < 2.1 and self.Tree.el_vidTight[0] > 0.:
			self.LEP = ROOT.TLorentzVector() # we will be using this to store the kinematic information of the lepton
			self.LEP.SetPtEtaPhiE(self.Tree.el_Pt[0], self.Tree.el_Eta[0], self.Tree.el_Phi[0], self.Tree.el_E[0])
			return True
		if self.leptype == "mu":
			return False
		return False
	def doTrigger(self):
		if self.Tree.HLT_Mu45_eta2p1 > 0 and not self.Tree.HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50 > 0:
			return True
		if self.Tree.HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50 > 0 and not self.Tree.HLT_Mu45_eta2p1 > 0:
			return True
		return False
	def __del__(self):
	        print "done! " + str(self.sadnumber)


# NOTES:


