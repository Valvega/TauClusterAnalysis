import ROOT
from ROOT import TFile, TTree, TMath, TChain, TH2F, TH1F
import sys
import os
import argparse
import ast
import glob
import math
from  ConfigParser import *
import modules.Calc 
from modules.Calc import v,ctau, DeltaR, DeltaPhi, DeltaEta

def GetChain(files_path,treename,signal):
    chain = TChain(treename)
    filelist = glob.glob(files_path+"/%s/*.root"%signal)
    for file in filelist:
        chain.Add(file)
    return chain    

def Looper(input_dir,output_dir,signal,maxevents):
    #Open input, nevents, and setup output file
    t       = GetChain(input_dir,'MuonSystem',signal)
    tmp     = TFile.Open("%s/histos_%s.root"%(output_dir,signal),'RECREATE')
    Nevents = t.GetEntries()
    if maxnevents>0 and maxnevents<=Nevents: Nevents = maxnevents


    #Book Histograms
    #Gen VLLs
    h_vll_e         = TH1F("h_vll_e",  "h_vll_e"  ,40,0,3000)
    h_vll_pt        = TH1F("h_vll_pt", "h_vll_pt" ,24,0,1200)
    h_vll_phi       = TH1F("h_vll_phi","h_vll_phi",14,-3.5,3.5)
    h_vll_eta       = TH1F("h_vll_eta","h_vll_eta",20,-5.0,5.0)

    #Gen taus
    h_tau_e         = TH1F("h_tau_e",  "h_tau_e"  ,32,0,1600)
    h_tau_pt        = TH1F("h_tau_pt", "h_tau_pt" ,24,0,1200)
    h_tau_phi       = TH1F("h_tau_phi","h_tau_phi",14,-3.5,3.5)
    h_tau_eta       = TH1F("h_tau_eta","h_tau_eta",20,-5.0,5.0)

    #Gen llps
    h_llp_e         = TH1F("h_llp_e",     "h_llp_e"  ,32,0,1600)
    h_llp_pt        = TH1F("h_llp_pt",    "h_llp_pt" ,24,0,1200)
    h_llp_phi       = TH1F("h_llp_phi",   "h_llp_phi",14,-3.5,3.5)
    h_llp_eta       = TH1F("h_llp_eta",   "h_llp_eta",20,-5.0,5.0)

    #Loop over events
    for e in range(0, Nevents):

        #Get entries from tree
        t.GetEntry(e)

        #Generator level variables
        for k in range(0,2):
           h_vll_e.Fill(  t.gVLL_e[k]  )
           h_vll_pt.Fill( t.gVLL_pt[k] )
           h_vll_phi.Fill(t.gVLL_phi[k])
           h_vll_eta.Fill(t.gVLL_eta[k])
           h_tau_e.Fill(  t.gTau_e[k]  )
           h_tau_pt.Fill( t.gTau_pt[k] )
           h_tau_phi.Fill(t.gTau_phi[k])
           h_tau_eta.Fill(t.gTau_eta[k])
           h_llp_e.Fill(  t.gLLP_e[k]  )
           h_llp_pt.Fill( t.gLLP_pt[k] )
           h_llp_phi.Fill(t.gLLP_phi[k])
           h_llp_eta.Fill(t.gLLP_eta[k])

    tmp.Write()
    tmp.Close()

###########OPTIONS
parser = argparse.ArgumentParser(description='Command line parser of skim options')
parser.add_argument('--config' ,       dest='cfgfile'   ,  help='Name of config file'   , required = True)
parser.add_argument('--tag',           dest='tag'       ,  help='Name of output folder' , required = True)
parser.add_argument('--maxnevents',    dest='maxnevents',  help='Number of maxevents'   , type=int, default  = -1 )
args           = parser.parse_args()
configfilename = args.cfgfile
tagname        = args.tag
maxnevents     = args.maxnevents

#Reading configuration
print "[INFO] Reading skim configuration file . . ."
cfgparser = ConfigParser()
cfgparser.read('%s'%configfilename)
print "[INFO] Getting configuration parameters . . ."
year         = ast.literal_eval(cfgparser.get("histoconfiguration","year"))
print "    -The year:"
print "      *",year
input_dir    = ast.literal_eval(cfgparser.get("histoconfiguration","input_dir"))
print "    -The input  directory:"
print "      *",input_dir+tagname
output_dir   = ast.literal_eval(cfgparser.get("histoconfiguration","output_dir"))
print "    -The output directory:"
print "      *",output_dir+tagname+"/%i"%year
signals      = ast.literal_eval(cfgparser.get("histoconfiguration","signals"))
print "    -The list of MC signal samples:"
for x in range(len(signals)): print "      *",signals[x]


#Create output directories
os.system("mkdir histograms/%s"%tagname)
os.system("mkdir histograms/%s/%i"%(tagname,year))
output_dir =  output_dir+tagname+"/%i"%year

#Run looper over signal samples
for k in range(0,len(signals)):
    print "[INFO] Running looper over %s"%signals[k]
    Looper(input_dir,output_dir,signals[k],maxnevents)