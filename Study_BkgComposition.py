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
from modules.Calc import DeltaR, DeltaPhi, DeltaEta

def GetChain(files_path,treename,signal):
    chain         = TChain(treename)
    cutflow       = TH1F("CutFlow","CutFlow",5,0,5)
    #print files_path+"/%s/*.root"%signal
    filelist = glob.glob(files_path+"/%s/*.root"%signal)
    for file in filelist:
        #Check file has a tree and cutflow histo, else skip it
        ftmp       = TFile.Open(file)
        if ftmp.GetListOfKeys().Contains(treename)==False: 
            print "This file doesn't have a tree, have a look at it: ",file 
            continue        
        chain.Add(file)   
        if ftmp.GetListOfKeys().Contains('CutFlow')==True: 
           ctmp       = ftmp.Get("CutFlow")
           cutflow.Add(ctmp)
        #Add tree and cutflow histograms
    return chain,cutflow


def Looper(input_dir,output_dir,signal,xs,lumi,maxevents):
    #Open input, nevents, and setup output file
    t,cutflow  = GetChain(input_dir,'tree',signal)
    tmp     = TFile.Open("%s/histos_%s.root"%(output_dir,signal),'RECREATE')
    Nevents = t.GetEntries()
    if maxnevents>0 and maxnevents<=Nevents: Nevents = maxnevents

    #Variables
    crosssection  = xs
    luminosity    = lumi
    totalMCevents = cutflow.GetBinContent(1)

    #counter
    counter = 0
    #Loop over events
    for e in range(0, Nevents):

        #Get entries from tree
        t.GetEntry(e)
        tauexist = False
        for k in range(0, t.nTau):
            if t.tauPt[e]> 50: 
                tauexist = True 

        if tauexist == False: continue
        counter+=1

    print "Number of events with >=1 tau: ",counter

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
print "[INFO] Reading configuration file . . ."
cfgparser = ConfigParser()
cfgparser.read('%s'%configfilename)
print "[INFO] Getting configuration parameters . . ."
year         = ast.literal_eval(cfgparser.get("histobkgconfiguration","year"))
print "    -The year:"
print "      *",year
lumi         = ast.literal_eval(cfgparser.get("histobkgconfiguration","lumi"))
print "    -The lumi:"
print "      *",lumi
input_dir    = ast.literal_eval(cfgparser.get("histobkgconfiguration","input_dir"))
print "    -The input  directory:"
print "      *",input_dir+tagname
output_dir   = ast.literal_eval(cfgparser.get("histobkgconfiguration","output_dir"))
print "    -The output directory:"
print "      *",output_dir+tagname+"/%i"%year
bkgs      = ast.literal_eval(cfgparser.get("histobkgconfiguration","bkgs"))
print "    -The list of MC bkg samples:"
for x in range(len(bkgs)): print "      *",bkgs[x]
bkgxs      = ast.literal_eval(cfgparser.get("histobkgconfiguration","bkgxs"))
print "    -The list of MC bkg xs:"
for x in range(len(bkgxs)): print "      *",bkgxs[x]

#Create output directories
os.system("mkdir histograms/%s"%tagname)
os.system("mkdir histograms/%s/%i"%(tagname,year))
output_dir =  output_dir+"/%s/%i"%(tagname,year)

#Run looper over signal samples
for k in range(0,len(bkgs)):
    print "[INFO] Running looper over %s"%bkgs[k]
    Looper(input_dir,output_dir,bkgs[k],bkgxs[k],lumi,maxnevents)
