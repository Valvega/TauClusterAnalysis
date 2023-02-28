import ROOT
from ROOT import TFile, TTree, TMath, TH1F,TChain,TH2F
import sys
import os
import argparse
import ast
import glob
import math
import string 
from  ConfigParser import *

def reweightcalculator(lenghtrf,oldctau,newctau):
    k        = float(oldctau)/float(newctau)
    alpha    = (lenghtrf)*(  1/float(newctau) - 1/float(oldctau)   )
    reweight = k*k*math.exp( -alpha   )
    return reweight

def GetChainSingle(files_path,treename,signal,isData):
    chain   = TChain(treename)
    cutflow       = TH1F("CutFlow","CutFlow",10,0,10)
    acceptancecsc = TH1F("AcceptanceCsc","AcceptanceCsc",10,0,10)
    acceptancedt  = TH1F("AcceptanceDt","AcceptanceDt" ,10,0,10)
    print files_path+"/%s/*.root"%signal
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
        if ftmp.GetListOfKeys().Contains('AcceptanceCsc')==True: 
            acsctmp    = ftmp.Get("AcceptanceCsc")
            acceptancecsc.Add(acsctmp)
        if ftmp.GetListOfKeys().Contains('AcceptanceDt')==True: 
            adttmp = ftmp.Get("AcceptanceDt")
            acceptancedt.Add(adttmp)
        #Add tree and cutflow histograms
    return chain,cutflow,acceptancecsc,acceptancedt   

def GetChainMultiple(files_path,treename,signals,isData):
    chain         = TChain(treename)
    cutflow       = TH1F("CutFlow","CutFlow",10,0,10)
    acceptancecsc = TH1F("AcceptanceCsc","AcceptanceCsc",10,0,10)
    acceptancedt  = TH1F("AcceptanceDt","AcceptanceDt" ,10,0,10)
    for signal in signals:
        print files_path+"/%s/*.root"%signal
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
            if ftmp.GetListOfKeys().Contains('AcceptanceCsc')==True: 
                acsctmp    = ftmp.Get("AcceptanceCsc")
                acceptancecsc.Add(acsctmp)
            if ftmp.GetListOfKeys().Contains('AcceptanceDt')==True: 
                adttmp = ftmp.Get("AcceptanceDt")
                acceptancedt.Add(adttmp)

    return chain,cutflow,acceptancecsc,acceptancedt 
