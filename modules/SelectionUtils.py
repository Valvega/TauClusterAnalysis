import ROOT
from ROOT import TFile, TTree, TMath, TH1F,TChain,TH2F
import sys
import os
import argparse
import ast
import glob
import math
import string 
import modules.Calc 
import modules.DataProcess 
from modules.Calc import v,ctau,DeltaR, DeltaPhi, DeltaEta

def SelectDtCluster(t):
        maxsize=-1
        clsid  =-1
        ncls   = False 
        for i in range(0,t.nDtRechitClusters):
            matchedtotau = False
            for k in range(0,t.nTau):
                deltaR = DeltaR(t.tauEta[k],t.tauPhi[k],t.dtRechitClusterEta[i],t.dtRechitClusterPhi[i])
                if deltaR < 0.2: matchedtotau=True 
            if matchedtotau: continue
            ncls = True
            if t.dtRechitClusterSize[i]>=maxsize:
                maxsize = t.dtRechitClusterSize[i]
                clsid   = i
        return ncls,clsid

def SelectCSCCluster(t):
        maxsize=-1
        clsid  =-1
        ncls   = False  
        for i in range(0,t.nCscRechitClusters):
            matchedtotau = False
            for k in range(0,t.nTau):
                deltaR = DeltaR(t.tauEta[k],t.tauPhi[k],t.cscRechitClusterEta[i],t.cscRechitClusterPhi[i])
                if deltaR < 0.2: matchedtotau=True 
            if matchedtotau: continue
            ncls = True
            if t.cscRechitClusterSize[i]>=maxsize:
                maxsize = t.cscRechitClusterSize[i]
                clsid   = i
        return ncls,clsid

def INTimeCSC(t,clsid):
    ipass = False
    if t.cscRechitClusterTimeWeighted[clsid]> -5 and t.cscRechitClusterTimeWeighted[clsid] < 12.5: ipass=True
    return ipass

def OOTimeCSC(t,clsid):
    ipass = False
    if t.cscRechitClusterTimeWeighted[clsid] < -12.5: ipass=True
    return ipass

def INTimeDT(t,clsid):
    ipass = False
    if t.dtRechitCluster_match_RPCBx_dPhi0p5[clsid]==0: ipass=True
    return ipass

def OOTimeDT(t,clsid):
    ipass = False
    if t.dtRechitCluster_match_RPCBx_dPhi0p5[clsid] <0: ipass=True
    return ipass

def SelectTau(t):
        maxtaupt_m  =-1
        tauid_m     =-1
        ntau_m      = False
        for i in range(0,t.nTau):
            if t.tauPt[i] > 50 and abs(t.tauEta[i]) <= 2.3 and t.tau_IsMedium[i]==True:
                   ntau_m=True
                   if t.tauPt[i]>=maxtaupt_m:
                       maxtaupt_m = t.tauPt[i]
                       tauid_m    = i    
        return ntau_m,tauid_m

def SelectAntiTau(t):
        ntau        = False
        maxtaupt_a  =-1
        tauid_a     =-1
        ntau_a      = False
        #Check if there is any medium tau
        for i in range(0,t.nTau):
            if t.tauPt[i]> 50 and abs(t.tauEta[i]) <= 2.3 and t.tau_IsMedium[i]==True: ntau=True
        if ntau==True: return False,-1
        #If not, get the highest pt tau f
        for i in range(0,t.nTau):
            if t.tauPt[i]> 50 and abs(t.tauEta[i]) <= 2.3: # and t.tau_IsLoose[i]==True:
                   ntau_a=True
                   if t.tauPt[i]>=maxtaupt_a:
                       maxtaupt_a = t.tauPt[i]
                       tauid_a    = i                  
        return ntau_a,tauid_a

def CSCMuonJetVeto(t,clsid):
    ipass = True
    if t.cscRechitClusterMuonVetoPt[clsid] > 20: ipass = False
    if t.cscRechitClusterJetVetoPt[clsid]  > 10: ipass = False
    return ipass
