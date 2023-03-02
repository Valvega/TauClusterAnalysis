import ROOT
from ROOT import TFile, TTree, TMath, TH1F,TChain,TH2F, TF1
import sys
import os
import argparse
import ast
import glob
import math
import string 
from  ConfigParser import *
import modules.Calc 
import modules.DataProcess
import modules.SelectionUtils
from modules.Calc import v,ctau,DeltaR, DeltaPhi, DeltaEta
from prettytable import PrettyTable
from modules.SelectionUtils import SelectDtCluster, SelectTau, SelectAntiTau, OOTimeDT, INTimeDT
from modules.DataProcess import GetChainSingle, GetChainMultiple
ROOT.gROOT.SetBatch(True)

def MakeDTEfficiencyTable(info,signal,isData):
        os.system("mkdir studies")
        os.system("mkdir studies/taudtcutflow")
        x = PrettyTable( ['Selections','Cut. Eff.(%)','Cul. Eff.(%)','Exp. Events'])
        x.align["Selections"]    = "l"
        x.align["Cut. Eff.(%)"]  = "r"
        x.align["Cul. Eff.(%)"]  = "r"
        x.align["Exp. Events"]   = "r"
        entries = len(info)
        if isData:
                for k in range(2,entries):
                    x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k-1][1]),'%.5f'%(100*info[k][1]/info[k][2][0]),'%.1f'%info[k][1]])                                  
        else:
                for k in range(0,entries):
                      if k==0:
                        x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k][2][0]),'%.5f'%(100*info[k][1]/info[k][2][0]),'%.1f'%info[k][1]])                       
                      elif k==1: 
                        x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k-1][1]),'%.5f'%(100*info[k][1]/info[k][2][0]),'%.1f'%info[k][1]])
                        x.add_row(['--------------------','--------------','--------------','-------------'])
                      else: 
                        if info[k-1][1]==0 or info[k][2][1]==0: 
                            x.add_row(['%s'%info[k][0],'0.00','0.00000','0.0'])
                        else:
                            x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k-1][1]),'%.5f'%(100*info[k][1]/info[k][2][1]),'%.1f'%info[k][1]])
        print x
        with open('studies/taudtcutflow/%s.txt'%(signal), 'w') as the_file:
           the_file.write(str(x))

def DTLooper(input_dir,output_dir,sample,samplename,isData,isSignal,signorm):
    #Open input and output file
    if isData==False:
        t, cf, ac_csc, ac_dt = GetChainSingle(input_dir,'MuonSystem',sample,isData)
        tmp  = TFile.Open("%s/histos_tau_dtcutflow_%s.root"%(output_dir,samplename),'RECREATE')
    else:
        t, cf, ac_csc, ac_dt = GetChainMultiple(input_dir,'MuonSystem',sample,isData)
        tmp  = TFile.Open("%s/histos_tau_dtcutflow_%s.root"%(output_dir,samplename),'RECREATE')   

    cut6    = 0
    cut7    = 0
    cut8    = 0
    cut9    = 0
    cut10   = 0
    cut11   = 0
    cut12   = 0
    cut13   = 0

    #CUT (in-time)
    h_cut8_clustertime            = TH1F("h_cut8_clustertime",         "h_cut8_clustertime",        20,-10,10)
    h_cut8_cut_dphi_cls_met       = TH1F("h_cut8_cut_dphi_cls_met",    "h_cut8_cut_dphi_cls_met",   4, 0, 3.2)
    h_cut9_clustertime            = TH1F("h_cut9_clustertime",         "h_cut9_clustertime",        20,-10,10)
    h_cut9_cut_dphi_cls_met       = TH1F("h_cut9_cut_dphi_cls_met",    "h_cut9_cut_dphi_cls_met",   4, 0, 3.2)

    #Loop over events
    for e in range(0, t.GetEntries()):
        t.GetEntry(e)
        #acceptance flag for mc
        if isData==False and isSignal==True and t.gLLP_dt[0]==0 and t.gLLP_dt[1]==0: continue
        #CUT (Cluster)
        clstagged,clsid = SelectDtCluster(t)
        if clstagged==False: continue
        cut6+=1
        plot = False
        if (isData==False) or (isData==True and t.dtRechitClusterSize[clsid]<80): plot=True
        #CUT (Tau)
        tautagged,tauid = SelectTau(t)
        if tautagged==False: continue
        cut7+=1
        #CUT (jet veto)
        if t.dtRechitClusterJetVetoPt[clsid]  > 20: continue
        cut8+=1
        if plot==True:
          h_cut8_clustertime.Fill(t.dtRechitCluster_match_RPCBx_dPhi0p5[clsid])
          h_cut8_cut_dphi_cls_met.Fill(abs(t.dtRechitClusterMet_dPhi[clsid]))
        #CUT (muon veto)
        if t.dtRechitClusterMuonVetoPt[clsid] > 10: continue
        cut9+=1
        if plot==True:
          h_cut9_clustertime.Fill(t.dtRechitCluster_match_RPCBx_dPhi0p5[clsid])
          h_cut9_cut_dphi_cls_met.Fill(abs(t.dtRechitClusterMet_dPhi[clsid]))
        #CUT (MB1 Adjacent)
        if t.dtRechitCluster_match_MB1hits_cosmics_minus[clsid] > 8 or t.dtRechitCluster_match_MB1hits_cosmics_plus[clsid] > 8: continue
        cut10+=1
        #CUT (RPCmatching)
        if t.dtRechitCluster_match_RPChits_dPhi0p5[clsid] <=0: continue
        cut11+=1
        #Cut (in time) 
        if INTimeDT(t,clsid)==True:
            cut12+=1
            #CUT (dphi) 
            if t.dtRechitClusterMet_dPhi[clsid]>=(math.pi/2):  continue
            cut13+=1

    #Collect dt cutflow info for table
    nums = []
    if isSignal==True:
        w = signorm/cf.GetBinContent(1) 
    else:
        w = 1 
    if isSignal==True:
        nums.append(w*ac_dt.GetBinContent(1))   #Acceptance
        nums.append(w*ac_dt.GetBinContent(2))   #METtrigger+met>200
        nums.append(w*ac_dt.GetBinContent(3))   #MET filters
        nums.append(w*ac_dt.GetBinContent(4))   #nrings
        nums.append(w*ac_dt.GetBinContent(5))    #at least one cluster   
    else:
        nums.append(w*cf.GetBinContent(1))     #Total
        nums.append(w*cf.GetBinContent(2))     #METtrigger+met>200
        nums.append(w*cf.GetBinContent(3))     #MET filters
        nums.append(w*cf.GetBinContent(4))     #nrings
        nums.append(w*cf.GetBinContent(5))     #at least one cluster        
    nums.append( w*float(cut6 ) )  
    nums.append( w*float(cut7 ) ) 
    nums.append( w*float(cut8 ) ) 
    nums.append( w*float(cut9 ) ) 
    nums.append( w*float(cut10) ) 
    nums.append( w*float(cut11) ) 
    nums.append( w*float(cut12) ) 
    if isSignal==True: 
        tag = 'Acceptance'
        den = [w*cf.GetBinContent(1),w*ac_dt.GetBinContent(2)]
    else:
        den = [w*cf.GetBinContent(2),w*cf.GetBinContent(2)]
        tag = 'Total'
    if isData: 
        nums.append(float(cut12)/2)  #Blinded, so model
    else: 
        nums.append(float(cut13))  

    effinfo = [
             [tag                        ,nums[0] , den],
             ['MET Trigger and MET200'   ,nums[1] , den],
             ['MET filters'              ,nums[2] , den], 
             ['CSC+DT rings <= 10'       ,nums[3] , den],
             ['# clusters >= 1'          ,nums[4] , den],
             ['# DT clusters >= 1'       ,nums[5] , den],
             ['# tauhads >= 1'           ,nums[6] , den], 
             ['Jet Veto'                 ,nums[7] , den], 
             ['Muon Veto'                ,nums[8] , den],
             ['MB1 Adjacent'             ,nums[9] , den], 
             ['RPC Matching'             ,nums[10], den], 
             ['Time (BX=0)'              ,nums[11], den], 
             #['|dPhi(cls,MET)| < pi/2'   ,nums[12], den],
    ]
    tmp.Write()
    tmp.Close()
    return effinfo

###########OPTIONS
parser = argparse.ArgumentParser(description='Command line parser of skim options')
parser.add_argument('--config' ,  dest='cfgfile',  help='Name of config file',  required = True)
parser.add_argument('--tag',    dest='tag',  help='Name of tag', required = True)
args           = parser.parse_args()
configfilename = args.cfgfile
tagname        = args.tag

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
print "      *",input_dir+tagname+"/%i"%year
output_dir   = ast.literal_eval(cfgparser.get("histoconfiguration","output_dir"))
print "    -The output directory:"
print "      *",output_dir+tagname+"/%i"%year
signals      = ast.literal_eval(cfgparser.get("histoconfiguration","signals"))
print "    -The list of MC signal samples:"
for x in range(len(signals)): print "      *",signals[x]
datas        = ast.literal_eval(cfgparser.get("histoconfiguration","datas"))
print "    -The list of data files:"
for x in range(len(datas)): print "      *",datas[x]
lumi        = ast.literal_eval(cfgparser.get("histoconfiguration","lumi"))
print "    -The lumi:"
print "      *",lumi
sigxs      = ast.literal_eval(cfgparser.get("histoconfiguration","sigxs"))
print "    -The list of MC signal xs (in pb):"
for x in range(len(sigxs)): print "      *",sigxs[x]

#Create output directories
os.system("mkdir histograms")
os.system("mkdir histograms/%s"%tagname)
os.system("mkdir histograms/%s/%i"%(tagname,year))
output_dir =  output_dir+tagname+"/%i"%year

#Run looper over signal samples
isData=False
isSignal=True
for k in range (0,len(signals) ):
    print "[INFO] Running looper over %s"%signals[k]
    effinfo = DTLooper(input_dir+tagname,output_dir,signals[k],signals[k],isData,isSignal,sigxs[k]*lumi)
    MakeDTEfficiencyTable(effinfo,signals[k],isData)

#Run looper over data samples
isData=True
isSignal=False
print "[INFO] Running looper over data"
effinfo = DTLooper(input_dir+tagname,output_dir,datas,'Data',isData,isSignal,sigxs[k]*lumi)
MakeDTEfficiencyTable(effinfo,'Data',isData)

