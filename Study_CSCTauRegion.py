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
from modules.SelectionUtils import SelectCSCCluster, SelectTau, SelectAntiTau, INTimeCSC, OOTimeCSC
from modules.DataProcess import GetChainSingle, GetChainMultiple
from modules.Calc import v,ctau,DeltaR, DeltaPhi, DeltaEta
from prettytable import PrettyTable
ROOT.gROOT.SetBatch(True)

def MakeCSCEfficiencyTable(info,signal,isData):
        os.system("mkdir studies")
        os.system("mkdir studies/taucsccutflow")
        x = PrettyTable( ['Selections','Cut. Eff.(%)','Cul. Eff.(%)','Exp. Events +/- Uncertainty'])
        x.align["Selections"]    = "l"
        x.align["Cut. Eff.(%)"]  = "r"
        x.align["Cul. Eff.(%)"]  = "r"
        x.align["Exp. Events +/- Uncertainty"]   = "r"
        entries = len(info)
        if isData:
                for k in range(2,entries):
                    x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k-1][1]),'%.5f'%(100*info[k][1]/info[k][2][0]),'%.3f +/- %.4f'%(info[k][1], info[k][3])])                                  
        else:
                for k in range(0,entries):
                      if k==0:
                        x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k][2][0]),'%.5f'%(100*info[k][1]/info[k][2][0]),'%.3f +/- %.4f'%(info[k][1], info[k][3])])                       
                      elif k==1: 
                        x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k-1][1]),'%.5f'%(100*info[k][1]/info[k][2][0]), '%.3f +/- %.4f'%(info[k][1], info[k][3])])
                        x.add_row(['--------------------','--------------','--------------','-------------'])
                      else: 
                        if info[k-1][1]==0 or info[k][2][1]==0: 
                            x.add_row(['%s'%info[k][0],'0.00','0.00000','0.000'])
                        else:
                            x.add_row(['%s'%info[k][0],'%.2f'%(100*info[k][1]/info[k-1][1]),'%.5f'%(100*info[k][1]/info[k][2][1]), '%.3f +/- %.4f'%(info[k][1], info[k][3])])
                      
        print x
        with open('studies/taucsccutflow/%s.txt'%(signal), 'w') as the_file:
           the_file.write(str(x))

def MakeCSCTimeTable(info,signal,isData):
        os.system("mkdir studies")
        os.system("mkdir studies/taucsccutflow")
        x = PrettyTable( ['Time Selections','Exp. Events +/- Uncertainty'])
        x.align["Time Selections"]    = "l"
        x.align["Exp. Events +/- Uncertainty"]   = "r"
        entries = len(info)
        for k in range(0,entries):
            x.add_row(['%s'%info[k][0],'%.5f +/- %.5f'%(info[k][1], info[k][3])])                                  
        print x
        with open('studies/taucsccutflow/%s_time.txt'%(signal), 'w') as the_file:
           the_file.write(str(x))

def CSCLooper(input_dir,output_dir,sample,samplename,isData,isSignal,signorm,tauregion):
    #Open input and output file
    if isData==False:
        t, cf, ac_csc, ac_dt = GetChainSingle(input_dir,'MuonSystem',sample,isData)
    else:
        t, cf, ac_csc, ac_dt = GetChainMultiple(input_dir,'MuonSystem',sample,isData)      
    #Counters
    cut6     = 0
    cut7     = 0
    cut8     = 0
    cut9     = 0
    cut10    = 0
    cut11    = 0
    cut12    = 0
    cut13_itime_n80 = 0
    cut13_itime_p80 = 0
    cut13_otime_n80 = 0 
    cut13_otime_p80 = 0 
    cut14_itime_n80 = 0
    cut14_itime_p80 = 0
    cut14_otime_n80 = 0 
    cut14_otime_p80 = 0 

    #Loop over events
    for e in range(0, t.GetEntries()):
        t.GetEntry(e)

        #acceptance for mc
        if isData==False and isSignal==True and t.gLLP_csc[0]==0 and t.gLLP_csc[1]==0: continue
 
        #CUT (Cluster)
        clstagged,clsid = SelectCSCCluster(t)
        if clstagged==False: continue
        cut6+=1
        plot = False
        if (isData==False) or (isData==True and t.cscRechitClusterSize[clsid]<80): plot=True

        #CUT (Select Tau or Antitau)
        if tauregion==True: 
                tautagged,tauid = SelectTau(t)
                if tautagged==False: continue
                cut7+=1 
                #CUT (Jet Veto)
                if t.cscRechitClusterJetVetoPt[clsid] > 10: continue
                cut8+=1
                #CUT (Muon Veto)
                if t.cscRechitClusterMuonVetoPt[clsid] > 20: continue
                cut9+=1
                #CUT (Chamber Veto)
                if t.cscRechitClusterNRechitChamberPlus11[clsid] > 0 or t.cscRechitClusterNRechitChamberMinus11[clsid] > 0: continue        
                if t.cscRechitCluster_match_MB1Seg_0p4[clsid] >0: continue 
                if t.cscRechitCluster_match_RB1_0p4[clsid] >0: continue 
                cut10+=1
                #Cut (Eta)
                if abs(t.cscRechitClusterEta[clsid]) > 2.2: continue
                cut11+=1
                #CUT(time spread)
                if t.cscRechitClusterTimeSpreadWeightedAll[clsid] > 20: continue
                cut12+=1
                #CUT IN-Time 
                if INTimeCSC(t,clsid)==True:
                        if isData==False: #MC
                            if t.cscRechitClusterSize[clsid]<80: cut13_itime_n80+=1
                            else: cut13_itime_p80+=1 
                            #CUT (dphi)
                            if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                            if t.cscRechitClusterSize[clsid]<80: cut14_itime_n80+=1
                            else: cut14_itime_p80+=1 
                        else: #Data (blinded > 80)
                            if t.cscRechitClusterSize[clsid]<80: cut13_itime_n80+=1
                            #CUT (dphi)
                            if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                            if t.cscRechitClusterSize[clsid]<80: cut14_itime_n80+=1                            
                #CUT OOT
                if OOTimeCSC(t,clsid)==True:
                        if t.cscRechitClusterSize[clsid]<80: cut13_otime_n80+=1
                        else: cut13_otime_p80+=1 
                        #CUT (dphi)
                        if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                        if t.cscRechitClusterSize[clsid]<80: cut14_otime_n80+=1
                        else: cut14_otime_p80+=1 
        else:
                tautagged,tauid = SelectAntiTau(t)
                if tautagged==False: continue
                cut7+=1 
                #CUT (Jet Veto)
                if t.cscRechitClusterJetVetoPt[clsid] > 10: continue
                cut8+=1
                #CUT (Muon Veto)
                if t.cscRechitClusterMuonVetoPt[clsid] > 20: continue
                cut9+=1
                #CUT (Chamber Veto)
                if t.cscRechitClusterNRechitChamberPlus11[clsid] > 0 or t.cscRechitClusterNRechitChamberMinus11[clsid] > 0: continue        
                if t.cscRechitCluster_match_MB1Seg_0p4[clsid] >0: continue 
                if t.cscRechitCluster_match_RB1_0p4[clsid] >0: continue 
                cut10+=1
                #Cut (Eta)
                if abs(t.cscRechitClusterEta[clsid]) > 2.2: continue
                cut11+=1
                #CUT(time spread)
                if t.cscRechitClusterTimeSpreadWeightedAll[clsid] > 20: continue
                cut12+=1
                #CUT IN-Time 
                if INTimeCSC(t,clsid)==True:
                        if t.cscRechitClusterSize[clsid]<80: cut13_itime_n80+=1
                        else: cut13_itime_p80+=1 
                        #CUT (dphi)
                        if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                        if t.cscRechitClusterSize[clsid]<80: cut14_itime_n80+=1
                        else: cut14_itime_p80+=1 
                #CUT OOT
                if OOTimeCSC(t,clsid)==True:
                        if t.cscRechitClusterSize[clsid]<80: cut13_otime_n80+=1
                        else: cut13_otime_p80+=1 
                        #CUT (dphi)
                        if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                        if t.cscRechitClusterSize[clsid]<80: cut14_otime_n80+=1
                        else: cut14_otime_p80+=1 

    #Collect CSC cutflow to make the table
    nums = []
    if isSignal==True:
        w = signorm/cf.GetBinContent(1) 
    else:
        w = 1 
    if isSignal==True:
        nums.append(w*ac_csc.GetBinContent(1))   #Acceptance
        nums.append(w*ac_csc.GetBinContent(2))   #METtrigger+met>200
        nums.append(w*ac_csc.GetBinContent(3))   #MET filters
        nums.append(w*ac_csc.GetBinContent(4))   #nrings
        nums.append(w*ac_csc.GetBinContent(5))    #at least one cluster   
    else:
        nums.append(w*cf.GetBinContent(1))     #Total
        nums.append(w*cf.GetBinContent(2))     #METtrigger+met>200
        nums.append(w*cf.GetBinContent(3))     #MET filters
        nums.append(w*cf.GetBinContent(4))     #nrings
        nums.append(w*cf.GetBinContent(5))     #at least one cluster        
    nums.append(w*float(cut6 ) )  
    nums.append(w*float(cut7 ) ) 
    nums.append(w*float(cut8 ) ) 
    nums.append(w*float(cut9 ) ) 
    nums.append(w*float(cut10) ) 
    nums.append(w*float(cut11) ) 
    nums.append(w*float(cut12) ) 
    if isSignal==True: 
        tag = 'Acceptance'
        den = [w*cf.GetBinContent(1),w*ac_csc.GetBinContent(2)]
    else:
        den = [w*cf.GetBinContent(2),w*cf.GetBinContent(2)]
        tag = 'Total'
    #Fill array     
    effinfo = [
             [tag                        ,nums[0] , den, math.sqrt(w*nums[0])],
             ['MET Trigger and MET>200'  ,nums[1] , den, math.sqrt(w*nums[1])],
             ['MET filters'              ,nums[2] , den, math.sqrt(w*nums[2])], 
             ['CSC+DT rings <= 10'       ,nums[3] , den, math.sqrt(w*nums[3])],
             ['# clusters >= 1'          ,nums[4] , den, math.sqrt(w*nums[4])],
             ['# CSC clusters >= 1'      ,nums[5] , den, math.sqrt(w*nums[5])],
             ['# tauhads >= 1'           ,nums[6] , den, math.sqrt(w*nums[6])], 
             ['Jet Veto'                 ,nums[7] , den, math.sqrt(w*nums[7])], 
             ['Muon Veto'                ,nums[8] , den, math.sqrt(w*nums[8])],
             ['ME1/1, MB1, RB1 Vetos'    ,nums[9] , den, math.sqrt(w*nums[9])],  
             ['|Eta| < 2.2'              ,nums[10], den, math.sqrt(w*nums[10])], 
    ]
    effinfotime = []
    effinfotime.append(['IN-time (rechits<80)' ,w*cut13_itime_n80, den, math.sqrt(w*w*cut13_itime_n80)])
    effinfotime.append(['IN-time (rechits>=80)',w*cut13_itime_p80, den, math.sqrt(w*w*cut13_itime_p80)])
    effinfotime.append(['OO-time (rechits<80)' ,w*cut13_otime_n80, den, math.sqrt(w*w*cut13_otime_n80)])
    effinfotime.append(['OO-time (rechits>=80)',w*cut13_otime_p80, den, math.sqrt(w*w*cut13_otime_p80)])

    return effinfo,effinfotime

###########OPTIONS
parser = argparse.ArgumentParser(description='Command line parser of skim options')
parser.add_argument('--config' ,  dest='cfgfile',  help='Name of config file',  required = True)
parser.add_argument('--tag',    dest='tag',  help='Name of tag', required = True)
parser.add_argument('--tauregion',    dest='tauregion', action='store_true', help='Apply Tau region')
parser.add_argument('--no-tauregion', dest='tauregion', action='store_false',help='Do not apply Tau region')
parser.set_defaults(tauregion=False)

args           = parser.parse_args()
configfilename = args.cfgfile
tagname        = args.tag
tauregion      = args.tauregion

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

#Run looper over MCsignal samples
isData=False
isSignal=True
for k in range (0,len(signals) ):
    print "[INFO] Running looper over %s"%signals[k]
    effinfo,effinfotime = CSCLooper(input_dir+tagname,output_dir,signals[k],signals[k],isData,isSignal,sigxs[k]*lumi,tauregion)
    MakeCSCEfficiencyTable(effinfo,signals[k],isData)
    MakeCSCTimeTable(effinfotime,signals[k],isData)

#Run looper over Data samples
isData=True
isSignal=False
print "[INFO] Running looper over data"
effinfo,effinfotime = CSCLooper(input_dir+tagname,output_dir,datas,'Data',isData,isSignal,1,tauregion)
MakeCSCEfficiencyTable(effinfo,'Data',isData)
MakeCSCTimeTable(effinfotime,'Data',isData)
