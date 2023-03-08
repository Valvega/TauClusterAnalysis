import ROOT
from ROOT import TFile, TTree, TMath, TH1F,TChain,TH2F, TF1, TCanvas, TLegend, TPaveText
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
ROOT.gStyle.SetOptStat(0)

def PlotAntiTauModel(input_dir,output_dir,sel):
    #Opening files
    tfile_data1        = TFile.Open("%s/histos_cscmodel_Data.root"%(input_dir) )
    tfile_data2        = TFile.Open("%s/histos_cscmodel_%s_Data.root"%(input_dir,sel) )
    h_tau_data        =  tfile_data1.Get("h_model_clustersize")
    h_antitau_data    =  tfile_data2.Get("h_model_clustersize")
    #Set format
    h_tau_data.SetMarkerColor(ROOT.kBlack)
    h_tau_data.SetMarkerSize(2)
    h_tau_data.SetMarkerStyle(20)
    h_tau_data.SetLineColor(ROOT.kBlack)
    h_tau_data.SetLineWidth(2)
    h_antitau_data_line = h_antitau_data.Clone()
    h_antitau_data_line.SetLineColor(ROOT.kOrange+2)
    h_antitau_data_line.SetLineWidth(3)
    h_antitau_data.SetLineColor(ROOT.kOrange+1)
    h_antitau_data.SetFillColor(ROOT.kOrange+1)

    #Create Canvas
    c1 = TCanvas("c1", "c1", 1300, 1300)
    c1.SetFrameLineWidth(4)
    c1.SetBottomMargin (0.15)
    c1.SetRightMargin (0.05)
    c1.SetLeftMargin (0.15)
    ROOT.gPad.SetLogy()
    #First just a frame to give the format to the plot
    histoframe = TH1F("","",h_tau_data.GetNbinsX(),h_tau_data.GetXaxis().GetXmin(), h_tau_data.GetXaxis().GetXmax())
    histoframe.GetYaxis().SetTitleSize(0.050)
    histoframe.GetXaxis().SetTitleSize(0.050)
    histoframe.GetYaxis().SetLabelSize(0.05)
    histoframe.GetXaxis().SetLabelSize(0.045)
    histoframe.GetXaxis().SetLabelOffset(0.010)
    histoframe.GetYaxis().SetTitleOffset(1.5)
    histoframe.GetXaxis().SetTitleOffset(1.1)
    histoframe.GetXaxis().SetTitle("Number of CSC RecHits")
    histoframe.GetYaxis().SetTitle("N. U.")
    histoframe.Draw()
    #Find maximumvalue in the frame
    maxs = [ h_tau_data.GetMaximum()/h_tau_data.Integral(), h_antitau_data.GetMaximum()/h_antitau_data.Integral()]
    maxvalue = max(maxs)
    histoframe.SetMaximum(10*maxvalue)
    h_antitau_data.DrawNormalized("E2 SAME")
    h_antitau_data_line.DrawNormalized("HIST SAME")
    h_tau_data.DrawNormalized("P SAME")

    #Draw a legend
    leg_1 = TLegend(0.25,0.70,0.95,0.90)
    leg_1.SetNColumns(1)
    leg_1.SetBorderSize(0)
    leg_1.SetTextSize(0.030)
    leg_1.SetTextFont(42)
    leg_1.SetLineColor(1)
    leg_1.SetLineWidth(10)
    leg_1.SetFillColor(0)
    leg_1.SetFillStyle(0)
    leg_1.Draw()
    leg_1.AddEntry(h_tau_data,"Data (Tau region)", "pe")
    leg_1.AddEntry(h_antitau_data_line,"Bkg model (AntiTau region)", "l")
    leg_1.AddEntry(h_antitau_data,"Bkg model stat. unc. (AntiTau region)", "f")
    #Add CMS labels
    pt1 = TPaveText(0.1463218,0.886316,0.3045977,0.978947,"brNDC")
    pt1.SetBorderSize(0)
    pt1.SetTextAlign(12)
    pt1.SetTextFont(62)
    pt1.SetTextSize(0.05)
    pt1.SetFillColor(0)
    pt1.SetFillStyle(0)
    pt1.AddText("CMS #font[52]{Work in progress}")
    pt1.Draw("SAME")
    #Redrawaxis
    histoframe.Draw("SAME AXIS")
    #Save plot
    c1.Update()
    c1.SaveAs("%s/csc_nrechits_%s.pdf"%(output_dir,sel))
    del c1

def CSCLooper(input_dir,output_dir,sample,samplename,isData,isSignal,norm,sel,tauregion):
    #Open input and output file
    if isData==False:
        t, cf, ac_csc, ac_dt = GetChainSingle(input_dir,'MuonSystem',sample,isData)          
    else:
        t, cf, ac_csc, ac_dt = GetChainMultiple(input_dir,'MuonSystem',sample,isData)
    if tauregion==True:
            tmp  = TFile.Open("%s/histos_cscmodel_%s.root"%(output_dir,samplename),'RECREATE')
    else:
            tmp  = TFile.Open("%s/histos_cscmodel_%s_%s.root"%(output_dir,sel,samplename),'RECREATE')         

    #CUT (in-time)
    h_model_clustersize     = TH1F("h_model_clustersize","h_model_clustersize"        ,6,50,80)
    #CUT (out-of-time)
    h_model_clustersize_oot = TH1F("h_model_clustersize_oot","h_model_clustersize_oot",30,50,200)
    
    #Loop over events
    for e in range(0, t.GetEntries()):
        t.GetEntry(e)

        #acceptance for mc
        if isData==False and isSignal==True and t.gLLP_csc[0]==0 and t.gLLP_csc[1]==0: continue
 
        #CUT (Cluster)
        clstagged,clsid = SelectCSCCluster(t)
        if clstagged==False: continue
        plot = False
        if (isData==False) or (isData==True and t.cscRechitClusterSize[clsid]<80): plot=True

        #CUT (Select Tau or Antitau)
        if tauregion==True: 
                tautagged,tauid = SelectTau(t)
                if tautagged==True:
                    #CUT (Jet Veto)
                    if t.cscRechitClusterJetVetoPt[clsid] > 10: continue
                    #CUT (Muon Veto)
                    if t.cscRechitClusterMuonVetoPt[clsid] > 20: continue
                    #CUT (Chamber Veto)
                    if t.cscRechitClusterNRechitChamberPlus11[clsid] > 0 or t.cscRechitClusterNRechitChamberMinus11[clsid] > 0: continue        
                    if t.cscRechitCluster_match_MB1Seg_0p4[clsid] >0: continue 
                    if t.cscRechitCluster_match_RB1_0p4[clsid] >0: continue 
                    #Cut (Eta)
                    if abs(t.cscRechitClusterEta[clsid]) > 2.2: continue
                    #CUT(time spread)
                    if t.cscRechitClusterTimeSpread[clsid] > 20: continue
                    #CUT IN-Time 
                    if INTimeCSC(t,clsid)==True:
                            #CUT (dphi)
                            #if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                            #Tau RecHits distribution (blinded)
                            if plot==True:
                               h_model_clustersize.Fill(t.cscRechitClusterSize[clsid])
                    if OOTimeCSC(t,clsid)==True:
                            #CUT (dphi)
                            #if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                            #Tau RecHits distribution (blinded)
                            h_model_clustersize_oot.Fill(t.cscRechitClusterSize[clsid])           
        else:
                antitautagged,antitauid = SelectAntiTau(t)
                if antitautagged==True:
                        #CUT (Jet Veto)
                        if t.cscRechitClusterJetVetoPt[clsid] > 10: continue
                        #CUT (Muon Veto)
                        if t.cscRechitClusterMuonVetoPt[clsid] > 20: continue
                        #CUT (Chamber Veto)
                        if t.cscRechitClusterNRechitChamberPlus11[clsid] > 0 or t.cscRechitClusterNRechitChamberMinus11[clsid] > 0: continue        
                        if t.cscRechitCluster_match_MB1Seg_0p4[clsid] >0: continue 
                        if t.cscRechitCluster_match_RB1_0p4[clsid] >0: continue 
                        #Cut (Eta)
                        if abs(t.cscRechitClusterEta[clsid]) > 2.2: continue
                        #CUT(time spread)
                        if t.cscRechitClusterTimeSpread[clsid] > 20: continue
                        #CUT IN-Time 
                        if INTimeCSC(t,clsid)==True:
                                ##CUT (dphi)
                                #if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                                #Anti-tau RecHits distribution
                                h_model_clustersize.Fill(t.cscRechitClusterSize[clsid])
                        if OOTimeCSC(t,clsid)==True:
                               ##CUT (dphi)
                               #if abs(t.cscRechitClusterMet_dPhi[clsid])> (math.pi/2): continue
                               #Anti-tau RecHits distribution
                               h_model_clustersize_oot.Fill(t.cscRechitClusterSize[clsid])       
    #weights
    if isSignal==True:
      w = norm/cf.GetBinContent(1) 
    else:
      w = 1 
    h_model_clustersize.Scale(w)
    tmp.Write()
    tmp.Close()

###########OPTIONS
parser = argparse.ArgumentParser(description='Command line parser of skim options')
parser.add_argument('--config' ,  dest='cfgfile',  help='Name of config file',  required = True)
parser.add_argument('--tag',    dest='tag',  help='Name of tag', required = True)
parser.add_argument('--sel',    dest='sel',  help='Name of anti selection', required = False)
parser.add_argument('--tauregion',    dest='tauregion', action='store_true', help='Apply Tau region')
parser.add_argument('--no-tauregion', dest='tauregion', action='store_false',help='Do not apply Tau region')
parser.add_argument('--draw',    dest='draw', action='store_true', help='Draw comparion')
parser.set_defaults(tauregion=False)
parser.set_defaults(draw=False)

args           = parser.parse_args()
configfilename = args.cfgfile
tagname        = args.tag
selname        = args.sel
tauregion      = args.tauregion
draw           = args.draw

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

##Run looper over MCsignal samples to make histograms
#isData=False
#isSignal=True
#for k in range (0,len(signals) ):
#    print "[INFO] Running looper over %s"%signals[k]
#    DTLooper(input_dir+tagname,output_dir,signals[k],signals[k],isData,isSignal,sigxs[k]*lumi,selname,tauregion)

#Run looper over Data samples to make histograms
isData=True
isSignal=False
print "[INFO] Running looper over data"
CSCLooper(input_dir+tagname,output_dir,datas,'Data',isData,isSignal,1,selname,tauregion)

if draw==True:
  print "[INFO] Plotting"
  os.system("mkdir plots")
  os.system("mkdir plots/bkgmodel")
  PlotAntiTauModel(output_dir,"plots/bkgmodel",selname)
