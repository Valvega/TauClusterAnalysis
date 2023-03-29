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

#Opening files
tfile_data1  = TFile.Open("%s/histos_cscmodel_Data.root"%(input_dir))
tfile_data2  = TFile.Open("%s/histos_cscmodel_0M_Data.root"%(input_dir))
tfile_data3  = TFile.Open("%s/histos_cscmodel_0L_Data.root"%(input_dir))
tfile_data4  = TFile.Open("%s/histos_cscmodel_0VL_Data.root"%(input_dir))
tfile_data5  = TFile.Open("%s/histos_cscmodel_0VVL_Data.root"%(input_dir))
tfile_data6  = TFile.Open("%s/histos_cscmodel_0VVVL_Data.root"%(input_dir))
##Out of time
h1           =  tfile_data1.Get("h_model_clustersize")
h2           =  tfile_data2.Get("h_model_clustersize")
h3           =  tfile_data3.Get("h_model_clustersize")
h4           =  tfile_data4.Get("h_model_clustersize")
h5           =  tfile_data5.Get("h_model_clustersize")
h6           =  tfile_data6.Get("h_model_clustersize")		
#In time
h7           =  tfile_data1.Get("h_model_clustersize_oot")
h8           =  tfile_data2.Get("h_model_clustersize_oot")
h9           =  tfile_data3.Get("h_model_clustersize_oot")
h10          =  tfile_data4.Get("h_model_clustersize_oot")
h11          =  tfile_data5.Get("h_model_clustersize_oot")
h12          =  tfile_data6.Get("h_model_clustersize_oot")

testp1 = h1.KolmogorovTest(h7)
testp2 = h2.KolmogorovTest(h8)
testp3 = h3.KolmogorovTest(h9)
testp4 = h4.KolmogorovTest(h10)
testp5 = h5.KolmogorovTest(h11)
testp6 = h6.KolmogorovTest(h12)

print "1MID difference between IT and OOT",testp1
print "0MID difference between IT and OOT",testp2
print "0LID difference between IT and OOT",testp3  
print "0VLID difference between IT and OOT",testp4 
print "0VVLID difference between IT and OOT",testp5
print "0VVVLID difference between IT and OOT",testp6

##############Histograms############################
#Set format
#Tau region
h1.SetMarkerColor(ROOT.kBlack)
h1.SetMarkerSize(2)
h1.SetMarkerStyle(20)
h1.SetLineColor(ROOT.kBlack)
h1.SetLineWidth(2)
h7_line = h7.Clone()
h7_line.SetLineColor(ROOT.kGray+2)
h7_line.SetLineWidth(3)
h7.SetLineColor(ROOT.kGray+1)
h7.SetFillColor(ROOT.kGray+1)
##0VVLID region
#h5.SetMarkerColor(ROOT.kRed+1)
#h5.SetMarkerSize(2)
#h5.SetMarkerStyle(20)
#h5.SetLineColor(ROOT.kRed+1)
#h5.SetLineWidth(2)
#h11_line = h11.Clone()
#h11_line.SetLineColor(ROOT.kRed+2)
#h11_line.SetLineWidth(3)
#h11.SetLineColor(ROOT.kRed-9)
#h11.SetFillColor(ROOT.kRed-9)
##0VVVLID region
#h6.SetMarkerColor(ROOT.kOrange+7)
#h6.SetMarkerSize(2)
#h6.SetMarkerStyle(20)
#h6.SetLineColor(ROOT.kOrange+7)
#h6.SetLineWidth(2)
#h12_line = h12.Clone()
#h12_line.SetLineColor(ROOT.kOrange+2)
#h12_line.SetLineWidth(3)
#h12.SetLineColor(ROOT.kOrange-9)
#h12.SetFillColor(ROOT.kOrange-9)

#Create Canvas
c1 = TCanvas("c1", "c1", 1300, 1300)
c1.SetFrameLineWidth(4)
c1.SetBottomMargin (0.15)
c1.SetRightMargin (0.05)
c1.SetLeftMargin (0.15)
ROOT.gPad.SetLogy()

#First just a frame to give the format to the plot
histoframe = TH1F("","", h1.GetNbinsX(), h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax())
#histoframe = TH1F("","", h5.GetNbinsX(), h5.GetXaxis().GetXmin(), h5.GetXaxis().GetXmax())
#histoframe = TH1F("","", h6.GetNbinsX(), h6.GetXaxis().GetXmin(), h6.GetXaxis().GetXmax())
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
maxs = [ h1.GetMaximum()/h1.Integral(), h7.GetMaximum()/h7.Integral()]
#maxs = [ h5.GetMaximum()/h5.Integral(), h11.GetMaximum()/h11.Integral()]
#maxs = [ h6.GetMaximum()/h6.Integral(), h12.GetMaximum()/h12.Integral()]
maxvalue = max(maxs)
histoframe.SetMaximum(10*maxvalue)
h7.DrawNormalized("E2 SAME")
h7_line.DrawNormalized("HIST SAME")
h1.DrawNormalized("P SAME")
#h11.DrawNormalized("E2 SAME")
#h11_line.DrawNormalized("HIST SAME")
#h5.DrawNormalized("P SAME")
#h12.DrawNormalized("E2 SAME")
#h12_line.DrawNormalized("HIST SAME")
#h6.DrawNormalized("P SAME")

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
leg_1.AddEntry(h1,"In Time, Data (Tau region)", "pe")
leg_1.AddEntry(h7_line,"Out of Time, Data (Tau region)", "l")
leg_1.AddEntry(h7,"Out of Time, Data unc. (Tau region)", "f")
#leg_1.AddEntry(h5,"In Time, Data (0VVLID region)", "pe")
#leg_1.AddEntry(h11_line,"Out of Time, Data (0VVLID region)", "l")
#leg_1.AddEntry(h11,"Out of Time, Data unc. (0VVLID region)", "f")
#leg_1.AddEntry(h6,"In Time, Data (0VVVLID region)", "pe")
#leg_1.AddEntry(h12_line,"Out of Time, Data (0VVLID region)", "l")
#leg_1.AddEntry(h12,"Out of Time, Data unc. (0VVVLID region)", "f")

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
c1.SaveAs("%s/csc_nrechits_Tau.pdf"%(output_dir))
#c1.SaveAs("%s/csc_nrechits_0VVL.pdf"%(output_dir))
#c1.SaveAs("%s/csc_nrechits_0VVVL.pdf"%(output_dir))

##############TABLES################################
####IN Time CSC
#a = PrettyTable(display_size=(3,8))
#a.field_names = ["Selection AntiTau", "P. KS test","S300(rechits<80)","S300(rechits>=80)", "S300/B(rechits<80)", "S300/B(rechits>=80)"]
#a.align["Selection AntiTau"] = "l"
#a.align["P. KS test"]        = "r"
#a.align["S300(rechits<80)"]  = "r"
#a.align["S300(rechits>=80)"]  = "r"
#a.align["S/B(rechits<80)"]   = "r"
#a.align["S/B(rechits>80)"]   = "r"
#a.add_row(["0 Medium"  ,"0.999","0.256 +/- 0.044","3.130 +/- 0.154","0.001","0.029"])
#a.add_row(["0 Loose"   ,"0.995","0.204 +/- 0.039","2.376 +/- 0.134","0.001","0.024"])
#a.add_row(["0 VLoose"  ,"0.963","0.173 +/- 0.036","1.621 +/- 0.111","0.001","0.018"])
#a.add_row(["0 VVLoose" ,"0.927","0.151 +/- 0.034","1.184 +/- 0.095","0.000","0.014"])
#a.add_row(["0 VVVLoose","0.977","0.128 +/- 0.031","0.965 +/- 0.085","0.001","0.016"])
#print(a.get_string())
####IN Time DT
#b = PrettyTable(display_size=(3,8))
#b.field_names = ["Selection AntiTau", "P. KS test","S300(rechits<80)","S300(rechits>=80)", "S300/B(rechits<80)", "S300/B(rechits>80)"]
#b.align["Selection AntiTau"] = "l"
#b.align["P. KS test"]        = "r"
#b.align["S300(rechits<80)"]  = "r"
#b.align["S300(rechits>=80)"]  = "r"
#b.align["S/B(rechits<80)"]   = "r"
#b.align["S/B(rechits>80)"]   = "r"
#b.add_row(["0 Medium"  ,"0.966","0.822 +/- 0.079","6.953 +/- 0.229","0.000","0.044"])
#b.add_row(["0 Loose"   ,"0.969","0.528 +/- 0.063","5.332 +/- 0.201","0.000","0.035"])
#b.add_row(["0 VLoose"  ,"0.963","0.332 +/- 0.050","3.711 +/- 0.167","0.000","0.025"])
#b.add_row(["0 VVLoose" ,"0.955","0.249 +/- 0.043","2.888 +/- 0.148","0.000","0.022"])
#b.add_row(["0 VVVLoose","0.933","0.204 +/- 0.039","2.330 +/- 0.133","0.000","0.025"])
#print(b.get_string())
####Out Of Time CSC
#c = PrettyTable(display_size=(3,8))
#c.field_names = ["Selection AntiTau", "P. KS test","S300(rechits<80)","S300(rechits>=80)", "S300/B(rechits<80)", "S300/B(rechits>80)"]
#c.align["Selection AntiTau"] = "l"
#c.align["P. KS test"]        = "r"
#c.align["S300(rechits<80)"]  = "r"
#c.align["S300(rechits>=80)"]  = "r"
#c.align["S/B(rechits<80)"]   = "r"
#c.align["S/B(rechits>80)"]   = "r"
#c.add_row(["0 Medium"  ,"0.880","0.023 +/- 0.013","0.000 +/- 0.000","0.000","0.000"])
#c.add_row(["0 Loose"   ,"0.869","0.015 +/- 0.011","0.000 +/- 0.000","0.000","0.000"])
#c.add_row(["0 VLoose"  ,"0.874","0.015 +/- 0.011","0.000 +/- 0.000","0.000","0.000"])
#c.add_row(["0 VVLoose" ,"0.871","0.008 +/- 0.008","0.000 +/- 0.000","0.000","0.000"])
#c.add_row(["0 VVVLoose","0.694","0.000 +/- 0.000","0.000 +/- 0.000","0.000","0.000"])
#print(c.get_string())
####Out Of Time DT
#d = PrettyTable(display_size=(3,8))
#d.field_names = ["Selection AntiTau", "P. KS test","S300(rechits<80)","S300(rechits>=80)", "S300/B(rechits<80)", "S300/B(rechits>80)"]
#d.align["Selection AntiTau"] = "l"
#d.align["P. KS test"]        = "r"
#d.align["S300(rechits<80)"]  = "r"
#d.align["S300(rechits>=80)"]  = "r"
#d.align["S/B(rechits<80)"]   = "r"
#d.align["S/B(rechits>80)"]   = "r"
#d.add_row(["0 Medium"  ,"0.999","0.008 +/- 0.008","0.008 +/- 0.008","0.000","0.000"])
#d.add_row(["0 Loose"   ,"0.999","0.008 +/- 0.008","0.008 +/- 0.008","0.000","0.000"])
#d.add_row(["0 VLoose"  ,"0.999","0.008 +/- 0.008","0.008 +/- 0.008","0.000","0.000"])
#d.add_row(["0 VVLoose" ,"0.999","0.008 +/- 0.008","0.008 +/- 0.008","0.000","0.000"])
#d.add_row(["0 VVVLoose","0.999","0.008 +/- 0.008","0.000 +/- 0.000","0.000","0.000"])
#print(d.get_string())
