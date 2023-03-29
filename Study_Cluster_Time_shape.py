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
parser.add_argument('--type',	dest='type', help='csc or dt', required = True)
#parser.add_argument('--tauregion',    dest='tauregion', action='store_true', help='Apply Tau region')
#parser.add_argument('--no-tauregion', dest='tauregion', action='store_false',help='Do not apply Tau region')
#parser.add_argument('--draw',    dest='draw', action='store_true', help='Draw comparion')
#parser.set_defaults(tauregion=False)
#parser.set_defaults(draw=False)

args           = parser.parse_args()
configfilename = args.cfgfile
tagname        = args.tag
selname        = args.sel
typename	   = args.type
#tauregion      = args.tauregion
#draw           = args.draw

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

def Histograms(histo1,histo2,regions,output_dir,typename):
	#Set format
	histo1.SetMarkerColor(ROOT.kBlack)
	histo1.SetMarkerSize(2)
	histo1.SetMarkerStyle(20)
	histo1.SetLineColor(ROOT.kBlack)
	histo1.SetLineWidth(2)
	histo2_line = histo2.Clone()
	histo2_line.SetLineColor(ROOT.kOrange+2)
	histo2_line.SetLineWidth(3)
	histo2.SetLineColor(ROOT.kOrange-9)
	histo2.SetFillColor(ROOT.kOrange-9)
	#Create Canvas
	c1 = TCanvas("c1", "c1", 1300, 1300)
	c1.SetFrameLineWidth(4)
	c1.SetBottomMargin (0.15)
	c1.SetRightMargin (0.05)
	c1.SetLeftMargin (0.15)
	ROOT.gPad.SetLogy()
	#First just a frame to give the format to the plot
	histoframe = TH1F("","", histo1.GetNbinsX(), histo1.GetXaxis().GetXmin(), histo1.GetXaxis().GetXmax())
	histoframe.GetYaxis().SetTitleSize(0.050)
	histoframe.GetXaxis().SetTitleSize(0.050)
	histoframe.GetYaxis().SetLabelSize(0.05)
	histoframe.GetXaxis().SetLabelSize(0.045)
	histoframe.GetXaxis().SetLabelOffset(0.010)
	histoframe.GetYaxis().SetTitleOffset(1.5)
	histoframe.GetXaxis().SetTitleOffset(1.1)
	histoframe.GetXaxis().SetTitle("Number of %s RecHits"%(typename))
	histoframe.GetYaxis().SetTitle("N. U.")
	histoframe.Draw()
	#Find maximumvalue in the frame
	maxs = [ histo1.GetMaximum()/histo1.Integral(), histo2.GetMaximum()/histo2.Integral()]
	maxvalue = max(maxs)
	histoframe.SetMaximum(10*maxvalue)
	histo2.DrawNormalized("E2 SAME")
	histo2_line.DrawNormalized("HIST SAME")
	histo1.DrawNormalized("P SAME")	
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
	leg_1.AddEntry(histo1,"In Time, Data (%s region)"%(regions), "pe")
	leg_1.AddEntry(histo2_line,"Out of Time, Data (%s region)"%(regions), "l")
	leg_1.AddEntry(histo2,"Out of Time, Data unc.(%s region)"%(regions), "f")	
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
	c1.SaveAs("%s/%s_nrechits_%s.pdf"%(output_dir,typename,regions))
	del c1

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
Histograms( h1, h7, "Tau",output_dir,typename)
Histograms( h5, h11,"0VVLID",output_dir,typename)
Histograms( h6, h12,"0VVVLID",output_dir,typename)
