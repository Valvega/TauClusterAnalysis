import ROOT
from ROOT import TFile, TTree, TMath, TChain, TH2F, TH1F, TCanvas,TPaveText,TLegend
import sys
import os
import argparse
import ast
import glob
import math
from  ConfigParser import *
import math
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

def Plot1DMultipleFiles(inputdirectory,outputdirectory,infiles,histo):    
        #Loop over booked histos
        allhistos = []
        for file  in infiles:
              #Get histograms 
              thefile = TFile.Open(inputdirectory+file[0])
              tmp  = thefile.Get(histo[0])
              htmp = tmp.Clone()
              htmp.SetDirectory(0)
              #Add format 
              htmp.SetMarkerColor(file[2])
              htmp.SetMarkerSize(2)
              htmp.SetMarkerStyle(20)
              htmp.SetLineColor(file[2])
              htmp.SetLineWidth(4)
              #Add histos for plotting  
              allhistos.append(htmp)    
              thefile.Close()
        #Create Canvas
        c1 = TCanvas("c1", "c1", 1300, 1300)
        c1.SetFrameLineWidth(4)
        c1.SetBottomMargin (0.15)
        c1.SetRightMargin (0.05)
        c1.SetLeftMargin (0.15)
        #First just a frame to give the format to the plot
        histoframe = TH1F("","",allhistos[0].GetNbinsX(),allhistos[0].GetXaxis().GetXmin(), allhistos[0].GetXaxis().GetXmax())
        histoframe.GetYaxis().SetTitleSize(0.050)
        histoframe.GetXaxis().SetTitleSize(0.050)
        histoframe.GetYaxis().SetLabelSize(0.05)
        histoframe.GetXaxis().SetLabelSize(0.045)
        histoframe.GetXaxis().SetLabelOffset(0.010)
        histoframe.GetYaxis().SetTitleOffset(1.5)
        histoframe.GetXaxis().SetTitleOffset(1.1)
        histoframe.GetXaxis().SetTitle(histo[1])
        histoframe.GetYaxis().SetTitle(histo[2])
        histoframe.Draw()
        maxs = []
        for j in range(0,len(allhistos)):  
          maxs.append(allhistos[j].GetMaximum()/allhistos[j].Integral())
        maxvalue = max(maxs)
        histoframe.SetMaximum(1.4*maxvalue)
        #Draw all histograms 
        for k in range(0,len(allhistos)):  allhistos[k].DrawNormalized("Histo SAME")
        #Draw a legend
        leg_1 = TLegend(0.15,0.75,0.85,0.90)
        leg_1.SetNColumns(1)
        leg_1.SetBorderSize(0)
        leg_1.SetTextSize(0.040)
        leg_1.SetTextFont(42)
        leg_1.SetLineColor(1)
        leg_1.SetLineWidth(10)
        leg_1.SetFillColor(0)
        leg_1.SetFillStyle(0)
        leg_1.Draw()
        for k in range(0,len(allhistos)): leg_1.AddEntry(allhistos[k],infiles[k][3], "l")
        #Add CMS labels
        pt1 = TPaveText(0.1463218,0.886316,0.3045977,0.978947,"brNDC")
        pt1.SetBorderSize(0)
        pt1.SetTextAlign(12)
        pt1.SetTextFont(62)
        pt1.SetTextSize(0.05)
        pt1.SetFillColor(0)
        pt1.SetFillStyle(0)
        #pt1.AddText("CMS #font[52]{Work in progress}")
        pt1.Draw("SAME")
        #Redrawaxis
        histoframe.Draw("SAME AXIS")
        #Save plot
        c1.Update()
        c1.SaveAs("%s/%s.pdf"%(outputdirectory,histo[3]))
        del c1



###########OPTIONS
parser = argparse.ArgumentParser(description='Command line parser of skim options')
parser.add_argument('--config' ,       dest='cfgfile'   ,  help='Name of config file'   , required = True)
parser.add_argument('--tag',           dest='tag'       ,  help='Name of output folder' , required = True)
args           = parser.parse_args()
configfilename = args.cfgfile
tagname        = args.tag
#Reading configuration
print "[INFO] Reading configuration file . . ."
cfgparser = ConfigParser()
cfgparser.read('%s'%configfilename)
print "[INFO] Getting configuration parameters . . ."
year         = ast.literal_eval(cfgparser.get("plotconfiguration","year"))
print "    -The year:"
print "      *",year
input_dir    = ast.literal_eval(cfgparser.get("plotconfiguration","input_dir"))
print "    -The input  directory:"
print "      *",input_dir+tagname+"/%i"%year
output_dir   = ast.literal_eval(cfgparser.get("plotconfiguration","output_dir"))
print "    -The output directory:"
print "      *",output_dir+tagname+"/%i"%year



#Files containing the histogram to plot ['Name of file','Legend','Color'] 
files_multiplefiles = [
            ['histos_VLLPair_VLLToTauS_MVLL300_MS10_ctau300.root' ,'VLLPair_VLLToTauS_MVLL300_MS10_ctau300' ,ROOT.kGreen+2 ,"m_{#tau^{/}}= 300 GeV, m_{A}=10 GeV, c#tau=0.3 m" ],
            ['histos_VLLPair_VLLToTauS_MVLL700_MS10_ctau100.root' ,'VLLPair_VLLToTauS_MVLL700_MS10_ctau100' ,ROOT.kRed+2   ,"m_{#tau^{/}}= 700 GeV, m_{A}=10 GeV, c#tau=0.1 m" ],
            ['histos_VLLPair_VLLToTauS_MVLL1000_MS10_ctau100.root','VLLPair_VLLToTauS_MVLL1000_MS10_ctau100',ROOT.kBlue+2  ,"m_{#tau^{/}}=1000 GeV, m_{A}=10 GeV, c#tau=0.1 m" ],
]

#Histograms to plot ['Histogram name in file','X axis label','Y axis label','file output name'] 
histos1d_multiplefiles = [
    ['h_vll_pt'    ,'Gen #tau^{/} p_{T} [GeV]'  ,'N.U.','plot_mult_vll_pt'],
    ['h_vll_eta'   ,'Gen #tau^{/} #eta'         ,'N.U.','plot_mult_vll_eta'],
    ['h_vll_e'     ,'Gen #tau^{/} Energy [GeV]' ,'N.U.','plot_mult_vll_e'],
    ['h_vll_phi'   ,'Gen #tau^{/} #phi'         ,'N.U.','plot_mult_vll_phi'],
    ['h_tau_pt'    ,'Gen #tau p_{T} [GeV]'      ,'N.U.','plot_mult_tau_pt'],
    ['h_tau_eta'   ,'Gen #tau #eta'             ,'N.U.','plot_mult_tau_eta'],
    ['h_tau_e'     ,'Gen #tau Energy [GeV]'     ,'N.U.','plot_mult_tau_e'],
    ['h_tau_phi'   ,'Gen #tau #phi'             ,'N.U.','plot_mult_tau_phi'],
    ['h_llp_pt'    ,'Gen A p_{T} [GeV]'         ,'N.U.','plot_mult_llp_pt'],
    ['h_llp_eta'   ,'Gen A #eta'                ,'N.U.','plot_mult_llp_eta'],
    ['h_llp_e'     ,'Gen A Energy [GeV]'        ,'N.U.','plot_mult_llp_e'],
    ['h_llp_phi'   ,'Gen A #phi'                ,'N.U.','plot_mult_llp_phi'],
] 

##Plot 1d distributions
os.system("mkdir plots/%s"%tagname)
os.system("mkdir plots/%s/%i"%(tagname,year))
output_dir =  output_dir+tagname+"/%i"%year
input_dir  =  input_dir+tagname+"/%i/"%year

#loop over the 1d histograms you want to compare
print "[INFO] Plotting distributions from multiple files . . ."
for histo in histos1d_multiplefiles:
   Plot1DMultipleFiles(input_dir,output_dir,files_multiplefiles,histo)
