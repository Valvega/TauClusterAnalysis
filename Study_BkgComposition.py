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
from prettytable import PrettyTable

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
    t,cutflow  = GetChain(input_dir,'tree',signal) #nombre del tree #MuonSystem
    tmp     = TFile.Open("%s/histos_%s.root"%(output_dir,signal),'RECREATE')
    Nevents = t.GetEntries()
    if maxnevents>0 and maxnevents<=Nevents: Nevents = maxnevents

    #Book Histograms
    #Reco Tau
    h_pt            = TH1F("h_pt",           "h_pt"           ,40,0,800)
    h_eta           = TH1F("h_eta",          "h_eta"          ,20,-4.0,4.0)
    h_changerdiso   = TH1F("h_changerdiso",  "h_changerdiso"  ,20,0,8)
    h_neutraliso    = TH1F("h_neutraliso",   "h_neutraliso"   ,20,0,8)
    h_met           = TH1F("h_met",          "h_met"          ,30,0,1000)
    h_Dphi_tau_met  = TH1F("h_Dphi_tau_met", "h_Dphi_tau_met" ,20,0,4.0)
    h_tauDeltaR     = TH1F("h_tauDeltaR",    "h_tauDeltaR"    ,40,0,4)

    #Variables
    crosssection  = xs
    luminosity    = lumi
    totalMCevents = cutflow.GetBinContent(1) #indice
    trigger_met200= cutflow.GetBinContent(2)
    met_filters   = cutflow.GetBinContent(3)
    nrings        = cutflow.GetBinContent(4)
    cluster       = cutflow.GetBinContent(5)

    #counter
    tau_counter = 0
    antitau_counter = 0
    #Loop over events
    for e in range(0, Nevents):

        #Get entries from tree
        t.GetEntry(e)
        
        #Tau region
        tauexist = False
        large = -99
        idx   = -99
        for k in range(0, t.nTau):
            if t.tauPt[k]>50 and abs(t.tauEta[k])<2.3 and t.tau_IsMedium[k]==1: 
                if t.tauPt[k]>large:
                    large=t.tauPt[k]
                    idx=k
                tauexist = True 
        if tauexist == True: tau_counter+=1
        
#        h_pt.Fill(           t.tauPt[idx],                            (crosssection*luminosity)/totalMCevents )
#        h_eta.Fill(          t.tauEta[idx],                           (crosssection*luminosity)/totalMCevents )
#        h_changerdiso.Fill(  t.tauChargedIso[idx],                    (crosssection*luminosity)/totalMCevents )
#        h_neutraliso.Fill(   t.tauNeutralIso[idx],                    (crosssection*luminosity)/totalMCevents )
#        h_met.Fill(          t.MET,                                   (crosssection*luminosity)/totalMCevents )
#        h_Dphi_tau_met.Fill( abs(DeltaPhi(t.tauPhi[idx], t.MET_Phi)), (crosssection*luminosity)/totalMCevents )
#        h_tauDeltaR.Fill(    t.tauDeltaR[idx],                        (crosssection*luminosity)/totalMCevents )

        #AntiTau region
        antitauexist = False
        antitau_large = -99
        antitau_idx   = -99
        for k in range(0, t.nTau):
            if t.tauPt[k]>50 and abs(t.tauEta[k])<2.3 and t.tau_IsVVLoose[k]==0: #change the antitau region
                if t.tauPt[k]>antitau_large:
                    antitau_large=t.tauPt[k]
                    antitau_idx=k
                antitauexist = True 
        if antitauexist == True: antitau_counter+=1

    tmp.Write()
    tmp.Close()
        
    v = []
    v.append(crosssection*luminosity) #No cuts
    v.append((crosssection*luminosity)/(totalMCevents**0.5)) #Uncertainty of No cuts
    v.append((trigger_met200/totalMCevents)*crosssection*luminosity) #MET trigger + MET200
    v.append((crosssection*luminosity*(trigger_met200**0.5))/totalMCevents ) #Uncertainty of MET trigger + MET200
    v.append(crosssection*luminosity*(tau_counter/totalMCevents)) #1 tau MediumID
    v.append((crosssection*luminosity*(tau_counter**0.5))/totalMCevents) #Uncertainty of 1 tau MID
    v.append(totalMCevents) #Data no cuts
    v.append(totalMCevents**0.5) #Data unc. no cuts
    v.append(trigger_met200) #Data Mets
    v.append(trigger_met200**0.5) #Data unc. mets
    v.append(tau_counter) #Data 1 tau MID
    v.append(tau_counter**0.5) #Data unc. 1 tau MID
    v.append(crosssection*luminosity*(antitau_counter/totalMCevents)) #0 tau VVLooseID
    v.append((crosssection*luminosity*(antitau_counter**0.5))/totalMCevents) #Uncertainty of 0 tau VVLooseID
    v.append(antitau_counter) #Data 0 tau VVLID
    v.append(antitau_counter**0.5) #Data unc. 0 tau VVLID
    return v

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

matrix = []

#Run looper over signal samples
for k in range(0,len(bkgs)):
   # print "[INFO] Running looper over %s"%bkgs[k]
    vector = Looper(input_dir,output_dir,bkgs[k],bkgxs[k],lumi,maxnevents)
    matrix.append(vector)
   # print vector

m =[]
w = []
#No cuts
w.append(matrix[0][0] + matrix[1][0] + matrix[2][0] + matrix[3][0] + matrix[4][0] + matrix[5][0] + matrix[6][0] + matrix[7][0] + matrix[8][0]) #QCD
w.append(matrix[9][0])#TT
w.append(matrix[10][0] + matrix[11][0] + matrix[12][0] + matrix[13][0] + matrix[14][0] + matrix[15][0] + matrix[16][0] + matrix[17][0]) #W+js
w.append(matrix[18][0] + matrix[19][0] + matrix[20][0] + matrix[21][0] + matrix[22][0] + matrix[23][0] + matrix[24][0]) #Z+js
#METs
w.append(matrix[0][2] + matrix[1][2] + matrix[2][2] + matrix[3][2] + matrix[4][2] + matrix[5][2] + matrix[6][2] + matrix[7][2] + matrix[8][2]) #QCD
w.append(matrix[9][2])#TT
w.append(matrix[10][2] + matrix[11][2] + matrix[12][2] + matrix[13][2] + matrix[14][2] + matrix[15][2] + matrix[16][2] + matrix[17][2]) #W+js
w.append(matrix[18][2] + matrix[19][2] + matrix[20][2] + matrix[21][2] + matrix[22][2] + matrix[23][2] + matrix[24][2]) #Z+js
#1 Tau MID
w.append(matrix[0][4] + matrix[1][4] + matrix[2][4] + matrix[3][4] + matrix[4][4] + matrix[5][4] + matrix[6][4] + matrix[7][4] + matrix[8][4]) #QCD
w.append(matrix[9][4])#TT
w.append(matrix[10][4] + matrix[11][4] + matrix[12][4] + matrix[13][4] + matrix[14][4] + matrix[15][4] + matrix[16][4] + matrix[17][4]) #W+js
w.append(matrix[18][4] + matrix[19][4] + matrix[20][4] + matrix[21][4] + matrix[22][4] + matrix[23][4] + matrix[24][4]) #Z+js
#0 Tau VVLID
w.append(matrix[0][12] + matrix[1][12] + matrix[2][12] + matrix[3][12] + matrix[4][12] + matrix[5][12] + matrix[6][12] + matrix[7][12] + matrix[8][12]) #QCD
w.append(matrix[9][12])#TT
w.append(matrix[10][12] + matrix[11][12] + matrix[12][12] + matrix[13][12] + matrix[14][12] + matrix[15][12] + matrix[16][12] + matrix[17][12]) #W+js
w.append(matrix[18][12] + matrix[19][12] + matrix[20][12] + matrix[21][12] + matrix[22][12] + matrix[23][12] + matrix[24][12]) #Z+js

j = []
#Uncert. No cuts
j.append((matrix[0][1]**2 + matrix[1][1]**2 + matrix[2][1]**2 + matrix[3][1]**2 + matrix[4][1]**2 + matrix[5][1]**2 + matrix[6][1]**2 + matrix[7][1]**2 + matrix[8][1]**2)**0.5) #QCD
j.append(matrix[9][1])#TT
j.append((matrix[10][1]**2 + matrix[11][1]**2 + matrix[12][1]**2 + matrix[13][1]**2 + matrix[14][1]**2 + matrix[15][1]**2 + matrix[16][1]**2 + matrix[17][1]**2)**0.5) #W+js
j.append((matrix[18][1]**2 + matrix[19][1]**2 + matrix[20][1]**2 + matrix[21][1]**2 + matrix[22][1]**2 + matrix[23][1]**2 + matrix[24][1]**2)**0.5) #Z+js
#Uncert. METs
j.append((matrix[0][3]**2 + matrix[1][3]**2 + matrix[2][3]**2 + matrix[3][3]**2 + matrix[4][3]**2 + matrix[5][3]**2 + matrix[6][3]**2 + matrix[7][3]**2 + matrix[8][3]**2)**0.5) #QCD
j.append(matrix[9][3])#TT
j.append((matrix[10][3]**2 + matrix[11][3]**2 + matrix[12][3]**2 + matrix[13][3]**2 + matrix[14][3]**2 + matrix[15][3]**2 + matrix[16][3]**2 + matrix[17][3]**2)**0.5) #W+js
j.append((matrix[18][3]**2 + matrix[19][3]**2 + matrix[20][3]**2 + matrix[21][3]**2 + matrix[22][3]**2 + matrix[23][3]**2 + matrix[24][3]**2)**0.5) #Z+js
#Uncert. 1 Tau MID
j.append((matrix[0][5]**2 + matrix[1][5]**2 + matrix[2][5]**2 + matrix[3][5]**2 + matrix[4][5]**2 + matrix[5][5]**2 + matrix[6][5]**2 + matrix[7][5]**2 + matrix[8][5]**2)**0.5) #QCD
j.append(matrix[9][5])#TT
j.append((matrix[10][5]**2 + matrix[11][5]**2 + matrix[12][5]**2 + matrix[13][5]**2 + matrix[14][5]**2 + matrix[15][5]**2 + matrix[16][5]**2 + matrix[17][5]**2)**0.5) #W+js
j.append((matrix[18][5]**2 + matrix[19][5]**2 + matrix[20][5]**2 + matrix[21][5]**2 + matrix[22][5]**2 + matrix[23][5]**2 + matrix[24][5]**2)**0.5) #Z+js
#Uncert. 0 Tau VVLID
j.append((matrix[0][13]**2 + matrix[1][13]**2 + matrix[2][13]**2 + matrix[3][13]**2 + matrix[4][13]**2 + matrix[5][13]**2 + matrix[6][13]**2 + matrix[7][13]**2 + matrix[8][13]**2)**0.5) #QCD
j.append(matrix[9][13])#TT
j.append((matrix[10][13]**2 + matrix[11][13]**2 + matrix[12][13]**2 + matrix[13][13]**2 + matrix[14][13]**2 + matrix[15][13]**2 + matrix[16][13]**2 + matrix[17][13]**2)**0.5) #W+js
j.append((matrix[18][13]**2 + matrix[19][13]**2 + matrix[20][13]**2 + matrix[21][13]**2 + matrix[22][13]**2 + matrix[23][13]**2 + matrix[24][13]**2)**0.5) #Z+js

m.append(w)
m.append(j)
m.append(w[0] + w[1] + w[2] + w[3]) #total no cuts
m.append(w[4] + w[5] + w[6] + w[7]) #total mets
m.append(w[8] + w[9] + w[10] + w[11]) #total 1 tau MID
m.append((j[0]**2 + j[1]**2 + j[2]**2 + j[3]**2)**0.5) #total uncer. no cuts
m.append((j[4]**2 + j[5]**2 + j[6]**2 + j[7]**2)**0.5) #total uncer. mets
m.append((j[8]**2 + j[9]**2 + j[10]**2 + j[11]**2)**0.5) #total uncer. 1 tau MID
m.append(matrix[25][6]) #Data no cuts
m.append(matrix[25][8]) #Data Mets
m.append(matrix[25][10]) #Data 1 tau MID
m.append(matrix[25][7]) #Data unc. no cuts
m.append(matrix[25][9]) #Data unc. mets
m.append(matrix[25][11]) #Data unc. 1 tau MID 
m.append(w[12] + w[13] + w[14] + w[15]) #total 0 tau VVLID
m.append((j[12]**2 + j[13]**2 + j[14]**2 + j[15]**2)**0.5) #total uncer. 0 tau VVLID
m.append(matrix[25][14]) #Data 1 tau MID
m.append(matrix[25][15]) #Data unc. 1 tau MID 

Table = PrettyTable(display_size=(3,8))
Table.field_names = ["Selections", "QCD", "TT", "W+jets", "Z+jets", "Total"]
Table.align["Selections"] = "l"
Table.align["QCD"]        = "r"
Table.align["TT"]         = "r"
Table.align["W+jets"]     = "r"
Table.align["Z+jets"]     = "r"
Table.align["Total"]      = "r"
Table.add_row(["No cuts", '%.1f +/- %.1f'%(m[0][0], m[1][0]), '%.1f +/- %.1f'%(m[0][1], m[1][1]), '%.1f +/- %.1f'%(m[0][2], m[1][2]), '%.1f +/- %.1f'%(m[0][3], m[1][3]), '%.1f +/- %.1f'%(m[2], m[5])])
Table.add_row(["MET trigger + MET200", '%.1f +/- %.1f'%(m[0][4], m[1][4]), '%.1f +/- %.1f'%(m[0][5], m[1][5]), '%.1f +/- %.1f'%(m[0][6], m[1][6]), '%.1f +/- %.1f'%(m[0][7], m[1][7]), '%.1f +/- %.1f'%(m[3], m[6])])
Table.add_row(["1 Tau MID", '%.1f +/- %.1f'%(m[0][8], m[1][8]), '%.1f +/- %.1f'%(m[0][9], m[1][9]), '%.1f +/- %.1f'%(m[0][10], m[1][10]), '%.1f +/- %.1f'%(m[0][11], m[1][11]), '%.1f +/- %.1f'%(m[4], m[7])])
Table.add_row(["0 Tau VVLID", '%.1f +/- %.1f'%(m[0][12], m[1][12]), '%.1f +/- %.1f'%(m[0][13], m[1][13]), '%.1f +/- %.1f'%(m[0][14], m[1][14]), '%.1f +/- %.1f'%(m[0][15], m[1][15]), '%.1f +/- %.1f'%(m[14], m[15])])
print(Table.get_string())

table =PrettyTable(display_size=(3,8))
table.field_names = ["Selections", "Total MC", "Data", "Data/Total"]
table.align["Selections"] = "l"
table.align["Total MC"]   = "r"
table.align["Data"]       = "r"
table.align["Data/Total"] = "r"
table.add_row(["No cuts",'%.1f +/- %.1f'%(m[2], m[5]), '%.1f +/- %.1f'%(m[8], m[11]), '%.2f'%(m[8]/m[2])])
table.add_row(["MET trigger + MET200",'%.1f +/- %.1f'%(m[3], m[6]), '%.1f +/- %.1f'%(m[9], m[12]), '%.2f'%(m[9]/m[3])])
table.add_row(["1 tau MID",'%.1f +/- %.1f'%(m[4], m[7]), '%.1f +/- %.1f'%(m[10], m[13]), '%.2f'%(m[10]/m[4])])
table.add_row(["0 tau VVLID",'%.1f +/- %.1f'%(m[14], m[15]), '%.1f +/- %.1f'%(m[16], m[17]), '%.2f'%(m[16]/m[14])])
print(table.get_string())
