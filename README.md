# VLL Run-2 Search
Code to analyze MC and data to search for vector-like lepton pairs, each decaying into a tau lepton and a long-lived pseudoscalar

## Install

Log in to the LPC machines
```
ssh -XY username@cmslpc-sl7.fnal.gov
```

Get code from git. You can use the nobackup area for more space
````
cd nobackup
mkdir VLLAnalysis
cd VLLAnalysis
cmsrel CMSSW_10_6_8
cd CMSSW_10_6_8/src/
cmsenv
git clone https://github.com/FNALLPC-VLL/TauClusterAnalysis.git
cd TauClusterAnalysis 
````

## Make histograms

A script template of how to read the ntuples and make histograms is Study_TauProperties.py. The command to run it is the following: 
````
python Study_TauProperties.py --config config/configuration_2018.cfg --tag taustudies --maxnevents 10000
````
where 
````
--config     The configuration file with the information about the inputs/output/year directories
--tag        the name of the output folder (defined by user)
--maxnevents the maximum number of events to be processed (all by default)
````

## Make bkg studies
A script template of how to read the ntuples and bkg studies is Study_BkgComposition.py. The command to run it is the following: 
````
python Study_BkgComposition.py --config config/configuration_2018.cfg --tag bkgtudies --maxnevents 10000
````
where 
````
--config     The configuration file with the information about the inputs/output/year directories,bkg samples, cross sections
--tag        the name of the output folder (defined by user)
--maxnevents the maximum number of events to be processed (all by default)
````

## Plot histograms

A script template of how to read the histograms and plot them in pyROOT is Plotter_TauProperties.py. The command to run it is the following: 
````
python Plotter_TauProperties.py --config config/configuration_2018.cfg --tag taustudies
````
where 
````
--config     configuration file with the information about the inputs and output directories
--tag        the name of the input/output folder (defined by user)
````

## Make CSC selection studies
A script template of how to read the ntuples and do CSC+Tau selection studies is Study_CSCTauRegion.py. The command to run it is the following:
````
python Study_CSCTauRegion.py  --config config/configuration_2018_table.cfg --tag metoptimizationbkgstudiesfull
````
where
````
--config     The configuration file with the information about the inputs/output/year directories, cross sections
--tag        the name of the input ntuples version folder (defined by user)
````
