# VLL Run-2 Analysis 
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
