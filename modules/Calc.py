import ROOT
from ROOT import TMath

def DeltaR(llp_eta,llp_phi,cluster_eta,cluster_phi):
    dphi = DeltaPhi(llp_phi,cluster_phi)
    deta = llp_eta-cluster_eta

    return TMath.Sqrt( dphi*dphi + deta*deta)

def DeltaPhi(llp_phi,cluster_phi):
    dphi = llp_phi-cluster_phi
    while dphi > TMath.Pi():
          dphi -= TMath.TwoPi()
    while dphi <= -TMath.Pi():
          dphi += TMath.TwoPi()
    return dphi

def DeltaEta(llp_eta,cluster_eta):
    deta = abs(llp_eta-cluster_eta)
    return deta