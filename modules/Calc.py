import ROOT
from ROOT import TMath

def ctau(pt,m,eta,vx,vy,vz):
    ct = TMath.Sqrt(vx*vx + vy*vy + vz*vz) * m /(pt*TMath.CosH(eta) )
    return ct
def v(vx,vy,vz):
    v  = TMath.Sqrt(vx*vx + vy*vy + vz*vz)
    return v

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
    dphi = abs(dphi)
    return dphi

def DeltaEta(llp_eta,cluster_eta):
    deta = abs(llp_eta-cluster_eta)
    return deta