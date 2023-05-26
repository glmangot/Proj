from ROOT import TFile, TTree, TH1F, TH2F, TObject
import numpy as np


ttTrees = {"2995", "5571", "5767", "9224", "9480"}

file1 = TFile("ttbar_write.root", "read")
file2 = TFile("ttbar_write.root", "read")
file3 = TFile("ttbar_write.root", "read")
file4 = TFile("ttbar_write.root", "read")

outFile = TFile("ttbar_py.root", "recreate")

# Go through trees, if they don't pass, -> bkgd, If they pass, ->sgl
for tree in range(5):
    id = np.zeros(1, dtype = np.int16)
    px = np.zeros(1, dtype = np.float32)
    py = np.zeros(1, dtype = np.float32)
    pz = np.zeros(1, dtype = np.float32)
    iso = np.zeros(1, dtype = np.float32)

    hist = TH1F("hist", "pT test", 200, 0, 200)
    h_id = TH1F("h_id", "id test", 100, -14, 14)
    h_iso = TH1F("h_iso", "iso test", 900, 0, 30)
    h_ispT = TH2F("h_ispT", "iso vs. pT", 200, 0, 200, 900, 0, 30)

    mu1 = file1.Get(ttTrees[tree])

    mu1.SetBranchAddress("id", id)
    mu1.SetBranchAddress("px", px)
    mu1.SetBranchAddress("py", py)
    mu1.SetBranchAddress("pz", pz)
    mu1.SetBranchAddress("iso", iso)

    for i in range(mu1.GetEntries()):
        mu1.GetEntry(i)
        pT = np.sqrt(px[0]*px[0] + py[0]*py[0])
        
        hist.Fill(pT)
        h_id.Fill(id)
        h_iso.Fill(iso)
        h_ispT.Fill(pT, iso)

hist.Write()
h_id.Write()
h_iso.Write()
h_ispT.Write()

file1.Close()
outFile.Close()