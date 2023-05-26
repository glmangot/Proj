from ROOT import TFile, TTree, TH1F, TH2F, TObject
import numpy as np
import scipy as sp

ttTrees = ["2995", "5571", "5767", "9224", "9480"]
zzTrees = ["4261", "5329"]#, "7046", "7643", "9875"]

file1 = TFile("ttbar_write.root", "read")
file2 = TFile("gghzz_write.root", "read")

outFile = TFile("Results.root", "recreate")

tree_ttOut = TTree("ttfrac", "ttfrac")
tree_zzOut = TTree("zzfrac", "zzfrac")

frac_tt = np.zeros(1, dtype = float)
cut_iso = np.zeros(1, dtype = float)
frac_zz = np.zeros(1, dtype = float)

tree_ttOut.Branch("frac", frac_tt, "frac/D")
tree_ttOut.Branch("cut_iso", cut_iso, "cut_iso/D")
tree_zzOut.Branch("frac", frac_zz, "frac/D")
tree_zzOut.Branch("cut_iso", cut_iso, "cut_iso/D")

x = np.linspace(0.0, 0.5, 100)

for run in range(len(x)):
    frac_tt[0] = 0.0
    net_tt = 0.0
    thru_tt = 0.0
    cut_iso[0] = x[run]

    frac_zz[0] = 0.0
    net_zz = 0.0
    thru_zz = 0.0
    
    for seed in ttTrees:
        tree = file1.Get(seed)
        for i in range(0, tree.GetEntries(), 4):
            tree.GetEntry(i)
            id1, px1, py1, pz1, iso1 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            tree.GetEntry(i + 1)
            id2, px2, py2, pz2, iso2 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            tree.GetEntry(i + 2)
            id3, px3, py3, pz3, iso3 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            tree.GetEntry(i + 3)
            id4, px4, py4, pz4, iso4 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            isos = [iso1, iso2, iso3, iso4]

            net_tt += 1.0

            pT1 = np.sqrt(px1**2 + py1**2)
            pT2 = np.sqrt(px2**2 + py2**2)
            pT3 = np.sqrt(px3**2 + py3**2)
            pT4 = np.sqrt(px4**2 + py4**2)
            pTs = [pT1, pT2, pT3, pT4]

            pTs = np.flip(np.sort(pTs))

            if all(iso > cut_iso for iso in isos) and pTs[0] > 20.0 and pTs[1] > 10.0 and pTs[2] > 10.0:
                
                thru_tt += 1.0

    frac_tt[0] = thru_tt / net_tt
    tree_ttOut.Fill()
    print(run, "tt: ", frac_tt, " ", cut_iso)

    for seed in zzTrees:
        tree = file2.Get(seed)
        for i in range(0, tree.GetEntries(), 4):
            tree.GetEntry(i)
            id1, px1, py1, pz1, iso1 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            tree.GetEntry(i + 1)
            id2, px2, py2, pz2, iso2 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            tree.GetEntry(i + 2)
            id3, px3, py3, pz3, iso3 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            tree.GetEntry(i + 3)
            id4, px4, py4, pz4, iso4 = tree.id, tree.px, tree.py, tree.pz, tree.iso

            isos = [iso1, iso2, iso3, iso4]

            net_zz += 1.0

            pT1 = np.sqrt(px1**2 + py1**2)
            pT2 = np.sqrt(px2**2 + py2**2)
            pT3 = np.sqrt(px3**2 + py3**2)
            pT4 = np.sqrt(px4**2 + py4**2)
            pTs = [pT1, pT2, pT3, pT4]

            pTs = np.flip(np.sort(pTs))

            if all(iso > cut_iso for iso in isos) and pTs[0] > 20.0 and pTs[1] > 10.0 and pTs[2] > 10.0:
                
                thru_zz += 1.0

    frac_zz[0] = thru_zz / net_zz
    tree_zzOut.Fill()
    print(run, "zz: ", frac_zz, " ", cut_iso, "\n")

tree_ttOut.Write()
tree_zzOut.Write()

outFile.Close()
file1.Close()