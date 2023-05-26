from ROOT import TFile, TTree, TH1F, TH2F, TObject
import numpy as np
import matplotlib.pyplot as plt
plt.close()
file1 = TFile("Results.root", "read")

tt_tree = file1.Get("ttfrac")
zz_tree = file1.Get("zzfrac")

tt = np.array([])
zz = np.array([])

for i in range(100):
    tt_tree.GetEntry(i)
    zz_tree.GetEntry(i)
    tt = np.append(tt, tt_tree.frac)
    zz = np.append(zz, zz_tree.frac)

plt.plot(tt, zz)
plt.Normalize()
plt.xlim([0.0, 1.1])
plt.ylim([0.0, 1.1])
plt.xlabel('tt')
plt.ylabel('zz')
plt.show()