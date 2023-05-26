#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH1F.h"

#include <vector>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <iostream>

//Plot iso from W muons
// ./ttbar_write > ttbar_out
using namespace std;

const double PI = 3.14159265358979;

inline double get_dphi(const double& p1, const double& p2) {
    if (p2 - p1 > PI) return p2 - p1 - PI;
    else if (p2 - p1 < -PI) return p2 - p1 + PI;
    else                    return p2 - p1;
}

inline double get_deta(const double& p1, const double& p2) {
    return p2 - p1;
}

static bool sortByPT(const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.second < b.second;
}

inline double get_dR(const double& a_eta, const double& a_phi, const double& b_eta, const double& b_phi) {
    double dphi = get_dphi(a_phi, b_phi);
    double deta = get_deta(a_eta, b_eta);

    return sqrt(deta * deta + dphi * dphi);
}

inline std::vector<std::pair<int, double>> get4mus(std::vector<pair<int, double>>& a) {
    std::vector<std::pair<int, double>> out;

    if (a.size() >= 4) {
        sort(a.begin(), a.end(), sortByPT);
        out.push_back(a[0]);
        out.push_back(a[1]);
        out.push_back(a[2]);
        out.push_back(a[3]);
        return out;
    }
    else
        return a;
}

int main(int argc, char* argv[]) {
    TApplication theApp("Hist", &argc, argv);

    //TFile* file = new TFile("ttbar_write.root", "update");
    TFile* file = new TFile("gghzz_write.root", "update");

    int nEvents, nfiles, id;

    float px, py, pz, e, iso;
    bool q;
    string proc;

    std::cout << "\n\nWhich process? \n\n   'tt' or 'zz'? ";
    std::cin >> proc;

    std::cout << "\n\nStarting process\n\n   How many pp collisions? ";
    std::cin >> nEvents;
    std::cout << "\n\n      Gonna work " << nEvents << " events!!";

    std::cout << "\n\n   How many runs? ";
    std::cin >> nfiles;
    std::cout << "\n\n      Gonna make " << nfiles << " files with " << nEvents << " events each!!!\n\n";

    sleep(2);

    for (int ifiles = 0; ifiles < nfiles; ifiles++) {
        Pythia8::Pythia pythia;
        Pythia8::Event& ev = pythia.event;
        
        std::srand(std::time(nullptr));
        int seed = std::rand() % 9000 + 1000;


        if (proc == "tt") {
            std::cout << "\n   ttbar run: \n";
            pythia.readString("Random:setSeed = on");
            pythia.readString(("Random:seed = " + to_string(seed)).c_str());

            pythia.readString("Top:gg2ttbar = on");
            pythia.readString("-24:onMode = off");
            pythia.readString("24:onMode = off");
            pythia.readString("24:onIfAny = 13 14");
            pythia.readString("-24:onIfAny = -13 -14");

            pythia.readString("Beams:eCM = 1.3e4");
            pythia.init();
        } 
        else if (proc == "zz") {
            std::cout << "\n   hzz run: \n";
            pythia.readString("Random:setSeed = on");
            pythia.readString(("Random:seed = " + to_string(seed)).c_str());

            pythia.readString("HiggsSM:gg2H = on");
            pythia.readString("25:onMode = off");
            pythia.readString("25:onIfMatch = 23 23");
            pythia.readString("23:onMode = off");
            pythia.readString("23:onIfAny = 13");
            pythia.readString("Beams:eCM = 1.3e4");
            pythia.init();
        }
        else {
            std::cout << "No process! abort";
            break;
        }

        TTree* data = new TTree(to_string(seed).c_str(), "data");
        // TH1F *h_Wiso = new TH1F(("h_Wiso" + std::to_string(seed)).c_str(), "Isolation for muons directly from W; Iso; num", 1e2, 0, 3);
        // TH1F *h_miso = new TH1F(("h_miso" + std::to_string(seed)).c_str(), "Isolation for muons w/o W requirement; Iso; num", 1e2, 0, 3);
        // TH1F *h_aiso = new TH1F(("h_aiso" + std::to_string(seed)).c_str(), "Isolation for muons w/ no requirement; Iso; num", 1e2, 0, 3);

        data->Branch("id", &id, "id/I");
        data->Branch("px", &px, "px/F");
        data->Branch("py", &py, "py/F");
        data->Branch("pz", &pz, "pz/F");
        data->Branch("iso", &iso, "iso/F");

        for (int iEvent = 0; iEvent < nEvents; iEvent++) {
            if (!pythia.next()) continue;
            
            if (iEvent % 10000 == 0) cout << " \n\n BingBong::next(): On event: " << iEvent << " in run: " << ifiles + 1 << "\n\n"; 

            std::vector<pair<int, double>> musa;

            for (int i = 0; i < ev.size(); i++) {
                if (ev[i].isFinal() && abs(ev[i].id()) == 13 && abs(ev[i].eta()) < 2.4 && ev[i].pT() > 3.)
                    musa.push_back(std::make_pair(i, ev[i].pT()));
                // if (ev[i].idAbs() == 13 && ev[ev[i].mother1()].idAbs() == 24) {
                //     double pT_sump(0);
                //     for (int j = 0; j < ev.size(); j++) {
                //         if (i != j && ev[j].idAbs() != 13 && abs(ev[j].eta()) < 2.4 && ev[j].isFinal()) {
                //             double dRs = get_dR(ev[i].eta(), ev[i].phi(), ev[j].eta(), ev[j].phi());
                //             if (dRs < 0.3)
                //                 pT_sump += ev[j].pT();
                //         }
                //     }
                //     // h_Wiso->Fill(pT_sump / ev[i].pT());
                // }
            }

            std::vector<pair<int, double>> mus = get4mus(musa);

            if (mus.size() == 4) {
                for (int i = 0; i < mus.size(); i++) {
                    id = ev[mus[i].first].id();
                    px = float(ev[mus[i].first].px());
                    py = float(ev[mus[i].first].py());
                    pz = float(ev[mus[i].first].pz());

                    double pT_sum(0);

                    for (int j = 0; j < ev.size(); j++) {
                        if (abs(ev[j].eta()) < 2.4 && ev[j].isFinal() && mus[i].first != j) {
                            double dR = get_dR(ev[mus[i].first].eta(), ev[mus[i].first].phi(), ev[j].eta(), ev[j].phi());
                            if (dR < 0.3)
                                pT_sum += ev[j].pT();
                        }
                    }
                    iso = float(pT_sum / mus[i].second);
                    //h_miso->Fill(iso);
                    data->Fill();
                }
            }
        }
        // h_Wiso->Write();
        // h_miso->Write();
        data->Write();
        pythia.stat();
        delete data;
    }
    file->Close();
}