#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include <vector>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <iostream>

using namespace std;

const double PI = 3.14159265358979;

inline double get_dphi(const double &p1, const double &p2) {
    if      (p2 - p1 >  PI) return p2 - p1 - PI;
    else if (p2 - p1 < -PI) return p2 - p1 + PI;
    else                    return p2 - p1;
}

inline double get_deta(const double &p1, const double &p2) {
    return p2 - p1;
}

static bool sortByPT(const std::pair<int, double> &a, const std::pair<int, double> &b) {
    return a.second < b.second;
}

inline double get_dR(const double &a_eta, const double &a_phi, const double &b_eta, const double &b_phi) {
    double dphi = get_dphi(a_phi, b_phi);
    double deta = get_deta(a_eta, b_eta);

    return sqrt(deta*deta + dphi*dphi);
}

inline std::vector<std::pair<int, double>> get4mus(std::vector<pair<int,double>> &a) {
    std::vector<std::pair<int,double>> out;

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

int main(int argc, char *argv[]) {
    TApplication theApp("Hist", &argc, argv);

    Pythia8::Pythia pythia;
    Pythia8::Event &ev = pythia.event;
    
    std::srand(std::time(nullptr));
    int seed = std::rand() % 9000 + 1000;


    pythia.readString("Random:setSeed = on");
    pythia.readString(("Random:seed = " + to_string(seed)).c_str());
    
    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("-24:onMode = off");
    pythia.readString("24:onMode = off");
    pythia.readString("24:onIfAny = 13 14");
    pythia.readString("-24:onIfAny = -13 -14");

    pythia.readString("Beams:eCM = 1.3e4");
    pythia.init();

    TFile *file = new TFile("gg2ttbar.root", "update");

    int nEvents, id;
    float px, py, pz, e, iso;
    bool q;

    std::cout << "\n\nHow many pp collisions? ";
    std::cin >> nEvents;
    std::cout << "\n\nOK! Generating " << nEvents << " events!!!";
    sleep(2);
        
    TTree *data = new TTree(to_string(seed).c_str(), "data");

    data->Branch("id", &id, "id/I");
    data->Branch("px", &px, "px/F");
    data->Branch("py", &py, "py/F");
    data->Branch("pz", &pz, "pz/F");
    data->Branch("e",  &e,  "e/F");
    data->Branch("q",  &q,  "q/B");
    data->Branch("iso", &iso, "iso/F");

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        if (!pythia.next()) continue;

        std::vector<pair<int, double>> musa;

        for (int i = 0; i < ev.size(); i++) {
            if (ev[i].isFinal() && abs(ev[i].id()) == 13 && ev[i].eta() < 2.4 && ev[i].pT() > 3)
                musa.push_back(std::make_pair(i,ev[i].pT()));
        }

        std::vector<pair<int,double>> mus = get4mus(musa);

        if (mus.size() == 4) {
            for (int i = 0; i < mus.size(); i++) {
                id = ev[mus[i].first].id();
                px = float(ev[mus[i].first].px());
                py = float(ev[mus[i].first].py());
                pz = float(ev[mus[i].first].pz());
                e  = float(ev[mus[i].first].e());
                q  = ev[mus[i].first].isCharged();

                double pT_sum(0);

                for (int j = 0; j < ev.size(); j++) {
                    if (i == j && j != ev.size() - 1) j++;
                    else if (abs(ev[j].eta()) < 2.4 && ev[j].isFinal()) {
                        double dR = get_dR(ev[mus[i].first].eta(), ev[mus[i].first].phi(), ev[j].eta(), ev[j].phi());
                        if (dR < 0.3)
                            pT_sum += ev[j].pT();                   
                    }
                }
                iso = float(pT_sum / mus[i].second);
                
                data->Fill();
            }
        }
        
    }
    data->Write();
    delete data;
}