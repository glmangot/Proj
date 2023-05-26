// Greg Mangot - 4/4/23
// Full implementation of most calcs

#include "Pythia8/Pythia.h"

#include "TFile.h"
#include "TApplication.h"
#include "TH1F.h"

#include <vector>
#include <chrono>
#include <tuple>
#include <algorithm>
// CHECK RADS VS DEG!!!!!
#include <cmath>
#include <iostream>

using namespace std;

// Architecture:
typedef tuple<int, double, double, double, double, Pythia8::Vec4, bool> MParticle;
typedef vector<MParticle> MEvent; 
typedef vector<MEvent> MRoot;

const double PI = 3.14159265358979323;

static bool sortByMass(const tuple<MParticle, MParticle, double> &a, const tuple<MParticle, MParticle, double> &b) {
    return get<2>(a) > get<2>(b);
}

static bool sortByPT(const double &a, const double &b) {
    return a > b;
}
inline double get_dphi(const double &p1, const double &p2) {
    if      (p2 - p1 >  PI) return p2 - p1 - PI;
    else if (p2 - p1 < -PI) return p2 - p1 + PI;
    else                    return p2 - p1;
}

inline double get_deta(const double &p1, const double &p2) {
    return p2 - p1;
}

inline double get_dR(const double &a_eta, const double &a_phi, const double &b_eta, const double &b_phi) {
    double dphi = get_dphi(a_phi, b_phi);
    double deta = get_deta(a_eta, b_eta);

    return sqrt(deta*deta + dphi*dphi);
}

inline float get_isolation(const double &p_pT, const vector<double> &pT_charged, const vector<double> &pT_neutral, const vector<double> &pT_gammas) {
    double _neutral(0), _charged(0), _gammas(0);

    for (int i = 0; i < pT_charged.size(); i++) _charged += pT_charged[i];
    for (int i = 0; i < pT_neutral.size(); i++) _neutral += pT_neutral[i];
    for (int i = 0; i < pT_gammas.size();  i++) _gammas  += pT_gammas[i];

    return (_charged + _neutral + _gammas) / p_pT;
}

int main(int argc, char* argv[]) {
    TApplication theApp("Hist", &argc, argv);

    Pythia8::Pythia pythia;
    Pythia8::Event &ev = pythia.event;
    
    pythia.readString("");
    pythia.readString("HiggsSM:gg2H = on");
    pythia.readString("25:onMode = off");
    pythia.readString("25:onIfMatch = 23 23");
    pythia.readString("23:onMode = off");

    pythia.readString("23:onIfAny = 11 13");
    pythia.readString("23:onIfAny = 11");
    pythia.readString("23:onIfAny = 13");
    pythia.readString("Beams:eCM = 1.3e4");
    pythia.init();
    
    TFile *file = new TFile("gg2H.root", "recreate");
    int nEvents;

    cout << "\n\nHow many pp collisions? ";
    cin >> nEvents;
    cout << "\n\nOK! Generating " << nEvents << " events!!!";
    sleep(2);    

    // Create fatty data matrix
    MRoot  _root;

    for (int iev = 0; iev < nEvents; iev++) {
        if (!pythia.next()) continue;
        // Create event vector
        MEvent _event;
        for (int i = 0; i < ev.size(); i++) {
            if (ev[i].isFinal()) {
                // Push a particle tuple into the event for every final particle
                _event.push_back(make_tuple(

                    ev[i].id(),        // 0

                    ev[i].eta(),       // 1

                    ev[i].phi(),       // 2

                    ev[i].m(),         // 3

                    ev[i].pT(),        // 4

                    ev[i].p(),         // 5

                    ev[i].isCharged()  // 6

                ));
            }
        }
        _root.push_back(_event);
    }
    

    auto start = chrono::high_resolution_clock::now();
    double eta_cut = 10;
    vector<int> unbuildable_events;
    for (int iruns = 1; iruns <= 19; iruns++) {
        cout << "\n    Starting run " << iruns << endl;
        auto count = chrono::high_resolution_clock::now();
        auto plchold = chrono::duration_cast<chrono::milliseconds>(count-start);
        cout << "    Elapsed time: " << plchold.count() * 1e-3 << endl;
        
        eta_cut -= 0.5;

        TH1F *h_dR = new TH1F(("h_dR" + to_string(iruns)).c_str(), ("#DeltaR test with #eta < |" + to_string(eta_cut) + "|; #DeltaR; num").c_str(), 1.3e2, -1, 12);
        TH1F *h_M2L = new TH1F(("h_M2L" + to_string(iruns)).c_str(), ("m_{2l} with #eta < |" + to_string(eta_cut) + "|; m_{2l}; num").c_str(), 1.3e2, -1, 1.5e2);
        TH1F *h_M4L = new TH1F(("h_M4L" + to_string(iruns)).c_str(), ("m_{4l} with #eta < |" + to_string(eta_cut) + "|; m_{4l}; num").c_str(), 1.3e2, 1e2, 1.4e2);
        TH1F *h_iso = new TH1F(("h_iso" + to_string(iruns)).c_str(), ("Isolation constant with #eta < |" + to_string(eta_cut) + "|; Isolation; num").c_str(), 1.3e2, -1, 25);

        int too_few_lep_events(0);
        for (int iev = 0; iev < nEvents; iev++) {
            MEvent curr_ev = _root[iev];

            vector<MParticle> mu_pos, mu_neg, be_pos, be_neg;
            for (int i = 0; i < curr_ev.size(); i++) {
                if ((abs(get<0>(curr_ev[i])) == 11 || abs(get<0>(curr_ev[i])) == 13) && abs(get<1>(curr_ev[i])) < eta_cut /*&& get<4>(curr_ev[i]) > 5*/) {
                    vector<double> charged, neutral, gammas;
                    for (int j = 0; j < curr_ev.size(); j++) {
                        if (i == j && j != curr_ev.size() - 1) j++;
                        else if (abs(get<1>(curr_ev[j])) < eta_cut && abs(get<0>(curr_ev[j])) != 11 && abs(get<0>(curr_ev[j])) != 13) {
                            double dR = get_dR(get<1>(curr_ev[i]), get<2>(curr_ev[i]), get<1>(curr_ev[j]), get<2>(curr_ev[j]));
                            if (dR > 0.02) {
                                h_dR->Fill(dR);
                            }
                            // dR < 0.3 for iso: P4
                            if (dR < 0.3) {
                                if (get<0>(curr_ev[j]) == 22)      gammas.push_back(get<4>(curr_ev[j]));
                                else if (get<6>(curr_ev[j]))      charged.push_back(get<4>(curr_ev[j]));
                                else                              neutral.push_back(get<4>(curr_ev[j]));
                            }
                        }
                    }
                    double iso = get_isolation(get<4>(curr_ev[i]), charged, neutral, gammas);
                    h_iso->Fill(iso);

                    // for muons (only IRL) but we'll pretend it works
                    if (iso < 0.35) {
                        if      (get<0>(curr_ev[i]) == 11)
                            be_neg.push_back(curr_ev[i]);
                        else if (get<0>(curr_ev[i]) == -11)
                            be_pos.push_back(curr_ev[i]);
                        else if (get<0>(curr_ev[i]) == 13)
                            mu_pos.push_back(curr_ev[i]);
                        else if (get<0>(curr_ev[i]) == -13)
                            mu_neg.push_back(curr_ev[i]);
                    }
                }
            }
            vector<tuple<MParticle, MParticle, double>> lepton_pairs;
            if (be_neg.size() >= 1 && be_pos.size() >= 1) {
                for (int i = 0; i < be_neg.size(); i++) {
                    for (int j = 0; j < be_pos.size(); j++) {
                        double MZ = (get<5>(be_neg[i]) + get<5>(be_pos[j])).mCalc();
                        if (MZ > 12 && MZ < 120)
                            lepton_pairs.push_back(make_tuple(be_neg[i], be_pos[j], MZ));
                    }
                }
            }
            if (mu_neg.size() >= 1 && mu_pos.size() >= 1) {
                for (int i = 0; i < mu_neg.size(); i++) {
                    for (int j = 0; j < mu_pos.size(); j++) {
                        double MZ = (get<5>(mu_neg[i]) + get<5>(mu_pos[j])).mCalc();
                        if (MZ > 12 && MZ < 120)
                            lepton_pairs.push_back(make_tuple(mu_neg[i], mu_pos[j], MZ));
                    }
                }
            }

            if (lepton_pairs.size() > 1) {
                sort(lepton_pairs.begin(), lepton_pairs.end(), sortByMass);
                while (lepton_pairs.size() > 2 && 
                ((get<0>(get<0>(lepton_pairs[0])) == get<0>(get<0>(lepton_pairs[1]))) || 
                 (get<0>(get<1>(lepton_pairs[0])) == get<0>(get<1>(lepton_pairs[1]))))) {
                    
                    lepton_pairs.erase(lepton_pairs.begin() + 1);
                    sort(lepton_pairs.begin(), lepton_pairs.end(), sortByMass);
                }

                double mass_4l = (get<5>(get<0>(lepton_pairs[0])) + get<5>(get<1>(lepton_pairs[0])) + get<5>(get<0>(lepton_pairs[1])) + get<5>(get<1>(lepton_pairs[1]))).mCalc();
                // Still need mixed flavors' inv mass > 4 GeV
                if ((get<0>(get<0>(lepton_pairs[0])) != get<0>(get<0>(lepton_pairs[1]))) && 
                    (get<0>(get<1>(lepton_pairs[0])) != get<0>(get<1>(lepton_pairs[1]))) &&
                    (mass_4l > 18 && mass_4l < 130) &&
                    (get<2>(lepton_pairs[0]) > 40)) {

                    vector<double> pT_check;
                    pT_check.push_back(get<4>(get<0>(lepton_pairs[0])));
                    pT_check.push_back(get<4>(get<1>(lepton_pairs[0])));
                    pT_check.push_back(get<4>(get<0>(lepton_pairs[1])));
                    pT_check.push_back(get<4>(get<1>(lepton_pairs[1])));
                    sort(pT_check.begin(), pT_check.end(), sortByPT);

                    if (pT_check[0] > 20 && pT_check[1] > 10 && pT_check[2] > 10) {
                        h_M2L->Fill(get<2>(lepton_pairs[0]));
                        h_M2L->Fill(get<2>(lepton_pairs[1]));
                        h_M4L->Fill(mass_4l);
                    }   
                }
                else too_few_lep_events++;
            }      
        }
        unbuildable_events.push_back(too_few_lep_events);
        h_dR->Write();
        h_M4L->Write();
        h_M2L->Write();
        h_iso->Write();

        delete h_dR;
        delete h_M4L;
        delete h_M2L;
        delete h_iso;
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "Analysis loop (best guess): " << duration.count() * 1e-3 << "\n";
    for (int i = 0; i < unbuildable_events.size(); i++) 
        cout << "\n Cut events per run: " << unbuildable_events[i] << endl;
    
    file->Write();
    delete file;
    return 0;
}
