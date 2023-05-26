#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TApplication.h"

#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <iostream>

const float PI = 3.14159265358979;

// inline float get_dphi(const float &p1, const float &p2) {
//     if      (p2 - p1 >  PI)    return p2 - p1 - PI;
//     else if (p2 - p1 < -PI)    return p2 - p1 + PI;
//     else                       return p2 - p1;
// }

// inline float get_dR(const float &a_eta, const float &a_phi, const float &b_eta, const float &b_phi) {
//     float dphi = get_dphi(a_phi, b_phi);
//     float deta = b_eta - a_eta;

//     return sqrt(deta*deta + dphi*dphi);
// }

// inline float get_isolation(const float &p_pT, const vector<float> &pT_charged, const vector<float> &pT_neutral, const vector<float> &pT_gammas) {
//     float _neutral(0), _charged(0), _gammas(0);

//     for (int i = 0; i < pT_charged.size(); i++) _charged += pT_charged[i];
//     for (int i = 0; i < pT_neutral.size(); i++) _neutral += pT_neutral[i];
//     for (int i = 0; i < pT_gammas.size();  i++) _gammas  += pT_gammas[i];

//     return (_charged + _neutral + _gammas) / p_pT;
// }

#define get_pT(px, py) sqrt(px*px + py*py)
#define get_p_scalar(px, py, pz) sqrt(px*px + py*py + pz*pz)
#define get_eta(px, py, pz) float(-log(tan(acos(pz/get_p_scalar(px,py,pz)) / 2.)))
#define get_phi(px, py, pz) float(atan2(py,px))

inline float get4M(const float &px1, const float &py1, const float &pz1, const float &e1,
                  const float &px2, const float &py2, const float &pz2, const float &e2,
                  const float &px3, const float &py3, const float &pz3, const float &e3,
                  const float &px4, const float &py4, const float &pz4, const float &e4) {

    return sqrt(pow(e1 + e2 + e3 + e4, 2) - pow(get_p_scalar(px1 + px2 + px3 + px4, py1 + py2 + py3 + py4, pz1 + pz2 + pz3 + pz4), 2));
}

int main(int argc, char* argv[]) {
    TApplication theApp("Hist", &argc, argv);

    TFile *file1 = TFile::Open("gg2ttbar.root", "read");
    TFile *file2 = TFile::Open("gg2ttbar.root", "read");
    TFile *file3 = TFile::Open("gg2ttbar.root", "read");
    TFile *file4 = TFile::Open("gg2ttbar.root", "read");

    TFile *outFile = new TFile("gg2ttbar_read.root", "recreate");

    TH1F *h_iso = new TH1F("h_iso", "Isolation constant for muons; num; Iso", 1.3e2, -1, 14);
    TH1F *h_m4l = new TH1F("h_m4l", "M_{4l} for ttbar muons; num; M", 1.3e2, 0, 1.3e2);
    

    int id1, id2, id3, id4;

    float px1, py1, pz1, e1, iso1;
    float px2, py2, pz2, e2, iso2;
    float px3, py3, pz3, e3, iso3;
    float px4, py4, pz4, e4, iso4;

    char q1, q2, q3, q4;
    std::vector<std::string> treeNames {"3209"};
    
    
    TTree *mu1 = (TTree*)file1->Get("3209");
    TTree *mu2 = (TTree*)file2->Get("3209");
    TTree *mu3 = (TTree*)file3->Get("3209");
    TTree *mu4 = (TTree*)file4->Get("3209");

    std::cout << mu1->GetEntries();
    // for (int i = 0; i < mu1->GetEntries() - 3; i += 4) {
        // if (!mu1 || !mu2 || !mu3 || !mu4) break;

        // mu1->SetBranchAddress("id", &id1);
        // mu1->SetBranchAddress("px", &px1);
        // mu1->SetBranchAddress("py", &py1);
        // mu1->SetBranchAddress("pz", &pz1);
        // mu1->SetBranchAddress("e", &e1);
        // mu1->SetBranchAddress("q", &q1);
        // mu1->SetBranchAddress("iso", &iso1);

        // mu2->SetBranchAddress("id", &id2);
        // mu2->SetBranchAddress("px", &px2);
        // mu2->SetBranchAddress("py", &py2);
        // mu2->SetBranchAddress("pz", &pz2);
        // mu2->SetBranchAddress("e", &e2);
        // mu2->SetBranchAddress("q", &q2);
        // mu2->SetBranchAddress("iso", &iso2);

        // mu3->SetBranchAddress("id", &id3);
        // mu3->SetBranchAddress("px", &px3);
        // mu3->SetBranchAddress("py", &py3);
        // mu3->SetBranchAddress("pz", &pz3);
        // mu3->SetBranchAddress("e", &e3);
        // mu3->SetBranchAddress("q", &q3);
        // mu3->SetBranchAddress("iso", &iso3);

        // mu4->SetBranchAddress("id", &id4);
        // mu4->SetBranchAddress("px", &px4);
        // mu4->SetBranchAddress("py", &py4);
        // mu4->SetBranchAddress("pz", &pz4);
        // mu4->SetBranchAddress("e", &e4);
        // mu4->SetBranchAddress("q", &q4);
        // mu4->SetBranchAddress("iso", &iso4);

        // mu1->GetEntry(i);
        // mu2->GetEntry(i + 1);
        // mu3->GetEntry(i + 2);
        // mu4->GetEntry(i + 3);

        // h_iso->Fill(iso1);
        // h_iso->Fill(iso2);
        // h_iso->Fill(iso3);
        // h_iso->Fill(iso4);

        // h_m4l->Fill(get4M(px1, py1, pz1, e1, px2, py2, pz2, e2, px3, py3, pz3, e3, px4, py4, pz4, e4));

        // mu1->ResetBranchAddresses();
        // mu2->ResetBranchAddresses();
        // mu3->ResetBranchAddresses();
        // mu4->ResetBranchAddresses();

        // if ((i + 1) % 16 == 0) {
        //     std::cout << "  Gone through " << i << " muons\n";
        // }
    // }
    delete mu1;
    delete mu2;
    delete mu3;
    delete mu4;

    h_iso->Write();
    h_m4l->Write();

    delete h_iso;
    delete h_m4l;
 

    file1->Close();
    file2->Close();
    file3->Close();
    file4->Close();
    
    outFile->Close();
}