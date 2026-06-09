
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TNamed.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <ctime>

static void StyleHist(TH1* h, int color, int marker){
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetLineWidth(2);
}

double Smearing(double pt, double C = 0.0246, double S = 1.2130, double N = 5.23) {
    double sigma = sqrt(
        C*C +
        pow(S/sqrt(pt), 2) +
        pow(N/pt, 2)
    );

    double smeared_pt = pt * (1.0 + gRandom->Gaus(0, sigma));
    return smeared_pt;
}

using namespace std;

int SaveHist(string inputFileName, string outFileName) {

    gStyle->SetOptStat(0); // Disable statistics box
    const int nPtBins_log = 35;
    
    const Double_t pTBins_log[nPtBins_log+1] = {
    8.5,  9.0,  9.5,    10.0,   11.,    12., 
    13.,  14.,  15.,    16.,    17.,    18.,
    19.,  20.,  21.,    22.6,   24.6,   26.6,
    28.6, 32.6, 36.6,   42.6,   48.6,   54.6,
    60.6, 74.0, 86.4,   103.6,  120.8,  140., 
    165., 250., 400.,   600.,   1000.,  2000.};

    TH1D* hAj = new TH1D("hAj", "Jet Asymmetry Distribution; A_{J}; Number of Jets", 12, 0, 0.7);    
    TH1D* hXj = new TH1D("hXj", "Jet Momentum Balance Distribution; X_{J}; Number of Jets", 100, -2, 2);
    TH1D* hdPhi = new TH1D("hdPhi", "Delta Phi Distribution; #Delta#phi; Number of Jets", 100, 0, TMath::Pi());

    TH1D* hjtpt = new TH1D("hjtpt", "Jet pT Distribution; Jet pT (GeV/c); Number of Jets", 100, 40, 2000);
    TH1D* hjtpt_log = new TH1D("hjtpt_log", "Jet pT Distribution; Jet pT (GeV/c); Number of Jets", nPtBins_log, pTBins_log);
    TH1D* hjteta = new TH1D("hjteta", "Jet eta Distribution; Jet eta; Number of Jets", 100, -5, 5);
    TH1D* hjtphi = new TH1D("hjtphi", "Jet phi Distribution; Jet phi; Number of Jets", 100, -TMath::Pi(), TMath::Pi());
    TH1D* hjtetmiss = new TH1D("hjtetmiss", "Jet Etmiss Distribution; Jet Etmiss (GeV); Number of Jets", 100, 0, 2000);

    TH1D* hAjsmeared = new TH1D("hAjsmeared", "Smeared Jet Asymmetry Distribution; A_{J}; Number of Jets", 12, 0, 0.7);
    TH1D* hjtptsmeared = new TH1D("hjtptsmeared", "Smeared Jet pT Distribution; Jet pT (GeV/c); Number of Jets", 100, 40, 2000);
    double jtpt1_min = 120.0;
    double jtpt2_min = 30.0;
    double dphicut_factor = 2.0/3.0; // Adjust this factor as needed
    double dphicut = dphicut_factor*TMath::Pi();
    double etacut = 2.0;

    TNamed* cutInfo = new TNamed(
        "Cuts",
        Form("Jet pT1 > %.1f GeV/c; "
            "Jet pT2 > %.1f GeV/c; "
            "#Delta#phi > 2/3 * #pi; "
            "|#eta| < %.1f",
            jtpt1_min,
            jtpt2_min,
            etacut)
    );
    
    hAj->Sumw2();
    hXj->Sumw2();
    hdPhi->Sumw2();
    hjtpt->Sumw2();
    hjtpt_log->Sumw2();
    hjteta->Sumw2();
    hjtphi->Sumw2();
    hjtetmiss->Sumw2();

    TFile* inFile = new TFile(inputFileName.c_str(),"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* JetTree = (TTree*)inFile->Get("JetTree");
    TTree* ParticleTree = (TTree*)inFile->Get("ParticleTree");
    vector<double>* jet_pt = nullptr;
    vector<double>* jet_eta = nullptr;
    vector<double>* jet_phi = nullptr;

    vector<double> jet_pt_smeared = {};
    double weights;
    
    Double_t jetpt1, jetpt2, jetphi1, jetphi2, ET_miss_jet, dPhi;
    Double_t A_J, X_J;
    ParticleTree->SetBranchAddress("weights", &weights);
    JetTree->SetBranchAddress("jet_pt", &jet_pt);
    JetTree->SetBranchAddress("jet_eta", &jet_eta);
    JetTree->SetBranchAddress("jet_phi", &jet_phi);
    JetTree->SetBranchAddress("ET_miss_jet", &ET_miss_jet);
    JetTree->SetBranchAddress("jetpt1", &jetpt1);
    JetTree->SetBranchAddress("jetpt2", &jetpt2);
    JetTree->SetBranchAddress("jetphi1", &jetphi1);
    JetTree->SetBranchAddress("jetphi2", &jetphi2);
    JetTree->SetBranchAddress("dPhi", &dPhi);  
    JetTree->SetBranchAddress("A_J", &A_J);
    JetTree->SetBranchAddress("X_J", &X_J);
    JetTree->SetBranchAddress("ET_miss_jet", &ET_miss_jet);
    Long64_t nEntries = JetTree->GetEntries();
    bool passEtaCuts = true;
    double newpt;

    int nEntriesPassSel = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nEntries; i++) {
        JetTree->GetEntry(i);
        ParticleTree->GetEntry(i);
        jet_pt_smeared.clear();
        for (size_t j = 0; j < jet_pt->size(); ++j) {
            newpt = Smearing(jet_pt->at(j));
            jet_pt_smeared.push_back(newpt);
        }

        double jetpt1_smeared = -1.0;
        double jetpt2_smeared = -1.0;
        double A_J_smeared = -1.0;

        int idx1 = -1, idx2 = -1;

        for (size_t j = 0; j < jet_pt->size(); ++j) {
            if (fabs(jet_eta->at(j)) > etacut) continue; // Apply eta cut
            if (jet_pt_smeared.at(j) > jetpt1_smeared) {
                jetpt2_smeared = jetpt1_smeared;
                jetpt1_smeared = jet_pt_smeared.at(j);
                idx2 = idx1;
                idx1 = j;
            } else if (jet_pt_smeared.at(j) > jetpt2_smeared) {
                jetpt2_smeared = jet_pt_smeared.at(j);
                idx2 = j;
            }
        }

        A_J_smeared = (jetpt1_smeared - jetpt2_smeared) / (jetpt1_smeared + jetpt2_smeared);
    
        if (i%10000 ==0){
            cout << "Processing entry " << i << "/" << nEntries << endl;
            cout << "A_J: " << A_J << endl;
            cout << "Smeared A_J: " << A_J_smeared << endl;
            cout << "Events passed selection: " << nEntriesPassSel << endl;
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
            cout << "Time elapsed: " << duration << " seconds" << endl;
        }

        if (jetpt1_smeared < 120 || jetpt2_smeared < 30) continue; // Apply jet pT cuts
        if (dPhi < dphicut) continue; // Apply delta phi cut
        if (fabs(jet_eta->at(idx1)) > etacut) continue;
        if (fabs(jet_eta->at(idx2)) > etacut) continue;

        nEntriesPassSel++;
        
        hAj->Fill(A_J, weights);
        hXj->Fill(X_J, weights);
        hdPhi->Fill(dPhi, weights);
        hAjsmeared->Fill(A_J_smeared, weights);
        for (size_t j = 0; j < jet_pt->size(); ++j) {
            hjtpt->Fill(jet_pt->at(j),weights);
            hjtptsmeared->Fill(jet_pt_smeared.at(j),weights);
            hjtpt_log->Fill(jet_pt->at(j),weights);
            hjteta->Fill(jet_eta->at(j),weights);
            hjtphi->Fill(jet_phi->at(j),weights);
        }
        hjtetmiss->Fill(ET_miss_jet,weights);

   }

    string outFolder = "/home/xirong/DijetAnalysis_2025_v3_svmit/HistRootFiles/";
    TFile* outFile = new TFile((outFolder + outFileName).c_str(), "RECREATE");
    cutInfo->Write();
    hAj->Write();
    hXj->Write();
    hdPhi->Write();
    hjtpt->Write();
    hjteta->Write();
    hjtphi->Write();
    hjtetmiss->Write();
    hjtpt_log->Write();
    hjtptsmeared->Write();
    hAjsmeared->Write();
    outFile->Close();
        
    return 0;
}
