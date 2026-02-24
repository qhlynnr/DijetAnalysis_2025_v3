#include <TFile.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <iostream>
#include <string>

// Use the standard namespace
using namespace std;

void Draw2D() {
    const char* filename = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021326_JewelBkgSubCheck/Jewelpbpb_5360GeV_ptm250_10Kevt_020426_C5_dum_med_jewelbkg_11050627_new.root";
    const char* filename2 = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021326_JewelBkgSubCheck/Jewelpbpb_5360GeV_ptm250_10Kevt_020126_C5_dum_med_11050627_nobkgsub.root";
    TFile* file = TFile::Open(filename, "READ");
    TFile* file2 = TFile::Open(filename2, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    else {
        cout << "File Opened" << endl;
    }

    gStyle->SetOptStat(0); // Disable statistics box

    TTree* jetTree = (TTree*)file->Get("JetTree");
    TTree* particleTree = (TTree*)file->Get("ParticleTree");
    TTree* jetTree_nobkg = (TTree*)file2->Get("JetTree");
    vector<float>* phi = nullptr;
    vector<float>* eta = nullptr;
    vector<float>* pt = nullptr;
    vector<float>* phi_nobkg = nullptr;
    vector<float>* eta_nobkg = nullptr;
    vector<float>* pt_nobkg = nullptr;
    Double_t jetpt1, jetpt1_nobkg;

    jetTree->SetBranchAddress("jet_phi", &phi);
    jetTree->SetBranchAddress("jet_eta", &eta);
    jetTree->SetBranchAddress("jet_pt", &pt);
    jetTree->SetBranchAddress("jetpt1", &jetpt1);
    jetTree_nobkg->SetBranchAddress("jet_phi", &phi_nobkg);
    jetTree_nobkg->SetBranchAddress("jet_eta", &eta_nobkg);
    jetTree_nobkg->SetBranchAddress("jet_pt", &pt_nobkg);
    jetTree_nobkg->SetBranchAddress("jetpt1", &jetpt1_nobkg);

    TCanvas* c1 = new TCanvas("c1", "Jet Distribution of 1 Event", 800, 600);
    TCanvas* c2 = new TCanvas("c2", "Jet Distribution of 1 Event without Bkg Subtraction", 800, 600);
    TH2D* graph = new TH2D ("h2","Jet Distribution of 1 Event;Phi;Eta", 50, -3.14, 3.14, 50,-2,2);
    TH2D* graph_nobkg = new TH2D ("h2_nobkg","Jet Distribution of 1 Event without Bkg Subtraction;Phi;Eta", 50, -3.14, 3.14, 50,-2,2);
    //    graph->GetXaxis()->SetRangeUser(-5, 5);
    // Change all trees
    Long64_t nEntries = jetTree->GetEntries();

    for (Long64_t i = 0; i < 10; ++i) {
        int entryNumber = i; // Change this to the desired entry number
        jetTree->GetEntry(entryNumber); // Get the first entry
        jetTree_nobkg->GetEntry(entryNumber); // Get the first entry
        for (size_t i = 0; i < phi->size(); ++i) {
            graph->Fill(phi->at(i), eta->at(i),pt->at(i));
            graph_nobkg->Fill(phi_nobkg->at(i), eta_nobkg->at(i),pt_nobkg->at(i));
        }

        c1->cd();
        graph->Draw("LEGO2Z");
        c1->SetTheta(45.);
        c1->SetPhi(60.);
        
        c2->cd();
        graph_nobkg->Draw("LEGO2Z");
        c2->SetTheta(45.);
        c2->SetPhi(60.);

        gStyle->SetPalette(kBird);
        graph->SetContour(1000);
        graph_nobkg->SetContour(1000);

        c1->SaveAs(Form("/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/021326Plot2DBkgCheck/Jewelpbpb_SingleEventDistribution_Jets_entry%d_afterBkgSubtraction.png", entryNumber));
        c2->SaveAs(Form("/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/021326Plot2DBkgCheck/Jewelpbpb_SingleEventDistribution_Jets_entry%d_beforeBkgSubtraction.png", entryNumber));

        graph->Reset();
        graph_nobkg->Reset();
    }
}

