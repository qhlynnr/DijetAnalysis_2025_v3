
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>

using namespace std;

vector<float> GetDPhi(const char* inFileName, string inTreeName, const char* branchName = "jtphi"){

    float phi1;
    float phi2;
    float dphi;
    float pi = TMath::Pi();
    int nref;
    Float_t jtphi[100];
    vector<float> jtdphi = {};

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* Tree = (TTree*)inFile->Get(inTreeName.c_str());
    Tree->SetBranchAddress(branchName,jtphi);
    Tree->SetBranchAddress("nref",&nref);

    for (int i = 0; i < Tree->GetEntries(); i++){
        Tree->GetEntry(i);
        phi1 = jtphi[0];
        phi2 = jtphi[1];
        dphi = phi1 - phi2;
        if (nref <= 2) {
            jtdphi.push_back(-999);  // or continue
            continue;
        }
        if (pi <= dphi && dphi <= 2 * pi){
            dphi = (dphi - 2*pi);
        }
        else if (-2*pi <= dphi && dphi <= -pi){
            dphi = dphi + 2*pi;
        }
        dphi = TMath::Abs(dphi);
        jtdphi.push_back(dphi);
    }
    return jtdphi;
}


int run() {

    string TreeName = "Tree";
    float pt1 = -1;
    float pt2 = -1;
    float A_J = 0;
    float ET_miss = 0;
   // float dPhi = -1;

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 600);
    TH1F* h_part_AJ = new TH1F("h_part_AJ", "Particle A_{J} of PbPb Jewel", 100, 0, 1);
    TH1F* h_part_EMiss = new TH1F("h_part_EMiss", "Particle ETMiss of PbPb Jewel", 100, 0, 200);
 
    TFile *f = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/PlottingMacros/0917PlotJewel/Jewel_pbpb_5360GeV_ptmin250_nevt26000_101325_med.root");            
    TTree *Tree = (TTree*)f->Get(TreeName.c_str());

    Tree->SetBranchAddress("pt1",&pt1);
    Tree->SetBranchAddress("pt2",&pt2);
    Tree->SetBranchAddress("A_J",&A_J);
    Tree->SetBranchAddress("ET_miss",&ET_miss);
    Long64_t nEntries = Tree->GetEntries();
    for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
            Tree->GetEntry(entrynum);
 //           if (pt1 < 50) continue;
  //          if (pt2 < 30) continue; // Apply jet pt cuts
 //           if (dPhi < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
            cout << "Event num: " << entrynum 
                    << ", ET_miss: " << ET_miss
    //               << ", dPhi: " << dPhi 
                    << ", A_J: " << A_J 
                    << ", pt1: " << pt1 
                    << ", pt2: " << pt2 << endl;
            cout << "Before filling A_J" << endl;
            h_part_AJ->Fill(A_J);
            cout << "After filling A_J" << endl;
            h_part_EMiss->Fill(ET_miss);
            cout << "After filling ETmiss" << endl;
        
        }
    f->Close();
    delete f;
    gStyle->SetOptStat(0); // Disable statistics box

    h_part_AJ->SetLineColor(kBlue+1);
    h_part_AJ->SetLineWidth(2);
    h_part_AJ->Draw("HIST");
    h_part_AJ->SetXTitle("A_{J}");

    TCanvas* c_Emiss = new TCanvas("c_Emiss", "c_Emiss", 800, 600);
    c_Emiss->cd();
    h_part_EMiss->SetLineColor(kRed+1);
    h_part_EMiss->SetLineWidth(2);
    h_part_EMiss->Draw("HIST");
    h_part_EMiss->SetXTitle("E_{T}^{miss} (GeV)");

    c_Emiss->SaveAs("Particle_ETMiss_PbPb_Jewel_115Evt.png");
    canvas->SaveAs("Particle_AJ_PbPb_Jewel_115Evt.png");
    return 0;
}