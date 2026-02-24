#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

int main(){

    const char* ppfilename = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021826_JewelWithWeights/Jewelpbpb_5360GeV_ptm250_10Kevt_021026_C5_vac_1105062_weightsT_nobkg.root";
    const char* pbpbfilename = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_dum_med_jewelbkg.root";

    TH1F *hRaa = new TH1F("hRaa", "Raa Distribution", 100, 0, 1.5);
    TH1F *h_jetpt = new TH1F("h_jetpt", " p_{T} of Jet in 10K PbPb Event using Jewel", 100, 0, 1000);
    TH1F *h_particlept = new TH1F("h_particlept", " p_{T} of Particle in 10K PbPb Event using Jewel", 100, 0, 500);
    TH1F *h_jetpt_vac = new TH1F("h_jetpt_vac", " p_{T} of Jet in 50K PbPb Vacuum Event using Jewel", 100, 0, 500);
    TH1F *h_particlept_vac = new TH1F("h_particlept_vac", " p_{T} of Particle in 50K PbPb Vacuum Event using Jewel", 100, 0, 500);

    TFile* f_pp = TFile::Open(ppfilename, "READ");
    TFile* f_pbpb = TFile::Open(pbpbfilename, "READ");
    if (!f_pp || f_pp->IsZombie()) {
        std::cerr << "Error: Could not open file " << ppfilename << std::endl;
        return 1;
    }
    if (!f_pbpb || f_pbpb->IsZombie()) {
        std::cerr << "Error: Could not open file " << pbpbfilename << std::endl;
        return 1;
    }

    TTree* t_pp_Particle = (TTree*)f_pp->Get("ParticleTree");
    TTree* t_pbpb_Particle = (TTree*)f_pbpb->Get("ParticleTree");
    TTree* t_pp_Jet = (TTree*)f_pp->Get("JetTree");
    TTree* t_pbpb_Jet = (TTree*)f_pbpb->Get("JetTree");

    if (!t_pp_Particle || !t_pbpb_Particle) {
        std::cerr << "Error: Could not find ParticleTree in one of the files." << std::endl;
        return 1;
    }
    if (!t_pp_Jet || !t_pbpb_Jet) {
        std::cerr << "Error: Could not find JetTree in one of the files." << std::endl;
        return 1;
    }

    float ppweight;
    vector<float>* pp_particle_pt = nullptr;
    vector<float>* pp_jet_pt = nullptr;
    vector<float>* pbpb_particle_pt = nullptr;
    vector<float>* pbpb_jet_pt = nullptr;

    t_pp_Particle->SetBranchAddress("weights", &ppweight);
    t_pp_Particle->SetBranchAddress("pt", &pp_particle_pt);
    t_pp_Jet->SetBranchAddress("jet_pt", &pp_jet_pt);
    t_pbpb_Particle->SetBranchAddress("pt", &pbpb_particle_pt);
    t_pbpb_Jet->SetBranchAddress("jet_pt", &pbpb_jet_pt);
    int nevents_pp = t_pp_Particle->GetEntries();
    int nevents_pbpb = t_pbpb_Particle->GetEntries();
    
    double sumW = 0.0;
    gStyle->SetOptStat(0);
    gPad->SetTicks();

    for (Long64_t i = 0; i < t_pp_Particle->GetEntries(); ++i) {
        t_pp_Particle->GetEntry(i);
            for (size_t j = 0; j < pp_particle_pt->size(); ++j) {
                h_particlept_vac->Fill(pp_particle_pt->at(j), ppweight);
            }
        sumW += ppweight;
    }

    for (Long64_t i = 0; i < t_pbpb_Particle->GetEntries(); ++i) {
        t_pbpb_Particle->GetEntry(i);
        for (size_t j = 0; j < pbpb_particle_pt->size(); ++j) {
            h_particlept->Fill(pbpb_particle_pt->at(j));
        }

    }

    h_particlept_vac->Scale(1.0 / sumW);
    h_particlept->Scale(1.0 / nevents_pbpb);
    TCanvas* c2 = new TCanvas("c2", "Particle pT Distribution in Vacuum", 800, 600);
    c2->cd();
    h_particlept_vac->SetTitle("Particle p_{T} Distribution in Vacuum");
    h_particlept_vac->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_particlept_vac->GetYaxis()->SetTitle("dN/dp_{T}");
    h_particlept_vac->SetMarkerStyle(20);
    h_particlept_vac->Draw("HIST");
    h_particlept->SetTitle("Particle p_{T} Distribution in Medium");
    h_particlept->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_particlept->GetYaxis()->SetTitle("dN/dp_{T}");
    h_particlept->SetMarkerStyle(20);
    h_particlept->SetLineColor(kRed);
    h_particlept->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(h_particlept_vac, "Vacuum", "l");
    legend->AddEntry(h_particlept, "Medium", "l");
    legend->Draw();
    
    c2->SetLogy();

    TCanvas* c1 = new TCanvas("c1", "Raa Distribution particle", 800, 800);
    TH1D *h_RAA = (TH1D*)h_particlept->Clone("h_RAA");
    h_RAA->Divide(h_particlept_vac);

    c1->cd();
    h_RAA->SetTitle("Raa Distribution of Particle p_{T} in PbPb vs pp");
    h_RAA->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_RAA->GetYaxis()->SetTitle("Raa");
    h_RAA->GetXaxis()->SetRangeUser(250, 100);
    h_RAA->GetYaxis()->SetRangeUser(0, 1.5);
    h_RAA->SetMarkerStyle(20);
    h_RAA->Draw("E1");
    c1->SaveAs("Raa_Distribution_50Kpbpb.png");
    c2->SaveAs("Particle_pT_Distribution_Comparison_50Kpbpb.png");



    return 0;
}
