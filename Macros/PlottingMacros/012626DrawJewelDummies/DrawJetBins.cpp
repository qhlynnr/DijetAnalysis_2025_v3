#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

using namespace std;

void plotHistograms(vector<TH1F*> histograms,
    vector<string> labels, 
    const string& varname, 
    const string& xtitle, 
    bool norm = false, 
    const string& outputfilename = "",
    const string& histtitle = "",
    bool logy = false,
    double xmax = 0.0,
    double ymax = 0.0
    ) {
    
    if (histograms.size() != labels.size() || histograms.empty()) {
        cout << "Error: histograms and labels vectors must have the same size and be non-empty!" << endl;
        return;
    }
    
    // Work on clones when normalizing so original histograms are not modified
    vector<TH1F*> hists = histograms;
    bool cloned = false;
    if (norm) {
        hists.clear();
        for (size_t i = 0; i < histograms.size(); i++) {
            TH1F* clone = (TH1F*)histograms[i]->Clone((string(histograms[i]->GetName()) + "_norm").c_str());
            double integral = clone->Integral();
            if (integral > 0) clone->Scale(1.0 / integral);
            hists.push_back(clone);
        }
        cloned = true;
    }

    // Determine y-axis maximum: use provided ymax if >0, otherwise compute from histograms
    double maxy = 0.0;
    if (ymax > 0.0) {
        maxy = ymax;
    } else {
        for (auto h : hists) {
            double max_val = h->GetMaximum();
            if (max_val > maxy) maxy = max_val;
        }
        maxy *= 1.2; // Add margin
        if (maxy == 0.0) maxy = 1.0; // avoid zero max
    }

    TCanvas* canvas_temp = new TCanvas(("c_" + varname).c_str(), "canvas", 800, 800);
    canvas_temp->SetLeftMargin(0.14);
    canvas_temp->SetBottomMargin(0.12);
    canvas_temp->SetTopMargin(0.08);
    canvas_temp->SetRightMargin(0.05);

    if (logy) canvas_temp->SetLogy(1);
    
    // Color palette for multiple histograms
    vector<int> colors = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+1, kOrange+1, kViolet+1, kSpring+1};
    
    // Optionally set x-axis upper limit if requested
    double xmin = hists[0]->GetXaxis()->GetXmin();
    if (xmax > 0.0) {
        for (auto h : hists) {
            h->GetXaxis()->SetRangeUser(xmin, xmax);
        }
    }

    // Draw first histogram to set up axes
    hists[0]->SetLineColor(colors[0]);
    hists[0]->SetLineWidth(2);
    if (!histtitle.empty()) hists[0]->SetTitle(histtitle.c_str());
    hists[0]->SetXTitle(xtitle.c_str());
    hists[0]->SetYTitle(norm ? "Normalized entries" : "Entries");
    hists[0]->SetMaximum(maxy);
    if (logy) {
        // ensure a positive non-zero minimum for log scale
        double miny = maxy > 0.0 ? maxy * 1e-6 : 1e-6;
        if (miny <= 0) miny = 1e-6;
        hists[0]->SetMinimum(miny);
    }
    hists[0]->Draw("HIST");

    // Draw remaining histograms
    for (size_t i = 1; i < hists.size(); i++) {
        int color = colors[i % colors.size()]; // Cycle through colors if more histograms than colors
        hists[i]->SetLineColor(color);
        hists[i]->SetLineWidth(2);
        hists[i]->Draw("HIST SAME");
    }

    // Create legend - adjust size based on number of entries
    double legend_height = 0.03 * hists.size() + 0.02; // Dynamic height
    double legend_y2 = 0.9;
    double legend_y1 = legend_y2 - legend_height;
    TLegend* legend = new TLegend(0.5, legend_y1, 0.85, legend_y2);
    
    for (size_t i = 0; i < hists.size(); i++) {
        legend->AddEntry(hists[i], labels[i].c_str(), "l");
    }
    legend->SetBorderSize(0);
    legend->SetTextSize(0.025);
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.025);

    // append suffix for log and normalization if requested
    string suffix = "";
    if (logy) suffix += "_log";
    if (norm)  suffix += "_norm";
    string filename = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/0126Dummies/" + outputfilename + varname + suffix + ".png";
    canvas_temp->SaveAs(filename.c_str());

    delete legend;
    delete canvas_temp;
    
    // Clean up cloned histograms
    if (cloned) {
        for (auto h : hists) {
            delete h;
        }
    }
};
int DrawJetBins(){
    gStyle->SetOptStat(0);

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 800);

    // Define histogram parameters
    
    //TH1F* h_part_EMiss = new TH1F("h_part_EMiss", "Particle ETMiss of PbPb Jewel", 100, 0, 200);
    TFile *f = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/Jewelpbpb_5360GeV_ptm250_10Kevt_020426_C5_dum_med_jewelbkg_new.root");         
    TFile *f_vac = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/JewelRootOutput/Jewel_pbpb_5360GeV_ptmin250_nevt10K_012026_vac_withdphi.root");
    TTree *ParticleTree= (TTree*)f->Get("ParticleTree");
    TTree *JetTree = (TTree*)f->Get("JetTree");
    TTree *JetTree_vac = (TTree*)f_vac->Get("JetTree");

    //Particle Level
    std::vector<float>* e = nullptr;
    std::vector<float>* px = nullptr;
    std::vector<float>* py = nullptr;
    std::vector<float>* pz = nullptr;
    std::vector<float>* pt = nullptr;
    std::vector<int>* pdgId= nullptr;
    double ET_miss_part = 0.0;
    ParticleTree->SetBranchAddress("e", &e);
    ParticleTree->SetBranchAddress("px", &px);
    ParticleTree->SetBranchAddress("py", &py);
    ParticleTree->SetBranchAddress("pz", &pz);
    ParticleTree->SetBranchAddress("pt", &pt);
    ParticleTree->SetBranchAddress("pdgId", &pdgId);


    //Jet Level
    double A_J=0, X_J=0, jetpt1=0, jetpt2=0, jetphi1=0, jetphi2=0, ET_miss_jet=0, dPhi=0;
    std::vector<float>* jet_pt = nullptr;
    std::vector<float>* jet_pt_thermal = nullptr;
    std::vector<float>* jet_pt_raw = nullptr;
    std::vector<float>* jet_eta = nullptr;
    std::vector<float>* jet_phi = nullptr;
    double ET_miss_jet_raw = 0;
    double ET_miss_jet_thermal = 0;

    std::vector<float>* jet_pt_vac = nullptr;

    int bins = 100;
    double xmax = 1000.0;
    TH1F* h_jetpt_raw = new TH1F("h_jetpt_raw", " Raw p_{T} of Jet in 10K PbPb Event using Jewel", bins, 0, xmax);
    TH1F* h_jetpt_thermal = new TH1F("h_jetpt_thermal", " Thermal p_{T} of Jet in 10K PbPb Event using Jewel", bins, 0, xmax);
    TH1F* h_jetpt = new TH1F("h_jetpt", " p_{T} of Jet in 10K PbPb Event using Jewel", bins, 0, xmax);
    TH1F* h_jetpt_vac = new TH1F("h_jetpt_vac", " p_{T} of Jet in 10K PbPb Vacuum Event using Jewel", bins, 0, xmax);

    TH1F* h_jetEMiss = new TH1F("h_jetETMiss", " ETMiss of Jet in 10K PbPb Event using Jewel", 120, -5, 100);
    TH1F* h_jetEMiss_raw = new TH1F("h_jetETMiss_raw", " Raw ETMiss of Jet in 10K PbPb Event using Jewel", 120, -5, 100);
    TH1F* h_jetEMiss_thermal = new TH1F("h_jetETMiss_thermal", " Thermal ETMiss of Jet in 10K PbPb Event using Jewel", 120, -5, 100);

    TH1F* h_jetptreconstructed_bin1 = new TH1F("h_jetptreconstructed_bin1", " p_{T} of Reconstructed Jet in 10K PbPb Event using Jewel (Bin 1: 0-50 GeV)", bins, 0, xmax);
    TH1F* h_jetptreconstructed_bin2 = new TH1F("h_jetptreconstructed_bin2", " p_{T} of Reconstructed Jet in 10K PbPb Event using Jewel (Bin 2: 50-100 GeV)", bins, 0, xmax);
    TH1F* h_jetptreconstructed_bin3 = new TH1F("h_jetptreconstructed_bin3", " p_{T} of Reconstructed Jet in 10K PbPb Event using Jewel (Bin 3: 100-200 GeV)", bins, 0, xmax);  
    TH1F* h_jetptreconstructed_bin4 = new TH1F("h_jetptreconstructed_bin4", " p_{T} of Reconstructed Jet in 10K PbPb Event using Jewel (Bin 4: 200-300 GeV)", bins, 0, xmax);  

    TH1F* h_jetptraw_bin1 = new TH1F("h_jetptraw_bin1", " p_{T} of Raw Jet in 10K PbPb Event using Jewel (Bin 1: 0-50 GeV)", bins, 0, xmax);
    TH1F* h_jetptraw_bin2 = new TH1F("h_jetptraw_bin2", " p_{T} of Raw Jet in 10K PbPb Event using Jewel (Bin 2: 50-100 GeV)", bins, 0, xmax);
    TH1F* h_jetptraw_bin3 = new TH1F("h_jetptraw_bin3", " p_{T} of Raw Jet in 10K PbPb Event using Jewel (Bin 3: 100-200 GeV)", bins, 0, xmax);  
    TH1F* h_jetptraw_bin4 = new TH1F("h_jetptraw_bin4", " p_{T} of Raw Jet in 10K PbPb Event using Jewel (Bin 4: 200-300 GeV)", bins, 0, xmax);

    TH1F* h_jetptthermal_bin1 = new TH1F("h_jetptthermal_bin1", " p_{T} of Thermal Jet in 10K PbPb Event using Jewel (Bin 1: 0-50 GeV)", bins, 0, xmax);
    TH1F* h_jetptthermal_bin2 = new TH1F("h_jetptthermal_bin2", " p_{T} of Thermal Jet in 10K PbPb Event using Jewel (Bin 2: 50-100 GeV)", bins, 0, xmax);
    TH1F* h_jetptthermal_bin3 = new TH1F("h_jetptthermal_bin3", " p_{T} of Thermal Jet in 10K PbPb Event using Jewel (Bin 3: 100-200 GeV)", bins, 0, xmax);  
    TH1F* h_jetptthermal_bin4 = new TH1F("h_jetptthermal_bin4", " p_{T} of Thermal Jet in 10K PbPb Event using Jewel (Bin 4: 200-300 GeV)", bins, 0, xmax);

    JetTree->SetBranchAddress("A_J",&A_J);
    JetTree->SetBranchAddress("X_J",&X_J);
    JetTree->SetBranchAddress("jetpt1",&jetpt1);
    JetTree->SetBranchAddress("jetpt2",&jetpt2);
    JetTree->SetBranchAddress("jetphi1",&jetphi1);
    JetTree->SetBranchAddress("jetphi2",&jetphi2);
    JetTree->SetBranchAddress("dPhi",&dPhi);
    JetTree->SetBranchAddress("jet_eta",&jet_eta);
    JetTree->SetBranchAddress("jet_phi",&jet_phi);

    JetTree->SetBranchAddress("ET_miss_jet",&ET_miss_jet);
    JetTree->SetBranchAddress("ET_miss_jet_raw",&ET_miss_jet_raw);
    JetTree->SetBranchAddress("ET_miss_jet_thermal",&ET_miss_jet_thermal);

    JetTree->SetBranchAddress("jet_pt",&jet_pt);
    JetTree->SetBranchAddress("jet_pt_thermal",&jet_pt_thermal);
    JetTree->SetBranchAddress("jet_pt_raw",&jet_pt_raw);

    JetTree_vac->SetBranchAddress("jet_pt",&jet_pt_vac);

    for(Long64_t i=0; i < JetTree->GetEntries(); i++){
        JetTree->GetEntry(i);
        ParticleTree->GetEntry(i);
        JetTree_vac->GetEntry(i);
        for (int j=0; j<jet_pt_raw->size(); j++){
            if (jet_pt->at(j) < 50){
                h_jetptreconstructed_bin1->Fill(jet_pt_raw->at(j) - jet_pt_thermal->at(j));
                h_jetptraw_bin1->Fill(jet_pt_raw->at(j));
                h_jetptthermal_bin1->Fill(jet_pt_thermal->at(j));
            }
            else if (jet_pt->at(j) >= 50 && jet_pt->at(j) < 100){
                h_jetptreconstructed_bin2->Fill(jet_pt_raw->at(j) - jet_pt_thermal->at(j));
                h_jetptraw_bin2->Fill(jet_pt_raw->at(j));
                h_jetptthermal_bin2->Fill(jet_pt_thermal->at(j));
            }
            else if (jet_pt->at(j) >= 100 && jet_pt->at(j) < 200){
                h_jetptreconstructed_bin3->Fill(jet_pt_raw->at(j) - jet_pt_thermal->at(j));
                h_jetptraw_bin3->Fill(jet_pt_raw->at(j));
                h_jetptthermal_bin3->Fill(jet_pt_thermal->at(j));
            }
            else if (jet_pt->at(j) >= 200 && jet_pt->at(j) < 300){
                h_jetptreconstructed_bin4->Fill(jet_pt_raw->at(j) - jet_pt_thermal->at(j));
                h_jetptraw_bin4->Fill(jet_pt_raw->at(j));
                h_jetptthermal_bin4->Fill(jet_pt_thermal->at(j));
        }
    }
}
    TCanvas *c = new TCanvas("c","",1000,800);
    c->Divide(2,2);   // 2 columns, 2 rows

    c->cd(1);
    h_jetptreconstructed_bin1->Draw("HIST");
    h_jetptreconstructed_bin1->SetLineColor(kRed);
    h_jetptreconstructed_bin1->SetLineWidth(2);
    h_jetptreconstructed_bin1->Draw("HIST");
    h_jetptraw_bin1->Draw("HIST SAME");
    h_jetptthermal_bin1->Draw("HIST SAME");
    TLegend *leg1 = new TLegend(0.6,0.6,0.88,0.88);
    leg1->AddEntry(h_jetptreconstructed_bin1,"Reconstructed Jet p_{T}","l");
    leg1->AddEntry(h_jetptraw_bin1,"Raw Jet p_{T}","l");
    leg1->AddEntry(h_jetptthermal_bin1,"Thermal Jet p_{T}","l");
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.025);
    leg1->Draw();

    c->cd(2);
    h_jetptreconstructed_bin2->Draw("HIST");
    h_jetptreconstructed_bin2->SetLineColor(kRed);
    h_jetptreconstructed_bin2->SetLineWidth(2);
    h_jetptreconstructed_bin2->Draw("HIST");
    h_jetptraw_bin2->Draw("HIST SAME");
    h_jetptthermal_bin2->Draw("HIST SAME");
    TLegend *leg2 = new TLegend(0.6,0.6,0.88,0.88);
    leg2->AddEntry(h_jetptreconstructed_bin2,"Reconstructed Jet p_{T}","l");
    leg2->AddEntry(h_jetptraw_bin2,"Raw Jet p_{T}","l");
    leg2->AddEntry(h_jetptthermal_bin2,"Thermal Jet p_{T}","l");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.025);
    leg2->Draw();

    c->cd(3);
    h_jetptreconstructed_bin3->Draw("HIST");
    h_jetptreconstructed_bin3->SetLineColor(kRed);
    h_jetptreconstructed_bin3->SetLineWidth(2);
    h_jetptreconstructed_bin3->Draw("HIST");
    h_jetptraw_bin3->Draw("HIST SAME");
    h_jetptthermal_bin3->Draw("HIST SAME");
    TLegend *leg3 = new TLegend(0.6,0.6,0.88,0.88);
    leg3->AddEntry(h_jetptreconstructed_bin3,"Reconstructed Jet p_{T}","l");
    leg3->AddEntry(h_jetptraw_bin3,"Raw Jet p_{T}","l");
    leg3->AddEntry(h_jetptthermal_bin3,"Thermal Jet p_{T}","l");
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.025);
    leg3->Draw();

    c->cd(4);
    h_jetptreconstructed_bin4->Draw("HIST");
    h_jetptreconstructed_bin4->SetLineColor(kRed);
    h_jetptreconstructed_bin4->SetLineWidth(2);
    h_jetptreconstructed_bin4->Draw("HIST");
    h_jetptraw_bin4->Draw("HIST SAME");
    h_jetptthermal_bin4->Draw("HIST SAME");
    TLegend *leg4 = new TLegend(0.6,0.6,0.88,0.88);
    leg4->AddEntry(h_jetptreconstructed_bin4,"Reconstructed Jet p_{T}","l");
    leg4->AddEntry(h_jetptraw_bin4,"Raw Jet p_{T}","l");
    leg4->AddEntry(h_jetptthermal_bin4,"Thermal Jet p_{T}","l");
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.025);
    leg4->Draw();

    c->SaveAs("/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/0126Dummies/020626JetPtComparison.png");
    f->Close();
    return 0;
}