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
    bool logx = false,
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
    if (logx) canvas_temp->SetLogx(1);
    
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
int DrawJet(){
    gStyle->SetOptStat(0);

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 800);

    // Define histogram parameters
    
    //TH1F* h_part_EMiss = new TH1F("h_part_EMiss", "Particle ETMiss of PbPb Jewel", 100, 0, 200);
    TFile *f = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/Jewelpbpb_5360GeV_ptm250_10Kevt_020426_C5_dum_med_jewelbkg_new.root");         
    TFile *f_vac = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/020126DebugDummies//Jewelpbpb_5360GeV_ptm250_10Kevt_020126_C5_vac_11050627.root");
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

        for (int j=0; j < jet_pt_raw->size(); j++){
            h_jetpt_raw->Fill(jet_pt_raw->at(j));
        }
        for (int j=0; j < jet_pt_thermal->size(); j++){
            h_jetpt_thermal->Fill(jet_pt_thermal->at(j));
        }
        for (int j=0; j < jet_pt->size(); j++){
            h_jetpt->Fill(jet_pt->at(j));
        }
        for (int j=0; j < jet_pt_vac->size(); j++){
            h_jetpt_vac->Fill(jet_pt_vac->at(j));
        }

        // Fill jet ETMiss histograms
        h_jetEMiss->Fill(ET_miss_jet);
        h_jetEMiss_raw->Fill(ET_miss_jet_raw);
        h_jetEMiss_thermal->Fill(ET_miss_jet_thermal);
    }

    vector<TH1F*> hists = {h_jetEMiss, h_jetEMiss_raw, h_jetEMiss_thermal
    };
    vector<string> labels = {
        "ETMiss (Corrected)", 
        "ETMiss (Raw Only)", 
        "ETMiss (Thermal Only)"
    };

    int nEvents = JetTree->GetEntries();
    int nEvents_vac = JetTree_vac->GetEntries();

    // create normalized clones for raw and vacuum pT
    TH1F* h_jetpt_raw_norm = (TH1F*)h_jetpt_raw->Clone("h_jetpt_raw_norm");
    if (h_jetpt_raw_norm->Integral() > 0) h_jetpt_raw_norm->Scale(1.0 / nEvents);

    TH1F* h_jetpt_vac_norm = (TH1F*)h_jetpt_vac->Clone("h_jetpt_vac_norm");
    if (h_jetpt_vac_norm->Integral() > 0) h_jetpt_vac_norm->Scale(1.0 / nEvents_vac);

    // also prepare a normalized ratio for the nominal jet_pt if desired
    TH1F* h_jetpt_norm = (TH1F*)h_jetpt->Clone("h_jetpt_norm");
    h_jetpt_norm->Scale(1.0 / nEvents);

    TH1F *h_ratio_raw = (TH1F*)h_jetpt_raw_norm->Clone("h_ratio_raw_norm");
    h_ratio_raw->Divide(h_jetpt_vac_norm);

    TH1F *h_ratio_pp = (TH1F*)h_jetpt_vac_norm->Clone("h_ratio_pp_norm");
    h_ratio_pp->Divide(h_jetpt_vac_norm); // should be all ones

    TH1F *h_ratio = (TH1F*)h_jetpt_norm->Clone("h_ratio");
    h_ratio->Divide(h_jetpt_vac_norm);

    // cleanup temporary clones you won't keep
    delete h_jetpt_raw_norm;
    delete h_jetpt_vac_norm;
    delete h_jetpt_norm; // only if you created h_ratio from it; otherwise keep as needed

    h_ratio_raw->GetXaxis()->SetRangeUser(25,1000);
    h_ratio->GetXaxis()->SetRangeUser(25,1000);
    h_ratio_pp->GetXaxis()->SetRangeUser(25,1000);
    
    plotHistograms({h_jetpt, h_jetpt_raw, h_jetpt_thermal,h_jetpt_vac}, 
                    {"Jet p_{T}", " Jet p_{T} Raw Only", " Jet p_{T} Thermal Only"," Jet p_{T} Vacuum Only"}, 
                    "jet_pt", "Jet p_{T} (GeV)", true, 
                    "0205_PbPb_Jewel_10KEvt_Dummies_JetPt_Comparison",
                    "Comparison of Jet p_{T} in PbPb Jewel After Clustering 10K Fixed Debug Seed",
                    false,true);
    plotHistograms({h_ratio, h_ratio_raw, h_ratio_pp}, 
                    {"R_{AA}", " R_{AA} Raw Only", " R_{AA} Vacuum Only"}, 
                    "jet_pt_ratio", "Ratio of Jet p_{T} to Vacuum Jet p_{T}", false, 
                    "0205_PbPb_Jewel_10KEvt_Dummies_JetPt_Ratio_Comparison",
                    "Comparison of Ratio of Jet p_{T} to Vacuum Jet p_{T} in PbPb Jewel After Clustering 10K Fixed Debug Seed",
                    true,false,
                    1000,2.0);
    /*plotHistograms({h_jetpt, h_jetpt_raw, h_jetpt_thermal}, 
                    {"Jet p_{T}", " Jet p_{T} Raw Only", " Jet p_{T} Thermal Only"}, 
                    "jet_pt", "Jet p_{T} (GeV)", false, 
                    "0204_PbPb_Jewel_10KEvt_Dummies_JetPt_Comparison",
                    "Comparison of Jet p_{T} in PbPb Jewel After Clustering 10K Fixed Debug Seed",
                    false);*/

   /* plotHistograms({h_jetEMiss, h_jetEMiss_raw, h_jetEMiss_thermal}, 
                    {"Jet ETMiss", "Jet ETMiss Raw Only", "Jet ETMiss Thermal Only"}, 
                    "jet_ETMiss", "Jet E_{T}^{Miss} (GeV)", true, 
                    "0204_PbPb_Jewel_10KEvt_Dummies_JetETMiss_Comparison",
                    "Comparison of Jet E_{T}^{Miss} in PbPb Jewel After Clustering 10K Fixed Debug Seed");*/

 /*   plotHistograms({h_part_mult_nodum,h_part_mult_dum, h_part_mult_norec}, 
                    {"Jewel PbPb No WriteDummies", "Jewel PbPb WriteDummies", "Jewel PbPb WriteDummies without Recoil"}, 
                    "particle_multiplicity", "Particle Multiplicity", 
                    true, 
                    "0203_PbPb_Jewel_10KEvt_Dummies_Multiplicity_Comparison",
                    "Comparison of Particle Multiplicity in PbPb Jewel 10K Fixed Debug Seed");
*/
/*    plotHistograms({h_part_ETMiss_nodum,h_part_ETMiss_dum}, 
                    {"Jewel PbPb No WriteDummies", "Jewel PbPb WriteDummies"}, 
                    "particle_ETMiss", "E_{T}^{Miss} (GeV)", 
                    true, 
                    "0201_PbPb_Jewel_10KEvt_Dummies_Comparison_twoplots");

    plotHistograms({h_part_status_nodum, h_part_status_dum}, 
                    {"Jewel PbPb No WriteDummies", "Jewel PbPb WriteDummies"}, 
                    "particle_status", "Particle Status", true, "0201_PbPb_Jewel_10KEvt_Status_Comparison_");
*/
    f->Close();
    return 0;
}