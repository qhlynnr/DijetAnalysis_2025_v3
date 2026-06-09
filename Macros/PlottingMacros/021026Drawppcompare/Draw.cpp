// ComparePythiaVsJewel.C
// Usage (ROOT):
//   root -l -q 'ComparePythiaVsJewel.C("pythia.root","jewel.root","out_png")'
//
// It will write PNGs into out_png/
//
// Notes:
// - Jet-level: uses JetTree branches: jetpt1, jetpt2, X_J, A_J, dPhi, jet_eta (if present) OR jet_eta vector (if present).
// - Particle-level: uses ParticleTree branches: pt (vector<float>), eta (vector<float>) if present.
// - All plots are normalized to unit area for shape comparison.

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TStyle.h"
#include <string>
#include <iostream>

static void StyleHist(TH1* h, int color, int marker){
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetLineWidth(2);
}

static void Normalize(TH1* h){
  double integral = h->Integral(0, h->GetNbinsX()+1);
  if (integral > 0) h->Scale(1.0 / integral);
}


#include <vector>

static void DrawAndSave(const std::vector<TH1*>& hs,
                        const std::vector<std::string>& labels,
                        const std::string& outPng,
                        const std::string& xTitle,
                        const std::string& yTitle = "Normalized entries",
                        const std::string& title = "",
                        float xMin = -999, float xMax = -999,
                        float yMin = -999, float yMax = -999,
                        bool logx = false,
                        bool logy = false) {
  if (hs.empty()) return;

  // basic sanity on labels
  std::vector<std::string> lab = labels;
  if (lab.size() != hs.size()) {
    lab.resize(hs.size());
    for (size_t i = 0; i < hs.size(); ++i) {
      if (i < labels.size()) lab[i] = labels[i];
      else lab[i] = std::string("hist") + std::to_string(i+1);
    }
  }

  TCanvas* c = new TCanvas("c", "canvas", 800, 800);
  c->SetTicks(1,1);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.08);
  c->SetRightMargin(0.05);

  if (logx) c->SetLogx();
  if (logy) c->SetLogy();

  // set titles on first histogram
  TH1* h0 = hs[0];
  if (!title.empty()) h0->SetTitle(title.c_str());
  h0->GetXaxis()->SetTitle(xTitle.c_str());
  h0->GetYaxis()->SetTitle(yTitle.c_str());

  // determine y range using all histograms
  double maxy = 0;
  for (auto h : hs) if (h) maxy = std::max(maxy, h->GetMaximum());
  if (maxy > 0) h0->SetMaximum(1.25 * maxy);

  if (xMin != -999 && xMax != -999) h0->GetXaxis()->SetRangeUser(xMin, xMax);
  if (yMin != -999 && yMax != -999) h0->GetYaxis()->SetRangeUser(yMin, yMax);
  // draw histograms (assume they are already styled)
  bool first = true;
  for (auto h : hs) {
    if (!h) continue;
    if (first) { h->Draw("HIST"); first = false; }
    else       { h->Draw("HIST SAME"); }
  }

  // legend placement: adapt height to number of entries
  const size_t n = hs.size();
  double ly2 = 0.9;
  double ly1 = ly2 - 0.05 * std::max<size_t>(n, 1);
  if (ly1 < 0.1) ly1 = 0.1;
  TLegend* leg = new TLegend(0.6, ly1, 0.9, ly2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (size_t i = 0; i < hs.size(); ++i) {
    if (!hs[i]) continue;
    leg->AddEntry(hs[i], lab[i].c_str(), "l");
  }
  leg->Draw();

  c->SaveAs(outPng.c_str());
}

void Draw() {
    // Hard-coded input files and output directory
    const char* pythiaFile = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021026_ppComparison/Pythia8Jets_pp_5360GeV_HardQCD_pTHatMin250_10000Evts_withjets.root";
    const char* jewelFile  = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021026_ppComparison/Jewelpbpb_5360GeV_ptm250_10Kevt_021026_C5_vac_11050627_withjets.root";
    const char* outDir     = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/021126ppComparePythiaVsJewel_noweight/";
    const char* jewelFile_weights = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021026_ppComparison/Jewelpbpb_5360GeV_ptm250_10Kevt_021026_C5_vac_1105062_weightsT.root";

    gStyle->SetOptStat(0);

    gSystem->mkdir(outDir, /*recursive=*/kTRUE);

    TFile fP(pythiaFile, "READ");
    TFile fJ(jewelFile,  "READ");
    TFile fJ_weights(jewelFile_weights, "READ");
    if (fP.IsZombie() || fJ.IsZombie()) {
        std::cerr << "ERROR: could not open one of the files.\n";
        return;
    }

    // --- Grab trees (try common names)
    TTree* pJet = (TTree*)fP.Get("JetTree");
    TTree* jJet = (TTree*)fJ.Get("JetTree");
    TTree* pPar = (TTree*)fP.Get("ParticleTree");
    TTree* jPar = (TTree*)fJ.Get("ParticleTree");
    TTree* jWeightPar = (TTree*)fJ_weights.Get("ParticleTree"); // weights stored in ParticleTree in this file

    if (!pJet || !jJet) std::cerr << "WARNING: JetTree missing in one of the files.\n";
    if (!pPar || !jPar) std::cerr << "WARNING: ParticleTree missing in one of the files.\n";

    // =========================
    // Jet-level comparisons
    // =========================
    if (pJet && jJet) {
        // jetpt1
        TH1D* h_pt1_P = new TH1D("h_pt1_P",";jet p_{T,1} [GeV];", 120, 0, 600);
        TH1D* h_pt1_J = new TH1D("h_pt1_J",";jet p_{T,1} [GeV];", 120, 0, 600);
        TH1D* h_pt1_J_weighted = new TH1D("h_pt1_J_weighted",";jet p_{T,1} [GeV];", 120, 0, 600);
        pJet->Draw("jetpt1>>h_pt1_P", "", "goff");
        jJet->Draw("jetpt1>>h_pt1_J", "", "goff");

        // jetpt2
        TH1D* h_pt2_P = new TH1D("h_pt2_P",";jet p_{T,2} [GeV];", 120, 0, 600);
        TH1D* h_pt2_J = new TH1D("h_pt2_J",";jet p_{T,2} [GeV];", 120, 0, 600);
        pJet->Draw("jetpt2>>h_pt2_P", "", "goff");
        jJet->Draw("jetpt2>>h_pt2_J", "", "goff");

        // X_J, A_J, dPhi
        TH1D* h_XJ_P = new TH1D("h_XJ_P",";X_{J} = p_{T,2}/p_{T,1};", 60, 0, 1.2);
        TH1D* h_XJ_J = new TH1D("h_XJ_J",";X_{J} = p_{T,2}/p_{T,1};", 60, 0, 1.2);
        pJet->Draw("X_J>>h_XJ_P", "", "goff");
        jJet->Draw("X_J>>h_XJ_J", "", "goff");

        TH1D* h_AJ_P = new TH1D("h_AJ_P",";A_{J} = (p_{T,1}-p_{T,2})/(p_{T,1}+p_{T,2});", 60, 0, 1.0);
        TH1D* h_AJ_J = new TH1D("h_AJ_J",";A_{J} = (p_{T,1}-p_{T,2})/(p_{T,1}+p_{T,2});", 60, 0, 1.0);
        pJet->Draw("A_J>>h_AJ_P", "", "goff");
        jJet->Draw("A_J>>h_AJ_J", "", "goff");

        TH1D* h_dphi_P = new TH1D("h_dphi_P",";#Delta#phi(j_{1},j_{2});", 64, 0, 3.2);
        TH1D* h_dphi_J = new TH1D("h_dphi_J",";#Delta#phi(j_{1},j_{2});", 64, 0, 3.2);
        pJet->Draw("dPhi>>h_dphi_P", "", "goff");
        jJet->Draw("dPhi>>h_dphi_J", "", "goff");


        // Style + normalize
        StyleHist(h_pt1_P, kBlue+1, 20); StyleHist(h_pt1_J, kRed+1, 24);
        StyleHist(h_pt2_P, kBlue+1, 20); StyleHist(h_pt2_J, kRed+1, 24);
        StyleHist(h_XJ_P,  kBlue+1, 20); StyleHist(h_XJ_J,  kRed+1, 24);
        StyleHist(h_AJ_P,  kBlue+1, 20); StyleHist(h_AJ_J,  kRed+1, 24);
        StyleHist(h_dphi_P,kBlue+1, 20); StyleHist(h_dphi_J,kRed+1, 24);
        Normalize(h_pt1_P); Normalize(h_pt1_J);
        Normalize(h_pt2_P); Normalize(h_pt2_J);
        Normalize(h_XJ_P);  Normalize(h_XJ_J);
        Normalize(h_AJ_P);  Normalize(h_AJ_J);
        Normalize(h_dphi_P);Normalize(h_dphi_J);

        // Save
        DrawAndSave({h_pt1_P, h_pt1_J},
                                {"Pythia8 pp", "JEWEL vacuum"},
                                std::string(outDir)+"/jet_pt1_pythia_vs_jewelVac.png",
                                "jet p_{T,1} [GeV]",
                                "dN/dp_{T,1}",
                                "Leading jet p_{T} distribution");
        DrawAndSave({h_pt2_P, h_pt2_J},
                                {"Pythia8 pp", "JEWEL vacuum"},
                                std::string(outDir)+"/jet_pt2_pythia_vs_jewelVac.png",
                                "jet p_{T,2} [GeV]",
                                "dN/dp_{T,2}",
                                "Subleading jet p_{T} distribution");
        DrawAndSave({h_XJ_P,  h_XJ_J},
                                {"Pythia8 pp", "JEWEL vacuum"},
                                std::string(outDir)+"/XJ_pythia_vs_jewelVac.png",
                                "X_{J}",
                                "dN/dX_{J}",
                                "Dijet momentum balance X_{J} distribution");
        DrawAndSave({h_AJ_P,  h_AJ_J},
                                {"Pythia8 pp", "JEWEL vacuum"},
                                std::string(outDir)+"/AJ_pythia_vs_jewelVac.png",
                                "A_{J}",
                                "dN/dA_{J}",
                                "Dijet momentum imbalance A_{J} distribution");
        DrawAndSave({h_dphi_P, h_dphi_J},
                                {"Pythia8 pp", "JEWEL vacuum"},
                                std::string(outDir)+"/dphi_pythia_vs_jewelVac.png",
                                "#Delta#phi(j_{1},j_{2})",
                                    "dN/d#Delta#phi",
                                "Dijet azimuthal angle difference distribution");
        
    }

    // =========================
    // Particle-level comparisons
    // =========================
    if (pPar && jPar) {
        // pt (vector<float>): draw all entries across all events
        bool havePt  = (pPar->GetBranch("pt")  && jPar->GetBranch("pt"));
        bool haveEta = (pPar->GetBranch("eta") && jPar->GetBranch("eta"));
        if (havePt) {
            TH1D* h_p_pt_P = new TH1D("h_p_pt_P",";particle p_{T} [GeV];", 80, 0, 500);
            TH1D* h_p_pt_J = new TH1D("h_p_pt_J",";particle p_{T} [GeV];", 80, 0, 500);
            TH1D* h_p_pt_J_weighted = new TH1D("h_p_pt_J_weighted",";particle p_{T} [GeV];", 80, 0, 500);
            TH1D* h_p_pt_J_weighted_Tnoapply = new TH1D("h_p_pt_J_weighted_Tnoapply",";particle p_{T} [GeV];", 80, 0, 500);
            pPar->Draw("pt>>h_p_pt_P", "", "goff");
            jPar->Draw("pt>>h_p_pt_J", "", "goff");
            jWeightPar->Draw("pt>>h_p_pt_J_weighted", "weights", "goff");
            jWeightPar->Draw("pt>>h_p_pt_J_weighted_Tnoapply", "", "goff");

            StyleHist(h_p_pt_P, kBlue+1, 20);
            StyleHist(h_p_pt_J, kRed+1, 24);
            StyleHist(h_p_pt_J_weighted, kGreen, 25);
            StyleHist(h_p_pt_J_weighted_Tnoapply, kMagenta, 26);
            Normalize(h_p_pt_P); Normalize(h_p_pt_J); Normalize(h_p_pt_J_weighted); Normalize(h_p_pt_J_weighted_Tnoapply);
            DrawAndSave({h_p_pt_P, h_p_pt_J, h_p_pt_J_weighted, h_p_pt_J_weighted_Tnoapply},
                                {"Pythia8 pp", "JEWEL vacuum", "JEWEL weighted", "JEWEL weighted (not applied)"},
                                    std::string(outDir)+"/particle_pt_pythia_vs_jewelVac.png",
                                    "particle p_{T} [GeV]", 
                                    "dN/dp_{T}",
                                    "Particle transverse momentum distribution",
                                    -999,-999,-999,-999,
                                    /*logx=*/false,
                                    /*logy=*/true);

            // make ratio of the already-normalized histograms
            TH1D *ptRatio = (TH1D*)h_p_pt_J->Clone("ptRatio");
            ptRatio->Divide(h_p_pt_P);
            TH1D *ptRatio_weighted = (TH1D*)h_p_pt_J_weighted->Clone("ptRatio_weighted");
            ptRatio_weighted->Divide(h_p_pt_P);
            TH1D *ptRatio_weighted_Tnoapply = (TH1D*)h_p_pt_J_weighted_Tnoapply->Clone("ptRatio_weighted_Tnoapply");
            ptRatio_weighted_Tnoapply->Divide(h_p_pt_P);

            ptRatio->GetYaxis()->SetRangeUser(0,5);

            DrawAndSave({ptRatio,ptRatio_weighted,ptRatio_weighted_Tnoapply}, {"JEWEL Unweighted/Pythia8 ratio", "JEWEL weighted / Pythia8 ratio", "JEWEL (weight not applied)) / Pythia8 ratio"},
                         std::string(outDir)+"/particle_pt_ratio_jewelVac_over_pythia.png",
                         "p_{T} [GeV]",
                         "Ratio of dN/dp_{T} (JEWEL/Pythia8)",
                         "Ratio of particle p_{T} distributions (JEWEL vacuum / Pythia8 pp)",
                          -999,-999,0,5,
                         /*logx=*/false,
                         /*logy=*/false);

        } else {
            std::cerr << "NOTE: ParticleTree missing pt branch in one of the files.\n";
        }

        if (haveEta) {
            TH1D* h_p_eta_P = new TH1D("h_p_eta_P",";particle #eta;", 80, -6.0, 6.0);
            TH1D* h_p_eta_J = new TH1D("h_p_eta_J",";particle #eta;", 80, -6.0, 6.0);
            TH1D* h_p_eta_J_weighted = new TH1D("h_p_eta_J_weighted",";particle #eta;", 80, -6.0, 6.0);
            TH1D* h_p_eta_J_weighted_Tnoapply = new TH1D("h_p_eta_J_weighted_Tnoapply",";particle #eta;", 80, -6.0, 6.0);
            pPar->Draw("eta>>h_p_eta_P", "", "goff");
            jPar->Draw("eta>>h_p_eta_J", "", "goff");
            jWeightPar->Draw("eta>>h_p_eta_J_weighted", "weights", "goff");
            jWeightPar->Draw("eta>>h_p_eta_J_weighted_Tnoapply", "", "goff");
            StyleHist(h_p_eta_P, kBlue+1, 20);
            StyleHist(h_p_eta_J, kRed+1, 24);
            StyleHist(h_p_eta_J_weighted, kGreen, 25);
            StyleHist(h_p_eta_J_weighted_Tnoapply, kMagenta, 26);
            Normalize(h_p_eta_P); Normalize(h_p_eta_J);Normalize(h_p_eta_J_weighted); Normalize(h_p_eta_J_weighted_Tnoapply);
            DrawAndSave({h_p_eta_P, h_p_eta_J, h_p_eta_J_weighted, h_p_eta_J_weighted_Tnoapply}, 
                {"Pythia8 pp", "JEWEL vacuum", "JEWEL weighted", "JEWEL weighted (not applied)"},
                                    std::string(outDir)+"/particl e_eta_pythia_vs_jewelVac.png",
                                    "particle #eta",
                                    "dN/d#eta",
                                    "Particle pseudorapidity distribution");


        } else {
            std::cerr << "NOTE: ParticleTree missing eta branch in one of the files.\n";
        }
    }

    std::cout << "Done. PNGs saved to: " << outDir << "\n";
}