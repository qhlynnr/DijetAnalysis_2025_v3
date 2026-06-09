// ReadVacTrees.C
// usage:
//   root -l -q 'ReadVacTrees.C("/path/to/file.root")'

#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <iostream>
// File: plot_embed_all.C
// Usage:
//   root -l -q 'plot_embed_all.C("embedded_angantyrjewel_withmedium.root")'
//
// Makes 800x800 plots, no stats box, and saves EACH plot as its own PDF+PNG.

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TMath.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <string>


struct BlockCaptionInfo{
  std::vector<std::string> text;
  double x;
  double y;
  double size;
  double spacing;
  bool useNDC;
};

static void DrawTLatexLines(const BlockCaptionInfo& info)
{
  TLatex lat;
  lat.SetNDC(info.useNDC);
  lat.SetTextSize(info.size);
  double yy = info.y;
  for (const auto& s : info.text) {
    lat.SetTextFont(42);  // Helvetica (normal, not bold)
    lat.DrawLatex(info.x, yy, s.c_str());
    yy -= info.spacing;
  }
}

static void StyleCanvas(TCanvas* c, bool logx = false, bool logy = false){
  c->SetTicks(1,1);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.08);
  c->SetRightMargin(0.05);

  if (logx) c->SetLogx();
  if (logy) c->SetLogy();

  c->SetGridx();        // vertical grid lines
  c->SetGridy();        // horizontal grid lines
}

static void StyleH(TH1* h){
  h->SetTitle("");
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitleOffset(1.25);
}

static void DrawAndSave(const std::vector<TH1*>& hs,
                        const std::vector<std::string>& labels,
                        const std::string& outPng,
                        const std::string& xTitle,
                        const std::string& yTitle = "Normalized entries",
                        const std::string& title = "",
                        float xMin = -999, float xMax = -999,
                        float yMin = -999, float yMax = -999,
                        bool logx = false,
                        bool logy = false,
                        BlockCaptionInfo captionInfo = {}) {
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
  StyleCanvas(c, logx, logy);

  TH1* h0 = hs[0];
  if (!title.empty()) h0->SetTitle(title.c_str());
  h0->GetXaxis()->SetTitle(xTitle.c_str());
  h0->GetYaxis()->SetTitle(yTitle.c_str());

  // determine y range using all histograms
  double maxy = 0;
  for (auto h : hs) if (h) maxy = std::max(maxy, h->GetMaximum());
  if (maxy > 0) h0->SetMaximum(1.25 * maxy);

  // X axis bounds: apply individually if provided
  double curXmin = h0->GetXaxis()->GetXmin();
  double curXmax = h0->GetXaxis()->GetXmax();
  double newXmin = curXmin;
  double newXmax = curXmax;
  if (xMin != -999) newXmin = xMin;
  if (xMax != -999) newXmax = xMax;
  if (newXmin != curXmin || newXmax != curXmax) h0->GetXaxis()->SetRangeUser(newXmin, newXmax);

  // Y axis bounds: set min/max independently if provided
  if (yMin != -999) h0->SetMinimum(yMin);
  if (yMax != -999) h0->SetMaximum(yMax);
  // draw histograms (assume they are already styled)
  bool first = true;
  for (auto h : hs) {
    if (!h) continue;
    if (first) { h->Draw("E1 HIST"); first = false; }
    else       { h->Draw("E1 HIST SAME"); }
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
    leg->AddEntry(hs[i], lab[i].c_str(), "lep");
  }
  leg->Draw();
  DrawTLatexLines(captionInfo);
  c->SaveAs(outPng.c_str());
}

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

void DrawHist()
{

  gStyle->SetOptStat(0); // Disable statistics box
  const char* f_med_bkg_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/HistRootFiles/0331CheckAJ2760TeV/H_Jewel_pbpb_100K_2760GeV_ptm250_C10_med_12345_4momsub_R0p3.root"; // change to your file path
  const char* f_vac_nobkg_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/HistRootFiles/0331CheckAJ2760TeV/H_Jewel_pbpb_100K_2760GeV_ptm250_C10_vac_12345_nobkg_R0p3.root"; // change to your file path

  TFile *f_vac_nobkg = TFile::Open(f_vac_nobkg_name, "READ");
  TFile *f_med_bkg = TFile::Open(f_med_bkg_name, "READ");

  TNamed* cutInfo = (TNamed*)f_vac_nobkg->Get("Cuts");
  if (cutInfo) {
    std::cout << "Cuts used in the analysis: " << cutInfo->GetTitle() << std::endl;
  } else {
    std::cerr << "Warning: Cuts information not found in the file." << std::endl;
  }
  if (!f_vac_nobkg || f_vac_nobkg->IsZombie()) {
    std::cerr << "Error: Cannot open file " << f_vac_nobkg_name << std::endl;
    return;
  }
  if (!f_med_bkg || f_med_bkg->IsZombie()) {
    std::cerr << "Error: Cannot open file " << f_med_bkg_name << std::endl;
    return;
  }

  // Parse the cuts information into a vector of strings
  std::stringstream ss(cutInfo->GetTitle());
  std::string item;
  std::vector<std::string> cutInfoVector;

  string collisionEnergy = "Collision Energy: 2.76 TeV";
  string ptHatMin = "pT-hat min: 40 GeV";
  string eventNum = "Events: 100K";
  string radius = "R = 0.3";

  cutInfoVector.push_back(ptHatMin);
  cutInfoVector.push_back(eventNum);
  cutInfoVector.push_back(radius);

  while (std::getline(ss, item, ';')) cutInfoVector.push_back(item);

  TH1D* hAJ_vac_nobkg = (TH1D*)f_vac_nobkg->Get("hAj");
  TH1D* hAJ_med_bkg = (TH1D*)f_med_bkg->Get("hAj");
  TH1D* hAJ_med_smeared = (TH1D*)f_med_bkg->Get("hAjsmeared");
  TH1D* hAJ_vac_smeared = (TH1D*)f_vac_nobkg->Get("hAjsmeared");

  TH1D* hXJ_vac_nobkg = (TH1D*)f_vac_nobkg->Get("hXj");
  TH1D* hXJ_med_bkg = (TH1D*)f_med_bkg->Get("hXj");

  vector<TH1D*> histVac = {hAJ_vac_nobkg, hXJ_vac_nobkg};
  vector<TH1D*> histMedBkg = {hAJ_med_bkg, hXJ_med_bkg};
 
  for (auto h : histVac) StyleHist(h,kBlue,20);
  for (auto h : histMedBkg) StyleHist(h,kRed,25);

  StyleHist(hAJ_med_smeared, kRed, 26);
  StyleHist(hAJ_vac_smeared, kBlue, 26);

  hAJ_med_bkg->Scale(1.0 / hAJ_med_bkg->GetSumOfWeights());
  hAJ_vac_nobkg->Scale(1.0 / hAJ_vac_nobkg->GetSumOfWeights());
  hXJ_med_bkg->Scale(1.0 / hXJ_med_bkg->GetSumOfWeights());
  hXJ_vac_nobkg->Scale(1.0 / hXJ_vac_nobkg->GetSumOfWeights());
  hAJ_med_smeared->Scale(1.0 / hAJ_med_smeared->GetSumOfWeights());
  hAJ_vac_smeared->Scale(1.0 / hAJ_vac_smeared->GetSumOfWeights());

  BlockCaptionInfo captionInfo = {cutInfoVector, 0.20, 0.54, 0.030, 0.04, true}; 
  const char* outdir = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/040126ValidationAJ/";
  
  captionInfo = {cutInfoVector, 0.30, 0.80, 0.030, 0.04, true}; 
  DrawAndSave({hAJ_vac_nobkg, hAJ_med_bkg},
                         {"Vacuum(NoBkgSub)", "Medium(Bkgsub)"},
                         std::string(outdir) + "AJComp_R0p3.png",
                         "A_{J}",
                         "dN/dA_{J}",
                         "A_{J} Distributions",
                         0, 1,
                         -999, -999,
                         false, false,
                         captionInfo);
  DrawAndSave({hAJ_vac_smeared, hAJ_med_smeared},
                          {"Vacuum(Smeared)", "Medium(Smeared)"},
                          std::string(outdir) + "AJSmearedComp_R0p3.png",
                          "A_{J}",
                          "dN/dA_{J}",
                          "Smeared A_{J} Distributions",
                          0, 1,
                          -999, -999,
                          false, false,
                          captionInfo);


/*
  DrawAndSave({hXJ_vac_nobkg, hXJ_med_bkg},
                         {"Vacuum(NoBkgSub)", "Medium(Bkgsub)"},
                         std::string(outdir) + "XJComp.png",
                         "X_{J}",
                         "dN/dX_{J}",
                         "X_{J} Distributions",
                         0, 1,
                         -999, -999,
                         false, false,
                         captionInfo);

  captionInfo = {cutInfoVector, 0.60, 0.54, 0.035, 0.04, true};


  DrawAndSave({hjtpt_vac_nobkg, hjtpt_med_bkg},
                         {"Vacuum(NoBkgSub)", "Medium(Bkgsub)"},
                         std::string(outdir) + "JetPtComp_log.png",
                         "Jet pT",
                         "dN/dpT",
                         "Jet pT Distributions",
                         0, 1600,
                         -999, -999,
                         false, true,
                         captionInfo);
  DrawAndSave({hjtpt_vac_nobkg, hjtpt_med_bkg},
                         {"Vacuum(NoBkgSub)", "Medium(Bkgsub)"},
                         std::string(outdir) + "JetPtComp.png",
                         "Jet pT",
                         "dN/dpT",
                         "Jet pT Distributions",
                         0, 1600,
                         -999, -999,
                         false, false,
                         captionInfo);

  DrawAndSave({hptratio},
                         {"Raa(withbkgsub)"},
                         std::string(outdir) + "JetPtRatio.png",
                         "Jet pT",
                         "dN/dpT_{pbpb} / dN/dpT_{pp}",
                         "Jet pT Ratio Distributions",
                         50, -999,
                         0, 2,
                         false, false,
                         captionInfo);
  captionInfo = {cutInfoVector, 0.60, 0.74, 0.035, 0.04, true};
*/

}
