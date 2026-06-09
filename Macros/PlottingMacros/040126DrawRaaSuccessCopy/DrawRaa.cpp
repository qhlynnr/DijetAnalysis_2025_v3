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
  //DrawTLatexLines(lines, latexX, latexY, spacing, textSize, useNDC);
  c->SetTicks(1,1);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.08);
  c->SetRightMargin(0.05);

  if (logx) c->SetLogx();
  if (logy) c->SetLogy();

  c->SetGridx();        // vertical grid lines
  c->SetGridy();        // horizontal grid lines
  // set titles on first histogram

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
  const char* f_med_bkg_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/HistRootFiles/031026PtHat40100KEvt/H_Jewelpbpb_5360GeV_ptm40_10Kevt_030426_C10_med_12345_jewel4momsub.root"; // change to your file path
  const char* f_vac_nobkg_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/HistRootFiles/031026PtHat40100KEvt/H_Jewelpbpb_5360GeV_ptm40_10Kevt_030426_C10_vac_12345_nobkg.root"; // change to your file path


  TFile *f_vac_nobkg = TFile::Open(f_vac_nobkg_name, "READ");
  TFile *f_med_bkg = TFile::Open(f_med_bkg_name, "READ");
  TFile *f_cmsdata = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/Data/HEPData-ins1848440-v1-Figure_7-3.root", "READ");
  if (!f_cmsdata || f_cmsdata->IsZombie()) {
        std::cerr << "Error: Could not open CMS data file." << std::endl;
        return;
  }
  TFile *f_atlasdata = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/Data/HEPData-ins1673184-v1-Table_20.root", "READ");
    if (!f_atlasdata || f_atlasdata->IsZombie()) {
          std::cerr << "Error: Could not open ATLAS data file." << std::endl;
          return;
    }
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

  string ptHatMin = "pT-hat min: 40 GeV";
  string eventNum = "Events: 100K";
  string radius = "R = 0.4";

  cutInfoVector.push_back(ptHatMin);
  cutInfoVector.push_back(eventNum);
  cutInfoVector.push_back(radius);

  while (std::getline(ss, item, ';')) cutInfoVector.push_back(item);

  TH1D* hjtpt_vac_nobkg = (TH1D*)f_vac_nobkg->Get("hjtpt");
  TH1D* hjtpt_med_bkg = (TH1D*)f_med_bkg->Get("hjtpt");

  TH1D* hjtpt_vac_nobkg_log = (TH1D*)f_vac_nobkg->Get("hjtpt_log");
  TH1D* hjtpt_med_bkg_log = (TH1D*)f_med_bkg->Get("hjtpt_log");

  TH1D* hAJ_vac_nobkg = (TH1D*)f_vac_nobkg->Get("hAj");
  TH1D* hAJ_med_bkg = (TH1D*)f_med_bkg->Get("hAj");

  TH1D* hXJ_vac_nobkg = (TH1D*)f_vac_nobkg->Get("hXj");
  TH1D* hXJ_med_bkg = (TH1D*)f_med_bkg->Get("hXj");

  vector<TH1D*> histVac = {hjtpt_vac_nobkg, hjtpt_vac_nobkg_log, hAJ_vac_nobkg, hXJ_vac_nobkg};
  vector<TH1D*> histMedBkg = {hjtpt_med_bkg, hjtpt_med_bkg_log, hAJ_med_bkg, hXJ_med_bkg};
 
  for (auto h : histVac) StyleHist(h,kBlue,20);
  for (auto h : histMedBkg) StyleHist(h,kRed,25);

  TH1D* hptratio_log = (TH1D*)hjtpt_med_bkg_log->Clone("hptratio_log");
  hptratio_log->Divide(hjtpt_vac_nobkg_log);
  TH1D* hptratio = (TH1D*)hjtpt_med_bkg->Clone("hptratio");
  hptratio->Divide(hjtpt_vac_nobkg);

  StyleHist(hptratio_log, kRed, 26);
  StyleHist(hptratio, kRed, 26);

  cout << "sum of weights vac: " << hjtpt_vac_nobkg->GetSumOfWeights() << endl;
  cout << "sum of weights med bkg: " << hjtpt_med_bkg->GetSumOfWeights() << endl; 
  hjtpt_vac_nobkg->Scale(1.0 / hjtpt_vac_nobkg->GetSumOfWeights());
  hjtpt_med_bkg->Scale(1.0 / hjtpt_med_bkg->GetSumOfWeights());

  hAJ_med_bkg->Scale(1.0 / hAJ_med_bkg->GetSumOfWeights());
  hAJ_vac_nobkg->Scale(1.0 / hAJ_vac_nobkg->GetSumOfWeights());

  TH1F* h_atlasdata_1 = (TH1F*)f_atlasdata->Get("Table 20/Hist1D_y1");
  TH1F* h_atlasdata_errorp = (TH1F*)f_atlasdata->Get("Table 20/Hist1D_y1_e1plus");
  TH1F* h_atlasdata_errorm = (TH1F*)f_atlasdata->Get("Table 20/Hist1D_y1_e1minus");
  TH1F* h_data_1 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y1");
  TH1F* h_data_errorp = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y1_e1plus");
  TH1F* h_data_errorm = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y1_e1minus");
  if (!h_data_1) {
        std::cerr << "Error: Could not find histogram 'Table 9/h1' in CMS data file." << std::endl;
        return;
  }

  TH1F* h_data_2 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y2");
  TH1F* h_data_errorp2 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y2_e1plus");
  TH1F* h_data_errorm2 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y2_e1minus");
  TH1F* h_data_3 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y3");
  TH1F* h_data_errorp3 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y3_e1plus");
  TH1F* h_data_errorm3 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y3_e1minus");
  TH1F* h_data_4 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y4");
  TH1F* h_data_errorp4 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y4_e1plus");
  TH1F* h_data_errorm4 = (TH1F*)f_cmsdata->Get("Figure 7-3/Hist1D_y4_e1minus");

  BlockCaptionInfo captionInfo;
  captionInfo = {cutInfoVector, 0.60, 0.54, 0.035, 0.04, true};
 
  const char* outdir = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/0331RAAwithData/";
 
  cout << "vac bin1 err " << hjtpt_vac_nobkg->GetBinError(1) << endl;
  cout << "med bkg bin1 err " << hjtpt_med_bkg->GetBinError(1) << endl;
  captionInfo = {cutInfoVector, 0.20, 0.54, 0.035, 0.04, true};

  StyleHist(h_data_1, kBlack, 20);
  StyleHist(h_data_2, kBlack, 20);
  StyleHist(h_data_3, kBlack, 20);
  StyleHist(h_data_4, kBlack, 20);
  StyleHist(h_atlasdata_1, kGreen+2, 21);
for (int i = 1; i <= h_atlasdata_1->GetNbinsX(); i++) {
    double err_up   = h_atlasdata_errorp->GetBinContent(i);
    double err_down = h_atlasdata_errorm->GetBinContent(i);
    double err = 0.5 * (fabs(err_up) + fabs(err_down)); // approximate

    h_atlasdata_1->SetBinError(i, err);
}
for (int i = 1; i <= h_data_1->GetNbinsX(); i++) {
    double err_up   = h_data_errorp->GetBinContent(i);
    double err_down = h_data_errorm->GetBinContent(i);
    double err = 0.5 * (fabs(err_up) + fabs(err_down)); // approximate

    h_data_1->SetBinError(i, err);
}
for (int i = 1; i <= h_data_2->GetNbinsX(); i++) {
    double err_up   = h_data_errorp2->GetBinContent(i);
    double err_down = h_data_errorm2->GetBinContent(i);
    double err = 0.5 * (fabs(err_up) + fabs(err_down)); // approximate

    h_data_2->SetBinError(i, err);
}
for (int i = 1; i <= h_data_3->GetNbinsX(); i++) {
    double err_up   = h_data_errorp3->GetBinContent(i);
    double err_down = h_data_errorm3->GetBinContent(i);
    double err = 0.5 * (fabs(err_up) + fabs(err_down)); // approximate

    h_data_3->SetBinError(i, err);
}
for (int i = 1; i <= h_data_4->GetNbinsX(); i++) {
    double err_up   = h_data_errorp4->GetBinContent(i);
    double err_down = h_data_errorm4->GetBinContent(i);
    double err = 0.5 * (fabs(err_up) + fabs(err_down));
    h_data_4->SetBinError(i, err);
  }

  DrawAndSave({h_data_1, h_data_2, h_data_3,h_data_4},
                         {"CMS Data 1", "CMS Data 2", "CMS Data 3"},
                         std::string(outdir) + "CMSData.png",
                         "Jet pT",
                         "Raa",
                         "CMS Data for Raa",
                         200, 1000,
                         -999, -999,
                         true, false,
                         captionInfo);


  captionInfo = {cutInfoVector, 0.20, 0.84, 0.030, 0.04, true};
  DrawAndSave({hptratio_log, h_data_1, h_atlasdata_1},
                         {"Jewel Pbpb Raa (5.37TeV)","CMS pbpb Jet Data Run2 (5.02 TeV)", "ATLAS pbpb Jet Data Run2 (5.02 TeV)"},
                         std::string(outdir) + "JetPtRatio_log.png",
                         "Jet pT",
                         "dN/dpT_{pbpb} / dN/dpT_{pp}",
                         "Jet pT Ratio Distributions (log scale)",
                         100, 1000,
                         0, 2,
                         true, false,
                         captionInfo);

}
