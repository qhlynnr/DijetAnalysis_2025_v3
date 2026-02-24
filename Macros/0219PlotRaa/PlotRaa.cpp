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


static void DrawTLatexLines(const std::vector<std::string>& lines,
              double x = 0.14, double y = 0.94,
              double spacing = 0.05, double textSize = 0.045,
              bool useNDC = true)
{
  TLatex lat;
  lat.SetNDC(useNDC);
  lat.SetTextSize(textSize);
  double yy = y;
  for (const auto& s : lines) {
  lat.SetTextFont(42);  // Helvetica (normal, not bold)
  lat.DrawLatex(x, yy, s.c_str());
  yy -= spacing;
  // prevent drawing off the canvas
  if (yy < 0.05) break;
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

static void SaveOne(TH1* h, const char* name, const char* title,
                    bool logy=false, const char* outdir="Plots/Compare/")
{
  gSystem->mkdir(outdir, kTRUE);

  TCanvas c("c","c",800,800);
  c.SetTicks(1,1);
  c.SetLeftMargin(0.13);
  c.SetBottomMargin(0.12);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.08);
  if (logy) c.SetLogy(1);

  StyleH(h);
  h->Draw("hist");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.050);
  lat.DrawLatex(0.15, 0.94, title);

 // c.SaveAs((std::string(outdir)+"/"+name+".pdf").c_str());
  c.SaveAs((std::string(outdir)+"/"+name+".png").c_str());
}


static void SaveTwo(TH1* h1, TH1* h2, const char* name, const char* title,
  const std::vector<std::string>& info,
    bool logy=false, bool normalize=false, const char* outdir="Plots")
{

  gSystem->mkdir(outdir, kTRUE);

  TCanvas c("c","c",800,800);
  c.SetTicks(1,1);
  c.SetLeftMargin(0.13);
  c.SetBottomMargin(0.12);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.08);
  if (logy) c.SetLogy(1);

  // Clone histograms so we don't modify originals
  std::string n1 = std::string(h1->GetName()) + "_tmp";
  std::string n2 = std::string(h2->GetName()) + "_tmp";
  TH1* h1c = (TH1*)h1->Clone(n1.c_str());
  TH1* h2c = (TH1*)h2->Clone(n2.c_str());

  // Optionally normalize to unit area (use sum of weights)
  if (normalize) {
    double s1 = h1c->GetSumOfWeights();
    double s2 = h2c->GetSumOfWeights();
    if (s1 > 0) h1c->Scale(1.0 / s1);
    if (s2 > 0) h2c->Scale(1.0 / s2);
  }

  // Apply styling to both histogram copies
  StyleH(h1c);
  StyleH(h2c);

  // Distinguish styles
  h1c->SetLineColor(kBlue);
  h1c->SetLineWidth(2);
  h2c->SetLineColor(kRed);
  h2c->SetLineWidth(2);
  h2c->SetLineStyle(2);

  // Make sure y-range accomodates both histograms
  const double m1 = h1c->GetMaximum();
  const double m2 = h2c->GetMaximum();
  double ymax = std::max(m1, m2);
  if (ymax <= 0) ymax = 1.0;
  if (logy) {
    // ensure positive minimum for log scale
    double ymin = std::min(h1c->GetMinimum() > 0 ? h1c->GetMinimum() : 1e-6,
       h2c->GetMinimum() > 0 ? h2c->GetMinimum() : 1e-6);
    if (ymin <= 0) ymin = 1e-6;
    h1c->SetMinimum(ymin);
    h1c->SetMaximum(ymax * 50.0);
  } else {
    h1c->SetMaximum(ymax * 1.3);
  }

  // Draw histograms (draw h1 first to set axis)
  h1c->Draw("hist");
  h2c->Draw("hist same");

  // Build simple legend using the histogram titles (text before first ';')
  TLegend leg(0.50, 0.75, 0.88, 0.88);
  std::string title1 = "med";
  std::string title2 = "vac";
  leg.AddEntry(h1c, title1.substr(0, title1.find(';')).c_str(), "l");
  leg.AddEntry(h2c, title2.substr(0, title2.find(';')).c_str(), "l");
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.030); // increase legend font size
  leg.Draw();

  // Draw the provided title at the top (NDC coordinates)
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.14, 0.94, title);

  // Draw info lines slightly below the title
  DrawTLatexLines(info, 0.20, 0.50, 0.03, 0.02);

 // c.SaveAs((std::string(outdir)+"/"+name+".pdf").c_str());
  c.SaveAs((std::string(outdir)+"/"+name+".png").c_str());

  // clean up clones
  delete h1c;
  delete h2c;
}
void PlotRaa()
{
  gStyle->SetOptStat(0);

  const char* f_vac_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021926_JewelPtHat40_10Kevt/Jewelpbpb_5360GeV_ptm40_10Kevt_021826_C5_vac_11050627_nobkg.root"; // change to your file path
  const char* f_med_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021926_JewelPtHat40_10Kevt/Jewelpbpb_5360GeV_ptm40_10Kevt_021826_C5_med_11050627_jewel4momsub.root"; // change to your file path
  TFile *f_vac = TFile::Open(f_vac_name, "READ");
  TFile *f_med = TFile::Open(f_med_name, "READ");
  if (!f_vac || f_vac->IsZombie()) {
    std::cerr << "Error: Cannot open file " << f_vac_name << std::endl;
    return;
  }
  if (!f_med || f_med->IsZombie()) {
    std::cerr << "Error: Cannot open file " << f_med_name << std::endl;
    return;
  }
  // input trees (in your file they are named ParticleTree and JetTree)
  TTree *ParticleTreeVac = (TTree*)f_vac->Get("ParticleTree");
  TTree *JetTreeVac      = (TTree*)f_vac->Get("JetTree");
  TTree *ParticleTreeMed = (TTree*)f_med->Get("ParticleTree");
  TTree *JetTreeMed      = (TTree*)f_med->Get("JetTree");

  if(!ParticleTreeVac){ std::cerr << "ERROR: missing ParticleTree\n"; return; }
  if(!JetTreeVac){      std::cerr << "ERROR: missing JetTree\n";      return; }
  if(!ParticleTreeMed){ std::cerr << "ERROR: missing ParticleTree\n"; return; }
  if(!JetTreeMed){      std::cerr << "ERROR: missing JetTree\n";      return; }

  // -------------------------
  // ParticleTreeVac branches
  // -------------------------
  std::vector<float> *pxVac  = nullptr;
  std::vector<float> *pyVac  = nullptr;
  std::vector<float> *pzVac  = nullptr;
  std::vector<float> *eVac   = nullptr;
  std::vector<float> *ptVac  = nullptr;
  std::vector<float> *etaVac = nullptr;
  std::vector<float> *phiVac = nullptr;
  std::vector<int>   *pdgIdVac = nullptr;
  std::vector<int>   *sourceVac = nullptr;
  double ET_miss_particleVac = 0.0;
  float  weightsVac = 1.0f;

  ParticleTreeVac->SetBranchAddress("px", &pxVac);
  ParticleTreeVac->SetBranchAddress("py", &pyVac);
  ParticleTreeVac->SetBranchAddress("pz", &pzVac);
  ParticleTreeVac->SetBranchAddress("e",  &eVac);
  ParticleTreeVac->SetBranchAddress("pt", &ptVac);
  ParticleTreeVac->SetBranchAddress("eta",&etaVac);
  ParticleTreeVac->SetBranchAddress("phi",&phiVac);
  ParticleTreeVac->SetBranchAddress("pdgId", &pdgIdVac);
  ParticleTreeVac->SetBranchAddress("ET_miss_particle", &ET_miss_particleVac);
  ParticleTreeVac->SetBranchAddress("source", &sourceVac);
  ParticleTreeVac->SetBranchAddress("weights", &weightsVac);

  // ---------------------
  // JetTreeVac branches
  // ---------------------
  std::vector<float> *jet_pxVac  = nullptr;
  std::vector<float> *jet_pyVac  = nullptr;
  std::vector<float> *jet_pzVac  = nullptr;
  std::vector<float> *jet_eVac   = nullptr;
  std::vector<float> *jet_ptVac  = nullptr;
  std::vector<float> *jet_etaVac = nullptr;
  std::vector<float> *jet_phiVac = nullptr;

  double A_JVac = -1.0;
  double X_JVac = -1.0;
  double jetpt1Vac = 0.0;
  double jetpt2Vac = 0.0;
  double jetphi1Vac = 0.0;
  double jetphi2Vac = 0.0;
  double ET_miss_jetVac = 0.0;
  double ET_miss_jet_rawVac = 0.0;
  double ET_miss_jet_thermalVac = 0.0;
  double dPhiVac = 0.0;

  std::vector<float> *jet_px_rawVac  = nullptr;
  std::vector<float> *jet_py_rawVac  = nullptr;
  std::vector<float> *jet_pz_rawVac  = nullptr;
  std::vector<float> *jet_e_rawVac   = nullptr;
  std::vector<float> *jet_pt_rawVac  = nullptr;
  std::vector<float> *jet_eta_rawVac = nullptr;
  std::vector<float> *jet_phi_rawVac = nullptr;

  std::vector<float> *jet_px_thermalVac  = nullptr;
  std::vector<float> *jet_py_thermalVac  = nullptr;
  std::vector<float> *jet_pz_thermalVac  = nullptr;
  std::vector<float> *jet_e_thermalVac   = nullptr;
  std::vector<float> *jet_pt_thermalVac  = nullptr;
  std::vector<float> *jet_eta_thermalVac = nullptr;
  std::vector<float> *jet_phi_thermalVac = nullptr;

  JetTreeVac->SetBranchAddress("jet_px",  &jet_pxVac);
  JetTreeVac->SetBranchAddress("jet_py",  &jet_pyVac);
  JetTreeVac->SetBranchAddress("jet_pz",  &jet_pzVac);
  JetTreeVac->SetBranchAddress("jet_e",   &jet_eVac);
  JetTreeVac->SetBranchAddress("jet_pt",  &jet_ptVac);
  JetTreeVac->SetBranchAddress("jet_eta", &jet_etaVac);
  JetTreeVac->SetBranchAddress("jet_phi", &jet_phiVac);

  JetTreeVac->SetBranchAddress("A_J",   &A_JVac);
  JetTreeVac->SetBranchAddress("X_J",   &X_JVac);
  JetTreeVac->SetBranchAddress("jetpt1",&jetpt1Vac);
  JetTreeVac->SetBranchAddress("jetpt2",&jetpt2Vac);
  JetTreeVac->SetBranchAddress("jetphi1",&jetphi1Vac);
  JetTreeVac->SetBranchAddress("jetphi2",&jetphi2Vac);
  JetTreeVac->SetBranchAddress("ET_miss_jet", &ET_miss_jetVac);
  JetTreeVac->SetBranchAddress("dPhi", &dPhiVac);

  // -------------------------
  // ParticleTreeMed branches (duplicated from Vac)
  // -------------------------
  std::vector<float> *pxMed  = nullptr;
  std::vector<float> *pyMed  = nullptr;
  std::vector<float> *pzMed  = nullptr;
  std::vector<float> *eMed   = nullptr;
  std::vector<float> *ptMed  = nullptr;
  std::vector<float> *etaMed = nullptr;
  std::vector<float> *phiMed = nullptr;
  std::vector<int>   *pdgIdMed = nullptr;
  std::vector<int>   *sourceMed = nullptr;
  double ET_miss_particleMed = 0.0;
  float  weightsMed = 1.0f;

  ParticleTreeMed->SetBranchAddress("px", &pxMed);
  ParticleTreeMed->SetBranchAddress("py", &pyMed);
  ParticleTreeMed->SetBranchAddress("pz", &pzMed);
  ParticleTreeMed->SetBranchAddress("e",  &eMed);
  ParticleTreeMed->SetBranchAddress("pt", &ptMed);
  ParticleTreeMed->SetBranchAddress("eta",&etaMed);
  ParticleTreeMed->SetBranchAddress("phi",&phiMed);
  ParticleTreeMed->SetBranchAddress("pdgId", &pdgIdMed);
  ParticleTreeMed->SetBranchAddress("ET_miss_particle", &ET_miss_particleMed);
  ParticleTreeMed->SetBranchAddress("source", &sourceMed);
  ParticleTreeMed->SetBranchAddress("weights", &weightsMed);

  // ---------------------
  // JetTreeMed branches (duplicated from Vac)
  // ---------------------
  std::vector<float> *jet_pxMed  = nullptr;
  std::vector<float> *jet_pyMed  = nullptr;
  std::vector<float> *jet_pzMed  = nullptr;
  std::vector<float> *jet_eMed   = nullptr;
  std::vector<float> *jet_ptMed  = nullptr;
  std::vector<float> *jet_etaMed = nullptr;
  std::vector<float> *jet_phiMed = nullptr;

  double A_JMed = -1.0;
  double X_JMed = -1.0;
  double jetpt1Med = 0.0;
  double jetpt2Med = 0.0;
  double jetphi1Med = 0.0;
  double jetphi2Med = 0.0;
  double ET_miss_jetMed = 0.0;
  double ET_miss_jet_rawMed = 0.0;
  double ET_miss_jet_thermalMed = 0.0;
  double dPhiMed = 0.0;

  std::vector<float> *jet_px_rawMed  = nullptr;
  std::vector<float> *jet_py_rawMed  = nullptr;
  std::vector<float> *jet_pz_rawMed  = nullptr;
  std::vector<float> *jet_e_rawMed   = nullptr;
  std::vector<float> *jet_pt_rawMed  = nullptr;
  std::vector<float> *jet_eta_rawMed = nullptr;
  std::vector<float> *jet_phi_rawMed = nullptr;

  std::vector<float> *jet_px_thermalMed  = nullptr;
  std::vector<float> *jet_py_thermalMed  = nullptr;
  std::vector<float> *jet_pz_thermalMed  = nullptr;
  std::vector<float> *jet_e_thermalMed   = nullptr;
  std::vector<float> *jet_pt_thermalMed  = nullptr;
  std::vector<float> *jet_eta_thermalMed = nullptr;
  std::vector<float> *jet_phi_thermalMed = nullptr;

  JetTreeMed->SetBranchAddress("jet_px",  &jet_pxMed);
  JetTreeMed->SetBranchAddress("jet_py",  &jet_pyMed);
  JetTreeMed->SetBranchAddress("jet_pz",  &jet_pzMed);
  JetTreeMed->SetBranchAddress("jet_e",   &jet_eMed);
  JetTreeMed->SetBranchAddress("jet_pt",  &jet_ptMed);
  JetTreeMed->SetBranchAddress("jet_eta", &jet_etaMed);
  JetTreeMed->SetBranchAddress("jet_phi", &jet_phiMed);

  JetTreeMed->SetBranchAddress("A_J",   &A_JMed);
  JetTreeMed->SetBranchAddress("X_J",   &X_JMed);
  JetTreeMed->SetBranchAddress("jetpt1",&jetpt1Med);
  JetTreeMed->SetBranchAddress("jetpt2",&jetpt2Med);
  JetTreeMed->SetBranchAddress("jetphi1",&jetphi1Med);
  JetTreeMed->SetBranchAddress("jetphi2",&jetphi2Med);
  JetTreeMed->SetBranchAddress("ET_miss_jet", &ET_miss_jetMed);
  JetTreeMed->SetBranchAddress("ET_miss_jet_raw", &ET_miss_jet_rawMed);
  JetTreeMed->SetBranchAddress("ET_miss_jet_thermal", &ET_miss_jet_thermalMed);
  JetTreeMed->SetBranchAddress("dPhi", &dPhiMed);

  JetTreeMed->SetBranchAddress("jet_px_raw",  &jet_px_rawMed);
  JetTreeMed->SetBranchAddress("jet_py_raw",  &jet_py_rawMed);
  JetTreeMed->SetBranchAddress("jet_pz_raw",  &jet_pz_rawMed);
  JetTreeMed->SetBranchAddress("jet_e_raw",   &jet_e_rawMed);
  JetTreeMed->SetBranchAddress("jet_pt_raw",  &jet_pt_rawMed);
  JetTreeMed->SetBranchAddress("jet_eta_raw", &jet_eta_rawMed);
  JetTreeMed->SetBranchAddress("jet_phi_raw", &jet_phi_rawMed);

  JetTreeMed->SetBranchAddress("jet_px_thermal",  &jet_px_thermalMed);
  JetTreeMed->SetBranchAddress("jet_py_thermal",  &jet_py_thermalMed);
  JetTreeMed->SetBranchAddress("jet_pz_thermal",  &jet_pz_thermalMed);
  JetTreeMed->SetBranchAddress("jet_e_thermal",   &jet_e_thermalMed);
  JetTreeMed->SetBranchAddress("jet_pt_thermal",  &jet_pt_thermalMed);
  JetTreeMed->SetBranchAddress("jet_eta_thermal", &jet_eta_thermalMed);
  JetTreeMed->SetBranchAddress("jet_phi_thermal", &jet_phi_thermalMed);

  // quick sanity read
  std::cout << "SANITY CHECK: reading first entries from trees\n";
  std::cout << "ParticleTreeVac entries = " << ParticleTreeVac->GetEntries() << "\n";
  std::cout << "JetTreeVac entries      = " << JetTreeVac->GetEntries() << "\n";
  std::cout << "ParticleTreeMed entries = " << ParticleTreeMed->GetEntries() << "\n";
  std::cout << "JetTreeMed entries      = " << JetTreeMed->GetEntries() << "\n";

  ParticleTreeVac->GetEntry(0);
  JetTreeVac->GetEntry(0);
  std::cout << "Event0: weightsVac=" << weightsVac
              << " nPart=" << (pxVac? (int)pxVac->size():-1)
              << " ETmissPart=" << ET_miss_particleVac << "\n";
  std::cout << "Weight check: weightsVac=" << weightsVac << "\n";
  std::cout << "Event0: nJets=" << (jet_ptVac? (int)jet_ptVac->size():-1)
              << " jetpt1Vac=" << jetpt1Vac
              << " A_JVac=" << A_JVac << "\n";
  for (size_t j = 0; j < (ptVac->size()); ++j) {
    std::cout << "pt:" << ptVac->at(j) << " eta:" << etaVac->at(j) << " phi:" << phiVac->at(j)
              << " pdgId:" << pdgIdVac->at(j) << "\n";
  }

  std::cout << "——————————————————" << std::endl;
  TH1F* hpt_vac = new TH1F("hpt_vac", "Particle p_{T} (Vacuum);p_{T} (GeV/c);Entries", 100, 40, 400);
  TH1F* hpt_med = new TH1F("hpt_med", "Particle p_{T} (Medium);p_{T} (GeV/c);Entries", 100, 40, 400);

  hpt_vac->Sumw2();
  hpt_med->Sumw2();

  for (int i = 0; i < ParticleTreeVac->GetEntries(); ++i) {
    ParticleTreeVac->GetEntry(i);
    for (size_t j = 0; j < pxVac->size(); ++j) {
        hpt_vac->Fill(ptVac->at(j), weightsVac);
    }
  }
  for (int i = 0; i < ParticleTreeMed->GetEntries(); ++i) {
    ParticleTreeMed->GetEntry(i);
    for (size_t j = 0; j < pxMed->size(); ++j) {
        hpt_med->Fill(ptMed->at(j), weightsMed);
    }
  }

  hpt_vac->Scale(1.0 / hpt_vac->GetSumOfWeights());
  hpt_med->Scale(1.0 / hpt_med->GetSumOfWeights());
  // -------------------------
  // Create histograms for MED and VAC (same binning, different names)
  // -------------------------
  std::vector<std::string> info = {
    "#sqrt(s)=5.36TeV",
    "JEWEL #hat{p}_{T}^{min} = 250",
    "Anti-k_{T}, R=0.4",
    "200 < p_{T} < 224",
    "#eta < 2.1",
    "#Delta #phi > 7/8 #pi"
  };
  
  const char* outdir = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/022026RaaPthat40/";
  SaveTwo(hpt_med, hpt_vac, "pt", "pt med vs vac", {""}, true, false, outdir);
}
