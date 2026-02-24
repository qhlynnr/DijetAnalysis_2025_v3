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
                    bool logy=false, const char* outdir="Plots/Compare/", double ymin = 999, double ymax = -1.0)
{
  gSystem->mkdir(outdir, kTRUE);

  TCanvas c("c","c",800,800);
  c.SetTicks(1,1);
  c.SetLeftMargin(0.13);
  c.SetBottomMargin(0.12);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.08);
  if (logy) c.SetLogy(1);

  std::string tmpName = std::string(h->GetName() ? h->GetName() : "h") + "_tmp";
  TH1* hc = (TH1*)h->Clone(tmpName.c_str());

  StyleH(hc);
  hc->SetLineColor(kBlack);
  hc->SetLineWidth(2);

  if (ymax > 0.0) {
      hc->SetMaximum(ymax);
  } else {
      double m = hc->GetMaximum();
      if (m <= 0) m = 1.0;
      hc->SetMaximum(m * 1.3);
  }

  if (ymin != 999) {
      hc->SetMinimum(ymin);
  } else if (logy) {
      double yminCandidate = hc->GetMinimum() > 0 ? hc->GetMinimum() : 1e-6;
      hc->SetMinimum(yminCandidate);
  }

  if(logy){
      double ymin = hc->GetMinimum() > 0 ? hc->GetMinimum() : 1e-6;
      hc->SetMinimum(ymin);
  }

  hc->Draw("hist E1");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.050);
  lat.DrawLatex(0.15, 0.94, title);

  c.SaveAs((std::string(outdir)+"/"+name+".png").c_str());

  delete hc;
}

static void SaveTwo(TH1* h1, TH1* h2, const char* name, const char* title,
    const std::vector<std::string>& info,
        bool logy=false, bool normalize=false, const char* outdir="Plots",
        const char* legendLabel1 = "med", const char* legendLabel2 = "vac",
        double userYmax = -1.0)
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

    // Make sure y-range accomodates both histograms
    const double m1 = h1c->GetMaximum();
    const double m2 = h2c->GetMaximum();
    double ymax = std::max(m1, m2);
    if (ymax <= 0) ymax = 1.0;

    bool useUserYmax = (userYmax > 0.0);

    if (logy) {
        // ensure positive minimum for log scale
        double yminCandidate1 = h1c->GetMinimum() > 0 ? h1c->GetMinimum() : 1e-6;
        double yminCandidate2 = h2c->GetMinimum() > 0 ? h2c->GetMinimum() : 1e-6;
        double ymin = std::min(yminCandidate1, yminCandidate2);
        if (ymin <= 0) ymin = 1e-6;
        h1c->SetMinimum(ymin);

        if (useUserYmax) {
            // If user provided an explicit ymax, use it directly
            h1c->SetMaximum(userYmax);
        } else {
            // otherwise scale relative to observed maximum
            h1c->SetMaximum(ymax * 50.0);
        }
    } else {
        if (useUserYmax) {
            h1c->SetMaximum(userYmax);
        } else {
            h1c->SetMaximum(ymax * 1.3);
        }
    }

    // Draw histograms (draw h1 first to set axis)
    h1c->Draw("hist E1");
    h2c->Draw("hist same E1");

    // Build simple legend using the provided labels (text before first ';')
    TLegend leg(0.50, 0.75, 0.88, 0.88);
    std::string title1 = std::string(legendLabel1);
    std::string title2 = std::string(legendLabel2);

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
void Draw()
{
  gStyle->SetOptStat(0);

  const char* f_vac_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/022426_JewelPtHat40_10KEvt_2/Jewelpbpb_5360GeV_ptm40_10Kevt_022026_C5_vac_11050627_precluster.root"; // change to your file path
  const char* f_med_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/022426_JewelPtHat40_10KEvt_2/Jewelpbpb_5360GeV_ptm40_10Kevt_022026_C5_med_11050627_precluster.root"; // change to your file path
  const char* f_vacref_name = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/Pythia8Jets_pp_5360GeV_HardQCD_pTHatMin40_10KEvts.root"; // change to your file path

  TFile *f_vac = TFile::Open(f_vac_name, "READ");
  TFile *f_med = TFile::Open(f_med_name, "READ");
  TFile *f_vacref = TFile::Open(f_vacref_name, "READ");
  if (!f_vac || f_vac->IsZombie()) {
    std::cerr << "Error: Cannot open file " << f_vac_name << std::endl;
    return;
  }
  if (!f_med || f_med->IsZombie()) {
    std::cerr << "Error: Cannot open file " << f_med_name << std::endl;
    return;
  }
  if (!f_vacref || f_vacref->IsZombie()) {
    std::cerr << "Error: Cannot open file " << f_vacref_name << std::endl;
    return;
  }
  // input trees (in your file they are named ParticleTree and JetTree)
  TTree *ParticleTreeVac = (TTree*)f_vac->Get("ParticleTree");
  TTree *ParticleTreeMed = (TTree*)f_med->Get("ParticleTree");
  TTree *ParticleTreeVacRef = (TTree*)f_vacref->Get("ParticleTree");

  if(!ParticleTreeVac){ std::cerr << "ERROR: missing ParticleTree\n"; return; }
  if(!ParticleTreeMed){ std::cerr << "ERROR: missing ParticleTree\n"; return; }
  if(!ParticleTreeVacRef){ std::cerr << "ERROR: missing ParticleTree\n"; return; }

  // -------------------------
  // ParticleTreeVac branches
  // -------------------------
  std::vector<double> *pxVac  = nullptr;
  std::vector<double> *pyVac  = nullptr;
  std::vector<double> *pzVac  = nullptr;
  std::vector<double> *eVac   = nullptr;
  std::vector<double> *ptVac  = nullptr;
  std::vector<double> *etaVac = nullptr;
  std::vector<double> *phiVac = nullptr;
  std::vector<int>   *pdgIdVac = nullptr;
  std::vector<int>   *statusVac = nullptr;
  double  weightsVac = 1.0f;

  ParticleTreeVac->SetBranchAddress("px", &pxVac);
  ParticleTreeVac->SetBranchAddress("py", &pyVac);
  ParticleTreeVac->SetBranchAddress("pz", &pzVac);
  ParticleTreeVac->SetBranchAddress("e",  &eVac);
  ParticleTreeVac->SetBranchAddress("pt", &ptVac);
  ParticleTreeVac->SetBranchAddress("eta",&etaVac);
  ParticleTreeVac->SetBranchAddress("phi",&phiVac);
  ParticleTreeVac->SetBranchAddress("pdgId", &pdgIdVac);
  ParticleTreeVac->SetBranchAddress("weights", &weightsVac);
  ParticleTreeVac->SetBranchAddress("status", &statusVac);
  // -------------------------
  // ParticleTreeMed branches (duplicated from Vac)
  // -------------------------
  std::vector<double> *pxMed  = nullptr;
  std::vector<double> *pyMed  = nullptr;
  std::vector<double> *pzMed  = nullptr;
  std::vector<double> *eMed   = nullptr;
  std::vector<double> *ptMed  = nullptr;
  std::vector<double> *etaMed = nullptr;
  std::vector<double> *phiMed = nullptr;
  std::vector<int>   *pdgIdMed = nullptr;
  std::vector<int>   *sourceMed = nullptr;
  std::vector<int>   *statusMed = nullptr;
  double ET_miss_particleMed = 0.0;
  double  weightsMed = 1.0f;

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
  ParticleTreeMed->SetBranchAddress("status", &statusMed);

  vector<double> *pxVacRef  = nullptr;
    vector<double> *pyVacRef  = nullptr;
    vector<double> *pzVacRef  = nullptr;
    vector<double> *eVacRef   = nullptr;
    vector<double> *ptVacRef  = nullptr;
    vector<double> *etaVacRef = nullptr;
    vector<double> *phiVacRef = nullptr;
    vector<int>   *pdgIdVacRef = nullptr;
    double ET_miss_particleVacRef = 0.0;
  ParticleTreeVacRef->SetBranchAddress("px", &pxVacRef);
  ParticleTreeVacRef->SetBranchAddress("py", &pyVacRef);
  ParticleTreeVacRef->SetBranchAddress("pz", &pzVacRef);
  ParticleTreeVacRef->SetBranchAddress("e",  &eVacRef);
  ParticleTreeVacRef->SetBranchAddress("pt", &ptVacRef);
  ParticleTreeVacRef->SetBranchAddress("eta",&etaVacRef);
  ParticleTreeVacRef->SetBranchAddress("phi",&phiVacRef);
  ParticleTreeVacRef->SetBranchAddress("pdgId", &pdgIdVacRef);
  ParticleTreeVacRef->SetBranchAddress("ET_miss_particle", &ET_miss_particleVacRef);

  // quick sanity read
  std::cout << "SANITY CHECK: reading first entries from trees\n";
  std::cout << "ParticleTreeVac entries = " << ParticleTreeVac->GetEntries() << "\n";
  std::cout << "ParticleTreeMed entries = " << ParticleTreeMed->GetEntries() << "\n";
  std::cout << "ParticleTreeVacRef entries = " << ParticleTreeVacRef->GetEntries() << "\n";

  ParticleTreeVac->GetEntry(0);
  std::cout << "Event0: weightsVac=" << weightsVac << std::endl;
  for (size_t j = 0; j < (ptVac->size()); ++j) {
    std::cout << "pt:" << ptVac->at(j) << " eta:" << etaVac->at(j) << " phi:" << phiVac->at(j)
              << " pdgId:" << pdgIdVac->at(j) << "\n";
  }

  std::cout << "——————————————————" << std::endl;
  TH1F* hpt_vac = new TH1F("hpt_vac", "Particle p_{T} (Vacuum);p_{T} (GeV/c);Entries", 25, 20, 350);
  TH1F* hpt_med = new TH1F("hpt_med", "Particle p_{T} (Medium);p_{T} (GeV/c);Entries", 25, 20, 350);
  TH1F* hpt_vacref = new TH1F("hpt_vacref", "Particle p_{T} (Vacuum Ref);p_{T} (GeV/c);Entries", 25, 20, 350);
  hpt_vac->Sumw2();
  hpt_med->Sumw2();

  for (int i = 0; i < ParticleTreeVac->GetEntries(); ++i) {
    ParticleTreeVac->GetEntry(i);
    for (size_t j = 0; j < pxVac->size(); ++j) {
        if (statusVac->at(j) != 1) continue;
        if (TMath::Abs(etaVac->at(j)) > 3.0) continue; 
        hpt_vac->Fill(ptVac->at(j), weightsVac);
    }
  }
  for (int i = 0; i < ParticleTreeMed->GetEntries(); ++i) {
    ParticleTreeMed->GetEntry(i);
    for (size_t j = 0; j < pxMed->size(); ++j) {
        if (statusMed->at(j) != 1) continue; 
        if (TMath::Abs(etaMed->at(j)) > 3.0) continue;
        hpt_med->Fill(ptMed->at(j), weightsMed);
    }
  }

 /* for (int i = 0; i < ParticleTreeVacRef->GetEntries(); ++i) {
    ParticleTreeVacRef->GetEntry(i);
    for (size_t j = 0; j < pxVacRef->size(); ++j) {
        hpt_vacref->Fill(ptVacRef->at(j),1);
    }
  }*/

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

  TH1F* hptratio = (TH1F*)hpt_med->Clone("hptratio");
  hptratio->Divide(hpt_vac);

  cout << "sum of weights vac: " << hpt_vac->GetSumOfWeights() << endl;
  cout << "sum of weights med: " << hpt_med->GetSumOfWeights() << endl;
  cout << "entries vacref: " << hpt_vacref->GetEntries() << endl;
  hpt_vac->Scale(1.0 / hpt_vac->GetSumOfWeights());
  hpt_med->Scale(1.0 / hpt_med->GetSumOfWeights());
  hpt_vacref->Scale(1.0 / hpt_vacref->GetSumOfWeights());


  const char* outdir = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/022326RaaPtHat40Attempt2/";
 // SaveOne(hpt_vac, "pt_vac", "Particle p_{T} (Vacuum)", false, outdir);
 // SaveOne(hpt_med, "pt_med", "Particle p_{T} (Medium)", false, outdir);
 // SaveOne(hpt_vacref, "pt_vacref", "Particle p_{T} (Vacuum Ref)", true, outdir);
  SaveOne(hptratio, "pt_raa", "Particle p_{T} ratio (Medium/Vacuum)", false, outdir,0,2);
  SaveTwo(hpt_med, hpt_vac, "ptcheck", "pt med vs vac", {""}, false, false, outdir,"med", "vac");
  SaveTwo(hpt_med, hpt_vac, "ptcheck_log", "pt med vs vac", {""}, true, false, outdir,"med", "vac");
 // SaveTwo(hpt_vac, hpt_vacref, "VacCompare", "pt pythia8 vs jewel debug", {""}, true, false, outdir,"jewel vac", "pythia8 pp");
}
