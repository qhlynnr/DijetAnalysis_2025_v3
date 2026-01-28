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
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

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
                    bool logy=false, const char* outdir="Plots")
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
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.14, 0.94, title);

  c.SaveAs((std::string(outdir)+"/"+name+".pdf").c_str());
  c.SaveAs((std::string(outdir)+"/"+name+".png").c_str());
}

static void SaveTwo(TH1* h1,TH1* h2, const char* name, const char* title,
                    bool logy=false, const char* outdir="Plots")
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
  h1->Draw("hist");
  h2->Draw("hist same");


  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.14, 0.94, title);

  c.SaveAs((std::string(outdir)+"/"+name+".pdf").c_str());
  c.SaveAs((std::string(outdir)+"/"+name+".png").c_str());
}


void Draw(const char* fname="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/JewelAngantyrEmbed/embedded_angantyrjewel_withmedium.root")
{
  gStyle->SetOptStat(0);

  TFile f(fname,"READ");
  auto* JetTree      = (TTree*)f.Get("JetTree");
  auto* ParticleTree = (TTree*)f.Get("ParticleTree");
  if (!JetTree || !ParticleTree) { printf("Missing JetTree or ParticleTree\n"); return; }

  // -------------------------
  // JetTree branches (vectors)
  // -------------------------
  std::vector<float>* jet_pt  = nullptr;
  std::vector<float>* jet_eta = nullptr;
  std::vector<float>* jet_phi = nullptr;
  JetTree->SetBranchAddress("jet_pt",  &jet_pt);
  JetTree->SetBranchAddress("jet_eta", &jet_eta);
  JetTree->SetBranchAddress("jet_phi", &jet_phi); // keep as reference

  // -------------------------
  // JetTree branches (scalars you listed)
  // NOTE: If any name differs in your file, ROOT will warn and the pointer stays unset.
  // -------------------------
  float A_J=0, X_J=0, jetpt1=0, jetpt2=0, jetphi1=0, jetphi2=0, ET_miss_jet=0, dPhi=0;

  // Only set if branch exists (safer when you're iterating)
  if (JetTree->GetBranch("A_J"))         JetTree->SetBranchAddress("A_J", &A_J);
  if (JetTree->GetBranch("X_J"))         JetTree->SetBranchAddress("X_J", &X_J);
  if (JetTree->GetBranch("jetpt1"))      JetTree->SetBranchAddress("jetpt1", &jetpt1);
  if (JetTree->GetBranch("jetpt2"))      JetTree->SetBranchAddress("jetpt2", &jetpt2);
  if (JetTree->GetBranch("jetphi1"))     JetTree->SetBranchAddress("jetphi1", &jetphi1);
  if (JetTree->GetBranch("jetphi2"))     JetTree->SetBranchAddress("jetphi2", &jetphi2);
  if (JetTree->GetBranch("ET_miss_jet")) JetTree->SetBranchAddress("ET_miss_jet", &ET_miss_jet);
  if (JetTree->GetBranch("dPhi"))        JetTree->SetBranchAddress("dPhi", &dPhi);

  // -------------------------
  // ParticleTree branches (vectors + scalars you listed)
  // -------------------------
  std::vector<float>* pt  = nullptr;
  std::vector<float>* eta = nullptr;
  std::vector<float>* phi = nullptr;
  std::vector<int>*   src = nullptr;

  ParticleTree->SetBranchAddress("pt",   &pt);
  ParticleTree->SetBranchAddress("eta",  &eta);
  ParticleTree->SetBranchAddress("phi",  &phi);
  ParticleTree->SetBranchAddress("source",&src);

  double ET_miss_particle = 0;
  if (ParticleTree->GetBranch("ET_miss_particle"))
    ParticleTree->SetBranchAddress("ET_miss_particle", &ET_miss_particle);

  // -------------------------
  // Histograms (reasonable defaults; adjust as needed)
  // -------------------------

  // New jet-level observables
  TH1F* hAJ     = new TH1F("hAJ",     "A_{J};A_{J};events / bin",                        50, 0, 1.0);
  TH1F* hXJ     = new TH1F("hXJ",     "x_{J};x_{J};events / bin",                        50, 0, 1.2);
  TH1F* hJet1Pt = new TH1F("hJet1Pt", "Leading jet;p_{T,1} [GeV];events / bin",          80, 0, 400);
  TH1F* hJet2Pt = new TH1F("hJet2Pt", "Subleading jet;p_{T,2} [GeV];events / bin",       80, 0, 400);
  TH1F* hDPhi   = new TH1F("hDPhi",   "#Delta#phi_{12};#Delta#phi;events / bin",         64, 0, 3.2);
  TH1F* hETJet  = new TH1F("hETJet",  "Missing E_{T} (jets);E_{T}^{miss} [GeV];events / bin", 80, 0, 400);

  // Jet phi1/phi2 (optional, but you asked to include the observables)
  TH1F* hPhi1   = new TH1F("hPhi1",   "Leading jet #phi;#phi_{1};events / bin",         64, -3.2, 3.2);
  TH1F* hPhi2   = new TH1F("hPhi2",   "Subleading jet #phi;#phi_{2};events / bin",      64, -3.2, 3.2);

  // Particle-level missing ET
  TH1F* hETPar  = new TH1F("hETPar",  "Missing E_{T} (particles);E_{T}^{miss} [GeV];events / bin", 80, 0, 400);

  for (auto* h : {hJetPt,hJetEta,hJetPhi,hParPt,hParEta,hParPhi,hSrc,
                  hAJ,hXJ,hJet1Pt,hJet2Pt,hDPhi,hETJet,hPhi1,hPhi2,hETPar}) {
    h->SetLineWidth(2);
  }

  // -------------------------
  // Fill from trees
  // -------------------------
  const Long64_t nE = std::min(JetTree->GetEntries(), ParticleTree->GetEntries());
  for (Long64_t i=0;i<nE;i++){
    JetTree->GetEntry(i);
    ParticleTree->GetEntry(i);

    // vector jets
    for (size_t j=0;j<jet_pt->size();j++){
      hJetPt ->Fill(jet_pt->at(j));
      hJetEta->Fill(jet_eta->at(j));
      hJetPhi->Fill(jet_phi->at(j));
    }

    // jet-level scalars (fill)
    hAJ->Fill(A_J);
    hXJ->Fill(X_J);
    hJet1Pt->Fill(jetpt1);
    hJet2Pt->Fill(jetpt2);
    hPhi1->Fill(jetphi1);
    hPhi2->Fill(jetphi2);
    hETJet->Fill(ET_miss_jet);
    hDPhi->Fill(dPhi);
    // particle-level scalar
    hETPar->Fill(ET_miss_particle);
  }

  // -------------------------
  // Save each histogram to its own file
  // -------------------------
  SaveOne(hJetPt,  "jet_pt",   "Jet p_{T}", true);
  SaveOne(hJetEta, "jet_eta",  "Jet #eta",  false);
  SaveOne(hJetPhi, "jet_phi",  "Jet #phi",  false);

  SaveOne(hParPt,  "par_pt",   "Particle p_{T}", true);
  SaveOne(hParEta, "par_eta",  "Particle #eta",  false);
  SaveOne(hParPhi, "par_phi",  "Particle #phi",  false);

  SaveOne(hSrc,    "par_source","Particle source", false);

  SaveOne(hAJ,     "A_J",      "A_{J}", false);
  SaveOne(hXJ,     "X_J",      "x_{J}", false);
  SaveOne(hJet1Pt, "jetpt1",   "Leading jet p_{T,1}", true);
  SaveOne(hJet2Pt, "jetpt2",   "Subleading jet p_{T,2}", true);
  SaveOne(hPhi1,   "jetphi1",  "Leading jet #phi_{1}", false);
  SaveOne(hPhi2,   "jetphi2",  "Subleading jet #phi_{2}", false);
  SaveOne(hDPhi,   "dPhi",     "#Delta#phi_{12}", false);
  SaveOne(hETJet,  "ET_miss_jet", "E_{T}^{miss} (jets)", true);

  SaveOne(hETPar,  "ET_miss_particle", "E_{T}^{miss} (particles)", true);

  printf("Saved plots to ./Plots/ (PDF+PNG)\n");
}