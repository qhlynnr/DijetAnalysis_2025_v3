// File: plot_embed_summary.C
// Usage:
//   root -l -q 'plot_embed_summary.C("embedded_angantyrjewel_withmedium.root")'

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLatex.h"
#include <vector>
#include <cmath>

void plot_embed_summary(const char* fname="embedded_angantyrjewel_withmedium.root")
{
  gStyle->SetOptStat(0);

  TFile f(fname, "READ");
  auto* JetTree      = (TTree*)f.Get("JetTree");
  auto* ParticleTree = (TTree*)f.Get("ParticleTree");
  if (!JetTree || !ParticleTree) { printf("Missing JetTree or ParticleTree\n"); return; }

  // -------------------------
  // JetTree branches
  // -------------------------
  std::vector<float>* jet_pt  = nullptr;
  std::vector<float>* jet_eta = nullptr;
  std::vector<float>* jet_phi = nullptr;
  JetTree->SetBranchAddress("jet_pt",  &jet_pt);
  JetTree->SetBranchAddress("jet_eta", &jet_eta);
  JetTree->SetBranchAddress("jet_phi", &jet_phi);

  // -------------------------
  // ParticleTree branches
  // -------------------------
  std::vector<float>* pt  = nullptr;
  std::vector<float>* eta = nullptr;
  std::vector<float>* phi = nullptr;
  std::vector<int>*   src = nullptr;
  ParticleTree->SetBranchAddress("pt",   &pt);
  ParticleTree->SetBranchAddress("eta",  &eta);
  ParticleTree->SetBranchAddress("phi",  &phi);
  ParticleTree->SetBranchAddress("source",&src);

  // -------------------------
  // Histograms (reasonable defaults)
  // -------------------------
  TH1F* hJetPt  = new TH1F("hJetPt",  "Jets;jet p_{T} [GeV];jets / bin",                80, 0, 400);
  TH1F* hJetEta = new TH1F("hJetEta", "Jets;jet #eta;jets / bin",                      60,-3, 3);
  TH1F* hJetPhi = new TH1F("hJetPhi", "Jets;jet #phi;jets / bin",                      64,-3.2, 3.2);

  TH1F* hParPt  = new TH1F("hParPt",  "Particles;p_{T} [GeV];particles / bin",         80, 0, 50);
  TH1F* hParEta = new TH1F("hParEta", "Particles;#eta;particles / bin",                60,-6, 6);
  TH1F* hParPhi = new TH1F("hParPhi", "Particles;#phi;particles / bin",                64,-3.2, 3.2);

  TH1F* hSrc    = new TH1F("hSrc",    "Particle source;source label;particles / bin",   3,-0.5, 2.5);
  hSrc->GetXaxis()->SetBinLabel(1,"A (Angantyr)");
  hSrc->GetXaxis()->SetBinLabel(2,"J (JEWEL)");
  hSrc->GetXaxis()->SetBinLabel(3,"Other");

  // nicer line width
  for (auto* h : {hJetPt,hJetEta,hJetPhi,hParPt,hParEta,hParPhi,hSrc}) h->SetLineWidth(2);

  // -------------------------
  // Fill from trees
  // -------------------------
  const Long64_t nE = std::min(JetTree->GetEntries(), ParticleTree->GetEntries());
  for (Long64_t i=0;i<nE;i++){
    JetTree->GetEntry(i);
    ParticleTree->GetEntry(i);

    // jets
    for (size_t j=0;j<jet_pt->size();j++){
      hJetPt ->Fill(jet_pt->at(j));
      hJetEta->Fill(jet_eta->at(j));
      hJetPhi->Fill(jet_phi->at(j));
    }

    // particles
    for (size_t k=0;k<pt->size();k++){
      hParPt ->Fill(pt->at(k));
      hParEta->Fill(eta->at(k));
      hParPhi->Fill(phi->at(k));

      int s = (src ? src->at(k) : -1);
      if (s==0) hSrc->Fill(0);
      else if (s==1) hSrc->Fill(1);
      else hSrc->Fill(2);
    }
  }

  // -------------------------
  // Canvas: 800x800, 2x3 pads
  // -------------------------
  TCanvas* c = new TCanvas("c","Embedding summary",800,800);
  c->Divide(2,3,0.01,0.01);

  auto prepPad = [&](int ipad){
    c->cd(ipad);
    gPad->SetTicks(1,1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.07);
  };

  auto drawH = [&](TH1F* h, int ipad, bool logy=false){
    prepPad(ipad);
    if (logy) gPad->SetLogy(1);
    h->SetTitle("");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->Draw("hist");
  };

  // reasonable: logy for pT spectra, linear for angles/source
  drawH(hJetPt,  1, true);
  drawH(hJetEta, 2, false);
  drawH(hJetPhi, 3, false);
  drawH(hParPt,  4, true);
  drawH(hParEta, 5, false);
  drawH(hSrc,    6, false);

  // Titles per pad
  TLatex lat; lat.SetNDC(true); lat.SetTextSize(0.055);

  c->cd(1); lat.DrawLatex(0.14,0.93,"Jet p_{T}");
  c->cd(2); lat.DrawLatex(0.14,0.93,"Jet #eta");
  c->cd(3); lat.DrawLatex(0.14,0.93,"Jet #phi");
  c->cd(4); lat.DrawLatex(0.14,0.93,"Particle p_{T}");
  c->cd(5); lat.DrawLatex(0.14,0.93,"Particle #eta");
  c->cd(6); lat.DrawLatex(0.14,0.93,"Particle source");

  c->SaveAs("embedding_summary_800x800.pdf");
  c->SaveAs("embedding_summary_800x800.png");
}