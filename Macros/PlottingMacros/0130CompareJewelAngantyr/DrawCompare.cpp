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
    bool logy=false, bool normalize=false, const char* outdir="Plots",
    float latx=0.20, float laty=0.50)
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
// Set a meaningful Y-axis title (change when normalizing)
  // Apply styling to both histogram copies
  StyleH(h1c);
  StyleH(h2c);
std::string ytitle = normalize ? "Normalized entries" : "Entries";
h1c->GetYaxis()->SetTitle(ytitle.c_str());
h2c->GetYaxis()->SetTitle(ytitle.c_str());

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
TLegend leg(0.14, 0.74, 0.56, 0.88);
  std::string title1 = "Embed + Jewel with Medium + FJ Bkg";
  std::string title2 = "Jewel Only with Medium + No Bkg";
  leg.AddEntry(h1c, title1.substr(0, title1.find(';')).c_str(), "l");
  leg.AddEntry(h2c, title2.substr(0, title2.find(';')).c_str(), "l");
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.025); // increase legend font size
  leg.Draw();

  // Draw the provided title at the top (NDC coordinates)
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.14, 0.94, title);

  // Draw info lines slightly below the title
  DrawTLatexLines(info, latx, laty, 0.03, 0.02);

 // c.SaveAs((std::string(outdir)+"/"+name+".pdf").c_str());
  c.SaveAs((std::string(outdir)+"/"+name+".png").c_str());

  // clean up clones
  delete h1c;
  delete h2c;
}


void Draw_Compareangjew()
{
  gStyle->SetOptStat(0);

  // Open two files (change filenames to your actual files)
  TFile fEmbed("/Users/lynnrong/Desktop/Research/012125_JewelEmbedAngantyr/Embed_AngJew_5360GeV_ptm250_50Kevt_C5B3p5_nodum_med_012926.root","READ");
  TFile fJew("/Users/lynnrong/Desktop/Research/MonteCarlo/0129Compare/Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_nodum_med_103428_nojetbkg.root","READ");

  // Get the identically-named trees from each file
  auto* JetTreeEmbed      = (TTree*)fEmbed.Get("JetTree");
  auto* ParticleTreeEmbed = (TTree*)fEmbed.Get("ParticleTree");
  auto* JetTreeJew      = (TTree*)fJew.Get("JetTree");
  auto* ParticleTreeJew = (TTree*)fJew.Get("ParticleTree");

  if (!JetTreeEmbed || !ParticleTreeEmbed) { printf("Missing JetTree or ParticleTree in embed file\n"); return; }
  if (!JetTreeJew || !ParticleTreeJew) { printf("Missing JetTree or ParticleTree in jew file\n"); return; }

  // Process BOTH embed and jew trees in parallel (duplicate histograms & fillings)

  printf("Loaded JetTreeEmbed/ParticleTreeEmbed and JetTreeJew/ParticleTreeJew. Processing both.\n");

  // -------------------------
  // JetTree branches (vectors) - EMBED
  // -------------------------
  std::vector<float>* jet_pt_embed  = nullptr;
  std::vector<float>* jet_eta_embed = nullptr;
  std::vector<float>* jet_phi_embed = nullptr;
  JetTreeEmbed->SetBranchAddress("jet_pt",  &jet_pt_embed);
  JetTreeEmbed->SetBranchAddress("jet_eta", &jet_eta_embed);
  JetTreeEmbed->SetBranchAddress("jet_phi", &jet_phi_embed);

  // -------------------------
  // JetTree branches (scalars) - EMBED
  // -------------------------
  double A_J_embed=0, X_J_embed=0, jetpt1_embed=0, jetpt2_embed=0, jetphi1_embed=0, jetphi2_embed=0, ET_miss_jet_embed=0, dPhi_embed=0;
  JetTreeEmbed->SetBranchAddress("A_J",         &A_J_embed);
  JetTreeEmbed->SetBranchAddress("X_J",         &X_J_embed);
  JetTreeEmbed->SetBranchAddress("jetpt1",      &jetpt1_embed);
  JetTreeEmbed->SetBranchAddress("jetpt2",      &jetpt2_embed);
  JetTreeEmbed->SetBranchAddress("jetphi1",     &jetphi1_embed);
  JetTreeEmbed->SetBranchAddress("jetphi2",     &jetphi2_embed);
  JetTreeEmbed->SetBranchAddress("ET_miss_jet", &ET_miss_jet_embed);
  JetTreeEmbed->SetBranchAddress("dPhi",        &dPhi_embed);

  // -------------------------
  // ParticleTree branches - EMBED
  // -------------------------
  std::vector<float>* pt_embed  = nullptr;
  std::vector<float>* eta_embed = nullptr;
  std::vector<float>* phi_embed = nullptr;
  std::vector<int>*   src_embed = nullptr;
  ParticleTreeEmbed->SetBranchAddress("pt",    &pt_embed);
  ParticleTreeEmbed->SetBranchAddress("eta",   &eta_embed);
  ParticleTreeEmbed->SetBranchAddress("phi",   &phi_embed);
  ParticleTreeEmbed->SetBranchAddress("source",&src_embed);
  double ET_miss_particle_embed = 0;
  ParticleTreeEmbed->SetBranchAddress("ET_miss_particle", &ET_miss_particle_embed);

  // -------------------------
  // JetTree branches (vectors) - JEW
  // -------------------------
  std::vector<float>* jet_pt_jew  = nullptr;
  std::vector<float>* jet_eta_jew = nullptr;
  std::vector<float>* jet_phi_jew = nullptr;
  JetTreeJew->SetBranchAddress("jet_pt",  &jet_pt_jew);
  JetTreeJew->SetBranchAddress("jet_eta", &jet_eta_jew);
  JetTreeJew->SetBranchAddress("jet_phi", &jet_phi_jew);

  // -------------------------
  // JetTree branches (scalars) - JEW
  // -------------------------
  double A_J_jew=0, X_J_jew=0, jetpt1_jew=0, jetpt2_jew=0, jetphi1_jew=0, jetphi2_jew=0, ET_miss_jet_jew=0, dPhi_jew=0;
  JetTreeJew->SetBranchAddress("A_J",         &A_J_jew);
  JetTreeJew->SetBranchAddress("X_J",         &X_J_jew);
  JetTreeJew->SetBranchAddress("jetpt1",      &jetpt1_jew);
  JetTreeJew->SetBranchAddress("jetpt2",      &jetpt2_jew);
  JetTreeJew->SetBranchAddress("jetphi1",     &jetphi1_jew);
  JetTreeJew->SetBranchAddress("jetphi2",     &jetphi2_jew);
  JetTreeJew->SetBranchAddress("ET_miss_jet", &ET_miss_jet_jew);
  JetTreeJew->SetBranchAddress("dPhi",        &dPhi_jew);

  // -------------------------
  // ParticleTree branches - JEW
  // -------------------------
  std::vector<float>* pt_jew  = nullptr;
  std::vector<float>* eta_jew = nullptr;
  std::vector<float>* phi_jew = nullptr;
  std::vector<int>*   src_jew = nullptr;
  ParticleTreeJew->SetBranchAddress("pt",    &pt_jew);
  ParticleTreeJew->SetBranchAddress("eta",   &eta_jew);
  ParticleTreeJew->SetBranchAddress("phi",   &phi_jew);
  ParticleTreeJew->SetBranchAddress("source",&src_jew);
  double ET_miss_particle_jew = 0;
  ParticleTreeJew->SetBranchAddress("ET_miss_particle", &ET_miss_particle_jew);

  // -------------------------
  // Create histograms for EMBED and JEW (same binning, different names)
  // -------------------------
  auto mk1 = [&](const char* name, const char* title, int n, double a, double b){
    return new TH1F(name, title, n, a, b);
  };

TH1F* hJetPt_embed  = mk1("hJetPt_embed",  "Jets (embed);jet p_{T} [GeV]",                80, 0, 1100);
TH1F* hJetEta_embed = mk1("hJetEta_embed", "Jets (embed);jet #eta",                      60,-3, 3);
TH1F* hJetPhi_embed = mk1("hJetPhi_embed", "Jets (embed);jet #phi",                      64,-3.2, 3.2);
TH1F* hParPt_embed  = mk1("hParPt_embed",  "Particles (embed);p_{T} [GeV]",         80, 0, 1100);
TH1F* hParEta_embed = mk1("hParEta_embed", "Particles (embed);#eta",                60,-6, 6);
TH1F* hParPhi_embed = mk1("hParPhi_embed", "Particles (embed);#phi",                64,-3.2, 3.2);
TH1F* hSrc_embed    = mk1("hSrc_embed",    "Particle source (embed);source label",   3,-0.5, 2.5);

TH1F* hAJ_embed     = mk1("hAJ_embed",     "A_{J} (embed);A_{J}",                        50, 0, 1.0);
TH1F* hXJ_embed     = mk1("hXJ_embed",     "x_{J} (embed);x_{J}",                        50, 0, 1.0);
TH1F* hJet1Pt_embed = mk1("hJet1Pt_embed", "Leading jet (embed);p_{T,1} [GeV]",          80, 0, 1100);
TH1F* hJet2Pt_embed = mk1("hJet2Pt_embed", "Subleading jet (embed);p_{T,2} [GeV]",       80, 0, 1100);
TH1F* hDPhi_embed   = mk1("hDPhi_embed",   "#Delta#phi_{12} (embed);#Delta#phi",         64, -3.2, 3.2);
TH1F* hETJet_embed  = mk1("hETJet_embed",  "Missing E_{T} (jets,embed);E_{T}^{miss} [GeV]", 80, 0, 150);
TH1F* hPhi1_embed   = mk1("hPhi1_embed",   "Leading jet #phi (embed);#phi_{1}",         64, -3.2, 3.2);
TH1F* hPhi2_embed   = mk1("hPhi2_embed",   "Subleading jet #phi (embed);#phi_{2}",      64, -3.2, 3.2);
TH1F* hETPar_embed  = mk1("hETPar_embed",  "Missing E_{T} (particles,embed);E_{T}^{miss} [GeV]", 80, 0, 150);

TH1F* hJetPt_jew  = mk1("hJetPt_jew",  "Jets (jew);jet p_{T} [GeV]",                80, 0, 1100);
TH1F* hJetEta_jew = mk1("hJetEta_jew", "Jets (jew);jet #eta",                      60,-3, 3);
TH1F* hJetPhi_jew = mk1("hJetPhi_jew", "Jets (jew);jet #phi",                      64,-3.2, 3.2);
TH1F* hParPt_jew  = mk1("hParPt_jew",  "Particles (jew);p_{T} [GeV]",         80, 0, 50);
TH1F* hParEta_jew = mk1("hParEta_jew", "Particles (jew);#eta",                60,-6, 6);
TH1F* hParPhi_jew = mk1("hParPhi_jew", "Particles (jew);#phi",                64,-3.2, 3.2);
TH1F* hSrc_jew    = mk1("hSrc_jew",    "Particle source (jew);source label",   3,-0.5, 2.5);

TH1F* hAJ_jew     = mk1("hAJ_jew",     "A_{J} (jew);A_{J}",                        50, 0, 1.0);
TH1F* hXJ_jew     = mk1("hXJ_jew",     "x_{J} (jew);x_{J}",                        50, 0, 1.0);
TH1F* hJet1Pt_jew = mk1("hJet1Pt_jew", "Leading jet (jew);p_{T,1} [GeV]",          80, 0, 550);
TH1F* hJet2Pt_jew = mk1("hJet2Pt_jew", "Subleading jet (jew);p_{T,2} [GeV]",       80, 0, 550);
TH1F* hDPhi_jew   = mk1("hDPhi_jew",   "#Delta#phi_{12} (jew);#Delta#phi",         64, 0, 3.2);
TH1F* hETJet_jew  = mk1("hETJet_jew",  "Missing E_{T} (jets,jew);E_{T}^{miss} [GeV]", 80, 0, 150);
TH1F* hPhi1_jew   = mk1("hPhi1_jew",   "Leading jet #phi (jew);#phi_{1}",         64, -3.2, 3.2);
TH1F* hPhi2_jew   = mk1("hPhi2_jew",   "Subleading jet #phi (jew);#phi_{2}",      64, -3.2, 3.2);
TH1F* hETPar_jew  = mk1("hETPar_jew",  "Missing E_{T} (particles,jew);E_{T}^{miss} [GeV]", 80, 0, 150);
  for (auto* h : {hJetPt_embed,hJetEta_embed,hJetPhi_embed,hParPt_embed,hParEta_embed,hParPhi_embed,hSrc_embed,
                  hAJ_embed,hXJ_embed,hJet1Pt_embed,hJet2Pt_embed,hDPhi_embed,hETJet_embed,hPhi1_embed,hPhi2_embed,hETPar_embed,
                  hJetPt_jew,hJetEta_jew,hJetPhi_jew,hParPt_jew,hParEta_jew,hParPhi_jew,hSrc_jew,
                  hAJ_jew,hXJ_jew,hJet1Pt_jew,hJet2Pt_jew,hDPhi_jew,hETJet_jew,hPhi1_jew,hPhi2_jew,hETPar_jew}) {
    h->SetLineWidth(2);
  }
  hSrc_embed->GetXaxis()->SetBinLabel(1,"E (Embed)");
  hSrc_embed->GetXaxis()->SetBinLabel(2,"J (JEWEL)");
  hSrc_embed->GetXaxis()->SetBinLabel(3,"Other");
  hSrc_jew->GetXaxis()->SetBinLabel(1,"A (Angantyr)");
  hSrc_jew->GetXaxis()->SetBinLabel(2,"J (JEWEL)");
  hSrc_jew->GetXaxis()->SetBinLabel(3,"Other");

  // -------------------------
  // Fill EMBED
  // -------------------------

    const Long64_t nE = std::min(JetTreeEmbed->GetEntries(), ParticleTreeEmbed->GetEntries());
    int JewelJetCount = 0;
    int JewelParticleCount = 0;
    int JewelTotalCount = 0;

    int EmbedJetCount = 0;
    int EmbedParticleCount = 0;
    int EmbedTotalCount = 0;

    for (Long64_t i=0;i<nE;i++){
      JetTreeEmbed->GetEntry(i);
      ParticleTreeEmbed->GetEntry(i);
      bool passEtaSel = true;
      // determine leading-jet eta safely
      double lead_eta_embed = 0;
      hJet1Pt_embed->Fill(jetpt1_embed);
      hJet2Pt_embed->Fill(jetpt2_embed);
      hPhi1_embed->Fill(jetphi1_embed);
      hPhi2_embed->Fill(jetphi2_embed);
      hDPhi_embed->Fill(abs(dPhi_embed));

      if (jetpt1_embed < 200) continue;
      if (jetpt2_embed < 30) continue;
      if (abs(dPhi_embed) < (7.0/8.0)*TMath::Pi()) continue;
      if (passEtaSel == false) continue;
      // vector jets
      if (jet_pt_embed) for (size_t j=0;j<jet_pt_embed->size();j++){
        hJetEta_embed->Fill(jet_eta_embed->at(j));
        if (abs(jet_eta_embed->at(j)) > 2.1) continue;
        hJetPt_embed ->Fill(jet_pt_embed->at(j));
        hJetPhi_embed->Fill(jet_phi_embed->at(j));
        EmbedJetCount++;
      }

    for (size_t j=0;j<jet_eta_embed->size();j++){
        if (abs(jet_eta_embed->at(j)) > 2.1){
            passEtaSel = false;
            break;
        }
      }
      hAJ_embed->Fill(A_J_embed);
      hXJ_embed->Fill(X_J_embed);
      hETJet_embed->Fill(ET_miss_jet_embed);

      EmbedTotalCount++;
      if (pt_embed) for (size_t k=0;k<pt_embed->size();k++){
        hParPt_embed ->Fill(pt_embed->at(k));
        hParEta_embed->Fill(eta_embed->at(k));
        hParPhi_embed->Fill(phi_embed->at(k));
        int s = (src_embed ? src_embed->at(k) : -1);
        if (s==0) hSrc_embed->Fill(0);
        else if (s==1) hSrc_embed->Fill(1);
        else hSrc_embed->Fill(2);
        EmbedParticleCount++;
      }

      if (ParticleTreeEmbed->GetBranch("ET_miss_particle")) hETPar_embed->Fill(ET_miss_particle_embed);
    }
    cout << "Embed Jet Count: " << EmbedJetCount << endl;
    cout << "Embed Particle Count: " << EmbedParticleCount << endl;
    cout << "Embed Event Count: " << EmbedTotalCount << endl;
  
  // -------------------------
  // Fill JEW
  // -------------------------
    const Long64_t nE_jew = std::min(JetTreeJew->GetEntries(), ParticleTreeJew->GetEntries());
    for (Long64_t i=0;i<nE_jew;i++){
      JetTreeJew->GetEntry(i);
      ParticleTreeJew->GetEntry(i);

      // vector jets
      bool passEtaSel = true;
      hJet1Pt_jew->Fill(jetpt1_jew);
      hJet2Pt_jew->Fill(jetpt2_jew);
      hPhi1_jew->Fill(jetphi1_jew);
      hPhi2_jew->Fill(jetphi2_jew);
      hDPhi_jew->Fill(abs(dPhi_jew));
      if (jetpt1_jew < 200) continue;
      if (jetpt2_jew < 30) continue;
      if (abs(dPhi_jew) < (7.0/8.0)*TMath::Pi()) continue;

    if (jet_pt_jew) for (size_t j=0;j<jet_pt_jew->size();j++){
        hJetEta_jew->Fill(jet_eta_jew->at(j));
        if (abs(jet_eta_jew->at(j)) > 2.1) continue;
        hJetPt_jew ->Fill(jet_pt_jew->at(j));
        hJetPhi_jew->Fill(jet_phi_jew->at(j));
        JewelJetCount++;
      }

      for (size_t j=0;j<jet_eta_jew->size();j++){
        if (abs(jet_eta_jew->at(j)) > 2.1){
            passEtaSel = false;
            break;
        }
      }
      if (passEtaSel == false) continue;

      JewelTotalCount++;

      hAJ_jew->Fill(A_J_jew);
      hXJ_jew->Fill(X_J_jew);
      hETJet_jew->Fill(ET_miss_jet_jew);


      if (pt_jew) for (size_t k=0;k<pt_jew->size();k++){
        hParPt_jew ->Fill(pt_jew->at(k));
        hParEta_jew->Fill(eta_jew->at(k));
        hParPhi_jew->Fill(phi_jew->at(k));
        int s = (src_jew ? src_jew->at(k) : -1);
        if (s==0) hSrc_jew->Fill(0);
        else if (s==1) hSrc_jew->Fill(1);
        else hSrc_jew->Fill(2);
        JewelParticleCount++;
      }

      if (ParticleTreeJew->GetBranch("ET_miss_particle")) hETPar_jew->Fill(ET_miss_particle_jew);
    }
    cout << "Jewel Jet Count: " << JewelJetCount << endl;
    cout << "Jewel Particle Count: " << JewelParticleCount << endl;
    cout << "Jewel Event Count: " << JewelTotalCount << endl;
    cout << "Number of JEWEL events processed: " << JetTreeJew->GetEntries() << endl;
    cout << "Number of embedded events processed: " << JetTreeEmbed->GetEntries() << endl;
    
std::vector<std::string> info = {
    "#sqrt(s)=5.36TeV",
    "JEWEL #hat{p}_{T}^{min} = 250GEV",
    "Anti-k_{T}, R=0.4",
    "200 < p_{T}",
    "#eta < 2.1",
    "#Delta #phi > 7/8 #pi",
    "0~5% centrality",
    "Jet Count (Embed): " + std::to_string(EmbedJetCount),
    "Jet Count (JEWEL): " + std::to_string(JewelJetCount)

};
  const char* outdir = "Plots/0129NewCompare";

  // -------------------------
  // Save comparisons (embed vs jew) using SaveTwo
  // -------------------------
 
SaveTwo(hAJ_embed, hAJ_jew, "AJ_embed_vs_jew", "PbPb 5360 GeV: A_{J}", info, false, true, outdir, 0.5,0.6);
SaveTwo(hXJ_embed, hXJ_jew, "XJ_embed_vs_jew", "PbPb 5360 GeV: x_{J}", info, false, false, outdir,0.2,0.6);
}