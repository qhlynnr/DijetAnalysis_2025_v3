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
    double ymax = 0.0) {
    
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
    
    // Color palette for multiple histograms
    vector<int> colors = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+1, kOrange+1, kViolet+1, kSpring+1};
    
    // Draw first histogram to set up axes
    hists[0]->SetLineColor(colors[0]);
    hists[0]->SetLineWidth(2);
    if (!histtitle.empty()) hists[0]->SetTitle(histtitle.c_str());
    hists[0]->SetXTitle(xtitle.c_str());
    hists[0]->SetYTitle(norm ? "Normalized entries" : "Entries");
    hists[0]->SetMaximum(maxy);
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

    string filename = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/0126Dummies/" + outputfilename + varname + ".png";
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
int DrawDummies(){
    gStyle->SetOptStat(0);

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 800);

    // Define histogram parameters
    
    //TH1F* h_part_EMiss = new TH1F("h_part_EMiss", "Particle ETMiss of PbPb Jewel", 100, 0, 200);
    TFile *f_nodum = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/020126DebugDummies/Jewelpbpb_5360GeV_ptm250_10Kevt_020126_C5_nodum_med_11050627_precluster.root");         
    TFile *f_dum = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/020126DebugDummies/Jewelpbpb_5360GeV_ptm250_10Kevt_020126_C5_dum_med_11050627_precluster.root");         
    TFile *f_norec = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/020126DebugDummies//Jewelpbpb_5360GeV_ptm250_10Kevt_020126_C5_recoff_med_11050627.root");
    //Old Files:
    // TFile *f_dum = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012626_JewelTestDummies_Preclustering/Jewelpbpb_5360GeV_ptm250_1Kevt_012226_Dum_C5_med.root");         
  //  TFile *f_nodum = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012626_JewelTestDummies_Preclustering/Jewelpbpb_5360GeV_ptm250_1Kevt_012226_noDum_C5_med.root");         

  //  TFile *f_nodum = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/020126DebugDummies//Jewelpbpb_5360GeV_ptm250_1Kevt_012226_noDum_C5_med.root");         
  //  TFile *f_dum = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/020126DebugDummies//Jewelpbpb_5360GeV_ptm250_1Kevt_012226_Dum_C5_med.root");         

    TTree *ParticleTree_dum = (TTree*)f_dum->Get("ParticleTree");
    TTree *ParticleTree_nodum = (TTree*)f_nodum->Get("ParticleTree");
    TTree *ParticleTree_norec = (TTree*)f_norec->Get("ParticleTree");

    TH1F* h_part_ETMiss_dum = new TH1F("h_ETMiss_dum", " ETMiss of PbPb Event using Jewel", 120, -5, 100);
    TH1F* h_part_ETMiss_nodum = new TH1F("h_ETMiss_nodum", " ETMiss of PbPb Event using Jewel", 120, -5, 100);
    TH1F* h_part_ETMiss_norec = new TH1F("h_ETMiss_norec", " ETMiss of PbPb Event using Jewel without Reconstruction", 120, -5, 100);
    TH1F* h_part_status_dum = new TH1F("h_status_dum", " Particle status using Jewel with Dummies", 10, -0.5, 9.5);
    TH1F* h_part_status_nodum = new TH1F("h_status_nodum", " Particle status using Jewel without Dummies", 10, -0.5, 9.5);
    TH1F* h_part_mult_dum = new TH1F("h_mult_dum", " Particle multiplicity using Jewel with Dummies", 100, 0, 5000);
    TH1F* h_part_mult_nodum = new TH1F("h_mult_nodum", " Particle multiplicity using Jewel without Dummies", 100, 0,  5000);
    TH1F* h_part_mult_norec = new TH1F("h_mult_norec", " Particle multiplicity using Jewel without Reconstruction", 100, 0,  5000);

    TH1F* h_part_ETMiss_dumwithstatus1 = new TH1F("h_ETMiss_dumwithstatus1", " ETMiss of PbPb Event using Jewel with status 1", 120, -5, 100);
    TH1F* h_part_ETMiss_nodumwithstatus1= new TH1F("h_ETMiss_nodumwithstatus1", " ETMiss of PbPb Event using Jewel with status 1", 120, -5, 100);
    TH1F* h_part_ETMiss_dumwithstatus3 = new TH1F("h_ETMiss_dumwithstatus3", " ETMiss of PbPb Event using Jewel with status 3", 120, -5, 100);
    TH1F* h_part_ETMiss_nodumwithstatus3 = new TH1F("h_ETMiss_nodumwithstatus3", " ETMiss of PbPb Event using Jewel with status 3", 120, -5, 100);
    TH1F* h_part_ETMiss_dumwithcuts = new TH1F("h_ETMiss_dumwithcuts", " ETMiss of PbPb Event using Jewel with cuts", 120, -5, 100);
    
    std::vector<float>* e_dum = nullptr;
    std::vector<float>* px_dum = nullptr;
    std::vector<float>* py_dum = nullptr;
    std::vector<float>* pz_dum = nullptr;
    std::vector<float>* pt_dum = nullptr;
    std::vector<float>* status_dum = nullptr;
    std::vector<int>* pdgId_dum = nullptr;
    ParticleTree_dum->SetBranchAddress("e", &e_dum);
    ParticleTree_dum->SetBranchAddress("px", &px_dum);
    ParticleTree_dum->SetBranchAddress("py", &py_dum);
    ParticleTree_dum->SetBranchAddress("pz", &pz_dum);
    ParticleTree_dum->SetBranchAddress("status", &status_dum);
    ParticleTree_dum->SetBranchAddress("pt", &pt_dum);
    ParticleTree_dum->SetBranchAddress("pdgId", &pdgId_dum);
    std::vector<float>* e_nodum = nullptr;
    std::vector<float>* px_nodum = nullptr;
    std::vector<float>* py_nodum = nullptr;
    std::vector<float>* pz_nodum = nullptr;
    std::vector<float>* pt_nodum = nullptr;
    std::vector<float>* status_nodum = nullptr;
    ParticleTree_nodum->SetBranchAddress("e", &e_nodum);
    ParticleTree_nodum->SetBranchAddress("px", &px_nodum);
    ParticleTree_nodum->SetBranchAddress("py", &py_nodum);
    ParticleTree_nodum->SetBranchAddress("pz", &pz_nodum);
    ParticleTree_nodum->SetBranchAddress("status", &status_nodum);
    ParticleTree_nodum->SetBranchAddress("pt", &pt_nodum);
    std::vector<float>* e_norec = nullptr;
    std::vector<float>* px_norec = nullptr;
    std::vector<float>* py_norec = nullptr;
    std::vector<float>* pz_norec = nullptr;
    std::vector<float>* pt_norec = nullptr;
    std::vector<float>* status_norec = nullptr;
    ParticleTree_norec->SetBranchAddress("e", &e_norec);
    ParticleTree_norec->SetBranchAddress("px", &px_norec);
    ParticleTree_norec->SetBranchAddress("py", &py_norec);
    ParticleTree_norec->SetBranchAddress("pz", &pz_norec);
    ParticleTree_norec->SetBranchAddress("status", &status_norec);
    ParticleTree_norec->SetBranchAddress("pt", &pt_norec);

    int nParticles_dum = 0;
    int nParticles_nodum = 0;
    int nParticles_norec = 0;
    ParticleTree_dum->SetBranchAddress("nParticles", &nParticles_dum);
    ParticleTree_nodum->SetBranchAddress("nParticles", &nParticles_nodum);  
    ParticleTree_norec->SetBranchAddress("nParticles", &nParticles_norec);
    cout << "ParticleTree_dum entries: " << ParticleTree_dum->GetEntries() << endl;
    cout << "ParticleTree_nodum entries: " << ParticleTree_nodum->GetEntries() << endl;

    Long64_t nentries_dum = ParticleTree_dum->GetEntries();
    Long64_t nentries_nodum = ParticleTree_nodum->GetEntries();
    // Condensed declarations for all per-event sums/ETMiss variables
    float ETMiss_dum=0.0, ETMiss_nodum=0.0, ETMiss_dumwithcuts=0.0,
            ETMiss_thermal=0.0,
          ETMiss_dumwithstatus1=0.0, ETMiss_nodumwithstatus1=0.0,
          ETMiss_dumwithstatus3=0.0, ETMiss_nodumwithstatus3=0.0,
          ETMiss_norec=0.0;

    float px_sum_nodum=0.0, px_sum_dum=0.0, py_sum_nodum=0.0, py_sum_dum=0.0,
          px_sum_cuts=0.0, py_sum_cuts=0.0,
          px_sum_dumstatus1=0.0, py_sum_dumstatus1=0.0,
          px_sum_nodumstatus1=0.0, py_sum_nodumstatus1=0.0,
          px_sum_dumstatus3=0.0, py_sum_dumstatus3=0.0,
          px_sum_nodumstatus3=0.0, py_sum_nodumstatus3=0.0,
          px_sum_norec=0.0, py_sum_norec=0.0;

    float px_thermal = 0.0, py_thermal = 0.0;

    for (Long64_t i=0; i<nentries_dum; i++) {
        ParticleTree_dum->GetEntry(i);
        ParticleTree_nodum->GetEntry(i); // assuming both trees have the same number of entries
        ParticleTree_norec->GetEntry(i);
        h_part_mult_dum->Fill(nParticles_dum);
        h_part_mult_nodum->Fill(nParticles_nodum);
        h_part_mult_norec->Fill(nParticles_norec);
        px_sum_dum = py_sum_dum = px_sum_nodum = py_sum_nodum =
        px_sum_cuts = py_sum_cuts =
        px_sum_dumstatus1 = py_sum_dumstatus1 = px_sum_nodumstatus1 = py_sum_nodumstatus1 =
        px_sum_dumstatus3 = py_sum_dumstatus3 = px_sum_nodumstatus3 = py_sum_nodumstatus3 =
        px_sum_norec = py_sum_norec = 0.0;

        for (size_t j=0; j<px_dum->size(); j++){
            px_sum_dum += px_dum->at(j);
            py_sum_dum += py_dum->at(j);
            h_part_status_dum->Fill(status_dum->at(j));
            if (status_dum->at(j) == 1){
                px_sum_dumstatus1 += px_dum->at(j);
                py_sum_dumstatus1 += py_dum->at(j);
            }
            if (status_dum->at(j) < 3.9 && status_dum->at(j) > 2.9){
                px_sum_dumstatus3 += px_dum->at(j);
                py_sum_dumstatus3 += py_dum->at(j);
            }
            if (e_dum->at(j) < 1e-4 && pdgId_dum->at(j) == 111){
                px_thermal += px_dum->at(j);
                py_thermal += py_dum->at(j);
            }
        }
        for (size_t j=0; j<px_nodum->size(); j++){
            px_sum_nodum += px_nodum->at(j);
            py_sum_nodum += py_nodum->at(j);
            h_part_status_nodum->Fill(status_nodum->at(j));
            if (status_nodum->at(j) == 1){
                px_sum_nodumstatus1 += px_nodum->at(j);
                py_sum_nodumstatus1 += py_nodum->at(j);
            }
        }

        for (size_t j=0; j<px_norec->size(); j++){
            px_sum_norec += px_norec->at(j);
            py_sum_norec += py_norec->at(j);
        }

        ETMiss_dum = sqrt(px_sum_dum*px_sum_dum + py_sum_dum*py_sum_dum);
        ETMiss_nodum = sqrt(px_sum_nodum*px_sum_nodum + py_sum_nodum*py_sum_nodum);
        ETMiss_norec = sqrt(px_sum_norec*px_sum_norec + py_sum_norec*py_sum_norec);
        
        ETMiss_thermal = sqrt(px_thermal*px_thermal + py_thermal*py_thermal);

        ETMiss_dumwithstatus1 = sqrt(px_sum_dumstatus1*px_sum_dumstatus1 + py_sum_dumstatus1*py_sum_dumstatus1);
        ETMiss_nodumwithstatus1 = sqrt(px_sum_nodumstatus1*px_sum_nodumstatus1 + py_sum_nodumstatus1*py_sum_nodumstatus1);
        ETMiss_dumwithstatus3 = sqrt(px_sum_dumstatus3*px_sum_dumstatus3 + py_sum_dumstatus3*py_sum_dumstatus3);

        ETMiss_norec = sqrt(px_sum_norec*px_sum_norec + py_sum_norec*py_sum_norec);

        ETMiss_dumwithcuts = ETMiss_dumwithstatus1 - ETMiss_dumwithstatus3;

        h_part_ETMiss_dumwithstatus1->Fill(ETMiss_dumwithstatus1);
        h_part_ETMiss_nodumwithstatus1->Fill(ETMiss_nodumwithstatus1);
        h_part_ETMiss_dum->Fill(ETMiss_dum);
        h_part_ETMiss_nodum->Fill(ETMiss_nodum);
        h_part_ETMiss_dumwithstatus3->Fill(ETMiss_dumwithstatus3);
        h_part_ETMiss_norec->Fill(ETMiss_norec);
        h_part_ETMiss_dumwithcuts->Fill(ETMiss_dumwithcuts);
    }

    //vector<TH1F*> hists = {h_part_ETMiss_dum, h_part_ETMiss_nodum,h_part_ETMiss_dumwithcuts};
    //vector<string> labels = {"Jewel PbPb WriteDummies", "Jewel PbPb no WriteDummies", "Jewel PbPb Write Dummies with cuts"};

    vector<TH1F*> hists = {h_part_ETMiss_nodum,
                            h_part_ETMiss_dum,
                            h_part_ETMiss_dumwithcuts};
    vector<string> labels = {"WriteDummies with Recoil",
                                "No Dummies with Recoil",
                                "Write Dummies (status 1 - status 3)"};

 /*   plotHistograms({h_part_mult_nodum,h_part_mult_dum, h_part_mult_norec}, 
                    {"Jewel PbPb No WriteDummies", "Jewel PbPb WriteDummies", "Jewel PbPb WriteDummies without Recoil"}, 
                    "particle_multiplicity", "Particle Multiplicity", 
                    true, 
                    "0203_PbPb_Jewel_10KEvt_Dummies_Multiplicity_Comparison",
                    "Comparison of Particle Multiplicity in PbPb Jewel 10K Fixed Debug Seed");
*/
    plotHistograms(hists, labels, "particle_ETMiss", "E_{T}^{Miss} (GeV)", true, 
        "0203_PbPb_Jewel_10KEvt_Dummies_Comparison_Bkgsub",
        "Comparison of E_{T}^{Miss} in PbPb Jewel 10K BkgSubtraction");

/*    plotHistograms({h_part_ETMiss_nodum,h_part_ETMiss_dum}, 
                    {"Jewel PbPb No WriteDummies", "Jewel PbPb WriteDummies"}, 
                    "particle_ETMiss", "E_{T}^{Miss} (GeV)", 
                    true, 
                    "0201_PbPb_Jewel_10KEvt_Dummies_Comparison_twoplots");

    plotHistograms({h_part_status_nodum, h_part_status_dum}, 
                    {"Jewel PbPb No WriteDummies", "Jewel PbPb WriteDummies"}, 
                    "particle_status", "Particle Status", true, "0201_PbPb_Jewel_10KEvt_Status_Comparison_");
*/
    f_nodum->Close();
    f_dum->Close();

    return 0;
}