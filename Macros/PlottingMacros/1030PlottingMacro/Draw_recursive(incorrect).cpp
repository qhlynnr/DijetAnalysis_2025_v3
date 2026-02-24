auto plotHistograms(TH1F* h_vac, TH1F* h_med, const string& varname, const string& xtitle) {
        TCanvas* canvas_temp = new TCanvas(("c_" + varname).c_str(), "canvas", 800, 800);

        h_vac->SetLineColor(kBlue+1);
        h_vac->SetLineWidth(2);
        h_vac->Draw("HIST");
        h_vac->SetXTitle(xtitle.c_str());
        h_med->Draw("HIST SAME");
        h_med->SetLineColor(kRed+1);
        h_med->SetLineWidth(2);

        TLegend* legend = new TLegend(0.55, 0.7, 0.75, 0.8);
        legend->AddEntry(h_vac, "Jewel PbPb without quenching", "l");
        legend->AddEntry(h_med, "Jewel PbPb with quenching", "l");
        legend->SetBorderSize(0);
        legend->SetTextSize(0.03);
        legend->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.15, 0.85, "NEvents #= 10K, #sqrt{s_{NN}} = 5.36 TeV");
        latex.DrawLatex(0.15, 0.80, "Anti-k_{T}, R=0.4, p_{T1} > 50 GeV/c, p_{T2} > 30 GeV/c");
        
        string filename = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/1030Jetlevel/" + varname + "_PbPb_Jewel_10KEvt_comparemedvsvac_011926_withcuts.png";
        canvas_temp->SaveAs(filename.c_str());
        delete canvas_temp;
    };

int Draw(){

    string TreeName = "Tree";
    float pt1 = -1;
    float pt2 = -1;
    float A_J = 0;
    float ET_miss = 0;

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 800);

    // Define histogram parameters
    struct HistConfig {
        string name;
        string title;
        string xtitle;
        int nbins;
        double xmin;
        double xmax;
    };
    
    vector<HistConfig> histConfigs = {
        {"AJ", "A_{J} of PbPb Jewel", "A_{J}", 100, 0, 1},
        {"jetpt1", "Leading Jet p_{T} of PbPb Jewel", "Leading Jet p_{T} [GeV/c]", 100, 0, 300},
        {"jetpt2", "Subleading Jet p_{T} of PbPb Jewel", "Subleading Jet p_{T} [GeV/c]", 100, 0, 200},
        {"jet_eta1", "Leading Jet #eta of PbPb Jewel", "Leading Jet #eta", 100, -3, 3},
        {"jet_eta2", "Subleading Jet #eta of PbPb Jewel", "Subleading Jet #eta", 100, -3, 3},
        {"jet_phi1", "Leading Jet #phi of PbPb Jewel", "Leading Jet #phi", 100, -TMath::Pi(), TMath::Pi()},
        {"jet_phi2", "Subleading Jet #phi of PbPb Jewel", "Subleading Jet #phi", 100, -TMath::Pi(), TMath::Pi()},
        {"ETmiss", "E_{T}^{miss} of PbPb Jewel", "E_{T}^{miss} [GeV]", 100, 0, 100}
    };
    
    // Create histograms using loops
    map<string, TH1F*> histograms_vac;
    map<string, TH1F*> histograms_med;
    
    for (const auto& config : histConfigs) {
        string name_vac = "h_" + config.name + "_vac";
        string name_med = "h_" + config.name + "_med";
        
        histograms_vac[config.name] = new TH1F(name_vac.c_str(), config.title.c_str(), 
                                               config.nbins, config.xmin, config.xmax);
        histograms_med[config.name] = new TH1F(name_med.c_str(), config.title.c_str(), 
                                               config.nbins, config.xmin, config.xmax);
    }

    //TH1F* h_part_EMiss = new TH1F("h_part_EMiss", "Particle ETMiss of PbPb Jewel", 100, 0, 200);
    TFile *f_vac = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/PlottingMacros/1014JewelWithJetClustering/RootOutput/Jewel_pbpb_5360GeV_ptmin250_nevt10K_103125_vac.root");         
    TFile *f_med = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/PlottingMacros/1014JewelWithJetClustering/RootOutput/Jewel_Pbpb_5360GeV_250pthat_nevt10K_101625_med_withjets.root");         

    TTree *JetTree_vac = (TTree*)f_vac->Get("JetTree");
    TTree *JetTree_med = (TTree*)f_med->Get("JetTree");

    double A_J_vac = 0;
    double A_J_med = 0;
    double jet_pt1_vac = 0;
    double jet_pt2_vac = 0;
    double jet_pt1_med = 0;
    double jet_pt2_med = 0;
    
    // Variables for jet eta and phi vectors
    vector<double>* jet_eta_vac = 0;
    vector<double>* jet_phi_vac = 0;
    vector<double>* jet_eta_med = 0;
    vector<double>* jet_phi_med = 0;
    
    // Variable for ET_miss_jet
    float ET_miss_jet_vac = 0;
    float ET_miss_jet_med = 0;

    JetTree_vac->SetBranchAddress("A_J",&A_J_vac);
    JetTree_med->SetBranchAddress("A_J",&A_J_med);
    JetTree_vac->SetBranchAddress("jetpt1",&jet_pt1_vac);
    JetTree_vac->SetBranchAddress("jetpt2",&jet_pt2_vac);
    JetTree_med->SetBranchAddress("jetpt1",&jet_pt1_med);
    JetTree_med->SetBranchAddress("jetpt2",&jet_pt2_med);
    JetTree_vac->SetBranchAddress("jet_eta",&jet_eta_vac);
    JetTree_vac->SetBranchAddress("jet_phi",&jet_phi_vac);
    JetTree_med->SetBranchAddress("jet_eta",&jet_eta_med);
    JetTree_med->SetBranchAddress("jet_phi",&jet_phi_med);
    JetTree_vac->SetBranchAddress("ET_miss_jet",&ET_miss_jet_vac);
    JetTree_med->SetBranchAddress("ET_miss_jet",&ET_miss_jet_med);

    Long64_t nEntries = JetTree_vac->GetEntries();
    for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
            JetTree_vac->GetEntry(entrynum);
            if (jet_pt1_vac < 50) continue;
            if (jet_pt2_vac < 30) continue; // Apply jet pt cuts
           // if (dPhi < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
            
            // Fill histograms
            histograms_vac["AJ"]->Fill(A_J_vac);
            histograms_vac["jetpt1"]->Fill(jet_pt1_vac);
            histograms_vac["jetpt2"]->Fill(jet_pt2_vac);
            histograms_vac["ETmiss"]->Fill(ET_miss_jet_vac);
    }

    nEntries = JetTree_med->GetEntries();
    for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
            JetTree_med->GetEntry(entrynum);
            if (jet_pt1_med < 50) continue;
            if (jet_pt2_med < 30) continue;
           // if (dPhi < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
            // Fill histograms
            histograms_med["AJ"]->Fill(A_J_med);
            histograms_med["jetpt1"]->Fill(jet_pt1_med);
            histograms_med["jetpt2"]->Fill(jet_pt2_med);
            histograms_med["ETmiss"]->Fill(ET_miss_jet_med);
    }
        f_vac->Close();
        f_med->Close();

    gStyle->SetOptStat(0); // Disable statistics box

    // Plot all variables using loop
    for (const auto& config : histConfigs) {
        plotHistograms(histograms_vac[config.name], histograms_med[config.name], 
                      config.name, config.xtitle);
    }
    return 0;
}