int Draw() {

    string TreeName = "Tree";
    float pt1 = -1;
    float pt2 = -1;
    float A_J = 0;
    float ET_miss = 0;

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 800);

    TH1F* h_jet_AJ_vac = new TH1F("h_AJ_vac", " A_{J} of PbPb Jewel", 100, 0, 1);
    TH1F* h_jet_AJ_med = new TH1F("h_AJ_med", " A_{J} of PbPb Jewel", 100, 0, 1);
    TH1F* h_jet_pt1 = new TH1F("h_jet_pt1", " Leading Jet p_{T} of PbPb Jewel", 100, 0, 1000);

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

    JetTree_vac->SetBranchAddress("A_J",&A_J_vac);
    JetTree_med->SetBranchAddress("A_J",&A_J_med);
    JetTree_vac->SetBranchAddress("jetpt1",&jet_pt1_vac);
    JetTree_vac->SetBranchAddress("jetpt2",&jet_pt2_vac);
    JetTree_med->SetBranchAddress("jetpt1",&jet_pt1_med);
    JetTree_med->SetBranchAddress("jetpt2",&jet_pt2_med);

    Long64_t nEntries = JetTree_vac->GetEntries();
    for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
            JetTree_vac->GetEntry(entrynum);
            if (jet_pt1_vac < 50) continue;
            if (jet_pt2_vac < 30) continue; // Apply jet pt cuts
           // if (dPhi < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
            h_jet_AJ_vac->Fill(A_J_vac);
    }

    nEntries = JetTree_med->GetEntries();
    for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
            JetTree_med->GetEntry(entrynum);
            if (jet_pt1_med < 50) continue;
            if (jet_pt2_med < 30) continue;
           // if (pt2 < 30) continue; // Apply jet pt cuts
           // if (dPhi < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
            h_jet_AJ_med->Fill(A_J_med);
    }
        f_vac->Close();
        f_med->Close();

    gStyle->SetOptStat(0); // Disable statistics box

    h_jet_AJ_vac->SetLineColor(kBlue+1);
    h_jet_AJ_vac->SetLineWidth(2);
    h_jet_AJ_vac->Draw("HIST");
    h_jet_AJ_vac->SetXTitle("A_{J}");
    h_jet_AJ_med->Draw("HIST SAME");
    h_jet_AJ_med->SetLineColor(kRed+1);
    h_jet_AJ_med->SetLineWidth(2);

    TLegend* legend = new TLegend(0.55, 0.7, 0.75, 0.8);
    legend->AddEntry(h_jet_AJ_vac, "Jewel PbPb without  quenching", "l");
    legend->AddEntry(h_jet_AJ_med, "Jewel PbPb with quenching", "l");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.15, 0.85, "NEvents #= 10K, #sqrt{s_{NN}} = 5.36 TeV");
    latex.DrawLatex(0.15, 0.80, "Anti-k_{T}, R=0.4, p_{T1} > 50 GeV/c, p_{T2} > 30 GeV/c");
    
    canvas->SaveAs("/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/1030Jetlevel/AJ_PbPb_Jewel_10KEvt_comparemedvsvac_011926_withcuts.png");
    return 0;
}