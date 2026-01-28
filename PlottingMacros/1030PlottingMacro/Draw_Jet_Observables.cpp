auto plotHistograms(TH1F* h_vac, TH1F* h_med, const string& varname, const string& xtitle, bool norm = false, double ymax = 0.0) {
    // Work on clones when normalizing so original histograms are not modified
    TH1F* hv = h_vac;
    TH1F* hm = h_med;
    bool cloned = false;
    if (norm) {
        hv = (TH1F*)h_vac->Clone((string(h_vac->GetName()) + "_norm").c_str());
        hm = (TH1F*)h_med->Clone((string(h_med->GetName()) + "_norm").c_str());
        cloned = true;
        double iv = hv->Integral();
        double im = hm->Integral();
        if (iv > 0) hv->Scale(1.0 / iv);
        if (im > 0) hm->Scale(1.0 / im);
    }

    // Determine y-axis maximum: use provided ymax if >0, otherwise compute from histograms
    double maxy = 0.0;
    if (ymax > 0.0) {
        maxy = ymax;
    } else {
        double mv = hv->GetMaximum();
        double mm = hm->GetMaximum();
        maxy = (mv > mm ? mv : mm);
        // Add a small margin so the curves are not flush to the top
        maxy *= 1.2;
        if (maxy == 0.0) maxy = 1.0; // avoid zero max
    }

    TCanvas* canvas_temp = new TCanvas(("c_" + varname).c_str(), "canvas", 800, 800);
    canvas_temp->SetLeftMargin(0.14);  // Increase left margin for y-axis label
    canvas_temp->SetBottomMargin(0.12); // Ensure bottom margin is adequate
    canvas_temp->SetTopMargin(0.08);    // Set top margin
    canvas_temp->SetRightMargin(0.05);  // Set right margin
    
    hv->SetLineColor(kBlue+1);
    hv->SetLineWidth(2);
    hv->SetXTitle(xtitle.c_str());
    hv->SetYTitle(norm ? "Normalized entries" : "Entries");
    hv->SetMaximum(maxy);
    hv->Draw("HIST");

    hm->SetLineColor(kRed+1);
    hm->SetLineWidth(2);
    // hm will be drawn on same axes; no need to call SetMaximum again
    hm->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.6, 0.8, 0.75, 0.9);
    legend->AddEntry(hv, "Jewel PbPb without quenching", "l");
    legend->AddEntry(hm, "Jewel PbPb with quenching", "l");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.025);
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.025);
    latex.DrawLatex(0.20, 0.85, "N_{Gen} = 10K, #sqrt{s_{NN}} = 5.36 TeV");
    latex.DrawLatex(0.20, 0.80, "Anti-k_{T}, R=0.4");
    latex.DrawLatex(0.20, 0.75, "p_{T1} > 50 GeV/c, p_{T2} > 30 GeV/c");
    latex.DrawLatex(0.20, 0.70, "#Delta #Phi > (5/6) #pi");

    string filename = "/home/xirong/DijetAnalysis_2025_v3_svmit/Plots/1030Jetlevel/" + varname + "_PbPb_Jewel_10KEvt_comparemedvsvac_011926_withcuts.png";
    canvas_temp->SaveAs(filename.c_str());

    delete legend;
    delete canvas_temp;
    if (cloned) {
        delete hv;
        delete hm;
    }
};

int Draw_Jet_Observables(){

    string TreeName = "Tree";
    float pt1 = -1;
    float pt2 = -1;
    float A_J = 0;
    float ET_miss = 0;

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 800);

    // Define histogram parameters
    
    //TH1F* h_part_EMiss = new TH1F("h_part_EMiss", "Particle ETMiss of PbPb Jewel", 100, 0, 200);
    TFile *f_vac = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/PlottingMacros/1014JewelWithJetClustering/RootOutput/Jewel_pbpb_5360GeV_ptmin250_nevt10K_012026_vac_withdphi.root");         
    TFile *f_med = TFile::Open("/home/xirong/DijetAnalysis_2025_v3_svmit/PlottingMacros/1014JewelWithJetClustering/RootOutput/Jewel_pbpb_5360GeV_ptmin250_nevt10K_012026_med_withdphi.root");         
   
    TTree *JetTree_vac = (TTree*)f_vac->Get("JetTree");
    TTree *JetTree_med = (TTree*)f_med->Get("JetTree");

    TH1F* h_jet_pt_med = new TH1F("h_pt_med", " pt of 10K PbPb Event  using Jewel", 100, 0, 1000);
    TH1F* h_jet_pt_vac = new TH1F("h_pt_vac", " pt of 10K PbPb Event  using Jewel", 100, 0, 1000);
    TH1F* h_jet_eta_med = new TH1F("h_eta_med", " #eta 10K PbPb Event  using Jewel", 100, -3, 3);
    TH1F* h_jet_eta_vac = new TH1F("h_eta_vac", " #eta 10K PbPb Event  using Jewel", 100, -3, 3);
    TH1F* h_jet_phi_med = new TH1F("h_phi_med", " phi of 10K PbPb Event  using Jewel", 100, -TMath::Pi(), TMath::Pi());
    TH1F* h_jet_phi_vac = new TH1F("h_phi_vac", " phi of 10K PbPb Event  using Jewel", 100, -TMath::Pi(), TMath::Pi());

    TH1F* h_pt1_vac = new TH1F("h_pt1_vac", "leading jet p_{T} 10K PbPb Event using Jewel", 100, 0, 1000);
    TH1F* h_pt1_med = new TH1F("h_pt1_med", "leading jet p_{T} 10K PbPb Event using Jewel", 100, 0, 1000);
    TH1F* h_pt2_vac = new TH1F("h_pt2_vac", "subleading jet p_{T} 10K PbPb Event using Jewel", 100, 0, 1000);
    TH1F* h_pt2_med = new TH1F("h_pt2_med", "subleading jet p_{T} 10K PbPb Event using Jewel", 100, 0, 1000);
    TH1F* h_AJ_vac = new TH1F("h_AJ_vac", "A_{J} (vac)", 100, 0.0, 1.0);
    TH1F* h_AJ_med = new TH1F("h_AJ_med", "A_{J} (med)", 100, 0.0, 1.0);

    TH1F* h_EMiss_vac = new TH1F("h_EMiss_vac", "Jet E_{T}^{miss} 10K PbPb Event using Jewel", 100, 0, 200);
    TH1F* h_EMiss_med = new TH1F("h_EMiss_med", "JetE_{T}^{miss} 10K PbPb Event using Jewel", 100, 0, 200);
    TH1F* h_dPhi_vac = new TH1F("h_dPhi_vac", "#Delta#phi 10K PbPb Event using Jewel", 100, -4, 4);
    TH1F* h_dPhi_med = new TH1F("h_dPhi_med", "#Delta#phi 10K PbPb Event using Jewel", 100, -4, 4);

 
    // Variables for jet vectors - IMPORTANT: Initialize to nullptr to prevent segfaults
    vector<double>* jet_pt_vac = nullptr;
    vector<double>* jet_eta_vac = nullptr;
    vector<double>* jet_phi_vac = nullptr;
    vector<double>* jet_pt_med = nullptr;
    vector<double>* jet_eta_med = nullptr;
    vector<double>* jet_phi_med = nullptr;

    double A_J_vac = 0;
    double A_J_med = 0;
    double jet_pt1_vac = 0;
    double jet_pt2_vac = 0;
    double jet_pt1_med = 0;
    double jet_pt2_med = 0;
    double dPhi_vac = 0;
    double dPhi_med = 0;
    float ET_miss_jet_vac = 0;
    float ET_miss_jet_med = 0;
    JetTree_vac->SetBranchAddress("ET_miss_jet",&ET_miss_jet_vac);
    JetTree_med->SetBranchAddress("ET_miss_jet",&ET_miss_jet_med);
    JetTree_vac->SetBranchAddress("A_J",&A_J_vac);
    JetTree_med->SetBranchAddress("A_J",&A_J_med);
    JetTree_vac->SetBranchAddress("dPhi",&dPhi_vac);
    JetTree_med->SetBranchAddress("dPhi",&dPhi_med);
    JetTree_vac->SetBranchAddress("jetpt1",&jet_pt1_vac);
    JetTree_vac->SetBranchAddress("jetpt2",&jet_pt2_vac);
    JetTree_med->SetBranchAddress("jetpt1",&jet_pt1_med);
    JetTree_med->SetBranchAddress("jetpt2",&jet_pt2_med);
    JetTree_vac->SetBranchAddress("jet_pt",&jet_pt_vac);
    JetTree_med->SetBranchAddress("jet_pt",&jet_pt_med);
    JetTree_vac->SetBranchAddress("jet_eta",&jet_eta_vac);
    JetTree_med->SetBranchAddress("jet_eta",&jet_eta_med);
    JetTree_med->SetBranchAddress("jet_phi",&jet_phi_med);
    JetTree_vac->SetBranchAddress("jet_phi",&jet_phi_vac);

    Long64_t nEntries = JetTree_vac->GetEntries();
    for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
            JetTree_vac->GetEntry(entrynum);
            if (jet_pt1_vac < 50) continue;
            if (jet_pt2_vac < 30) continue; // Apply jet pt cuts
            if (dPhi_vac < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
           //  cout << jet_pt_vac->size() << endl;
           // Add null pointer checks BEFORE accessing size() to prevent segfaults
           h_dPhi_vac->Fill(dPhi_vac);
           h_EMiss_vac->Fill(ET_miss_jet_vac);
        }

    Long64_t nEntries_med = JetTree_med->GetEntries();
    for (Long64_t entrynum = 0; entrynum < nEntries_med; entrynum++){
            JetTree_med->GetEntry(entrynum);
            if (jet_pt1_med < 50) continue;
            if (jet_pt2_med < 30) continue;
            if (dPhi_med < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
            // Fill histograms - Add null pointer checks to prevent segfaults
           h_dPhi_med->Fill(dPhi_med);
           h_EMiss_med->Fill(ET_miss_jet_med);
    }
    cout << h_dPhi_med->GetEntries() << endl;
    cout << h_dPhi_vac->GetEntries() << endl;
    cout << h_EMiss_med->GetEntries() << endl;
    cout << h_EMiss_vac->GetEntries() << endl;

    gStyle->SetOptStat(0); // Disable statistics box
    plotHistograms(h_dPhi_vac, h_dPhi_med, "dPhi", "#Delta#phi",true);
    plotHistograms(h_EMiss_vac, h_EMiss_med, "ETMiss", "Jet E_{T}^{miss}",true);

    f_vac->Close();
    f_med->Close();

    return 0;
}