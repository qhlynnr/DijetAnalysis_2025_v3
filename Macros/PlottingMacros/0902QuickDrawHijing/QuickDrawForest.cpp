#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include "CommandLine.h"
using namespace std;

int main(int argc, char* argv[]) {
    CommandLine CL(argc, argv);
    string inputFile = CL.Get("Input", "default.txt");
    bool debug = CL.GetBool("Debug", false);

    cout << "Input file: " << inputFile << endl;

    TFile* inFile = new TFile(inputFile.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return 1;
    }
    string JetAnalyserTreeString = "ak4PFJetAnalyzer/t";

    TTree* JetTree = (TTree*)inFile->Get(JetAnalyserTreeString.c_str());
    if (!JetTree) {
        cout << "Error: One or more TTrees not found in file " << inputFile << endl;
        inFile->Close();
        delete inFile;
        return 1;
    }
    
    int nref;
    float jtpt[500], jteta[500], jtphi[500];

    JetTree->SetBranchStatus("*", 0);
    JetTree->SetBranchStatus("nref", 1);
    JetTree->SetBranchStatus("jtpt", 1);
    JetTree->SetBranchStatus("jteta", 1);
    JetTree->SetBranchStatus("jtphi", 1);

    JetTree->SetBranchAddress("nref",&nref);
    JetTree->SetBranchAddress("jtpt",jtpt);
    JetTree->SetBranchAddress("jteta",jteta);
    JetTree->SetBranchAddress("jtphi",jtphi);

    Long64_t nEntries = JetTree->GetEntries();
    cout << "Total entries in the tree: " << nEntries << endl;
    for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
        JetTree->GetEntry(entrynum);
        if (debug && entrynum >= 10000) break; // Limit to first 10 entries if debug is true
        cout << "Entry " << entrynum << ": nref = " << nref << endl;
        for (int i = 0; i < nref; i++) {
            cout << "  Jet " << i << ": pt = " << jtpt[i] 
                 << ", eta = " << jteta[i] 
                 << ", phi = " << jtphi[i] << endl;
        }
    }
    
    return 0;
}