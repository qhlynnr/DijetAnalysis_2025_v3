#include "TFile.h"
#include "TTree.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char** argv) {
    if (argc < 3){
        std::cerr << "Usage: " << argv[0] << " <input.root> <output.root>" << std::endl;
        return 1;
    }

    std::string inFileName = argv[1];
    std::string outFileName = argv[2];

    TFile *fJ = TFile::Open(inFileName.c_str());
    if (!fJ || fJ->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inFileName << std::endl;
        return 1;
    }
    TTree* tJ = (TTree*)fJ->Get("ParticleTree");
    if (!tJ) {
        std::cerr << "Error: Cannot find ParticleTree in " << inFileName << std::endl;
        return 1;
    }

    TFile *outFile = TFile::Open(outFileName.c_str(), "RECREATE");
    TTree *jetTree = new TTree("JetTree", "Embedded particle-level events");
    TTree *particleTree = new TTree("ParticleTree", "Particle-level information from embedded events");

    Long64_t nEntriesJ = tJ->GetEntries();
    std::cout << "Number of entries in " << inFileName << ": " << nEntriesJ << std::endl;

    //Input Variables
    
    int nJ = 0;
    vector<double>* pxJ = nullptr;
    vector<double>* pyJ = nullptr;
    vector<double>* pzJ = nullptr;
    vector<double>* eJ = nullptr;
    vector<double>* ptJ = nullptr;
    vector<double>* etaJ = nullptr;
    vector<double>* phiJ = nullptr;
    vector<int>* pdgIdJ = nullptr;

    // Jet-level output vectors
    double jetpt1 = 0;
    double jetpt2 = 0;
    double jetphi1 = 0;
    double jetphi2 = 0;
    double dPhi = 0;
    double A_J = -1;
    double X_J = -1;
    double ET_miss_jet = 0;
    double ET_miss_particle = 0;
    const double pi = 3.14159265358979323846;
    double weights;

    vector<double>* jetPx = nullptr;
    vector<double>* jetPy = nullptr;
    vector<double>* jetPz = nullptr;
    vector<double>* jetE = nullptr;
    vector<double>* jetPt = nullptr;
    vector<double>* jetEta = nullptr;
    vector<double>* jetPhi = nullptr;
    vector<double>* partstatus = nullptr;

    vector<double>* partPx = nullptr;
    vector<double>* partPy = nullptr;
    vector<double>* partPz = nullptr;
    vector<double>* partE = nullptr;
    vector<double>* partPt = nullptr;
    vector<double>* partEta = nullptr;
    vector<double>* partPhi = nullptr;
    vector<int>* partPdgId = nullptr; // Combined PDG IDs
    vector<int>* partSource = nullptr; // 0 = Angantyr, 1 = JEWEL
    double partWeights;

    // jet def and bkg subtraction
    const double R = 0.3;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  
    double jetPtMin = 5.0;
    tJ->SetBranchAddress("nParticles", &nJ);
    tJ->SetBranchAddress("px", &pxJ);
    tJ->SetBranchAddress("py", &pyJ);
    tJ->SetBranchAddress("pz", &pzJ);
    tJ->SetBranchAddress("e",  &eJ);  // assuming JEWEL also uses lowercase 'e'
    tJ->SetBranchAddress("pt", &ptJ);
    tJ->SetBranchAddress("eta", &etaJ);
    tJ->SetBranchAddress("phi", &phiJ);
    tJ->SetBranchAddress("pdgId", &pdgIdJ);
    tJ->SetBranchAddress("status", &partstatus);
    tJ->SetBranchAddress("weights", &weights);

    jetTree->Branch("jet_px", &jetPx);
    jetTree->Branch("jet_py", &jetPy);
    jetTree->Branch("jet_pz", &jetPz);
    jetTree->Branch("jet_e",  &jetE);
    jetTree->Branch("jet_pt", &jetPt);
    jetTree->Branch("jet_eta", &jetEta);
    jetTree->Branch("jet_phi", &jetPhi);
    jetTree->Branch("A_J", &A_J);
    jetTree->Branch("X_J", &X_J);
    jetTree->Branch("jetpt1", &jetpt1);
    jetTree->Branch("jetpt2", &jetpt2);
    jetTree->Branch("jetphi1", &jetphi1);
    jetTree->Branch("jetphi2", &jetphi2);
    jetTree->Branch("ET_miss_jet", &ET_miss_jet);
    jetTree->Branch("dPhi", &dPhi);

    // Particle-level branches
    particleTree->Branch("px", &partPx);
    particleTree->Branch("py", &partPy);
    particleTree->Branch("pz", &partPz);
    particleTree->Branch("e",  &partE);
    particleTree->Branch("pt", &partPt);
    particleTree->Branch("eta", &partEta);
    particleTree->Branch("phi", &partPhi);
    particleTree->Branch("pdgId", &partPdgId);
    particleTree->Branch("ET_miss_particle", &ET_miss_particle);
    particleTree->Branch("source", &partSource); // 0 = Angantyr, 1 = JEWEL
    particleTree->Branch("weights", &partWeights);

    std::vector<fastjet::PseudoJet> particles;

    for (int i = 0; i < nEntriesJ; ++i) {
        if (i%1000 == 0) {
            std::cout << "Processing event " << i << " / " << nEntriesJ << std::endl;
        }
        particles.clear();
        tJ->GetEntry(i);

        partWeights = weights; // Store event-level weight for all particles in this event
        
        // Initialize momentum sums for this event
        double partpx_sum = 0, partpy_sum = 0, partpz_sum = 0;
        double jetpx_sum = 0, jetpy_sum = 0, jetpz_sum = 0;
        
        // Clear particle-level output vectors
        partPx->clear(); partPy->clear(); partPz->clear(); partE->clear(); 
        partPt->clear(); partEta->clear(); partPhi->clear(); 
        partPdgId->clear(); partSource->clear();
        
        if (i%1000 == 0) {
            std::cout << "  Number of JEWEL particles: " << nJ << std::endl;
        }

        // Add JEWEL particles
        for (int j = 0; j < nJ; ++j) {
            int pid = pdgIdJ->at(j);
            if (std::abs(pid) == 12 || std::abs(pid) == 14 || std::abs(pid) == 16) continue; // skip neutrinos (12,14,16)
            if (partstatus->at(j) != 1) continue; // only final state particles
            particles.push_back(fastjet::PseudoJet(pxJ->at(j), pyJ->at(j), pzJ->at(j), eJ->at(j)));
            // Store particle-level info
            partPx->push_back(pxJ->at(j));
            partPy->push_back(pyJ->at(j));
            partPz->push_back(pzJ->at(j));
            partE->push_back(eJ->at(j));
            partPt->push_back(ptJ->at(j));
            partEta->push_back(etaJ->at(j));
            partPhi->push_back(phiJ->at(j));
            partPdgId->push_back(pdgIdJ->at(j));
            partSource->push_back(1); // 1 = JEWEL

            partpx_sum += pxJ->at(j);
            partpy_sum += pyJ->at(j);
            partpz_sum += pzJ->at(j);
        }
        if (i%1000 == 0) {
            std::cout << "  Total number of particles for clustering: " << particles.size() << std::endl;
        }

        // Perform jet clustering for this event (subtract first, then sort by pT)
        fastjet::ClusterSequence cs(particles, jet_def);
        std::vector<fastjet::PseudoJet> jets_raw = cs.inclusive_jets(); // get all raw jets
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(jets_raw);

        // Clear output vectors and fill with jet 4-momenta (apply jetPtMin after subtraction)
        jetPx->clear(); jetPy->clear(); jetPz->clear(); jetE->clear();
        jetPt->clear(); jetEta->clear(); jetPhi->clear();
        for (const auto& J : jets) {
            if (J.pt() < jetPtMin) continue; // enforce pt cut after subtraction
            jetPx->push_back(J.px());
            jetPy->push_back(J.py());
            jetPz->push_back(J.pz());
            jetE->push_back(J.E());
            jetPt->push_back(J.pt());
            jetEta->push_back(J.eta());
            jetPhi->push_back(J.phi_std());
            jetpx_sum += J.px();
            jetpy_sum += J.py();
            jetpz_sum += J.pz();
        }

        jetpt1 = (jetPt->size() > 0) ? jetPt->at(0) : -1.0;
        jetpt2 = (jetPt->size() > 1) ? jetPt->at(1) : -1.0;
        jetphi1 = (jetPhi->size() > 0) ? jetPhi->at(0) : -999.0;
        jetphi2 = (jetPhi->size() > 1) ? jetPhi->at(1) : -999.0;

        if (jetpt1 > 0 && jetpt2 > 0) {
            if (jetpt1 < jetpt2){
                cout << "Warning: jetpt1 < jetpt2" << endl;
                cout << "jetpt1: " << jetpt1 << ", jetpt2: " << jetpt2 << endl;
            }
            A_J = (jetpt1 - jetpt2) / (jetpt1 + jetpt2);
            X_J = jetpt2 / jetpt1;
        } else {
            A_J = -1; // Undefined if less than 2 jets
            X_J = -1;
        }

        dPhi = jetphi1 - jetphi2;

        if (pi <= dPhi && dPhi <= 2 * pi){
            dPhi = (dPhi - 2*pi);
        }
        else if (-2*pi <= dPhi && dPhi <= -pi){
            dPhi = dPhi + 2*pi;
        }

        ET_miss_jet = sqrt(jetpx_sum*jetpx_sum + jetpy_sum*jetpy_sum);
        ET_miss_particle = sqrt(partpx_sum*partpx_sum + partpy_sum*partpy_sum);

        jetTree->Fill();
        particleTree->Fill();

        if (i % 100 == 0) {
            std::cout << "Processed " << i << " events." << std::endl;
            std::cout << "Event weight: " << partWeights << std::endl;
        }
        if (i<5){
            std::cout << "Event " << i << " summary:" << std::endl;
            std::cout << "  Number of particles: " << nJ << std::endl;
            std::cout << "  Number of jets found: " << jetPt->size() << std::endl;
            std::cout << "  Leading jet pt: " << jetpt1 << std::endl;
            std::cout << "  Subleading jet pt: " << jetpt2 << std::endl;
            std::cout << "  A_J: " << A_J << std::endl;
            std::cout << "  X_J: " << X_J << std::endl;
            std::cout << "  dPhi: " << dPhi << std::endl;
            std::cout << "  ET_miss_jet: " << ET_miss_jet << std::endl;
            std::cout << "  ET_miss_particle: " << ET_miss_particle << std::endl;
        }       
    }


    // Write and close files
    outFile->Write();
    outFile->Close();
    fJ->Close();


    return 0;
}
