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
inline double deltaR2(double eta1, double phi1, double eta2, double phi2){
    double dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    double deta = eta1 - eta2;
    return deta*deta + dphi*dphi;
}
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
    double ET_miss_jet_raw = 0;
    double ET_miss_jet_thermal = 0;
    const double pi = 3.14159265358979323846;
    double weights;

    vector<double>* jetPx = nullptr;
    vector<double>* jetPy = nullptr;
    vector<double>* jetPz = nullptr;
    vector<double>* jetE = nullptr;
    vector<double>* jetPt = nullptr;
    vector<double>* jetEta = nullptr;
    vector<double>* jetPhi = nullptr;

    vector<double>* jetPxRaw = nullptr;
    vector<double>* jetPyRaw = nullptr; 
    vector<double>* jetPzRaw = nullptr;
    vector<double>* jetERaw = nullptr;
    vector<double>* jetPtRaw = nullptr;
    vector<double>* jetEtaRaw = nullptr;
    vector<double>* jetPhiRaw = nullptr;

    vector<double>* jetPxThermal = nullptr;
    vector<double>* jetPyThermal = nullptr;
    vector<double>* jetPzThermal = nullptr;
    vector<double>* jetEThermal = nullptr;
    vector<double>* jetPtThermal = nullptr;
    vector<double>* jetEtaThermal = nullptr;
    vector<double>* jetPhiThermal = nullptr;

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
    double partWeights = 1.0;

    // jet def and bkg subtraction
    const double R = 03;
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
    jetTree->Branch("ET_miss_jet_raw", &ET_miss_jet_raw);
    jetTree->Branch("ET_miss_jet_thermal", &ET_miss_jet_thermal);
    jetTree->Branch("dPhi", &dPhi);

    jetTree->Branch("jet_px_raw", &jetPxRaw);
    jetTree->Branch("jet_py_raw", &jetPyRaw);
    jetTree->Branch("jet_pz_raw", &jetPzRaw);
    jetTree->Branch("jet_e_raw",  &jetERaw);
    jetTree->Branch("jet_pt_raw", &jetPtRaw);
    jetTree->Branch("jet_eta_raw", &jetEtaRaw);
    jetTree->Branch("jet_phi_raw", &jetPhiRaw);

    jetTree->Branch("jet_px_thermal", &jetPxThermal);
    jetTree->Branch("jet_py_thermal", &jetPyThermal);
    jetTree->Branch("jet_pz_thermal", &jetPzThermal);
    jetTree->Branch("jet_e_thermal",  &jetEThermal);
    jetTree->Branch("jet_pt_thermal", &jetPtThermal);
    jetTree->Branch("jet_eta_thermal", &jetEtaThermal);
    jetTree->Branch("jet_phi_thermal", &jetPhiThermal);


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
    std::vector<fastjet::PseudoJet> thermalmomenta;

    for (int i = 0; i < nEntriesJ; ++i) {
        if (i%1000 == 0) {
            std::cout << "Processing event " << i << " / " << nEntriesJ << std::endl;
        }
        particles.clear();
        thermalmomenta.clear();
        tJ->GetEntry(i);

        partWeights = weights; // Store event-level weight for all particles in this event
        
        // Initialize momentum sums for this event
        double partpx_sum = 0, partpy_sum = 0, partpz_sum = 0;
        double jetpx_sum = 0, jetpy_sum = 0, jetpz_sum = 0;
        double jetPxRaw_sum = 0, jetPyRaw_sum = 0, jetPzRaw_sum = 0;
        double jetPxThermal_sum = 0, jetPyThermal_sum = 0, jetPzThermal_sum = 0;

        ET_miss_jet = 0;
        ET_miss_particle = 0;
        ET_miss_jet_raw = 0;
        ET_miss_jet_thermal = 0;

        A_J = -1;
        X_J = -1;

        // Clear jet-level output vectors
        jetPx->clear(); jetPy->clear(); jetPz->clear(); jetE->clear();
        jetPt->clear(); jetEta->clear(); jetPhi->clear();
        jetPxRaw->clear(); jetPyRaw->clear(); jetPzRaw->clear(); jetERaw->clear();
        jetPtRaw->clear(); jetEtaRaw->clear(); jetPhiRaw->clear();
        jetPxThermal->clear(); jetPyThermal->clear(); jetPzThermal->clear(); jetEThermal->clear();
        jetPtThermal->clear(); jetEtaThermal->clear(); jetPhiThermal->clear();

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
            if (partstatus->at(j) == 1){
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
            //Thermal Particles
            if (partstatus->at(j) == 3){
                thermalmomenta.push_back(fastjet::PseudoJet(pxJ->at(j), pyJ->at(j), pzJ->at(j), eJ->at(j)));
            }
        }
        if (i%10 == 0) {
            std::cout << "Event " << i << " details:" << std::endl;
            std::cout << "  Total number of particles for clustering: " << particles.size() << std::endl;
            std::cout << "  Number of thermal particles: " << thermalmomenta.size() << std::endl;
            std::cout << "  Event weight: " << partWeights << std::endl;
        }

        // Perform jet clustering for this event (subtract first, then sort by pT)
        fastjet::ClusterSequence cs(particles, jet_def);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

        // Clear output vectors and fill with jet 4-momenta (apply jetPtMin after subtraction)
        jetPx->clear(); jetPy->clear(); jetPz->clear(); jetE->clear();
        jetPt->clear(); jetEta->clear(); jetPhi->clear();
        std::vector<fastjet::PseudoJet> jets_corrected;

        for (const auto& J : jets) {
            fastjet::PseudoJet ThermalJets(0, 0, 0, 0);

            for (const auto& c : J.constituents()){
                for (const auto& t : thermalmomenta) {
                    if (deltaR2(c.eta(), c.phi_std(), t.eta(), t.phi_std()) < 1e-5) {
                        ThermalJets += fastjet::PseudoJet(t.px(), t.py(), t.pz(), t.e());
                    }
                }
            }

            jetPxRaw->push_back(J.px());
            jetPyRaw->push_back(J.py());
            jetPzRaw->push_back(J.pz());
            jetERaw->push_back(J.E());
            jetPtRaw->push_back(J.pt());
            jetEtaRaw->push_back(J.eta());
            jetPhiRaw->push_back(J.phi_std());

            jetPxThermal->push_back(ThermalJets.px());
            jetPyThermal->push_back(ThermalJets.py());
            jetPzThermal->push_back(ThermalJets.pz());
            jetEThermal->push_back(ThermalJets.E());
            jetPtThermal->push_back(ThermalJets.pt());
            jetEtaThermal->push_back(ThermalJets.eta());
            jetPhiThermal->push_back(ThermalJets.phi());
        
            fastjet::PseudoJet correctedJet = J - ThermalJets;
            jetPx->push_back(correctedJet.px());
            jetPy->push_back(correctedJet.py());
            jetPz->push_back(correctedJet.pz());
            jetE->push_back(correctedJet.E());
            jetPt->push_back(correctedJet.pt());
            jetEta->push_back(correctedJet.eta());
            jetPhi->push_back(correctedJet.phi_std());

            jetpx_sum += correctedJet.px();
            jetpy_sum += correctedJet.py();
            jetpz_sum += correctedJet.pz();

            jetPxRaw_sum += J.px();
            jetPyRaw_sum += J.py();
            jetPzRaw_sum += J.pz();

            jetPxThermal_sum += ThermalJets.px();
            jetPyThermal_sum += ThermalJets.py(); 
            jetPzThermal_sum += ThermalJets.pz();    

            jets_corrected.push_back(correctedJet);
        }

        jetpt1 = 0;
        jetpt2 = 0;
        for (size_t j = 0; j < jets_corrected.size(); ++j) {
            if (jets_corrected[j].pt() > jetpt1) {
                jetpt2 = jetpt1;
                jetphi2 = jetphi1;
                jetpt1 = jets_corrected[j].pt();
                jetphi1 = jets_corrected[j].phi_std();
            } else if (jets_corrected[j].pt() > jetpt2) {
                jetpt2 = jets_corrected[j].pt();
                jetphi2 = jets_corrected[j].phi_std();
            }
        }

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
        ET_miss_jet_raw = sqrt(jetPxRaw_sum*jetPxRaw_sum + jetPyRaw_sum*jetPyRaw_sum);
        ET_miss_jet_thermal = sqrt(jetPxThermal_sum*jetPxThermal_sum + jetPyThermal_sum*jetPyThermal_sum);
        jetTree->Fill();
        particleTree->Fill();

        if (i % 1000 == 0) {
            std::cout << "Processed " << i << " events." << std::endl;
        }
        if (i<5){
            std::cout << "Event " << i << " summary:" << std::endl;
            std::cout << "  Number of JEWEL particles: " << nJ << std::endl;
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
