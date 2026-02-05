// JetClusteringHepMC.cpp
// Usage: ./cluster_from_hepmc in.hepmc out.root
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: cluster_from_hepmc in.hepmc out.root\n";
        return 1;
    }

    const char* infileHepMC = argv[1];
    const char* outfileRoot = argv[2];

    std::ifstream infile(infileHepMC);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << infileHepMC << std::endl;
        return 1;
    }

    HepMC::IO_GenEvent input(infile);
    HepMC::GenEvent* evt = nullptr;
    size_t ievt = 0;

    // Prepare ROOT output
    TFile* f = TFile::Open(outfileRoot, "RECREATE");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: could not open output file " << outfileRoot << std::endl;
        return 1;
    }

    TTree* particleTree = new TTree("ParticleTree", "An event tree straight from hep");

    int evtNumber;
    int nParticles, nVertices;

    std::vector<float>* px_v = new std::vector<float>();
    std::vector<float>* py_v = new std::vector<float>();
    std::vector<float>* pz_v = new std::vector<float>();
    std::vector<float>* e_v = new std::vector<float>();
    std::vector<int>* pdgId_v = new std::vector<int>();
    std::vector<float>* pt_v = new std::vector<float>();
    std::vector<float>* p_v = new std::vector<float>();
    std::vector<float>* eta_v = new std::vector<float>();
    std::vector<float>* phi_v = new std::vector<float>();
    std::vector<int>* status_v = new std::vector<int>();
    particleTree->Branch("nEvts", &evtNumber, "nEvts/I");
    particleTree->Branch("nParticles", &nParticles, "nParticles/I");
    particleTree->Branch("nVertices", &nVertices, "nVertices/I");
    particleTree->Branch("px", &px_v);
    particleTree->Branch("py", &py_v);
    particleTree->Branch("pz", &pz_v);
    particleTree->Branch("e", &e_v);
    particleTree->Branch("pdgId", &pdgId_v);
    particleTree->Branch("pt", &pt_v);
    particleTree->Branch("p", &p_v);
    particleTree->Branch("eta", &eta_v);
    particleTree->Branch("phi", &phi_v);
    particleTree->Branch("status", &status_v);

    while ((evt = input.read_next_event())) {
        px_v->clear(); py_v->clear(); pz_v->clear(); e_v->clear();
        pdgId_v->clear(); pt_v->clear(); p_v->clear(); eta_v->clear(); phi_v->clear(); status_v->clear();

        evtNumber = ievt;
        nVertices = evt->vertices_size();
        nParticles = 0;

        for (auto p = evt->particles_begin(); p != evt->particles_end(); ++p) {
            HepMC::GenParticle* part = *p;
            if (!part) continue;
            //if (part->status() != 1) continue;
            nParticles++;
            const HepMC::FourVector& mom = part->momentum();
            float pt = sqrt(mom.px()*mom.px() + mom.py()*mom.py());
            float part_p = sqrt(mom.px()*mom.px() + mom.py()*mom.py() + mom.pz()*mom.pz());
            float eta = -999;
            if ((fabs(part_p-fabs(mom.pz()))) > 1e-6){ // avoid division by zero
                eta = 0.5 * log((part_p + mom.pz()) / (part_p - mom.pz()));
            }
            if (ievt < 20 && nParticles <=5) {
                cout << "Status:" << part->status() << ", PdgID:" << part->pdg_id() << endl;
                cout << "E:" << mom.e() << endl;
            }
            float phi = atan2(mom.py(), mom.px());
            px_v->push_back(mom.px());
            py_v->push_back(mom.py());
            pz_v->push_back(mom.pz());
            e_v->push_back(mom.e());
            pdgId_v->push_back(part->pdg_id());
            pt_v->push_back(pt);
            p_v->push_back(part_p);
            eta_v->push_back(eta);
            phi_v->push_back(phi);
            status_v->push_back(part->status());
        }

        if (ievt < 10) {
            std::cout << "Event " << evtNumber << ": nParticles=" << nParticles 
                      << ", nVertices=" << nVertices << std::endl;
            std::cout << "PdgID\tpt\teta\tphi\n";
            for (size_t i = 0; i < 5 && i < pdgId_v->size(); ++i) {
                std::cout << pdgId_v->at(i) << "\t"
                          << pt_v->at(i) << "\t"
                          << eta_v->at(i) << "\t"
                          << phi_v->at(i) << "\n";
            }
        }
        if (ievt % 100 == 0) {
                std::cout << "Processed " << ievt << " events." << std::endl;
        }

        // Fill tree and free event
        particleTree->Fill();
        ievt++;
        delete evt;
    }
    delete px_v; delete py_v; delete pz_v; delete e_v; delete pdgId_v;
    delete pt_v; delete p_v; delete eta_v; delete phi_v;
    // Write and close output
    f->cd();
    particleTree->Write();
    f->Close();
    infile.close();

    std::cout << "Done. Processed " << ievt << " events.\n";
    return 0;
}