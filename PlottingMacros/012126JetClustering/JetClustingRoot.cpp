// JetClusteringHepMC.cpp
// Usage: ./cluster_from_hepmc in.hepmc out.root
#include <fastjet/ClusterSequence.hh>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>

using namespace fastjet;
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
    TTree* jettree = new TTree("JetTree", "AK4 jets using antikt alg with R=0.4 pt>20GeV");

    double jetpt1 = 0;
    double jetpt2 = 0;
    double jetphi1 = 0;
    double jetphi2 = 0;
    double dPhi = 0;
    double A_J = -1;
    int evtNumber;
    int nParticles, nVertices;
    float ET_miss, ET_miss_jet;
    double jetpx_sum = 0;
    double jetpy_sum = 0;
    double px_sum = 0, py_sum = 0, pz_sum = 0;

    double pi = 3.14159265358979323846;

    std::vector<double>* jet_px = new std::vector<double>();
    std::vector<double>* jet_py = new std::vector<double>();
    std::vector<double>* jet_pz = new std::vector<double>();
    std::vector<double>* jet_e = new std::vector<double>();
    std::vector<double>* jet_pt = new std::vector<double>();
    std::vector<double>* jet_eta = new std::vector<double>();
    std::vector<double>* jet_phi = new std::vector<double>();
    std::vector<double>* jet_m = new std::vector<double>();
    std::vector<float>* px_v = new std::vector<float>();
    std::vector<float>* py_v = new std::vector<float>();
    std::vector<float>* pz_v = new std::vector<float>();
    std::vector<float>* e_v = new std::vector<float>();
    std::vector<int>* pdgId_v = new std::vector<int>();
    std::vector<float>* pt_v = new std::vector<float>();
    std::vector<float>* p_v = new std::vector<float>();
    std::vector<float>* eta_v = new std::vector<float>();
    std::vector<float>* phi_v = new std::vector<float>();
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
    particleTree->Branch("ET_miss", &ET_miss, "ET_miss/F");
    jettree->Branch("jet_px", &jet_px);
    jettree->Branch("jet_py", &jet_py);
    jettree->Branch("jet_pz", &jet_pz);
    jettree->Branch("jet_e", &jet_e);
    jettree->Branch("jet_pt", &jet_pt);
    jettree->Branch("jet_eta", &jet_eta);
    jettree->Branch("jet_phi", &jet_phi);
    jettree->Branch("jet_m", &jet_m);
    jettree->Branch("A_J", &A_J);
    jettree->Branch("jetpt1", &jetpt1);
    jettree->Branch("jetpt2", &jetpt2);
    jettree->Branch("jetphi1", &jetphi1);
    jettree->Branch("jetphi2", &jetphi2);
    jettree->Branch("ET_miss_jet", &ET_miss_jet);
    jettree->Branch("dPhi", &dPhi);

    // Jet definition
    const double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);

    while ((evt = input.read_next_event())) {
        jetpt1 = jetpt2 = -1.0;
        A_J = -1.0;
        ET_miss = 0.0f;
        ET_miss_jet = 0.0f;
        nParticles = 0;
        px_sum = py_sum = pz_sum = 0.0;
        px_v->clear(); py_v->clear(); pz_v->clear(); e_v->clear();
        pdgId_v->clear(); pt_v->clear(); p_v->clear(); eta_v->clear(); phi_v->clear();
        jet_px->clear(); jet_py->clear(); jet_pz->clear(); jet_e->clear();
        jet_pt->clear(); jet_eta->clear(); jet_phi->clear(); jet_m->clear();
        
        std::vector<PseudoJet> particles;
        particles.reserve(evt->particles_size());
        evtNumber = ievt;
        nVertices = evt->vertices_size();
        nParticles = 0;

        for (auto p = evt->particles_begin(); p != evt->particles_end(); ++p) {
            HepMC::GenParticle* part = *p;
            if (!part) continue;
            if (part->status() != 1) continue;
            nParticles++;
            const HepMC::FourVector& mom = part->momentum();
            particles.emplace_back(mom.px(), mom.py(), mom.pz(), mom.e());

            float pt = sqrt(mom.px()*mom.px() + mom.py()*mom.py());
            float part_p = sqrt(mom.px()*mom.px() + mom.py()*mom.py() + mom.pz()*mom.pz());
            float eta = -999;
            if ((fabs(part_p-fabs(mom.pz()))) > 1e-6){ // avoid division by zero
                eta = 0.5 * log((part_p + mom.pz()) / (part_p - mom.pz()));
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
            px_sum += mom.px();
            py_sum += mom.py();
            pz_sum += mom.pz();
        }

        ET_miss = sqrt(px_sum*px_sum + py_sum*py_sum);

        // If no particles, skip
        if (particles.empty()) {
            delete evt;
            continue;
        }
        // Cluster jets
        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(20.0));

        // Fill ROOT vectors for this event
        jetpx_sum = 0.0;
        jetpy_sum = 0.0;
        for (const auto& j : jets) {
            jet_px->push_back(j.px());
            jet_py->push_back(j.py());
            jet_pz->push_back(j.pz());
            jet_e->push_back(j.e());
            jet_pt->push_back(j.pt());
            jet_eta->push_back(j.eta());
            jet_phi->push_back(j.phi_std()); // phi in [-pi,pi]
            jet_m->push_back(j.m());
            jetpx_sum += j.px();
            jetpy_sum += j.py();
        }
        
        jetpt1 = (jet_pt->size() > 0) ? jet_pt->at(0) : -1.0;
        jetpt2 = (jet_pt->size() > 1) ? jet_pt->at(1) : -1.0;

        jetphi1 = (jet_phi->size() > 0) ? jet_phi->at(0) : -999.0;
        jetphi2 = (jet_phi->size() > 1) ? jet_phi->at(1) : -999.0;


        if (jetpt1 > 0 && jetpt2 > 0) {
            if (jetpt1 < jetpt2){
                cout << "Warning: jetpt1 < jetpt2" << endl;
                cout << "jetpt1: " << jetpt1 << ", jetpt2: " << jetpt2 << endl;
            }
            A_J = (jetpt1 - jetpt2) / (jetpt1 + jetpt2);
        } else {
            A_J = -1; // Undefined if less than 2 jets
        }

        dPhi = jetphi1 - jetphi2;

        if (pi <= dPhi && dPhi <= 2 * pi){
            dPhi = (dPhi - 2*pi);
        }
        else if (-2*pi <= dPhi && dPhi <= -pi){
            dPhi = dPhi + 2*pi;
        }

        ET_miss_jet = sqrt(jetpx_sum*jetpx_sum + jetpy_sum*jetpy_sum);

        if (evtNumber < 5) {
            cout << "Event " << evtNumber << ": nParticles=" << nParticles 
                 << ", nJets=" << jets.size() 
                 << ", jetpt1=" << jetpt1 
                 << ", jetpt2=" << jetpt2 
                 << ", A_J=" << A_J 
                 << ", ET_miss=" << ET_miss 
                 << ", ET_miss_jet=" << ET_miss_jet
                 << ", jetphi1=" << jetphi1
                 << ", jetphi2=" << jetphi2
                 << ", dPhi=" << dPhi
                 << endl;
        }

        // Fill tree and free event
        particleTree->Fill();
        jettree->Fill();
        ievt++;
        delete evt;
    }
    delete px_v; delete py_v; delete pz_v; delete e_v; delete pdgId_v;
    delete pt_v; delete p_v; delete eta_v; delete phi_v;
    delete jet_px; delete jet_py; delete jet_pz; delete jet_e;
    delete jet_pt; delete jet_eta; delete jet_phi; delete jet_m;
    // Write and close output
    f->cd();
    jettree->Write();
    particleTree->Write();
    f->Close();
    infile.close();

    std::cout << "Done. Processed " << ievt << " events.\n";
    return 0;
}