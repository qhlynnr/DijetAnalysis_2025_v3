#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <fastjet/ClusterSequence.hh>
#include <vector>

using namespace fastjet;

int HepMCtoRootConverter() {

    //Opening HEPMC file
    const char* infilename = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/eventfiles/Jewel_pbpb_5360GeV_ptmin250_nevt10K_091125_med.hepmc";
    std::ifstream infile(infilename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << infilename << std::endl;
        return 1;
    }
    infile.clear(); // clear any eof/fail bits 

    HepMC::IO_GenEvent input(infile);
    HepMC::GenEvent* evt = nullptr;
    size_t ievt = 0;

    evt = input.read_next_event();

    //Opening root file

    TFile *f = new TFile("Jewel_pbpb_5360GeV_ptmin250_nevt10K_101325_med.root", "RECREATE");
    TTree *particleTree = new TTree("ParticleTree", "An event tree <3");

    int evtNumber;
    int nParticles;
    int nVertices;
    std::vector<float>* px = new std::vector<float>();
    std::vector<float>* py = new std::vector<float>();
    std::vector<float>* pz = new std::vector<float>();
    std::vector<float>* e = new std::vector<float>();
    std::vector<int>* pdgId = new std::vector<int>();

    std::vector<float>* pt_v = new std::vector<float>();
    std::vector<float>* p_v = new std::vector<float>();
    std::vector<float>* eta_v = new std::vector<float>();
    std::vector<float>* phi_v = new std::vector<float>();

    float pt1;
    float pt2;
    float A_J;
    float ET_miss;
    float phi;

    particleTree->Branch("nEvts", &evtNumber, "nEvts/I");
    particleTree->Branch("nParticles", &nParticles, "nParticles/I");
    particleTree->Branch("nVertices", &nVertices, "nVertices/I");
    particleTree->Branch("px", &px);
    particleTree->Branch("py", &py);
    particleTree->Branch("pz", &pz);
    particleTree->Branch("e", &e);
    particleTree->Branch("pdgId", &pdgId);
    particleTree->Branch("pt", &pt_v);
    particleTree->Branch("p", &p_v);
    particleTree->Branch("eta", &eta_v);
    particleTree->Branch("phi", &phi_v);
    particleTree->Branch("pt1", &pt1, "pt1/F");
    particleTree->Branch("pt2", &pt2, "pt2/F");
    particleTree->Branch("A_J", &A_J, "A_J/F");
    particleTree->Branch("ET_miss", &ET_miss, "ET_miss/F");

    //Event Loop
    while (evt) {
        std::cout << "Event " << ievt << " has " << evt->particles_size() << " particles and "
                  << evt->vertices_size() << " vertices\n";
        delete evt; // free memory

        evt = input.read_next_event();

    if (!evt){
        std::cout << "No more events found.\n";
        break;
    }

        evtNumber = ievt;
        nParticles = evt->particles_size();
        pt1 = 0;
        pt2 = 0;
        float px_sum = 0;
        float py_sum = 0;
        float pz_sum = 0;

        if (!evt){
            std::cout << "No more events found.\n";
            continue;
        }
    
        //Particle Loop
        for (HepMC::GenEvent::particle_iterator part = evt->particles_begin(); part != evt->particles_end(); ++part) {
            HepMC::GenParticle* particle = *part;
            const HepMC::FourVector& m = particle->momentum();

            if (!particle) continue;                          
            if (particle->status() != 1) continue;  //final state

          /* std::cout << " Particle ID: " << particle->pdg_id() 
                      << " Momentum: (" 
                      << m.px() << ", "
                      << m.py() << ", "
                      << m.pz() << ", "
                      << m.e() << ")\n";*/

            float pt = sqrt(m.px()*m.px() + m.py()*m.py());
            float p = sqrt(m.px()*m.px() + m.py()*m.py() + m.pz()*m.pz());
            float eta = 0.5 * log((p + m.pz()) / (p - m.pz()));
            float phi = atan2(m.py(), m.px());
            
            px->push_back(m.px());
            py->push_back(m.py());
            pz->push_back(m.pz());
            e->push_back(m.e());
            pdgId->push_back(particle->pdg_id());
            pt_v->push_back(pt);
            p_v->push_back(p);
            eta_v->push_back(eta);
            phi_v->push_back(phi);
            if (pt > pt1 && pt > pt2) {
                pt2 = pt1;
                pt1 = pt;
            } else if (pt > pt2 && pt < pt1) {
                pt2 = pt;
            }
            px_sum += m.px();
            py_sum += m.py();
            pz_sum += m.pz();
        }
        A_J = (pt1 - pt2) / (pt1 + pt2);
        ET_miss = sqrt(px_sum*px_sum + py_sum*py_sum);
        t->Fill();
        ievt++;
       // if (ievt >= 10) break; // limit to first 10 events for this test
    }
    t->Write();
    delete f;

    return 0;
}