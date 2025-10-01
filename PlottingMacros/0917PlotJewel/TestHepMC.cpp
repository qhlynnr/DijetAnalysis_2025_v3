#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include <iostream>
#include <fstream>


int TestHepMC() {

    const char* fn = "/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/eventfiles/Jewel_pbpb_5360GeV_ptmin250_nevt10_091125_med.hepmc";
    std::ifstream fin(fn);
    if (!fin.is_open()) {
        std::cerr << "Error: Could not open file " << fn << std::endl;
        return 1;
    }
    fin.clear(); // clear any eof/fail bits 

    HepMC::IO_GenEvent input(fin);
    HepMC::GenEvent* evt = input.read_next_event();
    if (!evt) {
        std::cerr << "Error: Could not read first event. File may be empty or corrupted." << std::endl;
        return 1;
    }
    size_t ievt = 0;
    while (evt) {
        std::cout << "Event " << ievt << " has " << evt->particles_size() << " particles and "
                  << evt->vertices_size() << " vertices\n";
        delete evt; // free memory
        evt = input.read_next_event();

        for (HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p) {
            HepMC::GenParticle* part = *p;
            if (!part) continue;                          
            if (part->status() != 1) continue;  //final state

            std::cout << "  Particle ID: " << part->pdg_id() 
                      << " Momentum: (" << part->momentum().px() << ", "
                      << part->momentum().py() << ", "
                      << part->momentum().pz() << ", "
                      << part->momentum().e() << ")\n";
        }

        ievt++;
        if (ievt >= 1) break; // limit to first 10 events for this test
    }
    return 0;
}