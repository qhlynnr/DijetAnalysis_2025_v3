#!/bin/bash
# run_cluster_simple.sh
# Minimal compile-and-run wrapper for ROOT + FastJet

# ======== USER SETTINGS ========
INPUT="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/eventfiles/Jewel_pbpb_5360GeV_ptmin250_nevt26000_091125_med.hepmc"          # input HepMC file
OUTPUT="/home/xirong/DijetAnalysis_2025_v3_svmit/PlottingMacros/1014JewelWithJetClustering/Jewel_pbpb_5360GeV_ptmin250_nevt26000_091125_med_withjets.root"   # output ROOT file
SRC="JetClusteringHepMC.cpp"  # your C++ source
EXE="JetClusteringHepMC"     # compiled binary name
# ===============================

set -e  # stop if any command fails
HEPMC_DIR=/home/data/public/xirong/Packages/HepMC2/HepMC-2.06.10_install
HEPMC_INC="$HEPMC_DIR/include"
HEPMC_LIB="$HEPMC_DIR/lib"

export LD_LIBRARY_PATH="$HEPMC_LIB:$LD_LIBRARY_PATH"

echo "[build] Compiling $SRC ..."
g++ -O2 -std=c++17 "$SRC" -o "$EXE" -I"$HEPMC_INC" \
 -L"$HEPMC_LIB" -lHepMC \
  `root-config --cflags --libs` \
  `fastjet-config --cxxflags --libs`
echo "[build] Done."


echo "[run] Running $EXE on $INPUT → $OUTPUT ..."
./"$EXE" "$INPUT" "$OUTPUT"
echo "[run] Finished. Output written to $OUTPUT"