#!/bin/bash
# run_cluster_simple.sh
# Minimal compile-and-run wrapper for ROOT + FastJet

# ======== USER SETTINGS ========
#/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/Jewel_pbpb_HepMC_5360GeV/Jewel_pbpb_5360GeV_ptmin250_nevt10K_103125_vac.hepmc
INPUT_DIR="/home/data/public/xirong/Packages/jewel/jewel-2.4.0/eventfiles/0127Events100K/"
INPUT="$INPUT_DIR/Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_nodum_med_103428.hepmc" # input HepMC file
OUTPUT_DIR="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_JewelAndAngantyr/"
OUTPUT="$OUTPUT_DIR/Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_nodum_med_103428_precluster.root"   # output ROOT file
SRC="HepMCtoRoot.cpp"  # your C++ source
EXE="HepMCtoRoot"     # compiled binary name
# ===============================

set -e  # stop if any command fails
HEPMC_DIR=/home/data/public/xirong/Packages/HepMC2/HepMC-2.06.10_install
HEPMC_INC="$HEPMC_DIR/include"
HEPMC_LIB="$HEPMC_DIR/lib"

export LD_LIBRARY_PATH="$HEPMC_LIB:$LD_LIBRARY_PATH"

echo "[build] Compiling $SRC ..."
g++ -O2 -std=c++17 "$SRC" -o "$EXE" -I"$HEPMC_INC" \
 -L"$HEPMC_LIB" -lHepMC \
  `root-config --cflags --libs`
echo "[build] Done."


echo "[run] Running $EXE on $INPUT → $OUTPUT ..."
./"$EXE" "$INPUT" "$OUTPUT"
echo "[run] Finished. Output written to $OUTPUT"