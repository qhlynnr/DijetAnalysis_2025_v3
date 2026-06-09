#!/bin/bash

SRC="JetClusteringWithNoBkg.cpp"
EXE="JetClusteringWithNoBkg"

INFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/033126ValidationAJ/Jewel_pbpb_100K_2760GeV_ptm250_C10_vac_12345.root"
OUTFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/033126ValidationAJ/Jewel_pbpb_100K_2760GeV_ptm250_C10_vac_12345_nobkg_0p3R.root"

echo "Compiling $SRC ..."
g++ $SRC -O2 -std=c++17 -o $EXE \
  $(root-config --cflags --libs) \
  $(fastjet-config --cxxflags --libs) \
  -lfastjettools

if [ $? -ne 0 ]; then
  echo "Compilation failed."
  exit 1
fi

echo "Running $EXE ..."
./$EXE "$INFILE" "$OUTFILE"