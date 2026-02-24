#!/bin/bash

SRC="JetClusteringWithNoBkg.cpp"
EXE="JetClusteringWithNoBkg"

INFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021926_JewelPtHat40_10Kevt/Jewelpbpb_5360GeV_ptm40_10Kevt_021826_C5_vac_11050627_precluster.root"
OUTFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/021926_JewelPtHat40_10Kevt/Jewelpbpb_5360GeV_ptm40_10Kevt_021826_C5_vac_11050627_nobkg.root"

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