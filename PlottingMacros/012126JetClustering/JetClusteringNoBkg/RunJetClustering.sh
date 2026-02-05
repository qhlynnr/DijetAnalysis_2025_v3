#!/bin/bash

SRC="JetClusteringWithNoBkg.cpp"
EXE="JetClusteringWithNoBkg"

INFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_JewelAndAngantyr//Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_nodum_med_103428_precluster.root"
OUTFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_JewelAndAngantyr//Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_nodum_med_103428_nojetbkg.root"

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