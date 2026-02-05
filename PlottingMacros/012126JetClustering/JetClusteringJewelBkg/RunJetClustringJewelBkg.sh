#!/bin/bash

SRC="JetClusteringWithJewelBkg.cpp"
EXE="JetClusteringWithJewelBkg"
#INFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/Precluster/Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_dum_med_103428_precluster.root"
INFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/020126DebugDummies/Jewelpbpb_5360GeV_ptm250_10Kevt_020126_C5_dum_med_11050627_precluster.root"
OUTFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/Jewelpbpb_5360GeV_ptm250_10Kevt_020426_C5_dum_med_jewelbkg_new.root"

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