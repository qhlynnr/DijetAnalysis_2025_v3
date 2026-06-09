#!/bin/bash

SRC="JetClusteringWithJewelBkg.cpp"
EXE="JetClusteringWithJewelBkg"
#INFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/Precluster/Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_dum_med_103428_precluster.root"
INFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/033126ValidationAJ/Jewel_pbpb_100K_2760GeV_ptm250_C10_med_12345.root"
OUTFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/033126ValidationAJ/Jewel_pbpb_100K_2760GeV_ptm250_C10_med_12345_4momsub_R0p3.root"

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