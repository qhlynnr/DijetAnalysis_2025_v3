#!/bin/bash

SRC="EmbedClusterWithFJBkg.cpp"
EXE="EmbedClusterWithFJBkg"

ANGANTYRFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/pythia8Jets_pbpb_5200GeV_50KEvents.root"
JEWELFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_Jewel50KEventsWithandWithoutWriting/Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_dum_med_103428_precluster.root"
OUTFILE="/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/012726_JewelAndAngantyr/Embedded_Jewelpbpb_5360GeV_ptm250_50Kevt_012726_C5_dum_med_withFJBkg.root"
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
./$EXE "$OUTFILE" "$ANGANTYRFILE" "$JEWELFILE"