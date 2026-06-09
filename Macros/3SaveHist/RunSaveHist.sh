#!/bin/bash

inputs=("/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/033126ValidationAJ/Jewel_pbpb_100K_2760GeV_ptm250_C10_med_12345_4momsub_R0p3.root" \
"/home/xirong/DijetAnalysis_2025_v3_svmit/MonteCarlo/033126ValidationAJ/Jewel_pbpb_100K_2760GeV_ptm250_C10_vac_12345_nobkg_0p3R.root")
outputs=("/0331CheckAJ2760TeV/H_Jewel_pbpb_100K_2760GeV_ptm250_C10_med_12345_4momsub_R0p3.root" \
"/0331CheckAJ2760TeV/H_Jewel_pbpb_100K_2760GeV_ptm250_C10_vac_12345_nobkg_R0p3.root")

for i in ${!inputs[@]}
do
  root -l -q "SaveHist.cpp(\"${inputs[$i]}\",\"${outputs[$i]}\")"
done