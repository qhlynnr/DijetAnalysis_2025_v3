#!/bin/bash

make

rm *.o
rm src/*.o
rm ../../UtilityCode/Source/*.o

INPUT="/afs/cern.ch/user/x/xirong/CMSSW_14_1_9/src/HijingRootFiles_PbPb2024/HiForestMiniAOD.root"
./QuickDrawForest \
  --Input $INPUT \
  --Debug true
