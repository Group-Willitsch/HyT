#!/bin/bash

./fly-fm.exe -i1 ../input/inputOH.dat -i2 ../timeseq/sequencesFM_vel/dec_450_328/T2jump.out
cp ../output/ff.dat fly_fm/dec_450_328.dat

./fly-fm.exe -i1 ../input/inputOH.dat -i2 ../timeseq/sequencesFM_vel/dec_450_273/T2jump.out
cp ../output/ff.dat fly_fm/dec_450_273.dat

./fly-fm.exe -i1 ../input/inputOH.dat -i2 ../timeseq/sequencesFM_vel/dec_450_202/T2jump.out
cp ../output/ff.dat fly_fm/dec_450_202.dat

./fly-fm.exe -i1 ../input/inputOH.dat -i2 ../timeseq/sequencesFM_vel/dec_450_92/T2jump.out
cp ../output/ff.dat fly_fm/dec_450_92.dat

./fly-fm.exe -i1 ../input/inputOH.dat -i2 ../timeseq/sequencesFM_vel/dec_450_49/T2jump.out
cp ../output/ff.dat fly_fm/dec_450_49.dat

./fly-fm.exe -i1 ../input/inputOH.dat -i2 ../timeseq/sequencesFM_vel/dec_450_30/T2jump.out
cp ../output/ff.dat fly_fm/dec_450_30.dat
