#!/bin/bash


./fly-fixed.exe -i1 ../input/test_L4/inputOH.dat -i2 ../timeseq/seq_yin/dec_450_30deg_328/T2jump.out
cp ../output/ff.dat test_L4/dec_450_30deg_382.dat

./fly-fixed.exe -i1 ../input/test_L4/inputOH.dat -i2 ../timeseq/seq_yin/dec_450_40deg_273/T2jump.out
cp ../output/ff.dat test_L4/dec_450_40deg_273.dat

./fly-fixed.exe -i1 ../input/test_L4/inputOH.dat -i2 ../timeseq/seq_yin/dec_450_50deg_202/T2jump.out
cp ../output/ff.dat test_L4/dec_450_50deg_202.dat

./fly-fixed.exe -i1 ../input/test_L4/inputOH.dat -i2 ../timeseq/seq_yin/dec_450_60deg_92/T2jump.out
cp ../output/ff.dat test_L4/dec_450_60deg_92.dat

./fly-fixed.exe -i1 ../input/test_L4/inputOH.dat -i2 ../timeseq/seq_yin/dec_450_62deg_49/T2jump.out
cp ../output/ff.dat test_L4/dec_450_62deg_49.dat

./fly-fixed.exe -i1 ../input/test_L4/inputOH.dat -i2 ../timeseq/seq_yin/dec_450_62p5_30/T2jump.out
cp ../output/ff.dat test_L4/dec_450_62p5_30.dat













