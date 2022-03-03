#!/bin/bash

for ((i = 80; i <= 110; i = i+5))
do
./fly-fixed.exe -i1 ../input/test_LB_new/inputOH_LB$i.dat -i2 ../timeseq/seq_yin/dec_450_60deg_92/T2jump.out
cp ../output/ff.dat test_LB_data_new/dec_450_60deg_LB$i.dat
done
