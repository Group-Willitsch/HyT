#!/bin/bash


for ((n=400;n<=480;n=n+10))
do
./fly-fixed.exe -i1 ../input/inputOH_velocity/inputOH_$n.dat -i2 ../timeseq/seq_yin_to35/from$n/T2jump.out
cp ../output/ff.dat compare_velocity_data_to35/from$n.dat
done






