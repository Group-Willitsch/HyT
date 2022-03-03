#!/bin/bash

for ((i = 0; i<= 50; i=i+10))
do
./timeseq.exe -i ../input/inputOH_470/inputOH_$i.dat
cp ../output/T2jump.out ../output/T2jump.dat ./seq_yin_470/dec_470_$i
done
