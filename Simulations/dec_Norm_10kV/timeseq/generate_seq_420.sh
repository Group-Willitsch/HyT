#!/bin/bash

for ((i = 1; i<= 9; i=i+1))
do
./timeseq.exe -i ../input/inputOH_420/inputOH_54p$i.dat
cp ../output/T2jump.out ../output/T2jump.dat ./seq_yin_420/dec_420_54p$i
done
