#!/bin/bash
#cd ../pot2acc
#./pot2acc.exe -i1 ../input/decelrods-firsthoriz-normal-12kv.PA -i2 ../input/OH_M32_K32.dat
#cd ../timeseq
./timeseq.exe -i ../input/inputOH.dat
cd ../fly
./fly-fixed.exe -i1 ../input/inputOH.dat -i2 ../output/T2jump.out
cd ../bin
gfortran -o2 bin_tof_NV.for -o bin_tof_NV.exe
./bin_tof_NV.exe
xmgrace fort.205&
./saveall.sh
#xmgrace fort.205&
#xmgrace fort.505&
#xmgrace fort.405&
#xmgrace fort.305&

