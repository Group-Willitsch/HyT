#!/bin/bash
echo "Give number for simulation folder:"
read foldernumber
echo $foldernumber
mkdir $foldernumber
cp fort* $foldernumber
#cp ../output/ff.dat $foldernumber
cp ../output/T2jump.dat $foldernumber
cp ../output/T2jump.out $foldernumber
cp ../output/outcx.dat $foldernumber
#cp ../output/outcx-baby.dat $foldernumber
cp ../input/inputOH.dat $foldernumber
cp ../input/OH* $foldernumber
cp $foldernumber/fort.205 $foldernumber/fort-$foldernumber

