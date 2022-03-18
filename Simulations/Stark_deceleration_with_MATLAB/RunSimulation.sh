#!/bin/bash

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks=6 # make sure this matches the number of cores in your matlab input script (cf. LASTN=maxNumCompThreads(N))
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=short

HomeDir="/swdata/$USER/TrajectorySimulations/" # make sure we are in right directory
cd $HomeDir

INFILE=$1  # name of matlab script that shall be run taken from input
hostname  # in case the job crashes we know where it was

# merging all functions from FunctionList.txt with INFILE into a single script
cat $INFILE > temp_stark.m # create temp_trajectoryBASH.m and write content of INFILE into it
FILE=temp_trajectoryBASH.m # define variable FILE
while IFS= read -r line; do
         cat $line >> $FILE # append the content of all functions in FunctionList.txt
	 echo $'\r' >> $FILE # add a line break
done < "FunctionList.txt"

module load matlab/matlab-R2018a # load matlab shell environment

# HomeDir=$(pwd) # save current folder in a variable for later
SlurmID=$SLURM_JOBID # save slurm job ID 
ScratchDir="/scratch/$USER/matlab.$SLURM_JOBID" # create a name for the directory
mkdir -p ${ScratchDir} # make the directory
cp $FILE $ScratchDir # copy matlab script to scratch folder
cd $ScratchDir # change directory to scratch folder
stub=$(basename $FILE ".m") # get filename without ".m" extension for matlab
stubINFILE=$(basename $INFILE ".m") # get filename without ".m" extension for matlab
srun matlab -r "run $stub;exit;" > ${stubINFILE}.log # run matlab (srun for slurm only, otherwise without srun) and save terminal output in a log file

cp -r ${ScratchDir}/* $HomeDir/ && rm -r ${ScratchDir} # copy all inputs and outputs back to original folder and remove scratch folder

cd $HomeDir
rm $FILE # remove temporary script


