#!/bin/bash

#SBATCH --job-name=NAME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=wcpu

export COMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#HomeDir="/swdata/$USER/Simulations/Stark_deceleration_with_MATLAB/" # make sure we are in right directory in studix where your matlab file is stored
#cd $HomeDir

#INFILE=$1  # name of matlab script that shall be run taken from input
hostname  # in case the job crashes we know where it was

# merging all functions from FunctionList.txt with INFILE into a single script
# cat $INFILE > temp_trajectoryBASH.m # create temp_trajectoryBASH.m and write contenst of INFILE into it

FILE=test_CMAES.m # define variable FILEs

module load matlab/matlab-R2021b # load matlab shell environment

# HomeDir=$(pwd) # save current folder in a variable for later
SlurmID=$SLURM_JOBID # save slurm job IDs
#ScratchDir="/scratch/$USER/matlab.$SlurmID" # create a name for the directory
#mkdir -p ${ScratchDir} # make the directory
#cp $FILE $ScratchDir # copy matlab script to scratch folder
#cd $ScratchDir # change directory to scratch folder
stub=$(basename $FILE ".m") # get filename without ".m" extension for matlab
stubFILE=$(basename $FILE ".m") # get filename without ".m" extension for matlab
echo ${stubFILE}_${SlurmID}.log
matlab -r "run $stub;exit;" > ${stubFILE}_${SlurmID}.log # run matlab (srun for slurm only, otherwise without srun) and save terminal output in a log file

#cp -r ${ScratchDir}/* $HomeDir/ && rm -r ${ScratchDir} # copy all inputs and outputs back to original folder and remove scratch folder

#cd $HomeDir
#rm $FILE # remove temporary script
