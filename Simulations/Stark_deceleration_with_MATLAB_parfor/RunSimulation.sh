#!/bin/bash

#SBATCH --job-name=opti_atan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --index_funct=1
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=wcpu

export COMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#HomeDir=$(pwd) # save current folder in a variable for later


#HomeDir="/swdata/$USER/test_scratch" # make sure we are in right directory in studix where your matlab file is stored
#cd $HomeDir

#INFILE=$1  # name of matlab script that shall be run taken from input
hostname  # in case the job crashes we know where it was

#FileFolder='/swdata/brunner/Optmization_CMAES_scratch/Optimizations/' # Folder to copy
File=run_CMAES_atan.m # intialize the script which is used to run the optimization
#FileFull=$FileFolder$File # full path

# target folder



#SlurmID=$SLURM_JOBID # save slurm job IDs
#ScratchDir="/scratch/$USER/matlab_$SlurmID" # create a name for the directory
#mkdir -p ${ScratchDir} # make the directory

#cp -r $FileFolder/* $ScratchDir # copy every level in FileFolder to ScarthDir

#cd $ScratchDir # change directory to scratch folder



stub=$(basename $File ".m") # get filename without ".m" extension for matlab

# echo ${stub}_${SlurmID}.log

# load matlab shell environment
module load matlab/matlab-R2021b

matlab -r "run $stub;exit;" > ${stub}_${SlurmID}.log # run matlab (srun for slurm only, otherwise without srun) and save terminal output in a log file

# cp -r ${ScratchDir}/* $FileFolder/ # copy back home

# rm -r ${ScratchDir} # copy all inputs and outputs back to original folder and remove scratch folder

# cd $FileFolder
