#!/bin/bash -l

#SBATCH
#SBATCH --job-name=unzip folder
#SBATCH --time=4:0:0
#SBATCH --error=./unzip folder.err
#SBATCH --output=./unzip folder.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


gunzip *