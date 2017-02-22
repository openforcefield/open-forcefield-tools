#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration multistate reweighting analysis r51 k'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --time 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python manipulateparameters.py AlkEthOH_r51
