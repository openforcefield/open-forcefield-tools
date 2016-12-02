#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration multistate reweighting analysis r48 [#6X4:1]-[6X4:2] k'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --time 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python manipulateparameters.py AlkEthOH_r48
