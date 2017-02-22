#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration AlkEthOH_r48'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --time 23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python run_molecule.py AlkEthOH_r48 1.0700 1.0800 0.01 [#6X4:1]-[#1:2]
