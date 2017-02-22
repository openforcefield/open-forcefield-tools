#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration AlkEthOH_c100'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --time 19:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python run_molecule.py AlkEthOH_c100 200 300 50 [#6X4:1]-[#1-2]
