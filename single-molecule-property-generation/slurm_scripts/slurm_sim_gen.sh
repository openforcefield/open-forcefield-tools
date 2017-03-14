#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration general'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 12
#SBATCH --time 20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python run_molecule.py AlkEthOH_c100 630 680 5 [#6X4:1]-[#1:2] &
python run_molecule.py AlkEthOH_c100 685 735 5 [#6X4:1]-[#1:2] &
python run_molecule.py AlkEthOH_c100 580 630 5 [#6X4:1]-[#6X4:2] &
python run_molecule.py AlkEthOH_c100 635 685 5 [#6X4:1]-[#6X4:2] &
wait
