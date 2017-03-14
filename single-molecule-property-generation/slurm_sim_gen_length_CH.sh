#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration general C-H'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 12
#SBATCH --time 20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python run_molecule.py AlkEthOH_c100 0.99 1.10 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_c1266 0.99 1.10 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_c581 0.99 1.10 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_r0 0.99 1.10 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_r48 0.99 1.10 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_r51 0.99 1.10 0.01 [#6X4:1]-[#1:2] &
python run_molecule.py AlkEthOH_c100 1.10 1.20 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_c1266 1.10 1.20 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_c581 1.10 1.20 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_r0 1.10 1.20 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_r48 1.10 1.20 0.01 [#6X4:1]-[#1:2] &
#python run_molecule.py AlkEthOH_r51 1.10 1.20 0.01 [#6X4:1]-[#1:2] &
wait
