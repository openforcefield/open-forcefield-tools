#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration general'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 12
#SBATCH --time 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python run_molecule.py AlkEthOH_c100 690 695 5 [#6X4:1]-[#6X4:2] &
python run_molecule.py AlkEthOH_c1266 690 695 5 [#6X4:1]-[#6X4:2] &
python run_molecule.py AlkEthOH_c581 690 695 5 [#6X4:1]-[#6X4:2] &
python run_molecule.py AlkEthOH_r0 690 695 5 [#6X4:1]-[#6X4:2] &
python run_molecule.py AlkEthOH_r48 685 695 5 [#6X4:1]-[#6X4:2] &
wait
