#!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration multistate reweighting analysis k2'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 12
#SBATCH --time 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python manipulateparameters.py AlkEthOH_c100 &
#python manipulateparameters.py AlkEthOH_c1266 &
#python manipulateparameters.py AlkEthOH_c581 &
#python manipulateparameters.py AlkEthOH_r0 &
#python manipulateparameters.py AlkEthOH_r48 &
#python manipulateparameters.py AlkEthOH_r51 &
wait
