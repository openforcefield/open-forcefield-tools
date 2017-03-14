!/bin/bash

#SBATCH --job-name 'Single Molecule Parameter Exploration general C-C'
#SBATCH --qos janus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 12
#SBATCH --time 20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

#some modules I use
ml slurm

#commands
python run_molecule.py AlkEthOH_c100 1.3734 1.526 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_c1266 1.3734 1.526 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_c581 1.3734 1.526 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_r0 1.3734 1.526 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_r48 1.3734 1.526 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_r51 1.3734 1.526 0.01526 [#6X4:1]-[#6X4:2] &
python run_molecule.py AlkEthOH_c100 1.54126 1.6786 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_c1266 1.54126 1.6786 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_c581 1.54126 1.6786 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_r0 1.54126 1.6786 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_r48 1.54126 1.6786 0.01526 [#6X4:1]-[#6X4:2] &
#python run_molecule.py AlkEthOH_r51 1.54126 1.6786 0.01526 [#6X4:1]-[#6X4:2] &
wait
