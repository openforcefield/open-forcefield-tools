#!/bin/bash

rename [#6X4:2] [#1:2] slurm_script_mbar_analysis_*.sh

#sed -i 's/k/length/g' slurm_script_mbar_analysis_*.sh

for i in slurm_script_mbar_analysis_*.sh
do
    sbatch "$i"
done

