#!/bin/bash
for i in slurm_sim_c100_C-H_length*.sh
do
    sbatch "$i"
done

rename c100 c1266 slurm_sim_c100_C-H_length*.sh

sed -i 's/c100/c1266/g' slurm_sim_c1266_C-H_length*.sh

for i in slurm_sim_c1266_C-H_length*.sh
do
    sbatch "$i"
done

rename c1266 c581 slurm_sim_c1266_C-H_length*.sh

sed -i 's/c1266/c581/g' slurm_sim_c581_C-H_length*.sh

for i in slurm_sim_c581_C-H_length*.sh
do
    sbatch "$i"
done

rename c581 r0 slurm_sim_c581_C-H_length*.sh

sed -i 's/c581/r0/g' slurm_sim_r0_C-H_length*.sh

for i in slurm_sim_r0_C-H_length*.sh
do
    sbatch "$i"
done

rename r0 r48 slurm_sim_r0_C-H_length*.sh

sed -i 's/r0/r48/g' slurm_sim_r48_C-H_length*.sh

for i in slurm_sim_r48_C-H_length*.sh
do
    sbatch "$i"
done

rename r48 r51 slurm_sim_r48_C-H_length*.sh

sed -i 's/r48/r51/g' slurm_sim_r51_C-H_length*.sh

for i in slurm_sim_r51_C-H_length*.sh
do
    sbatch "$i"
done
