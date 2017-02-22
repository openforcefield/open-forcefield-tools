#!/bin/bash
for i in slurm_sim_c100_C-H*.sh
do
    sbatch "$i"
done

rename c100 c1266 slurm_sim_c100_C-H*.sh

sed -i 's/c100/c1266/g' slurm_sim_c1266_C-H*.sh

for i in slurm_sim_c1266_C-H*.sh
do
    sbatch "$i"
done

rename c1266 c581 slurm_sim_c1266_C-H*.sh

sed -i 's/c1266/c581/g' slurm_sim_c581_C-H*.sh

for i in slurm_sim_c581_C-H*.sh
do
    sbatch "$i"
done

rename c581 r0 slurm_sim_c581_C-H*.sh

sed -i 's/c581/r0/g' slurm_sim_r0_C-H*.sh

for i in slurm_sim_r0_C-H*.sh
do
    sbatch "$i"
done

rename r0 r48 slurm_sim_r0_C-H*.sh

sed -i 's/r0/r48/g' slurm_sim_r48_C-H*.sh

for i in slurm_sim_r48_C-H*.sh
do
    sbatch "$i"
done

rename r48 r51 slurm_sim_r48_C-H*.sh

sed -i 's/r48/r51/g' slurm_sim_r51_C-H*.sh

for i in slurm_sim_r51_C-H*.sh
do
    sbatch "$i"
done

rename r51 c100 slurm_sim_r51_C-H*.sh

sed -i 's/r51/c100/g' slurm_sim_c100_C-H*.sh


# Start the C-C simulations

for i in slurm_sim_c100_C-C*.sh
do
    sbatch "$i"
done

rename c100 c1266 slurm_sim_c100_C-C*.sh

sed -i 's/c100/c1266/g' slurm_sim_c1266_C-C*.sh

for i in slurm_sim_c1266_C-C*.sh
do
    sbatch "$i"
done

rename c1266 c581 slurm_sim_c1266_C-C*.sh

sed -i 's/c1266/c581/g' slurm_sim_c581_C-C*.sh

for i in slurm_sim_c581_C-C*.sh
do
    sbatch "$i"
done

rename c581 r0 slurm_sim_c581_C-C*.sh

sed -i 's/c581/r0/g' slurm_sim_r0_C-C*.sh

for i in slurm_sim_r0_C-C*.sh
do
    sbatch "$i"
done

rename r0 r48 slurm_sim_r0_C-C*.sh

sed -i 's/r0/r48/g' slurm_sim_r48_C-C*.sh

for i in slurm_sim_r48_C-C*.sh
do
    sbatch "$i"
done

rename r48 r51 slurm_sim_r48_C-C*.sh

sed -i 's/r48/r51/g' slurm_sim_r51_C-C*.sh

for i in slurm_sim_r51_C-C*.sh
do
    sbatch "$i"
done

rename r51 c100 slurm_sim_r51_C-C*.sh

sed -i 's/r51/c100/g' slurm_sim_c100_C-C*.sh

