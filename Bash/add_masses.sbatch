#!/bin/bash

#SBATCH --mail-user=vrodriguez-gomez@cfa.harvard.edu
#SBATCH --mail-type=ALL
##SBATCH --array=0-135%16
#SBATCH -J add_L75n1820FP
#SBATCH -o add_L75n1820FP.out
#SBATCH -e add_L75n1820FP.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 7-00:00:00
#SBATCH -p hernquist
#SBATCH --mem=60000

# OMP_NUM_THREADS controls the number of threads your application uses:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

SIMULATION=L75n1820FP
IN_RAD=0
SNAPNUM_FIRST=0
SNAPNUM_LAST=135

BASEDIR=/n/hernquistfs1/Illustris/Runs/$SIMULATION/output
ASSEMBLYDIR=/n/hernquistfs3/vrodrigu/StellarAssembly/output/Illustris/$SIMULATION

#~ # Single snapshot (also uncomment --array instruction):
#~ python add_masses.py $BASEDIR $ASSEMBLYDIR ${IN_RAD} ${SLURM_ARRAY_TASK_ID}

# Loop over snapshots:
python ./SubLink/Python/add_masses.py $BASEDIR $ASSEMBLYDIR ${IN_RAD} ${SNAPNUM_FIRST} ${SNAPNUM_LAST}
