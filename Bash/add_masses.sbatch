#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
##SBATCH --array=0-99%16
#SBATCH -J add_L35n2160TNG
#SBATCH -o add_L35n2160TNG.out
#SBATCH -e add_L35n2160TNG.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 1-00:00:00
#SBATCH -p hernquist
#SBATCH --mem=60000

# OMP_NUM_THREADS controls the number of threads your application uses:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# Recommended memory: 30 GB for L75n1820 runs; 120 GB for L205n2500 runs?

# IllustrisTNG
SIMULATION=L35n2160TNG
SIMDIR=/n/hernquistfs3/IllustrisTNG/Runs/${SIMULATION}
IN_RAD=0
SNAPNUM_FIRST=0
SNAPNUM_LAST=99

BASEDIR=${SIMDIR}/output
ASSEMBLYDIR=${SIMDIR}/postprocessing/StellarAssembly

#~ # Single snapshot (also uncomment --array instruction):
#~ python add_masses.py $BASEDIR $ASSEMBLYDIR ${IN_RAD} ${SLURM_ARRAY_TASK_ID}

# Loop over snapshots:
python add_masses.py $BASEDIR $ASSEMBLYDIR ${IN_RAD} ${SNAPNUM_FIRST} ${SNAPNUM_LAST}
