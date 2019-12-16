#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J off_L35n2160TNG
#SBATCH -o off_L35n2160TNG.out
#SBATCH -e off_L35n2160TNG.err
#SBATCH --ntasks-per-node=1
#SBATCH -t 2-00:00:00
#SBATCH -p hernquist
##SBATCH --exclusive
#SBATCH --mem=60000

## OMP_NUM_THREADS controls the number of threads your application uses
## This variable cat be set by the following command:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# Recommended memory: 240 GB for L205n2500 runs; 60 GB for L75n1820 runs

# For IllustrisTNG only
SUITE=IllustrisTNG
SIMULATION=L35n2160TNG
SIMDIR=/n/hernquistfs3/${SUITE}/Runs/${SIMULATION}
SNAPNUM_FIRST=0
SNAPNUM_LAST=99

BASEDIR=${SIMDIR}/output

# Tracking DM
TREEDIR=${SIMDIR}/postprocessing/trees/SubLink
python ./SubLink/Python/compute_offsets.py $BASEDIR $TREEDIR $SNAPNUM_FIRST $SNAPNUM_LAST
echo Finished for $SIMULATION

#~ # Tracking baryons
#~ TREEDIR=${SIMDIR}/postprocessing/trees/SubLink_gal
#~ python ./SubLink/Python/compute_offsets.py $BASEDIR $TREEDIR $SNAPNUM_FIRST $SNAPNUM_LAST
#~ echo Finished for $SIMULATION

