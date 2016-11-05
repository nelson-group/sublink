#!/bin/bash

#SBATCH --mail-user=vrg@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -J stars_L75n1820FP
#SBATCH -o stars_L75n1820FP.out
#SBATCH -e stars_L75n1820FP.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 7-00:00:00
#SBATCH -p hernquist
#SBATCH --exclusive
#SBATCH --mem=240000

# OMP_NUM_THREADS controls the number of threads your application uses:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

## Note: For L75n455SPT simulation, use snap-groupordered_xxx.hdf5
## instead of snap_xxx.hdf5 (hard-coded in SubLink code).
TRACKING=Galaxies
SIMULATION=L75n1820FP

SNAPNUM_FIRST=0
SNAPNUM_LAST=135
SNAPNUM_RESTART=-1

BASEDIR=/n/hernquistfs1/Illustris/Runs/$SIMULATION/output
TREEDIR=/n/ghernquist/vrodrigu/MergerTrees/output/$TRACKING/Illustris/$SIMULATION
WRITEDIR=/n/hernquistfs3/vrodrigu/StellarAssembly/output/Illustris/$SIMULATION
# Make sure that write directory exists
if [ ! -d $WRITEDIR ]; then
  mkdir -p $WRITEDIR
fi
WRITEPATH=$WRITEDIR/stars

./SubLink/Examples/stellar_assembly $BASEDIR $TREEDIR $WRITEPATH ${SNAPNUM_FIRST} ${SNAPNUM_LAST} ${SNAPNUM_RESTART}
