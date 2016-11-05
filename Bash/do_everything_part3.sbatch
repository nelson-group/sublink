#!/bin/bash

#SBATCH --mail-user=vrg@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -J gal_L75n910FP_part3
#SBATCH -o gal_L75n910FP_part3.out
#SBATCH -e gal_L75n910FP_part3.err
#SBATCH --ntasks-per-node=8
#SBATCH -t 7-00:00:00
#SBATCH -p hernquist
##SBATCH --exclusive
#SBATCH --mem=30000

## OMP_NUM_THREADS controls the number of threads your application uses
## This variable cat be set by the following command:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

TRACKING=Galaxies
SIMULATION=L75n910FP

# ----------------------- FIND DESCENDANTS

BASEDIR=/n/hernquistfs1/Illustris/Runs/${SIMULATION}/output
SKIP_SUFFIX=
if [ $SIMULATION = L75n1820FP ]; then
  SKIP_FILENAME=skip_snapshots/skip_snapshots${SKIP_SUFFIX}_${SIMULATION}.txt
else
  SKIP_FILENAME=skip_snapshots/skip_snapshots${SKIP_SUFFIX}.txt
fi
WRITEDIR=/n/ghernquist/vrodrigu/MergerTrees/output/$TRACKING/Illustris/${SIMULATION}_test/descendants
if [ ! -d $WRITEDIR ]; then
  mkdir -p $WRITEDIR
fi
WRITEPATH=$WRITEDIR/desc${SKIP_SUFFIX}

SNAPNUM_FIRST=0
SNAPNUM_LAST=135
SNAPNUM_START=0
SNAPNUM_END=135

#~ # First pass
#~ PASS=first
#~ ./SubLink/Descendants/find_descendants $BASEDIR ${WRITEPATH}_${PASS} $SNAPNUM_FIRST $SNAPNUM_LAST $SNAPNUM_START $SNAPNUM_END $TRACKING ${PASS} ${SKIP_FILENAME}
#~ 
#~ # Second pass
#~ PASS=second
#~ ./SubLink/Descendants/find_descendants $BASEDIR ${WRITEPATH}_${PASS} $SNAPNUM_FIRST $SNAPNUM_LAST $SNAPNUM_START $SNAPNUM_END $TRACKING ${PASS} ${SKIP_FILENAME}

# Compare descendants
./SubLink/Descendants/compare_descendants ${WRITEPATH} $SNAPNUM_FIRST $SNAPNUM_LAST ${SKIP_FILENAME}

# ---------------------------- BUILD TREES

TREEDIR=/n/ghernquist/vrodrigu/MergerTrees/output/$TRACKING/Illustris/${SIMULATION}_test
INPUT_PATH=$TREEDIR/descendants/desc${SKIP_SUFFIX}
OUTPUT_PATH=$TREEDIR/tree${SKIP_SUFFIX}

./SubLink/SubhaloTrees/build_trees ${INPUT_PATH} ${OUTPUT_PATH} ${SNAPNUM_FIRST} ${SNAPNUM_LAST} ${SKIP_FILENAME}

