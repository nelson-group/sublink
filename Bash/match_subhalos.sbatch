#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J dm_match_L75n1820TNG_0.0
#SBATCH -o dm_match_L75n1820TNG_0.0.out
#SBATCH -e dm_match_L75n1820TNG_0.0.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 2-00:00:00
#SBATCH -p hernquist
#SBATCH --exclusive
#SBATCH --mem=240000

## OMP_NUM_THREADS controls the number of threads your application use
## This variable cat be set by the following command:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

TRACKING=Subhalos
SKIP_SUFFIX=
ALPHA_WEIGHT=0.0

#~ # Illustris
#~ FILESYSTEM=hernquistfs1
#~ SUITE=Illustris
#~ SIMULATION1=L75n1820FP
#~ SIMULATION2=L75n1820DM
#~ SNAPNUM_FIRST=0
#~ SNAPNUM_LAST=135
#~ SNAPNUM_START=0
#~ SNAPNUM_END=135

# IllustrisTNG
FILESYSTEM=hernquistfs3
SUITE=IllustrisTNG
SIMULATION1=L75n1820TNG
SIMULATION2=L75n1820TNG_DM
SNAPNUM_FIRST=0
SNAPNUM_LAST=99
SNAPNUM_START=0
SNAPNUM_END=99

BASEDIR1=/n/${FILESYSTEM}/${SUITE}/Runs/${SIMULATION1}/output
BASEDIR2=/n/${FILESYSTEM}/${SUITE}/Runs/${SIMULATION2}/output
TREEDIR=/n/hernquistfs3/vrodrigu/MergerTrees/output/${TRACKING}/${SUITE}/${SIMULATION1}/
if [ $SIMULATION = L75n1820FP ]; then
  SKIP_FILENAME=skip_snapshots/skip_snapshots${SKIP_SUFFIX}_${SIMULATION}.txt
else
  SKIP_FILENAME=skip_snapshots/skip_snapshots${SKIP_SUFFIX}.txt
fi
WRITEDIR=${TREEDIR}/dm_matches_${ALPHA_WEIGHT}
if [ ! -d $WRITEDIR ]; then
  mkdir -p $WRITEDIR
fi
WRITEPATH=${WRITEDIR}/dm_match

./SubLink/Descendants/match_subhalos $BASEDIR1 $BASEDIR2 $WRITEPATH $SNAPNUM_FIRST $SNAPNUM_LAST $SNAPNUM_START $SNAPNUM_END $TRACKING ${SKIP_FILENAME} ${ALPHA_WEIGHT}
./SubLink/Descendants/match_subhalos $BASEDIR2 $BASEDIR1 ${WRITEPATH}_reverse $SNAPNUM_FIRST $SNAPNUM_LAST $SNAPNUM_START $SNAPNUM_END $TRACKING ${SKIP_FILENAME} ${ALPHA_WEIGHT}
