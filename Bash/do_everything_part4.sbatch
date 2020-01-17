#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=FAIL
#SBATCH --array=0-7
#SBATCH -J gal_L75n910FP_part4
#SBATCH -o gal_L75n910FP_part4.out
#SBATCH -e gal_L75n910FP_part4.err
#SBATCH --ntasks-per-node=1
#SBATCH -t 7-00:00:00
#SBATCH -p hernquist
#SBATCH --mem=3750

## OMP_NUM_THREADS controls the number of threads your application uses
## This variable cat be set by the following command:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# Recommended memory:
# L75n1820FP, Galaxies: 30000 MB
# L75n910FP,  Subhalos: 7500 MB
# L75n1820FP, Subhalos: 60000 MB

# ----------- General-purpose variables --------------

NUM_JOBS=8
TRACKING=Galaxies
SIMULATION=L75n910FP

BASEDIR=/n/hernquistfs1/Illustris/Runs/$SIMULATION/output
TREEDIR=/n/ghernquist/vrodrigu/MergerTrees/output/$TRACKING/Illustris/$SIMULATION
SNAPNUM_FIRST=0
SNAPNUM_LAST=135

# Create extra columns
echo ""; echo "CREATING EXTRA COLUMNS..."
python ./SubLink/Python/create_columns.py $BASEDIR $TREEDIR $SNAPNUM_FIRST $SNAPNUM_LAST $NUM_JOBS ${SLURM_ARRAY_TASK_ID}
