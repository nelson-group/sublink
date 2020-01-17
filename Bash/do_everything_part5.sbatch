#!/bin/bash

#SBATCH --mail-user=v.rodriguez@irya.unam.mx
#SBATCH --mail-type=ALL
#SBATCH -J gal_L75n910FP_part5
#SBATCH -o gal_L75n910FP_part5.out
#SBATCH -e gal_L75n910FP_part5.err
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

TRACKING=Galaxies
SIMULATION=L75n910FP

BASEDIR=/n/hernquistfs1/Illustris/Runs/$SIMULATION/output
TREEDIR=/n/ghernquist/vrodrigu/MergerTrees/output/$TRACKING/Illustris/$SIMULATION
SNAPNUM_FIRST=0
SNAPNUM_LAST=135

# Create extended trees
echo ""; echo "CREATING EXTENDED TREES..."
python ./SubLink/Python/create_extended_trees.py $BASEDIR $TREEDIR $SNAPNUM_LAST

# Concatenate minimal trees
echo ""; echo "CONCATENATING MINIMAL TREES..."
python ./SubLink/Python/concatenate_trees.py ${TREEDIR}/tree

# Concatenate extended trees
echo ""; echo "CONCATENATING EXTENDED TREES..."
python ./SubLink/Python/concatenate_trees.py ${TREEDIR}/tree_extended

# Compute offsets
echo ""; echo "COMPUTING OFFSETS..."
python ./SubLink/Python/compute_offsets.py $BASEDIR $TREEDIR $SNAPNUM_FIRST $SNAPNUM_LAST

