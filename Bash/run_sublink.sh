# run_sublink.sh: do full tree construction for a simulation ()

# config
export OMP_NUM_THREADS=16

TRACKING=Galaxies #Galaxies,Subhalos,Stars

SNAP_START=0
SNAP_END=100
BASEDIR=/n/sim/path/output/

# prepare
OUTPATH=$BASEDIR/sublink/
mkdir -p $OUTPATH
touch dummy

# (1) find descendants
./Descendants/find_descendants $BASEDIR ${OUTPATH}_first $SNAP_START $SNAP_END $SNAP_START $SNAP_END $TRACKING first dummy
./Descendants/find_descendants $BASEDIR ${OUTPATH}_second $SNAP_START $SNAP_END $SNAP_START $SNAP_END $TRACKING second dummy

# (2) compare descendants
./Descendants/compare_descendants $OUTPATH $SNAP_START $SNAP_END dummy

# (3) build trees
./SubhaloTrees/build_trees $OUTPATH ${OUTPATH}tree $SNAP_START $SNAP_END dummy

# (4) concatenate, create 'basic' and 'extended' trees, calculate offsets
python Python/create_columns.py $BASEDIR $OUTPATH $SNAP_START $SNAP_END 1 0
python Python/create_extended_trees.py $BASEDIR $OUTPATH $SNAP_END
python Python/concatenate_trees.py ${OUTPATH}tree
python Python/concatenate_trees.py ${OUTPATH}tree_extended
python Python/compute_offsets.py $BASEDIR ${OUTPATH} $SNAP_START $SNAP_END

# cleanup
rm $OUTPATH/_*.hdf5
rm -rf $OUTPATH/columns
rm $OUTPATH/tree.*.hdf5
rm $OUTPATH/tree_extended.*.hdf5
rm dummy
echo "All done."
