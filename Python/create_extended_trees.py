import numpy as np
import h5py
import sys
import os

"""
Read a tree in "minimal" format and add all the fields
from the subhalo catalogs, which should already be ordered
and found in the columns folder.
Store all fields as separate arrays (datasets) rather than
a single compound datatype.
"""

def create_extended_trees(treedir, basedir, snapnum_last):
    """
    Create extended trees using information in the columns folder.
    """

    # Do this for each independent tree file
    filenum = 0
    filename_minimal = '%s/tree.%d.hdf5' % (treedir, filenum)
    while os.path.exists(filename_minimal):

        # Open minimal tree
        f_minimal = h5py.File(filename_minimal, 'r')
        nrows = f_minimal['Tree'].shape[0]

        # Open new HDF5 file for extended format
        filename_extended = '%s/tree_extended.%d.hdf5' % (treedir, filenum)
        f_extended = h5py.File(filename_extended, 'w')

        # Add "minimal" fields, one by one
        for field_name in f_minimal['Tree'].dtype.names:
            tmp_array = f_minimal['Tree'][field_name][:]

            # Very ad hoc: snapnum type should be int16
            if field_name == 'SnapNum':
                tmp_array = tmp_array.astype(np.int16)

            f_extended.create_dataset(field_name, data=tmp_array)

        # Add additional fields (assume that info inside $TREEDIR/columns is complete).
        for field_name in os.listdir('%s/columns' % (treedir)):
            filename_column = '%s/columns/%s/column.%d.hdf5' % (treedir, field_name, filenum)
            f_column = h5py.File(filename_column, 'r')
            tmp_array = f_column[field_name][:]
            f_extended.create_dataset(field_name, data=tmp_array)
            f_column.close()

        # Close files
        f_minimal.close()
        f_extended.close()

        # Next file
        print('Finished for file %d.' % (filenum))
        filenum += 1
        filename_minimal = treedir + '/tree.' + str(filenum) + '.hdf5'

if __name__ == '__main__':

    # Get arguments
    try:
        basedir = sys.argv[1]
        treedir = sys.argv[2]
        snapnum_last = int(sys.argv[3])
    except:
        print('Arguments: basedir treedir snapnum_last')
        sys.exit()

    create_extended_trees(treedir, basedir, snapnum_last)
