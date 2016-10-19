import numpy as np
import h5py
import sys
import os

"""
Concatenate tree files into a single file.
Older version used pytables to create a column "index"
for faster searching. This is not supported anymore.
"""

def concatenate_datasets(treepath):
    # Open HDF5 file for the full tree
    h5filename_full = treepath + '.hdf5'
    h5file_full = h5py.File(h5filename_full, 'w')

    # Get dataset names
    filename = treepath + '.0.hdf5'
    f = h5py.File(filename, 'r')
    field_name_list = f.keys()
    f.close()

    # Add datasets one by one
    for field_name in field_name_list:
        first = True
        filenum = 0
        filename = '%s.%d.hdf5' % (treepath, filenum)
        while os.path.exists(filename):
            # Read data
            f = h5py.File(filename, 'r')
            if first:
                tmp_array = f[field_name][:]
                first = False
            else:
                tmp_array = np.concatenate((tmp_array, f[field_name][:]), axis=0)
            f.close()
            # Next iteration
            filenum += 1
            filename = '%s.%d.hdf5' % (treepath, filenum)

        # Add dataset to file
        h5file_full.create_dataset(field_name, data=tmp_array)
        print('Finished for %s.' % (field_name))

    # Close (and flush) the big file
    h5file_full.close()

if __name__ == '__main__':
    try:
        treepath = sys.argv[1]
    except:
        print('Arguments: treepath')
        sys.exit()
    concatenate_datasets(treepath)
