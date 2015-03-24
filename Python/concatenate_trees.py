import numpy as np
import tables
import sys
import os

"""
Concatenate tree files into a single file.
I use pytables for this because h5py probably
does not support creating search indices.

"""

def concatenate_tables(treepath, particle_resolution):
    # Rough estimate of the number of rows
    total_nparts = particle_resolution**3
    expected_total_nsubs = total_nparts/10

    # Have a peek at the first file to get the
    # data type (pytables description)
    filename = treepath + '.0.hdf5'
    f = tables.openFile(filename)
    tree_description = f.root.Tree.description
    f.close()
    
    # Open HDF5 file for the *full, indexed* tree and create table
    h5filename_full = treepath + '.hdf5'
    h5file_full = tables.openFile(h5filename_full, mode='w', title='Full tree')
    table_full = h5file_full.createTable("/", "Tree", tree_description,
                     expectedrows=expected_total_nsubs)

    # Deal with one file at a time
    filenum = 0
    filename = treepath + '.' + str(filenum) + '.hdf5'
    while os.path.exists(filename):
        # Read data
        f = tables.openFile(filename)
        data = f.root.Tree.read()
        # Add to the full table
        table_full.append(data)
        # Close file        
        f.close()
        # Next file
        print 'Finished for file %d.\n' % (filenum)
        filenum += 1
        filename = treepath + '.' + str(filenum) + '.hdf5'

    # Create indices for the big file
    print 'Creating indices for big file...'
    table_full.cols.SubhaloID.createCSIndex()
    table_full.cols.SubhaloIDRaw.createCSIndex()

    # Close (and flush) the big file
    h5file_full.close()

def concatenate_arrays(treepath, particle_resolution):
    # Open HDF5 file for the full tree
    h5filename_full = treepath + '.hdf5'
    h5file_full = tables.openFile(h5filename_full, mode='w', title='Full tree')

    # Repeat for all fields
    filename = treepath + '.0.hdf5'
    f = tables.openFile(filename)
    field_name_list = f.root.__members__
    f.close()

    for field_name in field_name_list:
        first = True
        filenum = 0
        filename = '%s.%d.hdf5' % (treepath, filenum)
        while os.path.exists(filename):
            f = tables.openFile(filename)
            if first:
                tmp_array = f.root._f_getChild(field_name)[:]
                first = False
            else:
                tmp_array = np.concatenate(
                        (tmp_array, f.root._f_getChild(field_name)[:]), axis=0)
            f.close()
            # Next iteration
            filenum += 1
            filename = '%s.%d.hdf5' % (treepath, filenum)
        # Add to file in one shot
        h5file_full.createArray("/", field_name, tmp_array)
        print 'Finished for %s.' % (field_name)
    # Close (and flush) the big file
    h5file_full.close()

if __name__ == '__main__':
    try:
        treedir = sys.argv[1]
        name = sys.argv[2]
        particle_resolution = int(sys.argv[3])
    except:
        print 'Arguments: treedir name particle_resolution'
        sys.exit()

    treepath = treedir + '/' + name
    if name == 'tree':
        concatenate_tables(treepath, particle_resolution)
    elif name == 'tree_extended':
        concatenate_arrays(treepath, particle_resolution)

