import numpy as np
import h5py
import sys
import os

"""
Read a tree in "minimal" format and add all Subhalo fields
from the Subfind catalogs, as well as a few FoF group quantities.
Store each new field in a separate file. Deal with one field at a time.
"""

def read_block(basedir, snapnum, field_name, group_name):
    """
    Read block for a certain snapshot (all files).
    We assume that if the first file is empty/inexistent,
    then all of them are.
    """
    assert group_name == 'Subhalo' or group_name == 'Group'
    
    # By default, the return value is None:
    retval = None

    # Iterate over files    
    filenum = 0
    skip = 0
    doneflag = False
    while not doneflag:
        #name = 'groups'
        name = 'fof_subhalo_tab'
        filename = (basedir + '/groups_%s/%s_%s.%d.hdf5' %
                (str(snapnum).zfill(3), name, str(snapnum).zfill(3), filenum))

        # Only proceed if file exists
        if not os.path.exists(filename):
            if filenum == 0:
                print 'Snapshot %d seems to be empty. Skipping...' % (snapnum)
                return None
            else:
                # Note that nfiles is only defined for filenum > 0
                filenum += 1
                if filenum == nfiles:
                    doneflag = True
                continue

        # Some files fail when opening. Therefore,
        # only proceed if file is not corrupt:
        try:
            # Open file
            f = h5py.File(filename, 'r')
        except:
            print 'Corrupt file:', filename
            print 'Skipping...'
            # Get next file
            filenum += 1
            continue

        # Only proceed if field exists
        if field_name not in f[group_name].keys():
            f.close()
            if filenum == 0:
                #print 'Snapshot %d seems to be empty. Skipping...' % (snapnum)
                return None
            else:
                # Note that nfiles is only defined for filenum > 0
                filenum += 1
                if filenum == nfiles:
                    doneflag = True
                continue

        # Array "handle"
        handle = f[group_name][field_name]

        # Allocate memory
        if filenum == 0:
            nfiles = f['Header'].attrs['NumFiles']
            if group_name == 'Subhalo':
                total_num_objects = f['Header'].attrs['Nsubgroups_Total']
            else:
                total_num_objects = f['Header'].attrs['Ngroups_Total']
            
            if len(handle.shape) > 1:
                retval = np.empty((total_num_objects, handle.shape[1]), dtype=handle.dtype)
            else:
                retval = np.empty(total_num_objects, dtype=handle.dtype)

        # Read data
        if group_name == 'Subhalo':
            cur_num_objects = f['Header'].attrs['Nsubgroups_ThisFile']
        else:
            cur_num_objects = f['Header'].attrs['Ngroups_ThisFile']
            
        retval[skip:(skip+cur_num_objects)] = handle[:]
        skip += cur_num_objects

        # Close file and get next
        f.close()
        filenum += 1
        if filenum == nfiles:
            doneflag = True

    # Sanity check:
    if skip != total_num_objects:
        print 'Bad number of objects when reading %s, snapshot %d!' % (field_name, snapnum)

    #print 'Finished loading %s (cur_dimension %d) from snapshot %d.' % (field_name, cur_dimension, snapnum)
    return retval

# --------------------------- MAIN FUNCTIONS -------------------------------

def create_extra_columns(treedir, basedir, snapnum_first, snapnum_last,
        num_jobs, job_position):
    """
    Create HDF5 files with the merger-tree-ordered data corresponding to
    all subhalo quantities, as well as a few select halo quantities.
    """

    # Load snapshot numbers and Subfind IDs, for once, into memory
    print 'Loading snapnum and subfind_id...'
    snapnums_dict = {}
    subfind_ids_dict = {}
    nrows_dict = {}
    # Repeat for every tree file
    filenum = 0
    filename = treedir + '/tree.' + str(filenum) + '.hdf5'
    while os.path.exists(filename):
        # Read data
        f = h5py.File(filename, 'r')
        snapnums_dict[filenum] = f['Tree']['SnapNum'][:]
        subfind_ids_dict[filenum] = f['Tree']['SubfindID'][:]
        nrows_dict[filenum] = f['Tree'].shape[0]
        f.close()
        # Next file
        print 'Finished for file %d.' % (filenum)
        filenum += 1
        filename = treedir + '/tree.' + str(filenum) + '.hdf5'

    # Have a peek at first file from last snapshot and
    # get info for all subhalo and FoF group quantities.
    # (Note that 'shape' is only used used to retrieve info about the second
    # dimension of the field, if any.)

    #name = 'groups'
    name = 'fof_subhalo_tab'
    groupfilename = (basedir + '/groups_%s/%s_%s.0.hdf5' %
                (str(snapnum_last).zfill(3), name, str(snapnum_last).zfill(3)))

    f = h5py.File(groupfilename, 'r')
    subhalo_quantities = f['Subhalo'].keys()
    #group_quantities = ['Group_M_Crit200', 'Group_M_Mean200', 'Group_M_TopHat200']
    group_quantities = f['Group'].keys()
    field_name_list = subhalo_quantities + group_quantities
    group_name_list = (len(subhalo_quantities) * ['Subhalo'] +
                       len(group_quantities) * ['Group'])
    field_dtype_list = []; field_shape_list = []
    for field_name, group_name in zip(field_name_list, group_name_list):
        field_dtype_list.append(f[group_name][field_name].dtype)
        field_shape_list.append(f[group_name][field_name].shape)
    f.close()

    # We only deal with a subset of the fields
    for k in range(job_position, len(field_name_list), num_jobs):
        field_name = field_name_list[k]
        group_name = group_name_list[k]
        field_dtype = field_dtype_list[k]
        field_shape = field_shape_list[k]

        # Small check
        assert group_name == 'Subhalo' or group_name == 'Group'

        # Load data from catalog
        data_catalog = {}
        sub_gr_nr = {}  # only needed if dealing with FoF group quantity
        print 'Loading %s from catalog...' % (field_name)
        for snapnum in range(snapnum_first, snapnum_last+1):
            # Note that some snapshots may be empty:
            retval = read_block(basedir, snapnum, field_name, group_name)
            if retval != None:
                data_catalog[snapnum] = retval

            if group_name == 'Group':
                retval = read_block(basedir, snapnum, "SubhaloGrNr", "Subhalo")
                if retval != None:
                    sub_gr_nr[snapnum] = retval

        # Create write folder if necessary
        dir_newcolumn = treedir + '/columns/' + field_name
        if not os.path.exists(dir_newcolumn):
            os.popen('mkdir -p ' + dir_newcolumn)

        print 'Ordering values and writing to files...'
        # Repeat for every tree file
        for filenum in nrows_dict.keys():
            # Create HDF5 file for the new column
            filename_newcolumn = dir_newcolumn + '/column.' + str(filenum) + '.hdf5'
            f_newcolumn = h5py.File(filename_newcolumn, 'w')

            # Iterate through the tree rows, adding info to an array in memory.
            if len(field_shape) == 1:
                data_for_table = np.empty(nrows_dict[filenum], dtype=field_dtype)
            elif len(field_shape) == 2:
                data_for_table = np.empty((nrows_dict[filenum], field_shape[1]), dtype=field_dtype)
            else:
                print 'No support for %d-dimensional data.' % (len(field_shape))
                raise

            if group_name == 'Subhalo':
                for i in xrange(nrows_dict[filenum]):
                    data_for_table[i] = data_catalog[snapnums_dict[filenum][i]][subfind_ids_dict[filenum][i]]
            else:  # FoF groups
                for i in xrange(nrows_dict[filenum]):
                    cur_snap = snapnums_dict[filenum][i]
                    cur_subfind_id = subfind_ids_dict[filenum][i]
                    data_for_table[i] = data_catalog[cur_snap][sub_gr_nr[cur_snap][cur_subfind_id]]

            # Save to HDF5 file
            f_newcolumn.create_dataset(field_name, data=data_for_table)
            f_newcolumn.close()
            
        print 'Finished writing %s to files.\n' % (field_name)

if __name__ == '__main__':
    # Get arguments
    try:
        basedir = sys.argv[1]
        treedir = sys.argv[2]
        snapnum_first = int(sys.argv[3])
        snapnum_last = int(sys.argv[4])
        num_jobs = int(sys.argv[5])
        job_position = int(sys.argv[6])
    except:
        print 'Arguments: basedir treedir snapnum_first snapnum_last use_mpi num_jobs job_position'
        sys.exit()

    # Note: must run a job array with, e.g., num_jobs=8 and job_position=0-7
    create_extra_columns(treedir, basedir, snapnum_first, snapnum_last, num_jobs, job_position)
