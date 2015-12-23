import numpy as np
import h5py
import sys
import os
import time

"""
Compute offsets for 'tree_extended.hdf5'
"""

def compute_offsets(basedir, treedir, snapnum_first, snapnum_last):
    # Read necessary fields from tree
    start = time.time()
    print 'Reading info from tree...'
    f_tree = h5py.File(treedir + '/tree_extended.hdf5', 'r')
    SubhaloID_tree = f_tree['SubhaloID'][:]
    LastProgenitorID_tree = f_tree['LastProgenitorID'][:]
    MainLeafProgenitorID_tree = f_tree['MainLeafProgenitorID'][:]
    SnapNum = f_tree['SnapNum'][:]
    SubfindID = f_tree['SubfindID'][:]
    f_tree.close()
    print 'Time: %f' % (time.time() - start)

    # Store offsets in lists of arrays
    start = time.time()
    print 'Creating empty arrays...'
    RowNum_offsets = []
    SubhaloID_offsets = []
    LastProgenitorID_offsets = []
    MainLeafProgenitorID_offsets = []

    for snapnum in range(snapnum_last+1):
        if snapnum >= snapnum_first:
            # Get number of subhalos
            filename = (basedir + '/groups_%s/fof_subhalo_tab_%s.0.hdf5' %
                        (str(snapnum).zfill(3), str(snapnum).zfill(3)))
            f = h5py.File(filename, 'r')
            nsubs = f['Header'].attrs['Nsubgroups_Total']
            f.close()
        else:
            nsubs = 0

        # Create arrays
        RowNum_offsets.append(-1*np.ones(nsubs, dtype=np.int64))
        SubhaloID_offsets.append(-1*np.ones(nsubs, dtype=np.int64))
        LastProgenitorID_offsets.append(-1*np.ones(nsubs, dtype=np.int64))
        MainLeafProgenitorID_offsets.append(-1*np.ones(nsubs, dtype=np.int64))
    print 'Time: %f' % (time.time() - start)

    start = time.time()
    print 'Computing offsets...'
    for i in xrange(len(SubhaloID_tree)):
        RowNum_offsets[SnapNum[i]][SubfindID[i]] = i
        SubhaloID_offsets[SnapNum[i]][SubfindID[i]] = SubhaloID_tree[i]
        LastProgenitorID_offsets[SnapNum[i]][SubfindID[i]] = LastProgenitorID_tree[i]
        MainLeafProgenitorID_offsets[SnapNum[i]][SubfindID[i]] = MainLeafProgenitorID_tree[i]

    print 'Time: %f' % (time.time() - start)

    # Write to files
    start = time.time()
    print 'Writing to files...'

    # Create folder if necessary
    dir_offsets = treedir + '/offsets'
    if not os.path.exists(dir_offsets):
        os.popen('mkdir -p ' + dir_offsets)

    for snapnum in range(snapnum_first, snapnum_last+1):
        filename_offsets = '%s/offsets_%s.hdf5' % (dir_offsets, str(snapnum).zfill(3))
        f_offsets = h5py.File(filename_offsets, 'w')
        f_offsets.create_dataset("RowNum", data=RowNum_offsets[snapnum])
        f_offsets.create_dataset("SubhaloID", data=SubhaloID_offsets[snapnum])
        f_offsets.create_dataset("LastProgenitorID", data=LastProgenitorID_offsets[snapnum])
        f_offsets.create_dataset("MainLeafProgenitorID", data=MainLeafProgenitorID_offsets[snapnum])
        f_offsets.close()

    print 'Time: %f' % (time.time() - start)


if __name__ == '__main__':

    try:
        basedir = sys.argv[1]
        treedir = sys.argv[2]
        snapnum_first = int(sys.argv[3])
        snapnum_last = int(sys.argv[4])
    except:
        print 'Arguments: basedir treedir snapnum_first snapnum_last'
        sys.exit()

    compute_offsets(basedir, treedir, snapnum_first, snapnum_last)
