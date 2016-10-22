import numpy as np
import h5py
import sys
import os
import time

def compute_offsets(basedir, treedir, snapnum_first, snapnum_last):
    # Get number of subhalos and "offset" of each tree file
    file_nsubs = []
    file_offsets = []
    filenum = 0
    while True:
        filename = treedir + '/tree_extended.%d.hdf5' % (filenum)
        # Only proceed if file exists
        if not os.path.exists(filename):
            break
        # Number of subhalos
        f = h5py.File(filename, 'r')
        file_nsubs.append(f['SubhaloID'].size)
        f.close()
        # Offset
        if filenum > 0:
            file_offsets.append(file_offsets[filenum-1] + file_nsubs[filenum-1])
        else:
            file_offsets.append(0)
        # Next iteration
        filenum += 1
    # Convert to arrays
    file_nsubs = np.array(file_nsubs, dtype=np.int64)
    file_offsets = np.array(file_offsets, dtype=np.int64)
    assert(filenum == len(file_nsubs) == len(file_offsets))

    # Store tree-ordered data in these arrays
    nsubs_trees = sum(file_nsubs)
    SubhaloID_trees = -1*np.ones(nsubs_trees, dtype=np.int64)
    LastProgenitorID_trees = -1*np.ones(nsubs_trees, dtype=np.int64)
    MainLeafProgenitorID_trees = -1*np.ones(nsubs_trees, dtype=np.int64)
    SnapNum_trees = -1*np.ones(nsubs_trees, dtype=np.int16)
    SubfindID_trees = -1*np.ones(nsubs_trees, dtype=np.int32)

    # Iterate over files to read tree-ordered data
    start = time.time()
    print('Reading info from trees...')
    filenum = 0
    while True:
        filename = treedir + '/tree_extended.%d.hdf5' % (filenum)
        # Only proceed if file exists
        if not os.path.exists(filename):
            break
        # Read data from current tree file
        f = h5py.File(filename, 'r')
        locs = slice(file_offsets[filenum], file_offsets[filenum] + file_nsubs[filenum])
        SubhaloID_trees[locs] = f['SubhaloID'][:]
        LastProgenitorID_trees[locs] = f['LastProgenitorID'][:]
        MainLeafProgenitorID_trees[locs] = f['MainLeafProgenitorID'][:]
        SnapNum_trees[locs] = f['SnapNum'][:]
        SubfindID_trees[locs] = f['SubfindID'][:]
        f.close()
        # Next iteration
        filenum += 1
    print('Time: %f' % (time.time() - start))

    # Store offsets in lists of arrays
    start = time.time()
    print('Creating null arrays...')
    RowNum_offsets = []
    SubhaloID_offsets = []
    LastProgenitorID_offsets = []
    MainLeafProgenitorID_offsets = []
    for snapnum in xrange(snapnum_last+1):
        # Get number of subhalos
        if snapnum >= snapnum_first:
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
    print('Time: %f' % (time.time() - start))

    start = time.time()
    print('Computing offsets...')
    for i in xrange(nsubs_trees):
        RowNum_offsets[SnapNum_trees[i]][SubfindID_trees[i]] = i
        SubhaloID_offsets[SnapNum_trees[i]][SubfindID_trees[i]] = SubhaloID_trees[i]
        LastProgenitorID_offsets[SnapNum_trees[i]][SubfindID_trees[i]] = LastProgenitorID_trees[i]
        MainLeafProgenitorID_offsets[SnapNum_trees[i]][SubfindID_trees[i]] = MainLeafProgenitorID_trees[i]
    print('Time: %f' % (time.time() - start))

    # Write to files
    start = time.time()
    print('Writing to files...')

    # Create folder if necessary
    dir_offsets = treedir + '/offsets'
    if not os.path.exists(dir_offsets):
        os.popen('mkdir -p ' + dir_offsets)

    for snapnum in range(snapnum_first, snapnum_last+1):
        filename_offsets = '%s/offsets_%s.hdf5' % (dir_offsets, str(snapnum).zfill(3))
        f_offsets = h5py.File(filename_offsets, 'w')
        f_offsets.create_dataset("FileOffsets", data=file_offsets)
        f_offsets.create_dataset("RowNum", data=RowNum_offsets[snapnum])
        f_offsets.create_dataset("SubhaloID", data=SubhaloID_offsets[snapnum])
        f_offsets.create_dataset("LastProgenitorID", data=LastProgenitorID_offsets[snapnum])
        f_offsets.create_dataset("MainLeafProgenitorID", data=MainLeafProgenitorID_offsets[snapnum])
        f_offsets.close()

    print('Time: %f' % (time.time() - start))


if __name__ == '__main__':

    try:
        basedir = sys.argv[1]
        treedir = sys.argv[2]
        snapnum_first = int(sys.argv[3])
        snapnum_last = int(sys.argv[4])
    except:
        print('Arguments: basedir treedir snapnum_first snapnum_last')
        sys.exit()

    compute_offsets(basedir, treedir, snapnum_first, snapnum_last)
