import numpy as np
import sys
import h5py
import illustris_python as il
import time

def calculate_subhalo_offsets(basedir, snapnum, parttype):
    """
    Parameters
    ----------
    basedir : string
    snapnum : int
    parttype : int

    Returns
    -------
    sub_offset : ndarray, dtype=int64
    sub_len : ndarray, dtype=int64
    """
    sub_len_type = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloLenType'])
    nsubs = len(sub_len_type)
    group_len_type = il.groupcat.loadHalos(basedir, snapnum, fields=['GroupLenType'])
    group_nsubs = il.groupcat.loadHalos(basedir, snapnum, fields=['GroupNsubs'])
    ngroups = len(group_len_type)

    if (type(sub_len_type) is dict) or (type(group_len_type) is dict):  # empty return
        return np.array([], dtype="int64"), np.array([], dtype="int64")

    group_offset = np.zeros(ngroups, dtype="int64")
    sub_offset = np.zeros(nsubs, dtype="int64")
    k=0
    for i in range(0, ngroups):
        if (i>0):
            group_offset[i] =  group_offset[i-1] + group_len_type[i-1, parttype]
        if (group_nsubs[i]>0):
            sub_offset[k] = group_offset[i]
            k+=1
            for j in range(1, group_nsubs[i]):
                sub_offset[k] =  sub_offset[k-1] + sub_len_type[k-1, parttype]
                k+=1
    if (k!=nsubs):
        raise Exception('Problem with subhalo offsets.')

    return sub_offset, sub_len_type[:, parttype]

def add_masses(snapnum, in_rad):
    """Add masses corresponding to different stellar components of a galaxy.
    """
    parttype_stars = 4
    
    # For performance checks
    start_all = time.time()

    # AD HOC
    if ('L75n1820FP' in basedir) and (snapnum in [53,55]):
        return

    # Write output to this file
    if in_rad:
        f_write = h5py.File('%s/galaxies_in_rad_%s.hdf5' % (assemblydir, str(snapnum).zfill(3)), 'w')
    else:
        f_write = h5py.File('%s/galaxies_%s.hdf5' % (assemblydir, str(snapnum).zfill(3)), 'w')

    # Load particle assembly info
    start = time.time()
    print('Loading info from assembly file...')
    f = h5py.File('%s/stars_%s.hdf5' % (assemblydir, str(snapnum).zfill(3)), 'r')
    # Only proceed if file is non-empty
    if len(list(f.keys())) == 0:
        print("No stellar particles in snapshot %d. Skipping..." % (snapnum))
        f.close()
        f_write.close()
        return
    in_situ = f['InSitu'][:]
    after_infall = f['AfterInfall'][:]
    subfind_id = f['SubfindID'][:]
    subfind_id_at_formation = f['SubfindIDAtFormation'][:]
    accretion_origin = f['AccretionOrigin'][:]
    merger_mass_ratio = f['MergerMassRatio'][:]
    f.close()
    print('Time: %f s.' % (time.time() - start))

    # Sanity check
    locs1 = np.logical_or(subfind_id_at_formation == -1, subfind_id == -1)
    locs2 = np.logical_and(after_infall == -1, in_situ != 1)
    assert np.sum(locs1) == np.sum(locs2)

    # Boolean arrays for different stellar particle classes
    in_situ_flag = in_situ == 1
    ex_situ_flag = in_situ == 0
    before_infall_flag = after_infall == 0
    after_infall_flag = after_infall == 1
    completed_merger_flag = accretion_origin == 0
    ongoing_merger_flag = accretion_origin == 1
    flyby_flag = accretion_origin == 2
    formed_outside_galaxy_flag = subfind_id_at_formation == -1
    merger_flag = merger_mass_ratio >= 0
    major_merger_flag = merger_mass_ratio >= 1.0/4.0
    majorminor_merger_flag = merger_mass_ratio >= 1.0/10.0

    # Load SUBFIND catalog
    start = time.time()
    print('Loading SUBFIND catalog...')
    cat_subs = il.groupcat.loadSubhalos(basedir, snapnum,
            fields=['SubhaloHalfmassRadType', 'SubhaloMassType'])
    print('Time: %f s.' % (time.time() - start))
    # Only proceed if there is at least one subhalo
    nsubs = cat_subs['count']
    if nsubs == 0:
        print("No subhalos in snapshot %d. Skipping..." % (snapnum))
        f_write.close()
        return
    # Get subhalo properties
    sub_halfmassrad = cat_subs['SubhaloHalfmassRadType'][:, parttype_stars]

    # Load some info from snapshot header
    f = h5py.File('%s/snapdir_%s/snap_%s.0.hdf5' % (basedir, str(snapnum).zfill(3), str(snapnum).zfill(3)), 'r')
    header = dict(f['Header'].attrs.items())
    h = header['HubbleParam']
    box_size = header['BoxSize']
    f.close()

    # Load stellar particle info from snapshot file
    start = time.time()
    print('Loading info from snapshot...')
    all_masses = il.snapshot.loadSubset(basedir, snapnum, parttype_stars,
            fields=['Masses'])
    all_formtimes = il.snapshot.loadSubset(basedir, snapnum, parttype_stars,
            fields=['GFM_StellarFormationTime'])
    if in_rad:
        all_pos = il.snapshot.loadSubset(basedir, snapnum, parttype_stars,
                fields=['Coordinates']).astype(np.float32)  # note conversion to float
    print('Time: %f s.' % (time.time() - start))

    # Calculate subhalo offsets
    print('Calculating subhalo offsets...')
    start = time.time()
    sub_offset, sub_len = calculate_subhalo_offsets(basedir, snapnum, parttype_stars)
    print('Time: %f s.' % (time.time() - start))

    # Calculate "ex situ" stellar masses of different components.
    print('Calculating "ex situ" stellar mass fractions...')
    start = time.time()
    StellarMassInSitu = np.zeros(nsubs, dtype=np.float32)
    StellarMassExSitu = np.zeros(nsubs, dtype=np.float32)
    StellarMassTotal = np.zeros(nsubs, dtype=np.float32)
    StellarMassBeforeInfall = np.zeros(nsubs, dtype=np.float32)
    StellarMassAfterInfall = np.zeros(nsubs, dtype=np.float32)
    StellarMassFormedOutsideGalaxies = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromCompletedMergers = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromCompletedMergersMajor = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromCompletedMergersMajorMinor = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromOngoingMergers = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromOngoingMergersMajor = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromOngoingMergersMajorMinor = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromFlybys = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromFlybysMajor = np.zeros(nsubs, dtype=np.float32)
    StellarMassFromFlybysMajorMinor = np.zeros(nsubs, dtype=np.float32)
    
    #~ StellarMassFromMergers = np.zeros(nsubs, dtype=np.float32)
    #~ StellarMassFromMajorMergers = np.zeros(nsubs, dtype=np.float32)
    #~ StellarMassFromMajorMinorMergers = np.zeros(nsubs, dtype=np.float32)
    for sub_index in range(nsubs):
        # Careful: note difference between slice and boolean array
        locs = slice(sub_offset[sub_index], sub_offset[sub_index] + sub_len[sub_index])
        assert np.sum(in_situ[locs] == -1) == 0

        # If zero stellar particles, nothing to do.
        if locs.stop - locs.start <= 0:
            continue

        # Take care of wind particles
        locs_notwind = all_formtimes[locs] >= 0
        locs_valid = locs_notwind[:]

        # Only consider particles within 2*rhalf.
        if in_rad:
            # Distances to the most bound particle (inclusive)
            pos = all_pos[locs]  # remove wind particles later
            dx = pos[:] - pos[0]  # including first particle (r=0)
            # Periodic boundary conditions 
            dx = dx - (np.abs(dx) > 0.5*box_size) * np.copysign(box_size, dx - 0.5*box_size)
            r = np.sqrt(np.sum(dx**2, axis=1))
            rhalf = sub_halfmassrad[sub_index]
            locs_valid = locs_notwind * (r < 2.0*rhalf)
        
        StellarMassInSitu[sub_index] = np.sum(all_masses[locs][in_situ_flag[locs]*locs_valid])
        StellarMassExSitu[sub_index] = np.sum(all_masses[locs][ex_situ_flag[locs]*locs_valid])
        StellarMassTotal[sub_index] = np.sum(all_masses[locs][locs_valid])
        StellarMassBeforeInfall[sub_index] = np.sum(all_masses[locs][before_infall_flag[locs]*locs_valid])
        StellarMassAfterInfall[sub_index] = np.sum(all_masses[locs][after_infall_flag[locs]*locs_valid])
        StellarMassFormedOutsideGalaxies[sub_index] = np.sum(all_masses[locs][formed_outside_galaxy_flag[locs]*locs_valid])
        StellarMassFromCompletedMergers[sub_index] = np.sum(all_masses[locs][completed_merger_flag[locs]*locs_valid])
        StellarMassFromCompletedMergersMajor[sub_index] = np.sum(all_masses[locs][completed_merger_flag[locs]*major_merger_flag[locs]*locs_valid])
        StellarMassFromCompletedMergersMajorMinor[sub_index] = np.sum(all_masses[locs][completed_merger_flag[locs]*majorminor_merger_flag[locs]*locs_valid])
        StellarMassFromOngoingMergers[sub_index] = np.sum(all_masses[locs][ongoing_merger_flag[locs]*locs_valid])
        StellarMassFromOngoingMergersMajor[sub_index] = np.sum(all_masses[locs][ongoing_merger_flag[locs]*major_merger_flag[locs]*locs_valid])
        StellarMassFromOngoingMergersMajorMinor[sub_index] = np.sum(all_masses[locs][ongoing_merger_flag[locs]*majorminor_merger_flag[locs]*locs_valid])
        StellarMassFromFlybys[sub_index] = np.sum(all_masses[locs][flyby_flag[locs]*locs_valid])
        StellarMassFromFlybysMajor[sub_index] = np.sum(all_masses[locs][flyby_flag[locs]*major_merger_flag[locs]*locs_valid])
        StellarMassFromFlybysMajorMinor[sub_index] = np.sum(all_masses[locs][flyby_flag[locs]*majorminor_merger_flag[locs]*locs_valid])

        #~ StellarMassFromMergers[sub_index] = np.sum(all_masses[locs][merger_flag[locs]*locs_valid])
        #~ StellarMassFromMajorMergers[sub_index] = np.sum(all_masses[locs][major_merger_flag[locs]*locs_valid])
        #~ StellarMassFromMajorMinorMergers[sub_index] = np.sum(all_masses[locs][majorminor_merger_flag[locs]*locs_valid])

    print('Time: %f s.' % (time.time() - start))

    # Write to file
    print('Writing to file...')
    start = time.time()
    f_write.create_dataset("StellarMassInSitu", data=StellarMassInSitu)
    f_write.create_dataset("StellarMassExSitu", data=StellarMassExSitu)
    f_write.create_dataset("StellarMassTotal", data=StellarMassTotal)
    f_write.create_dataset("StellarMassBeforeInfall", data=StellarMassBeforeInfall)
    f_write.create_dataset("StellarMassAfterInfall", data=StellarMassAfterInfall)
    f_write.create_dataset("StellarMassFormedOutsideGalaxies", data=StellarMassFormedOutsideGalaxies)
    f_write.create_dataset("StellarMassFromCompletedMergers", data=StellarMassFromCompletedMergers)
    f_write.create_dataset("StellarMassFromCompletedMergersMajor", data=StellarMassFromCompletedMergersMajor)
    f_write.create_dataset("StellarMassFromCompletedMergersMajorMinor", data=StellarMassFromCompletedMergersMajorMinor)
    f_write.create_dataset("StellarMassFromOngoingMergers", data=StellarMassFromOngoingMergers)
    f_write.create_dataset("StellarMassFromOngoingMergersMajor", data=StellarMassFromOngoingMergersMajor)
    f_write.create_dataset("StellarMassFromOngoingMergersMajorMinor", data=StellarMassFromOngoingMergersMajorMinor)
    f_write.create_dataset("StellarMassFromFlybys", data=StellarMassFromFlybys)
    f_write.create_dataset("StellarMassFromFlybysMajor", data=StellarMassFromFlybysMajor)
    f_write.create_dataset("StellarMassFromFlybysMajorMinor", data=StellarMassFromFlybysMajorMinor)

    #~ f_write.create_dataset("StellarMassFromMergers", data=StellarMassFromMergers)
    #~ f_write.create_dataset("StellarMassFromMajorMergers", data=StellarMassFromMajorMergers)
    #~ f_write.create_dataset("StellarMassFromMajorMinorMergers", data=StellarMassFromMajorMinorMergers)
    f_write.close()
    print('Time: %f s.' % (time.time() - start))
    
    print('Finished for snapshot %d.\n' % (snapnum))
    print('Total time: %f s.\n' % (time.time() - start_all))
    
    
if __name__ == '__main__':

    """
    if in_rad == True, only consider stellar particles within 2*rhalf.
    
    Some invariants (up to rounding errors):
    
    StellarMassInSitu + StellarMassExSitu == StellarMassTotal
    StellarMassBeforeInfall + StellarMassAfterInfall+ StellarMassFormedOutsideGalaxies == StellarMassExSitu
    StellarMassFromMergers + StellarMassFormedOutsideGalaxies == StellarMassExSitu
    StellarMassFromCompletedMergers + StellarMassFromOngoingMergers + StellarMassFromFlybys == StellarMassFromMergers
    """

    if len(sys.argv) == 6:
        basedir = sys.argv[1]
        assemblydir = sys.argv[2]
        in_rad = bool(int(sys.argv[3]))
        snapnum_first = int(sys.argv[4])
        snapnum_last = int(sys.argv[5])
        for snapnum in reversed(range(snapnum_first, snapnum_last + 1)):
            add_masses(snapnum, in_rad)
    elif len(sys.argv) == 5:
        basedir = sys.argv[1]
        assemblydir = sys.argv[2]
        in_rad = bool(int(sys.argv[3]))
        snapnum = int(sys.argv[4])
        add_masses(snapnum, in_rad)
    else:
        print('Arguments: basedir assemblydir in_rad snapnum [snapnum_last]')
        sys.exit()
