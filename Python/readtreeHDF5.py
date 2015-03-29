import numpy as np
import h5py
import sys
import os

"""
Simple Python script for reading merger tree HDF5 files
in "database mode," which is optimized for extracting
quantities along the main branch of a subhalo.

For convenience, there are some built-in functions to do
simple tasks, such as:

    get_main_branch
    get_all_progenitors
    get_direct_progenitors
    get_future_branch

The merger trees can also be loaded in "linked-list mode."
This allows for more flexibility and faster tree traversal,
but is only supported in C++ (this approach would be too
memory-expensive in Python).

Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)

--------------- USAGE EXAMPLE: PRINT STELLAR MASS HISTORY ---------------

import readtreeHDF5
treedir = '/n/ghernquist/vrodrigu/MergerTrees/output/Subhalos/Illustris/L75n1820FP'
tree = readtreeHDF5.TreeDB(treedir)
snapnum = 135; subfind_id = 0
branch = tree.get_main_branch(snapnum, subfind_id, keysel=['SubhaloMassType'])
print branch.SubhaloMassType[:, 4]

-------------------------------------------------------------------------

"""

class _Subset(object):
    """
    Used to represent a subset of _AdjacentRows.
    Can be initialized with an integer or boolean array.
    """
    def __init__(self, adj_rows, indices):
        # Copy fields
        for field_name in adj_rows._fields:
            setattr(self, field_name, getattr(adj_rows, field_name)[indices])

class _AdjacentRows(object):
    """
    Used by the TreeDB class. Consists of 
    a set of adjacent rows from the merger tree file.
    Since subhalo IDs are assigned in a depth-first fashion,
    a "chunk" of adjacent rows can represent, e.g., a main branch
    or a "subtree."
    For a given file number and range of row numbers,
    create arrays containing information from the merger tree
    for the specified rows.
    """
    def __init__(self, treefile, row_start, row_end, row_original=None, filenum=-1, keysel=None):
        # Public attributes
        self._row_start = row_start
        self._row_end = row_end
        if row_original == None:
            self._index_given_sub = 0
        else:
            self._index_given_sub = row_original - row_start

        # Only interested in these row numbers:
        self.nrows = row_end - row_start + 1
        locs = slice(row_start, row_end+1)

        # Find out which fields to add
        if keysel == None:
            self._fields = treefile.keys()
        else:
            self._fields = keysel
        
        # Add them
        for field_name in self._fields:
            setattr(self, field_name, treefile[field_name][locs])

    def _get_subset(self, indices):
        return _Subset(self, indices)

class TreeDB:
    """
    Python class to extract information from merger tree files
    in "database mode."
    
    --------------- USAGE EXAMPLE: PRINT STELLAR MASS HISTORY ---------------
    import readtreeHDF5
    treedir = '/n/ghernquist/vrodrigu/MergerTrees/output/Subhalos/Illustris/L75n1820FP'
    tree = readtreeHDF5.TreeDB(treedir)
    snapnum = 135; subfind_id = 0
    branch = tree.get_main_branch(snapnum, subfind_id, keysel=['SubhaloMassType'])
    print branch.SubhaloMassType[:, 4]
    -----------------------------------------------------------------------
    """

    def __init__(self, treedir, name='tree_extended', filenum=-1):
        """
        Create a TreeDB object.
        
        Parameters
        ----------
        treedir : string
                  Directory where the merger tree files are located.
        name : string, optional
               Base name of the HDF5 files, which by default is 'tree_extended'.
        filenum : int, optional
               File number of the HDF5 file of interest; -1 refers to the full,
               "concatenated" tree file.
        """

        # Check that a few files/paths exist
        for rel_path in ['tree_extended.hdf5', 'offsets']:
            if not os.path.exists(treedir + '/' + rel_path):
                print 'Path not found: ' + treedir + '/' + rel_path
                sys.exit()
        if filenum != -1:
            print 'Currently no support for individual tree files.'
            sys.exit()

        # Open tree file
        treefile = h5py.File(treedir + '/tree_extended.hdf5', 'r')

        # Set some attributes
        self._treedir = treedir
        self._name = name
        self._filenum = filenum
        self._treefile = treefile
        self._offset_files = {}  # one per snapshot

    def __del__(self):
        """
        Close files. All open objects become invalid.
        """
        self._treefile.close()
        for f in self._offset_files.values():
            f.close()

    def _get_offset_file(self, snapnum):
        """
        If necessary, add new entry to offset files dictionary.
        Otherwise, return existing one.
        """
        if snapnum not in self._offset_files.keys():
            self._offset_files[snapnum] = h5py.File(
                    '%s/offsets/offsets_%s.hdf5' % (
                    self._treedir, str(snapnum).zfill(3)), 'r')
        return self._offset_files[snapnum]

    def get_main_branch(self, snapnum, subfind_id, keysel=None):
        """
        For a subhalo specified by its snapshot number and Subfind ID,
        return the progenitors along its main branch, i.e. all subhalos
        with IDs between SubhaloID and MainLeafProgenitorID.
        
        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be loaded. By default, all fields are loaded, which
                can be very time- and memory-expensive.
        """
        # Get row number and other info from offset tables
        f = self._get_offset_file(snapnum)
        rownum = f['RowNum'][subfind_id]
        subhalo_id = f['SubhaloID'][subfind_id]
        main_leaf_progenitor_id = f['MainLeafProgenitorID'][subfind_id]
        if rownum == -1:
            print 'Subhalo not found: snapnum = %d, subfind_id = %d.' % (snapnum, subfind_id)
            return None

        # Create branch instance
        row_start = rownum
        row_end = rownum + (main_leaf_progenitor_id - subhalo_id)
        branch = _AdjacentRows(self._treefile, row_start, row_end, keysel=keysel)
        return branch

    def get_all_progenitors(self, snapnum, subfind_id, keysel=None):
        """
        For a subhalo specified by its snapshot number and Subfind ID,
        return all the objects in the subtree which is rooted on the
        subhalo of interest, i.e. all subhalos with IDs between SubhaloID
        and LastProgenitorID. Note that this includes the given subhalo itself.
        
        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be loaded. By default, all fields are loaded, which
                can be very time- and memory-expensive.
        """

        # Get row number and other info from offset tables
        f = self._get_offset_file(snapnum)
        rownum = f['RowNum'][subfind_id]
        subhalo_id = f['SubhaloID'][subfind_id]
        last_progenitor_id = f['LastProgenitorID'][subfind_id]
        if rownum == -1:
            print 'Subhalo not found: snapnum = %d, subfind_id = %d.' % (snapnum, subfind_id)
            return None

        # Create branch instance
        row_start = rownum
        row_end = rownum + (last_progenitor_id - subhalo_id)
        subtree = _AdjacentRows(self._treefile, row_start, row_end, keysel=keysel)
        
        return subtree

    def _get_subhalos_between_root_and_given(self, snapnum, subfind_id, keysel=None):
        """
        Return all subhalos with IDs between RootDescendantID and
        SubhaloID (of the given subhalo), in a depth-first fashion.
        This function is used by "get_forward_branch."
        
        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be loaded. By default, all fields are loaded, which
                can be very time- and memory-expensive.
        """
        # Get row number and other info from offset tables
        f = self._get_offset_file(snapnum)
        rownum = f['RowNum'][subfind_id]
        subhalo_id = f['SubhaloID'][subfind_id]
        if rownum == -1:
            print 'Subhalo not found: snapnum = %d, subfind_id = %d.' % (snapnum, subfind_id)
            return None

        # Get root_descendant_id from merger tree
        root_descendant_id = self._treefile['RootDescendantID'][rownum]

        # We know the row number of the root descendant without searching for it
        row_start = rownum - (subhalo_id - root_descendant_id)
        row_end = rownum

        # Create branch instance
        branch = _AdjacentRows(self._treefile, row_start, row_end, row_original=rownum, keysel=keysel)
        return branch

    def get_direct_progenitors(self, snapnum, subfind_id, **kwargs):
        """
        Return the subhalos for which DescendantID corresponds to the
        current subhalo.

        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be loaded. By default, all fields are loaded, which
                can be very time- and memory-expensive.
        """

        # Make sure that some fields are included.
        include_fields = ['SubhaloID', 'DescendantID']
        if 'keysel' in kwargs:
            tmp_list = kwargs['keysel'][:]  # make copy
            for field_name in include_fields:
                if field_name not in tmp_list:
                    tmp_list.append(field_name)
            kwargs['keysel'] = tmp_list

        subtree = self.get_all_progenitors(snapnum, subfind_id, **kwargs)
        subhalo_id = subtree.SubhaloID[0]  # unique ID of given subhalo
        indices = subtree.DescendantID == subhalo_id
        return subtree._get_subset(indices)

    def get_future_branch(self, snapnum, subfind_id, **kwargs):
        """
        Return the subhalos found in a sort of "forward" branch between
        SubhaloID and RootDescendantID. Note that these subhalos are not
        necessarily stored in adjacent rows, as is the case
        with a main branch (following FirstProgenitor links).

        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be loaded. By default, all fields are loaded, which
                can be very time- and memory-expensive.
        """

        # Make sure that some fields are included.
        include_fields = ['SubhaloID', 'DescendantID', 'RootDescendantID']
        if 'keysel' in kwargs:
            tmp_list = kwargs['keysel'][:]  # make copy
            for field_name in include_fields:
                if field_name not in tmp_list:
                    tmp_list.append(field_name)
            kwargs['keysel'] = tmp_list

        subtree = self._get_subhalos_between_root_and_given(snapnum, subfind_id, **kwargs)
        # Unfortunately, there are no shortcuts in this case and we must
        # proceed iteratively. This is almost at the limit of what one
        # can do when reading trees in "database mode."
        desc_id = subtree.DescendantID[subtree._index_given_sub]
        root_desc_id = subtree.RootDescendantID[subtree._index_given_sub]
        indices = [subtree._index_given_sub]
        while desc_id >= root_desc_id:
            cur_index = np.where(subtree.SubhaloID == desc_id)[0][0]
            indices.append(cur_index)
            desc_id = subtree.DescendantID[cur_index]
        indices = indices[::-1]  # reverse
        return subtree._get_subset(indices)
