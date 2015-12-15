# README #

SubLink is a C++ library used for constructing and analyzing merger trees in numerical simulations of galaxy formation.

### Brief description ###

SubLink constructs merger trees for *subhalos* (e.g. the gravitationally bound objects produced by the SUBFIND algorithm) rather than their parent *halos* (e.g. friends-of-friends groups). The algorithm proceeds in two main stages:

1. Finding subhalo descendants: each subhalo from a given snapshot is assigned a unique descendant from the next snapshot by comparing particle IDs in a weighted fashion, putting a higher priority on particles which are more tightly bound. In some cases, a small subhalo is allowed to "skip" a snapshot when finding a descendant, in order to account for situations in which a small subhalo is temporarily "lost" while traversing a larger structure.

2. Once all the subhalo descendants have been identified, this information is re-arranged into a "usable" form (i.e., the merger trees are created). In particular, the data is stored in a depth-first fashion ([Lemson & Springel 2006](http://adsabs.harvard.edu/abs/2006ASPC..351..212L)), which allows for efficient retrieval of merging histories. The first progenitor (or main progenitor) is defined as the one with the "most massive history" behind it ([De Lucia & Blaizot 2007](http://adsabs.harvard.edu/abs/2007MNRAS.375....2D)). Once the merger trees have been constructed, additional subhalo/galaxy properties (e.g., velocity dispersions, metallicities, etc.) can be added in post-processing.

Furthermore, merger trees can be constructed in two different modes: (1) tracking only DM particles or (2) tracking stellar particles and star-forming gas cells. These two types of merger tree are very similar, since they are constructed for the same objects. The baryonic trees have been used in most studies of galaxy formation with the Illustris simulation (http://www.illustris-project.org/).

### Further information ###

* The details of the algorithm can be found in [Rodriguez-Gomez et al. (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.449...49R), while a complete description of the data format is presented in [Nelson et al. (2015)](http://adsabs.harvard.edu/abs/2015A%26C....13...12N).

* Function and type specifications can be found here: http://vrodgom.bitbucket.org/sublink_docs/files.html.

* Further documentation and examples can be found here: https://www.cfa.harvard.edu/~vrodrigu/merger_trees/index.html.

### Code organization ###

Briefly, the code is organized into the following directories:

* **Descendants:** Code for finding subhalo/galaxy descendants.
* **SubhaloTrees:** Code for constructing the merger trees, once the descendants have been determined.
* **Python:** Some Python scripts used during the final stages of merger tree construction. This folder also contains the reading script readtreeHDF5.py, which is recommended for casual users of merger trees.
* **Examples:** Examples of how to use the merger trees using C++. In particular, count_mergers.cpp and stellar_assembly.cpp were used in publications about the merger rate and stellar mass assembly of Illustris galaxies.
* **InputOutput:** Some libraries used for manipulating HDF5 files.
* **Spatial:** A library for performing spatial queries, based on Z-order curves.
* **Util:** Miscellaneous functions and type definitions used by other parts of the code.

### How to use ###

* Run the Makefile.
* Edit and execute the bash scripts "do_everything_part{1-5}.sbatch".

### Dependencies ###

* g++ 4.7 or higher
* HDF5
* OpenMP (optional)

### Authors ###
* Vicente Rodriguez-Gomez (vrodriguez-gomez [at] cfa.harvard.edu)

### Contributors ###
* Dylan Nelson

### Acknowledgments ###

* Some ideas were taken from Shy Genel's halo merger tree code.
* Some code was developed during the course "Systems Development for Computational Science," taught at Harvard by Cris Cecka.