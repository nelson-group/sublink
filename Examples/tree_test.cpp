/** @file tree_test.cpp
 * @brief Example code that shows how to work with SubLink merger trees.
 *
 * This example finds the "root descendant" of all subhalos from
 * snapshot 60 iteratively and compares them with the RootDescendant
 * pointer included in the trees.
 *
 * @author Vicente Rodriguez-Gomez (v.rodriguez@irya.unam.mx)
 */

// Include last_progenitor, main_leaf_progenitor, and root_descendant
#define EXTRA_POINTERS

#include "../InputOutput/ReadTreeHDF5.hpp"

int main()
{
  // Load Tree object
  std::string treedir = "/n/ghernquist/vrodrigu/MergerTrees/output/Galaxies/Illustris/L75n455FP";
  int filenum = -1;  // read from all merger tree files
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);

  // Iterate over subhalos in snapshot of interest.
  std::size_t count = 0;
  auto snap = tree.snapshot(60);
  for (auto it = snap.begin(); it != snap.end(); ++it) {

    // Find root descendant iteratively and compare with tree value.
    auto orig_sub = *it;
    auto cur_sub = orig_sub;
    while (cur_sub.descendant().is_valid())
      cur_sub = cur_sub.descendant();
    if (cur_sub != orig_sub.root_descendant())
      std::cerr << "Root descendant does not match!\n";

    // Find main leaf progenitor iteratively and compare with tree value.
    cur_sub = orig_sub;
    while (cur_sub.first_progenitor().is_valid())
      cur_sub = cur_sub.first_progenitor();
    if (cur_sub != orig_sub.main_leaf_progenitor())
      std::cerr << "Main leaf progenitor does not match!\n";

    // Print some output occasionally
    ++count;
    if (count%100000 == 0)
      std::cout << "Already checked " << count << " subhalos.\n";
  }

  std::cout << "Checked " << count << " subhalos.\n";

  return 0;
}
