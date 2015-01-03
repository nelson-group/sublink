/** @file tree_test.cpp
 * @brief Simple test: find the "root descendant" of all subhalos from
 *        snapshot 60 iteratively and compare with the RootDescendant
 *        pointer included in the trees.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include "../InputOutput/ReadTreeHDF5.hpp"

int main()
{
  // Load tree
  std::string treedir = "/n/ghernquist/vrodrigu/MergerTrees/output/Subhalos/Illustris/L75n455FP";
  int filenum = -1;  // "concatenated" tree file
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);

  // Iterate over subhalos in snapshot of interest.
  auto snap = tree.snapshot(60);
  for (auto it = snap.begin(); it != snap.end(); ++it) {
    auto orig_sub = *it;
    auto root_desc = orig_sub.descendant();
    while (root_desc.is_valid())
      root_desc = root_desc.descendant();

    // Compare with root descendant stored in tree
    if (root_desc != orig_sub.root_descendant())
      std::cerr << "Root descendant does not match!\n";

  }

  return 0;
}
