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
  // Load Tree object
  std::string treedir = "/n/ghernquist/vrodrigu/MergerTrees/output/Subhalos/Illustris/L75n455FP";
  int filenum = -1;  // "concatenated" tree file
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);

  // Iterate over subhalos in snapshot of interest.
  std::size_t count = 0;
  auto snap = tree.snapshot(60);
  for (auto it = snap.begin(); it != snap.end(); ++it) {
    auto orig_sub = *it;
    auto cur_sub = orig_sub;

    // Find root descendant iteratively.
    while (cur_sub.descendant().is_valid())
      cur_sub = cur_sub.descendant();

    // Compare with root descendant as stored in tree
    if (cur_sub != orig_sub.root_descendant())
      std::cerr << "Root descendant does not match!\n";

    // Print some output occasionally
    ++count;
    if (count%100000 == 0)
      std::cout << "Already checked " << count << " subhalos.\n";
  }

  std::cout << "Checked " << count << " subhalos.\n";

  return 0;
}
