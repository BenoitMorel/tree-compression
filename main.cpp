#include <assert.h>
#include <stdarg.h>

#include "compress_functions.h"
#include "uncompress_functions.h"
#include "datastructure_compression_functions.h"

/* static functions */
static void fatal (const char * format, ...);

/**
 * Run the simple compression on the given newick tree file.
 * @param tree_file path to tree file
 */
void simpleCompressionEqual(const char * tree_file) {
  std::cout << "Simple compression: \n\n";

  // compress the given tree and store the structures to file
  simple_compression(tree_file, "succinct_tree.sdsl", "node_permutation.sdsl",
                      "branch_lengths_uncompressed.sdsl", PRINT_COMPRESSION + PRINT_COMPRESSION_STRUCTURES);

  pll_utree_t * tree = pll_utree_parse_newick(tree_file);

  pll_unode_t * root = searchRoot(tree);
  setTree(root);
  orderTree(root);

  // load the decompress structures from disc
  sdsl::bit_vector succinct_tree_loaded = uncompressSuccinctStructure("succinct_tree.sdsl");
  sdsl::int_vector<> node_permutation_loaded = uncompressSimplePermutation("node_permutation.sdsl");
  std::vector<double> branch_lengths = uncompressBranchLengths("branch_lengths_uncompressed.sdsl");

  // reconstruct the tree
  pll_unode_t * tree_loaded = simple_uncompression(succinct_tree_loaded, node_permutation_loaded, branch_lengths);

  // print the newick reconstruction of the loaded tree
  std::cout << toNewick(tree_loaded) << "\n";

  // print whether input tree and loaded tree are equal
  std::cout << "Trees equal: " << std::boolalpha << treesEqual(root->back, tree_loaded->back) << "\n-----------------------------------------------------------\n";
}

/**
 *  Run the RF compression between the two given newick tree files.
 * @param tree_file1 path to tree file 1
 * @param tree_file2 path to tree file 2
 */
void RFCompressionEqual(const char * tree_file1, const char * tree_file2) {
    std::cout << "RF-compression: \n\n";

    // compress the given trees and store the structures to file
    rf_distance_compression(tree_file1, tree_file2, "edges_to_contract.sdsl",
                      "subtrees_succinct.sdsl", "node_permutations.sdsl",
                      "consensus_branches.sdsl", "non_consensus_branches.sdsl",
                      PRINT_COMPRESSION + PRINT_COMPRESSION_STRUCTURES);

    // load the decompress structures from disc
    sdsl::int_vector<> edges_to_contract_loaded = uncompressRFEdgesToContract("edges_to_contract.sdsl");
    sdsl::bit_vector subtrees_succinct_loaded = uncompressSuccinctStructure("subtrees_succinct.sdsl");
    sdsl::int_vector<> permutations_loaded = uncompressRFSubtreePermutations("node_permutations.sdsl");

    std::vector<double> consensus_branches = uncompressBranchLengths("consensus_branches.sdsl");
    std::vector<double> non_consensus_branches = uncompressBranchLengths("non_consensus_branches.sdsl");

    // load the first tree to recontruct the second tree applying the topology changes
    pll_utree_t * tree1 = pll_utree_parse_newick (tree_file1);
    pll_unode_t * root1 = searchRoot(tree1);
    setTree(root1);
    orderTree(root1);

    // decompress the structures; recontruct the second tree
    pll_unode_t * tree_rf = rf_distance_uncompression(root1, edges_to_contract_loaded, subtrees_succinct_loaded,
                    permutations_loaded, consensus_branches, non_consensus_branches);

    pll_utree_t * tree2 = pll_utree_parse_newick (tree_file2);
    pll_unode_t * root2 = searchRoot(tree2);
    setTree(root2);
    orderTree(root2);

    // std::cout << toNewick(root2) << "\n\n\n";
    // std::cout << toNewick(tree_rf) << "\n";

    assert(treesEqual(tree_rf->back, root2->back));
    // check whether the loaded second tree is equal to the decompressed second tree.
    std::cout << "\n" << std::boolalpha << "trees equal: " << treesEqual(tree_rf->back, root2->back) << "\n";
}

/**
 * Run the compression.
 * Input are paths to two newick tree files.
 */
int main (int argc, const char * argv[])
{
  if (argc != 3)
    fatal (" syntax: %s [newick] [newick]", argv[0]);

  // simple compression on each of the trees
  simpleCompressionEqual(argv[1]);
  simpleCompressionEqual(argv[2]);

  // RF compression between the two tree files
  RFCompressionEqual(argv[1], argv[2]);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void fatal (const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf (stderr, format, argptr);
  va_end(argptr);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}
