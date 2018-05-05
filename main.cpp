#include <assert.h>
#include <stdarg.h>

#include "compress_functions.h"
#include "uncompress_functions.h"
#include "datastructure_compression_functions.h"

/* static functions */
static void fatal (const char * format, ...);

void simpleCompressionEqual(const char * tree_file) {
  std::cout << "Simple compression: \n\n";
  simple_compression(tree_file, "output_files/succinct_tree.sdsl", "output_files/node_permutation.sdsl",
                      "output_files/branch_lengths_uncompressed.sdsl", PRINT_COMPRESSION + PRINT_COMPRESSION_STRUCTURES);

  pll_utree_t * tree = pll_utree_parse_newick(tree_file);

  pll_unode_t * root = searchRoot(tree);
  setTree(root);
  orderTree(root);

  sdsl::bit_vector succinct_tree_loaded = uncompressSuccinctStructure("output_files/succinct_tree.sdsl");
  sdsl::int_vector<> node_permutation_loaded = uncompressSimplePermutation("output_files/node_permutation.sdsl");
  std::vector<double> branch_lengths = uncompressBranchLengths("output_files/branch_lengths_uncompressed.sdsl");

  pll_unode_t * tree_loaded = simple_uncompression(succinct_tree_loaded, node_permutation_loaded, branch_lengths);

  // std::cout << toNewick(root) << "\n\n\n";
  // std::cout << toNewick(tree_loaded) << "\n";

  std::cout << "Trees equal: " << std::boolalpha << treesEqual(root->back, tree_loaded->back) << "\n-----------------------------------------------------------\n";
}

void RFCompressionEqual(const char * tree_file1, const char * tree_file2) {
    std::cout << "RF-compression: \n\n";
    rf_distance_compression(tree_file1, tree_file2, "output_files/edges_to_contract.sdsl",
                      "output_files/subtrees_succinct.sdsl", "output_files/node_permutations.sdsl",
                      "output_files/consensus_branches.sdsl", "output_files/non_consensus_branches.sdsl",
                      PRINT_COMPRESSION + PRINT_COMPRESSION_STRUCTURES);

    sdsl::int_vector<> edges_to_contract_loaded = uncompressRFEdgesToContract("output_files/edges_to_contract.sdsl");
    sdsl::bit_vector subtrees_succinct_loaded = uncompressSuccinctStructure("output_files/subtrees_succinct.sdsl");
    sdsl::int_vector<> permutations_loaded = uncompressRFSubtreePermutations("output_files/node_permutations.sdsl");

    std::vector<double> consensus_branches = uncompressBranchLengths("output_files/consensus_branches.sdsl");
    std::vector<double> non_consensus_branches = uncompressBranchLengths("output_files/non_consensus_branches.sdsl");

    pll_utree_t * tree1 = pll_utree_parse_newick (tree_file1);
    pll_unode_t * root1 = searchRoot(tree1);
    setTree(root1);
    orderTree(root1);

    //std::cout << "\n\n\n-------------------------RF-UNCOMPRESSION-------------------------\n\n\n";

    pll_unode_t * tree_rf = rf_distance_uncompression(root1, edges_to_contract_loaded, subtrees_succinct_loaded,
                    permutations_loaded, consensus_branches, non_consensus_branches);

    pll_utree_t * tree2 = pll_utree_parse_newick (tree_file2);
    pll_unode_t * root2 = searchRoot(tree2);
    setTree(root2);
    orderTree(root2);

    // std::cout << toNewick(root2) << "\n\n\n";
    // std::cout << toNewick(tree_rf) << "\n";

    assert(treesEqual(tree_rf->back, root2->back));
    std::cout << "\n" << std::boolalpha << "trees equal: " << treesEqual(tree_rf->back, root2->back) << "\n";

}

int main (int argc, const char * argv[])
{
  if (argc != 3)
    fatal (" syntax: %s [newick] [newick]", argv[0]);

  simpleCompressionEqual(argv[1]);
  simpleCompressionEqual(argv[2]);

  // for (int i = 1; i < 1002; i++) {
  //     std::cout << i << "\n";
  //
  //     std::stringstream ss1;
  //     ss1 << "500\ \(500000\ gen\)/tree_" << i << ".nwk";
  //
  //     simpleCompressionEqual(strdup(ss1.str().c_str()));
  // }

  RFCompressionEqual(argv[1], argv[2]);

  //for (int i = 1; i < 1001; i++) {

      // for (int j = i+1; j < 1001; j+=20) {
      //     //std::cout << "RF-compression: \n";
      //
      //     std::cout << i << " & " << j << ":  ";
      //
      //     std::stringstream ss1;
      //     ss1 << "500\ \(500000\ gen\)/tree_" << i << ".nwk";
      //
      //     std::stringstream ss2;
      //     ss2 << "500\ \(500000\ gen\)/tree_" << j << ".nwk";
      //
      //     //std::cout << ss1.str() << " & " << ss2.str() << "\n";
      //
      //     //simpleCompressionEqual(strdup(ss1.str().c_str()));
      //     RFCompressionEqual(strdup(ss1.str().c_str()), strdup(ss2.str().c_str()));
      //     //rf_distance_compression(strdup(ss1.str().c_str()), strdup(ss2.str().c_str()));
      // }


      //std::cout << "RF-compression: \n";

      // std::cout << i << " & " << i+1 << ":  ";
      //
      // std::stringstream ss1;
      // ss1 << "500\ \(500000\ gen\)/tree_" << i << ".nwk";
      //
      // std::stringstream ss2;
      // ss2 << "500\ \(500000\ gen\)/tree_" << (i+1) << ".nwk";

      //std::cout << ss1.str() << " & " << ss2.str() << "\n";

      //simpleCompressionEqual(strdup(ss1.str().c_str()));
      //RFCompressionEqual(strdup(ss1.str().c_str()), strdup(ss2.str().c_str()));
      //rf_distance_compression(strdup(ss1.str().c_str()), strdup(ss2.str().c_str()));
  //}

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
