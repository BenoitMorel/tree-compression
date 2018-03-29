#include <assert.h>
#include <stdarg.h>

#include "compress_functions.h"
#include "uncompress_functions.h"

/* static functions */
static void fatal (const char * format, ...);

int main (int argc, const char * argv[])
{
  if (argc != 3)
    fatal (" syntax: %s [newick] [newick]", argv[0]);

  std::cout << "Simple compression of tree 1: \n";
  simple_compression(argv[1], "output_files/succinct_tree.sdsl", "output_files/node_permutation.sdsl", 0/*PRINT_COMPRESSION*/);
  std::cout << "\n";

  pll_utree_t * tree = pll_utree_parse_newick (argv[1]);
  pll_unode_t * root = searchRoot(tree);
  setTree(root);
  orderTree(root);

  //printTree(root);
  std::cout << "\n";

  sdsl::bit_vector succinct_tree_loaded;
  load_from_file(succinct_tree_loaded, "output_files/succinct_tree.sdsl");

  sdsl::int_vector<> node_permutation_loaded;
  load_from_file(node_permutation_loaded, "output_files/node_permutation.sdsl");

  pll_unode_t * tree_loaded = simple_uncompression(succinct_tree_loaded, node_permutation_loaded);
  //printTree(tree_loaded);

  std::cout << "Trees equal: " << std::boolalpha << treesEqual(root, tree_loaded) << "\n";

  std::cout << "Simple compression of tree 2: \n";
  simple_compression(argv[2], "output_files/succinct_tree.sdsl", "output_files/node_permutation.sdsl", PRINT_COMPRESSION);
  std::cout << "\n";

  std::cout << "RF-compression: \n";
  rf_distance_compression(argv[1], argv[2], "output_files/edges_to_contract.sdsl",
                    "output_files/subtrees_succinct.sdsl", "output_files/node_permutations.sdsl", PRINT_COMPRESSION);
  std::cout << "\n";

  sdsl::int_vector<> edges_to_contract_loaded;
  load_from_file(edges_to_contract_loaded, "output_files/edges_to_contract.sdsl");

  sdsl::bit_vector subtrees_succinct_loaded;
  load_from_file(subtrees_succinct_loaded, "output_files/subtrees_succinct.sdsl");

  sdsl::int_vector<> permutations_loaded;
  load_from_file(permutations_loaded, "output_files/node_permutations.sdsl");


  pll_utree_t * tree1 = pll_utree_parse_newick (argv[1]);
  pll_unode_t * root1 = searchRoot(tree1);
  setTree(root1);
  orderTree(root1);

  std::cout << "\n\n\n-------------------------RF-UNCOMPRESSION-------------------------\n\n\n";

  pll_unode_t * tree_rf = rf_distance_uncompression(root1, edges_to_contract_loaded, subtrees_succinct_loaded, permutations_loaded);

  pll_utree_t * tree2 = pll_utree_parse_newick (argv[2]);
  pll_unode_t * root2 = searchRoot(tree2);
  setTree(root2);
  orderTree(root2);

  // std::cout << "\n";
  // printTree(tree_rf);
  // std::cout << "\n";
  // printTree(root2);
  std::cout << "\n" << std::boolalpha << "trees equal: " << treesEqual(tree_rf, root2) << "\n";

  /*for (int i = 1; i < 1001; i++) {
      //std::cout << "RF-compression: \n";

      std::stringstream ss1;
      ss1 << "500\ \(500000\ gen\)/tree_" << i << ".nwk";

      std::stringstream ss2;
      ss2 << "500\ \(500000\ gen\)/tree_" << (i+1) << ".nwk";

      //std::cout << ss1.str() << " & " << ss2.str() << "\n";

      //simple_compression(strdup(ss1.str().c_str()));
      rf_distance_compression(strdup(ss1.str().c_str()), strdup(ss2.str().c_str()));
  }*/

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
