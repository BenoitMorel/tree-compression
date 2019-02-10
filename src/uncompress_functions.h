#include <assert.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <libpll/pll_tree.h>

#include "modified_library_functions.h"
#ifdef __cplusplus
}
#endif

#include <sdsl/bit_vectors.hpp>
#include <algorithm>

#include "util.h"

/**
 * Decompresses a tree stored with simple compression.
 * @param  succinct_structure vector containing succint structure
 * @param  node_permutation   vector containing node permutation
 * @param  branch_lengths     vector containing branch lengths
 * @return                    root of the decompressed tree
 */
pll_unode_t * simple_uncompression(sdsl::bit_vector &succinct_structure, sdsl::int_vector<> &node_permutation, std::vector<double> branch_lengths);

/**
 * Decompresses a tree stored with rf distance compression.
 * @param  predecessor_tree       the predecessor tree (compressed tree is stored
 *                                relative to this tree)
 * @param  edges_to_contract      vector containing the edges to contract
 * @param  subtrees_succinct      vector containing the succinct subtrees
 * @param  succinct_permutations  vector containing the permutations
 * @param  consensus_branches     vector containing the consensus branch lengths (diffs)
 * @param  non_consensus_branches vector containing the non consensus branch lengths
 * @return                        root of the decompressed tree
 */
pll_unode_t * rf_distance_uncompression(const pll_unode_t * predecessor_tree, sdsl::int_vector<> &edges_to_contract,
          sdsl::bit_vector &subtrees_succinct, sdsl::int_vector<> &succinct_permutations,
          std::vector<double> consensus_branches, std::vector<double> non_consensus_branches);
