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

pll_unode_t * simple_uncompression(sdsl::bit_vector &succinct_structure, sdsl::int_vector<> &node_permutation, std::vector<double> branch_lengths);

pll_unode_t * rf_distance_uncompression(const pll_unode_t * predecessor_tree, sdsl::int_vector<> &edges_to_contract,
          sdsl::bit_vector &subtrees_succinct, sdsl::int_vector<> &succinct_permutations,
          std::vector<double> consensus_branches, std::vector<double> non_consensus_branches);
