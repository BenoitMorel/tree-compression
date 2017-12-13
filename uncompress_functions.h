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

void simple_uncompression(sdsl::bit_vector &succinct_structure, sdsl::int_vector<> &node_permutation);

void rf_distance_uncompression(pll_unode_t * predecessor_tree, sdsl::int_vector<> &edges_to_contract, 
          sdsl::bit_vector &subtrees_succinct, sdsl::int_vector<> &succinct_permutations);
