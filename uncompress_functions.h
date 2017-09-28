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
