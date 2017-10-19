#include <assert.h>
#include <stdarg.h>
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <libpll/pll_tree.h>

#include "modified_library_functions.h"
#ifdef __cplusplus
}
#endif

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <algorithm>

#include "util.h"



#include "uncompress_functions.h"

/**
 * Takes a tree file and computes a simple compression of the tree.
 *
 * @param tree_file the tree file
 */
void simple_compression(char * tree_file);

/**
 * Takes two tree files and computes a compression between the two
 * trees using the rf distance.
 *
 * @param tree1_file first tree
 * @param tree2_file second tree
 */
void rf_distance_compression(char * tree1_file, char * tree2_file);
