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
#include <sdsl/wavelet_trees.hpp>
#include <algorithm>

#include "util.h"

#include "uncompress_functions.h"

enum Flags{
    // print out size that is needed to store the compression
    PRINT_COMPRESSION              = 0x01,

    // print out all the compression structures used
    PRINT_COMPRESSION_STRUCTURES   = 0x02,

    // just for rf-distance compression: print out splits
    PRINT_SPLITS                   = 0x04
};

/**
 * Takes a tree file and computes a simple compression of the tree.
 *
 * @param  tree_file               tree in newick format to compress
 * @param  succinct_structure_file file to store succinct structure
 * @param  node_permutation_file   file to store node permutation
 * @param  flags                   flags
 * @return                         value < 0 in case of an eŕror
 */
int simple_compression(const char * tree_file, const char * succinct_structure_file,
        const char * node_permutation_file, int flags);

/**
 * Takes two tree files and computes a compression between the two
 * trees using the rf distance.
 *
 * @param tree1_file             first tree in newick format
 * @param tree2_file             second tree in newick format
 * @param edges_to_contract_file file to store edges to contract
 * @param subtrees_succinct_file file to store subtrees succinct
 * @param node_permutations_file file to store node permutations
 * @param flags                  flags
 * @return                       value < 0 in case of an eŕror
 */
int rf_distance_compression(const char * tree1_file, const char * tree2_file,
        const char * edges_to_contract_file, const char * subtrees_succinct_file,
        const char * node_permutations_file, int flags);
