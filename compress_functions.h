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
#include <sdsl/int_vector.hpp>
#include <sdsl/vectors.hpp>

#include "util.h"

void simple_compression(char * tree_file);

void rf_distance_compression(char * tree1_file, char * tree2_file);
