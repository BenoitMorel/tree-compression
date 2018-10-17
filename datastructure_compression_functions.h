#include <assert.h>
#include <stdarg.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <algorithm>

/**
 * This class contains methods to compress and store the individual data structures
 * as well as methods to decompress them.
 */


size_t compressAndStoreSuccinctStructure(sdsl::bit_vector &succinct_structure, std::string filename);

size_t compressAndStoreSimplePermutation(sdsl::int_vector<> &permutation, std::string filename);

size_t compressAndStoreRFEdgesToContract(sdsl::int_vector<> &edges_to_contract, std::string filename);

size_t compressAndStoreRFSubtreePermutations(sdsl::int_vector<> &edges_to_contract, std::string filename);

size_t compressAndStoreBranchLengths(std::vector<double> branch_lengths, std::string filename);


sdsl::bit_vector uncompressSuccinctStructure(std::string filename);

sdsl::int_vector<> uncompressSimplePermutation(std::string filename);

sdsl::int_vector<> uncompressRFEdgesToContract(std::string filename);

sdsl::int_vector<> uncompressRFSubtreePermutations(std::string filename);

std::vector<double> uncompressBranchLengths(std::string filename);
