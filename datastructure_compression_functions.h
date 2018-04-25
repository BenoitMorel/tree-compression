#include <assert.h>
#include <stdarg.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <algorithm>

size_t writeDoubleVector(const std::vector<double>& myVector, std::string filename);

std::vector<double> readDoubleVector(std::string filename);

size_t compressBranchLengthsAndStore(std::vector<double> branch_lengths, const char * file);

std::vector<double> uncompressBranchLengths(const char * file);

bool storeBranchLengthsUncompressed(std::vector<double> branch_lengths, const std::string& file_path);

std::vector<double> loadBranchLengthsUncompressed(const std::string& file_path);
