#include "datastructure_compression_functions.h"

size_t writeDoubleVectorUncompressed(const std::vector<double>& myVector, std::string filename) {
  std::ofstream ofs(filename, std::ios::out | std::ofstream::binary);
  std::ostream_iterator<char> osi{ ofs };
  const char* beginByte = (char*)&myVector[0];

  const char* endByte = (char*)&myVector.back() + sizeof(double);
  std::copy(beginByte, endByte, osi);

  return myVector.size() * 8;
}

std::vector<double> readDoubleVectorUncompressed(std::string filename) {
  std::vector<char> buffer{};
  std::ifstream ifs(filename, std::ios::in | std::ifstream::binary);
  std::istreambuf_iterator<char> iter(ifs);
  std::istreambuf_iterator<char> end{};
  std::copy(iter, end, std::back_inserter(buffer));
  std::vector<double> newVector(buffer.size() / sizeof(double));
  memcpy(&newVector[0], &buffer[0], buffer.size());
  return newVector;
}

uint64_t enc(double x, size_t precision) {
    double expo = 1;
    for (size_t i=0; i<precision; ++i) {
      expo *= 10;
    }

    int64_t y;
    uint64_t z;
    y = llround(x*expo);
    if ( y < 0 ) {
        z = 2*static_cast<uint64_t>(-y)-1;
    } else {
        z = 2*static_cast<uint64_t>(y);
    }
    return z;
}

double dec(uint64_t z, size_t precision) {
    double expo = 1;
    for (size_t i=0; i<precision; ++i) {
      expo *= 10;
    }

    int64_t y;
    if ( z % 2 ) {
        y = -static_cast<int64_t>((z+1)/2);
    } else {
        y = static_cast<int64_t>(z/2);
    }
    return (static_cast<double>(y))/expo;
}

size_t compressBranchLengthsAndStore(std::vector<double> branch_lengths, std::string filename) {
  std::vector<uint64_t> vec;

  uint64_t z;
  for (size_t i = 0; i < branch_lengths.size(); i++) {
    z = enc(branch_lengths[i], PRECISION);
    vec.push_back(z);
  }

  uint64_t max_number = *max_element(vec.begin(), vec.end());
  uint32_t width = sdsl::bits::hi(max_number);

  sdsl::int_vector<0> seq(branch_lengths.size(), 0, width+1);

  for (size_t i=0; i<vec.size(); ++i) {
       seq[i] = vec[i];
  }

  sdsl::wt_int<sdsl::rrr_vector<63>> wt;
  sdsl::construct_im(wt, seq);

  sdsl::store_to_file(wt, filename);

  return sdsl::size_in_bytes(wt);
}

std::vector<double> uncompressBranchLengthsWV(std::string filename) {
  sdsl::wt_int<sdsl::rrr_vector<63>> loaded_wt;
  sdsl::load_from_file(loaded_wt, filename);

  std::vector<double> loaded(loaded_wt.size());
  for (size_t i = 0; i < loaded_wt.size(); i++) {
    loaded[i] = dec(loaded_wt[i], PRECISION);
  }

  return loaded;
}


size_t compressAndStoreSuccinctStructure(sdsl::bit_vector &succinct_structure, std::string filename) {

  auto size = sdsl::size_in_bytes(succinct_structure);
  store_to_file(succinct_structure, filename);

  return size;
}

size_t compressAndStoreSimplePermutation(sdsl::int_vector<> &permutation, std::string filename) {
    // TODO: compress

    sdsl::util::bit_compress(permutation);
    auto size = sdsl::size_in_bytes(permutation);

    store_to_file(permutation, filename);

    return size;
}

size_t compressAndStoreRFEdgesToContract(sdsl::int_vector<> &edges_to_contract, std::string filename) {
    std::vector<uint64_t> vec;
    vec.push_back(edges_to_contract[0]);
    for (size_t i = 0; i < edges_to_contract.size()-1; i++) {
        vec.push_back(edges_to_contract[i+1]-edges_to_contract[i]);
    }

    uint64_t max_number = *max_element(vec.begin(), vec.end());
    uint32_t width = sdsl::bits::hi(max_number);

    sdsl::int_vector<0> seq(edges_to_contract.size(), 0, width+1);

    for (size_t i=0; i<vec.size(); ++i) {
       seq[i] = vec[i];
    }

    sdsl::util::bit_compress(seq);

    auto size = sdsl::size_in_bytes(seq);    
    store_to_file(seq, filename);

    return size;
}

size_t compressAndStoreRFSubtreePermutations(sdsl::int_vector<> &subtree_permutations, std::string filename) {
    sdsl::util::bit_compress(subtree_permutations);
    auto size = sdsl::size_in_bytes(subtree_permutations);

    store_to_file(subtree_permutations, filename);

    return size;
}

size_t compressAndStoreBranchLengths(std::vector<double> branch_lengths, std::string filename) {
    return compressBranchLengthsAndStore(branch_lengths, filename);
}



sdsl::bit_vector uncompressSuccinctStructure(std::string filename) {
    sdsl::bit_vector succinct_tree_loaded;
    load_from_file(succinct_tree_loaded, filename);
    return succinct_tree_loaded;
}

sdsl::int_vector<> uncompressSimplePermutation(std::string filename) {
    sdsl::int_vector<> node_permutation_loaded;
    load_from_file(node_permutation_loaded, filename);
    return node_permutation_loaded;
}

sdsl::int_vector<> uncompressRFEdgesToContract(std::string filename) {
    sdsl::int_vector<> edges_to_contract_loaded;
    load_from_file(edges_to_contract_loaded, filename);

    std::vector<uint64_t> vec;
    vec.push_back(edges_to_contract_loaded[0]);
    for (size_t i = 1; i < edges_to_contract_loaded.size(); i++) {
        vec.push_back(vec[i-1] + edges_to_contract_loaded[i]);
    }

    uint64_t max_number = *max_element(vec.begin(), vec.end());
    uint32_t width = sdsl::bits::hi(max_number);

    sdsl::int_vector<0> seq(edges_to_contract_loaded.size(), 0, width+1);

    for (size_t i=0; i<vec.size(); ++i) {
       seq[i] = vec[i];
    }

    return seq;
}

sdsl::int_vector<> uncompressRFSubtreePermutations(std::string filename) {
    sdsl::int_vector<> permutations_loaded;
    load_from_file(permutations_loaded, filename);
    return permutations_loaded;
}

std::vector<double> uncompressBranchLengths(std::string filename) {
    return uncompressBranchLengthsWV(filename);
}
