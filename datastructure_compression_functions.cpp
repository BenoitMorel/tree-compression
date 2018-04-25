#include "datastructure_compression_functions.h"

size_t writeDoubleVectorSimple(const std::vector<double>& myVector, std::string filename)
{
  std::ofstream ofs(filename, std::ios::out | std::ofstream::binary);
  std::ostream_iterator<char> osi{ ofs };
  const char* beginByte = (char*)&myVector[0];

  const char* endByte = (char*)&myVector.back() + sizeof(double);
  std::copy(beginByte, endByte, osi);

  return myVector.size() * 8;
}

std::vector<double> readDoubleVectorSimple(std::string filename)
{
  std::vector<char> buffer{};
  std::ifstream ifs(filename, std::ios::in | std::ifstream::binary);
  std::istreambuf_iterator<char> iter(ifs);
  std::istreambuf_iterator<char> end{};
  std::copy(iter, end, std::back_inserter(buffer));
  std::vector<double> newVector(buffer.size() / sizeof(double));
  memcpy(&newVector[0], &buffer[0], buffer.size());
  return newVector;
}

size_t writeDoubleVector(const std::vector<double>& myVector, std::string filename) {
    return writeDoubleVectorSimple(myVector, filename);
    //return compressBranchLengthsAndStore(myVector, filename.c_str());
}

std::vector<double> readDoubleVector(std::string filename) {
    return readDoubleVectorSimple(filename);
    //return uncompressBranchLengths(filename.c_str());
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

size_t compressBranchLengthsAndStore(std::vector<double> branch_lengths, const char * file) {
  std::vector<uint64_t> vec;

  // std::cout << "stored: ";
  // for (size_t i = 0; i < branch_lengths.size(); i++) {
  //   std::cout << branch_lengths[i] << "  ";
  // }
  // std::cout << "\n";

  uint64_t z;
  for (size_t i = 0; i < branch_lengths.size(); i++) {
    z = enc(branch_lengths[i], 9);
    vec.push_back(z);
  }

  // std::cout << "vec: ";
  // for (size_t i = 0; i < vec.size(); i++) {
  //   std::cout << vec[i] << "  ";
  // }
  // std::cout << "\n";

  uint64_t max_number = *max_element(vec.begin(), vec.end());
  uint32_t width = sdsl::bits::hi(max_number);
  sdsl::int_vector<0> seq(branch_lengths.size(), 0, width);
  for (size_t i=0; i<vec.size(); ++i) {
       seq[i] = vec[i];
  }
  sdsl::wt_int<sdsl::rrr_vector<63>> wt;
  sdsl::construct_im(wt, seq);
  sdsl::store_to_file(wt, file);

  // for (size_t i = 0; i < wt.size(); i++) {
  //   std::cout << wt[i] << "  ";
  // }
  // std::cout << "\n";

  return sdsl::size_in_bytes(wt);
}

std::vector<double> uncompressBranchLengths(const char * file) {
  sdsl::wt_int<sdsl::rrr_vector<63>> loaded_wt;
  sdsl::load_from_file(loaded_wt, file);

  std::vector<double> loaded(loaded_wt.size());
  for (size_t i = 0; i < loaded_wt.size(); i++) {
    loaded[i] = dec(loaded_wt[i], 9);
  }

  // std::cout << "loaded: ";
  // for (size_t i = 0; i < loaded.size(); i++) {
  //   std::cout << loaded[i] << "  ";
  // }
  // std::cout << "\n";

  return loaded;
}

bool storeBranchLengthsUncompressed(std::vector<double> branch_lengths, const std::string& file_path ) {
    std::ofstream os(file_path.c_str(), std::ios::binary | std::ios::out);
    if ( !os.is_open() )
        return false;
    os.write(reinterpret_cast<const char*>(&branch_lengths[0]), std::streamsize(branch_lengths.size()*sizeof(double)));
    os.close();
    return true;
}

std::vector<double> loadBranchLengthsUncompressed(const std::string& file_path) {
    std::vector<double> branch_lengths;
    std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
    if ( !is.is_open() )
        return branch_lengths;
    is.read(reinterpret_cast<char*>(&branch_lengths[0]), std::streamsize(branch_lengths.size()*sizeof(double)));
    is.close();
    return branch_lengths;
}
