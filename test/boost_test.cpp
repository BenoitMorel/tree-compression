#define BOOST_TEST_MODULE test module name
#include <boost/test/included/unit_test.hpp>

#include "../compress_functions.cpp"
#include "../compress_functions.h"
#include "../uncompress_functions.cpp"
#include "../uncompress_functions.h"
#include "../util.cpp"
#include "../util.h"
#include "../modified_library_functions.c"
#include "../modified_library_functions.h"

bool vectorsEqual(sdsl::bit_vector a, std::vector<int> b) {
    if(a.size() != b.size()) {
        return false;
    }
    for(size_t i = 0; i < a.size(); i++) {
        if(a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

bool vectorsEqual(sdsl::int_vector<> a, std::vector<int> b) {
    if(a.size() != b.size()) {
        return false;
    }
    for(size_t i = 0; i < a.size(); i++) {
        if(a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

BOOST_AUTO_TEST_SUITE(simple_compression_test_suite)

BOOST_AUTO_TEST_CASE(simple_compression_easy_tree1) {
    simple_compression("test_inputs/easy_tree1.nwk", "simple_comp_succinct_tree.sdsl", "simple_comp_node_permutation.sdsl");

    sdsl::bit_vector subtrees_succinct_loaded;
    load_from_file(subtrees_succinct_loaded, "simple_comp_succinct_tree.sdsl");
    sdsl::int_vector<> permutations_loaded;
    load_from_file(permutations_loaded, "simple_comp_node_permutation.sdsl");

    std::vector<int> structure {0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1};
    std::vector<int> permutation {1, 2, 3, 4, 5};

    BOOST_CHECK(vectorsEqual(subtrees_succinct_loaded, structure));
    BOOST_CHECK(vectorsEqual(permutations_loaded, permutation));
}

BOOST_AUTO_TEST_CASE(simple_compression_small_tree1) {
    simple_compression("test_inputs/small_tree1.nwk", "simple_comp_succinct_tree.sdsl", "simple_comp_node_permutation.sdsl");

    sdsl::bit_vector subtrees_succinct_loaded;
    load_from_file(subtrees_succinct_loaded, "simple_comp_succinct_tree.sdsl");
    sdsl::int_vector<> permutations_loaded;
    load_from_file(permutations_loaded, "simple_comp_node_permutation.sdsl");

    std::vector<int> structure {0, 0, 1, 0, 0, 1, 0, 1, 1, 1};
    std::vector<int> permutation {1, 2, 3};

    BOOST_CHECK(vectorsEqual(subtrees_succinct_loaded, structure));
    BOOST_CHECK(vectorsEqual(permutations_loaded, permutation));
}

BOOST_AUTO_TEST_CASE(simple_compression_trivial_tree1) {
    int result = simple_compression("test_inputs/trivial_tree1.nwk", "simple_comp_succinct_tree.sdsl", "simple_comp_node_permutation.sdsl");

    BOOST_CHECK(result == -1);
}

BOOST_AUTO_TEST_CASE(simple_compression_trivial_tree2) {
    int result = simple_compression("test_inputs/trivial_tree2.nwk", "simple_comp_succinct_tree.sdsl", "simple_comp_node_permutation.sdsl");

    BOOST_CHECK(result == -1);
}

BOOST_AUTO_TEST_CASE(simple_compression_invalid_tree) {
    int result = simple_compression("test_inputs/invalid_tree.nwk", "simple_comp_succinct_tree.sdsl", "simple_comp_node_permutation.sdsl");

    BOOST_CHECK(result == -1);
}


BOOST_AUTO_TEST_SUITE_END()





BOOST_AUTO_TEST_SUITE(rf_distance_compression_test_suite)

BOOST_AUTO_TEST_CASE(rf_distance_compression_easy_tree1) {

    rf_distance_compression("test_inputs/rf_easy_1_1.nwk", "test_inputs/rf_easy_1_2.nwk", "rf_edges_to_contract.sdsl", "rf_subtrees_succinct.sdsl", "rf_permutations.sdsl");

    sdsl::int_vector<> edges_to_contract_loaded;
    load_from_file(edges_to_contract_loaded, "rf_edges_to_contract.sdsl");

    sdsl::bit_vector subtrees_succinct_loaded;
    load_from_file(subtrees_succinct_loaded, "rf_subtrees_succinct.sdsl");

    sdsl::int_vector<> permutations_loaded;
    load_from_file(permutations_loaded, "rf_permutations.sdsl");

    std::vector<int> edges_to_contract {4};
    std::vector<int> subtrees_succinct {0, 1};
    std::vector<int> permutations {0, 2, 1};

    BOOST_CHECK(vectorsEqual(edges_to_contract_loaded, edges_to_contract));
    BOOST_CHECK(vectorsEqual(subtrees_succinct_loaded, subtrees_succinct));
    BOOST_CHECK(vectorsEqual(permutations_loaded, permutations));
}

BOOST_AUTO_TEST_CASE(rf_distance_compression_easy_tree2) {

    rf_distance_compression("test_inputs/rf_easy_2_1.nwk", "test_inputs/rf_easy_2_2.nwk", "rf_edges_to_contract.sdsl", "rf_subtrees_succinct.sdsl", "rf_permutations.sdsl");

    sdsl::int_vector<> edges_to_contract_loaded;
    load_from_file(edges_to_contract_loaded, "rf_edges_to_contract.sdsl");

    sdsl::bit_vector subtrees_succinct_loaded;
    load_from_file(subtrees_succinct_loaded, "rf_subtrees_succinct.sdsl");

    sdsl::int_vector<> permutations_loaded;
    load_from_file(permutations_loaded, "rf_permutations.sdsl");

    std::vector<int> edges_to_contract {4, 10};
    std::vector<int> subtrees_succinct {0, 1, 0, 1};
    std::vector<int> permutations {0, 2, 1, 0, 2, 1};

    BOOST_CHECK(vectorsEqual(edges_to_contract_loaded, edges_to_contract));
    BOOST_CHECK(vectorsEqual(subtrees_succinct_loaded, subtrees_succinct));
    BOOST_CHECK(vectorsEqual(permutations_loaded, permutations));
}

BOOST_AUTO_TEST_CASE(rf_distance_compression_easy_tree3) {

    rf_distance_compression("test_inputs/rf_easy_3_1.nwk", "test_inputs/rf_easy_3_2.nwk", "rf_edges_to_contract.sdsl", "rf_subtrees_succinct.sdsl", "rf_permutations.sdsl");

    sdsl::int_vector<> edges_to_contract_loaded;
    load_from_file(edges_to_contract_loaded, "rf_edges_to_contract.sdsl");

    sdsl::bit_vector subtrees_succinct_loaded;
    load_from_file(subtrees_succinct_loaded, "rf_subtrees_succinct.sdsl");

    sdsl::int_vector<> permutations_loaded;
    load_from_file(permutations_loaded, "rf_permutations.sdsl");

    std::vector<int> edges_to_contract {4, 5, 8};
    std::vector<int> subtrees_succinct {0, 0, 1, 1, 0, 1};
    std::vector<int> permutations {0, 3, 2, 4, 1};

    BOOST_CHECK(vectorsEqual(edges_to_contract_loaded, edges_to_contract));
    BOOST_CHECK(vectorsEqual(subtrees_succinct_loaded, subtrees_succinct));
    BOOST_CHECK(vectorsEqual(permutations_loaded, permutations));
}

BOOST_AUTO_TEST_SUITE_END()
