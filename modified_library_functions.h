/*
 This header contains modified methods from the library for the purpose of
 working with the rf distance in different ways.
 */

#include <libpll/pll_tree.h>

// this method is just copied because "compare_splits" uses it internally and
// it is not part of pll_tree.h
unsigned int bitv_length(unsigned int bit_count);

// original method from the library.
// return the amount of same splits.
static int compare_splits (pll_split_t s1,
                           pll_split_t s2,
                           unsigned int split_len);

// original method: pllmod_utree_rf_distance(pll_unode_t * t1,
//                                                 pll_unode_t * t2,
//                                                 unsigned int tip_count)
// this method is additionally given two int arrays that will tell which
// splits were common in s1 and s2 (set to 1 if in common; set to 0 otherwise)
unsigned int pllmod_utree_split_rf_distance_extended(pll_split_t * s1,
                                                       pll_split_t * s2,
                                                       int * s1_present,
                                                       int * s2_present,
                                                       unsigned int tip_count);

// original method: pllmod_utree_rf_distance(pll_unode_t * t1,
//                                                 pll_unode_t * t2,
//                                                 unsigned int tip_count)
// this method is additionally given two splits_to_node arrays that allow the method
// to compare the branch lengths of the branches that are associated with equal splits.
// The method returns the amount of branches in the tree having the same length
// PROBLEM: trivial bipartitions are not part of s1 and s1, but have a branch length as well!
unsigned int same_branch_lengths(pll_split_t * s1,
                                                       pll_split_t * s2,
                                                       pll_unode_t ** splits_to_node1,
                                                       pll_unode_t ** splits_to_node2,
                                                       unsigned int tip_count);
