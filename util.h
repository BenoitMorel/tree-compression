#include <libpll/pll_tree.h>
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <fstream>

#include <assert.h>

/*
 * Prints out a given pll_unode_t on the console.
 */
void printNode(pll_unode_t * node);

/*
 * Searches the root (the node with the label "1") in the given tree
 * and returns a pointer to the root node.
 */
pll_unode_t * searchRoot(pll_utree_t * tree);

/*
 * Prints out a given tree on the console.
 */
void printTree(pll_unode_t * tree);

/*
 * Traverses the given tree, storing for each node in the (void*) data
 * field the number of the smallest leaf underneath this node.
 * After this method has been called, the tree can be traversed such a way
 * that the leafs can be visited in ascending order.
 */
void setTree(pll_unode_t * tree);

/*
 *  Takes a binary tree represented by its smallest leaf and orderes
 *  the tree such that when using depth-first search, always the sub-
 *  tree containing the smallest leaf label is visited first.
 *
 *  Precondition: setTree must have been called on the tree!
 */
 void orderTree(pll_unode_t * tree);

/*
 * Traverses the given tree ordered. setTree had to be called before.
 * Fill an array that maps node_id's to branch numbers (according to later
 * position in the succinct balanced parantheses representation)
 */
void assignBranchNumbers(pll_unode_t * tree, sdsl::bit_vector &bp, sdsl::int_vector<> &iv, std::vector<double> &branch_lengths, unsigned int* node_id_to_branch_id);

/*
 * Contracts the edge between node and node->back.
 */
void contractEdge(pll_unode_t * node);

/**
 * Stores a given array of doubles to a file.
 *
 * @param  pdata     the array of doubles
 * @param  length    length of the array
 * @param  file_path file to store the array in
 * @return           true in case of success, false otherwise
 */
bool saveArray( const double* pdata, size_t length, const std::string& file_path);

/**
 * Loads an array of doubles from file
 * @param  pdata     array to store the doubles in
 * @param  length    length of the array
 * @param  file_path file to read the doubles from
 * @return           true in case of success, false otherwise
 */
bool loadArray( double* pdata, size_t length, const std::string& file_path);

void pll_utree_destroy_consensus(pll_utree_t * tree);

/**
 * Traverses the given tree (represented by its root) in postorder and applies
 * the function leaf_func to the leaves and the function inner_node_func to the
 * inner nodes, after the nodes have been visited.
 *
 * @param root            root of the tree
 * @param leaf_func       function to apply to the leaves
 * @param inner_node_func function to apply to the inner nodes
 */
void traverseTree(pll_unode_t * root, void (*leaf_func)(pll_unode_t *),
            void (*inner_node_func)(pll_unode_t *));

uint64_t enc(double x, size_t precision);

double dec(uint64_t z, size_t precision);

/**
 * Traverses the consensus tree, searches for nodes with outdegree > 2 and appends
 * the order of the children to a vector.
 * @param tree  root of the consensus tree
 * @param perms vector to store permutations of children of nodes with outdegree > 2
 */
void traverseConsensus(pll_unode_t * tree, std::vector<std::vector<int>> &perms);
