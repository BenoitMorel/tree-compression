#include <libpll/pll_tree.h>
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <fstream>

#include <assert.h>

/**
 * Returns the given tree in newick format
 * @param  tree tree
 * @return      string containing newick format
 */
std::string toNewick(pll_unode_t * tree);

/**
 * Prints out a given pll_unode_t on the console.
 * @param node the node
 */
void printNode(pll_unode_t * node);

/**
 * Searches for the root of the tree, i.e. the node with the label "1" and returns it.
 * @param  tree the tree
 * @return      leaf with label "1"
 */
pll_unode_t * searchRoot(pll_utree_t * tree);

/**
 * Prints out a given tree on the console.
 * @param tree root of the tree
 */
void printTree(pll_unode_t * tree);

/**
 * Prints out a given tree descending from the given node on the console.
 * @param tree root of the tree
 */
void printTreeDescending(pll_unode_t * tree);

/**
 * Traverses the given tree, storing for each node in the (void*) data
 * field the number of the smallest leaf underneath this node.
 * After this method has been called, the tree can be traversed such a way
 * that the leafs can be visited in ascending order.
 * @param tree root of the tree
 */
void setTree(pll_unode_t * tree);

/**
 * Takes a binary tree represented by its smallest leaf and orderes
 * the tree such that when using depth-first search, always the sub-
 * tree containing the smallest leaf label is visited first.
 *
 * Precondition: setTree must have been called on the tree!
 * @param tree root of the tree
 */
 void orderTree(pll_unode_t * tree);

/**
 * Traverses the given tree in depth-first order.
 * setTree had to be called on the tree before.
 *
 * Fills an array that maps node_id's to branch numbers (according to later
 * position in the succinct balanced parantheses representation)
 * @param tree                 root of the tree
 * @param bp                   vector to store succinct representation (balaced parantheses, ( = 0, ) = 1)
 * @param iv                   vector to store labels of leafs in depth-first order
 * @param branch_lengths       vector to store all branch_lengths in dfs
 * @param node_id_to_branch_id mapping of node_id to branch_id
 */
void assignBranchNumbers(pll_unode_t * tree, sdsl::bit_vector &bp, sdsl::int_vector<> &iv, std::vector<double> &branch_lengths, unsigned int* node_id_to_branch_id);

/**
 * Contracts the edge between node and node->back.
 * @param node node with incident edge
 */
void contractEdge(pll_unode_t * node);

// TODO
void pll_utree_destroy_consensus(pll_utree_t * tree);

/**
 * Traverses the given tree (represented by its root) and applies
 * the function leaf_func to the leaves.
 * @param root      root of the tree
 * @param leaf_func function to apply to the leaves
 */
void traverseTree(pll_unode_t * root, void (*leaf_func)(pll_unode_t *),
            void (*inner_node_func)(pll_unode_t *));

/**
 * Traverses the consensus tree, searches for nodes with outdegree > 2 and appends
 * the order of the children to a vector.
 * @param tree  root of the consensus tree
 * @param perms vector to store permutations of children of nodes with outdegree > 2
 */
void traverseConsensus(pll_unode_t * tree, std::vector<std::vector<int>> &perms);

/**
 * Traverses the consensus tree, searches for nodes with outdegree > 2 and appends
 * the order of the children and the "root node" to a vector.
 * @param tree          root of the consensus tree
 * @param subtree_roots vector to store "root nodes" of nodes with outdegree > 2
 * @param perms         vector to store permutations of children of nodes with outdegree > 2
 */
void traverseConsensus(pll_unode_t * tree, std::vector<pll_unode_t *> &subtree_roots,
                      std::vector<std::vector<pll_unode_t *>> &perms);

/**
 * Check if vector a is a permutation of vector b, i.e. check if both vectors
 * contain the same elements but in a different order.
 * @param  a first vector
 * @param  b second vector
 * @return true iff a is permutation of b
 */
bool isPermutation(const std::vector<int> &a, const std::vector<int> &b);

/**
 * Traverses a given tree and returns a vector of all nodes with outdegree > 2,
 * i.e. all nodes that make the tree non-binary, in DFS order.
 * @param tree tree to traverse
 * @return vector with all non binary nodes
 */
std::vector<pll_unode_t *> getNonBinaryNodesDFS(pll_unode_t * tree);

/**
 * Checks whether the two trees are equal (same topology, same labels at the leafs
 * ans same branch lengths).
 * @param  tree1 first tree
 * @param  tree2 second tree
 * @return       true iff the trees are equal
 */
bool treesEqual(pll_unode_t * tree1, pll_unode_t * tree2);

/**
 * Returns the internal predecessor of the given node.
 * For binary trees: node->next->next
 * @param  node the node
 * @return      internal predecessor
 */
pll_unode_t * internalPredecessor(pll_unode_t * node);
