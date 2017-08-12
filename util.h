#include <libpll/pll_tree.h>
#include <sdsl/bit_vectors.hpp>

#include <assert.h>

/*
 * Prints out a given pll_unode_t on the console.
 */
void printNode(pll_unode_t * node);

/*
 * Searches the root (the node with the label "1") in the given tree
 * and return a pointer to the root node.
 */
pll_unode_t * searchRoot(pll_utree_t * tree);

/*
 * Prints out a given tree on the console.
 */
void printTree(pll_unode_t * tree);

/*
 * Traverses the given tree, storing for each node in the (void*) data
 * field the number of the smallest leaf underneath this node.
 * After this method has been performed, the tree can be traversed such
 * that the leafs can be visited in ascending order.
 */
void setTree(pll_unode_t * tree);

/*
 * Traverses the given tree ordered. setTree had to be called before.
 * Fill an array that maps node_id's to branch numbers (according to later
 * position in the succinct balanced parantheses representation)
 */
void assignBranchNumbers(pll_unode_t * tree, sdsl::bit_vector &bp, unsigned int* node_id_to_branch_id);
