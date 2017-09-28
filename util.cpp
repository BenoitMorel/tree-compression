#include <assert.h>

#include "util.h"


void printNode(pll_unode_t * node) {
  assert(node != NULL);
  printf("Index: %i\t\tLabel: %s\t\tLength: %f\n", node->node_index, node->label, node->length);
}

pll_unode_t * searchRoot(pll_utree_t * tree) {

  for (size_t i = 0; i < tree->tip_count; i++) {
    if(strcmp(tree->nodes[i]->label, "1") == 0) {
      //printf("\nROOT: ");
      //printNode(tree->nodes[i]);
      return tree->nodes[i];
    }
  }
  return NULL;
}

void printTreeRec(pll_unode_t * tree) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
    printNode(tree);
  } else {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next != NULL);
    assert(tree->next->back != NULL);
    assert(tree->next->next->back != NULL);
    printNode(tree);
    printNode(tree->next);
    printTreeRec(tree->next->back);
    printNode(tree->next->next);
    printTreeRec(tree->next->next->back);
  }
}

void printTree(pll_unode_t * tree) {
  assert(tree->back != NULL);

  printNode(tree);
  printTreeRec(tree->back);

}

int setTreeRec(pll_unode_t * tree) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
    int n = atoi(tree->label);
    tree->data = (void*)(intptr_t) n;
    tree->back->data = (void*)(intptr_t) n;
    return n;
  } else {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next != NULL);
    assert(tree->next->back != NULL);
    assert(tree->next->next->back != NULL);
    int n1 = setTreeRec(tree->next->back);
    int n2 = setTreeRec(tree->next->next->back);
    int n;
    n1 < n2 ? (n = n1) : (n = n2);
    tree->data = (void*)(intptr_t) n;
    tree->back->data = (void*)(intptr_t) n;
    return n;
  }
}

void setTree(pll_unode_t * tree) {
  assert(tree->next == NULL);
  setTreeRec(tree->back);
  //printf("\n\nSmallest: %i\n\n", setTreeRec(tree->back));

}

void orderTreeRec(pll_unode_t * tree) {
  assert(tree != NULL);
  if(tree->next != NULL) {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next->next == tree); // tree is binary

    if(((intptr_t) tree->next->data) > ((intptr_t) tree->next->next->data)) {
        // swap tree->next and tree->next->next
        pll_unode_t * temp = tree->next->next;
        tree->next->next = tree;
        temp->next = tree->next;
        tree->next = temp;
    }

    orderTreeRec(tree->next->back);
    orderTreeRec(tree->next->next->back);
  }
}

/*
 *  Takes a binary tree represented by its smallest leaf and orderes
 *  the tree such that when using depth-first search, always the sub-
 *  tree containing the smallest leaf label is visited first.
 *
 *  Precondition: setTree must have been called on the tree!
 */
void orderTree(pll_unode_t * tree) {
  assert(tree->next == NULL);
  assert(atoi(tree->label) == 1);

  orderTreeRec(tree->back);
}

bool innerNodeCompare(pll_unode_t * node1, pll_unode_t * node2) {
  return ((intptr_t) node1->data) < ((intptr_t) node2->data);
}

void assignBranchNumbersRec(pll_unode_t * tree, unsigned int * bp_idx, sdsl::bit_vector &bp,
              unsigned int * iv_idx, sdsl::int_vector<> &iv, unsigned int * bl_idx, std::vector<double> &branch_lengths,
              unsigned int * n, unsigned int* node_id_to_branch_id) {
  assert(tree != NULL);
  branch_lengths[*bl_idx] = tree->length;
  (*bl_idx)++;
  if(tree->next == NULL) {
    // leaf

    //printf("%i\t", *n);
    node_id_to_branch_id[tree->node_index] = *n;
    node_id_to_branch_id[tree->back->node_index] = *n;
    (*n)++;
    iv[*iv_idx] = atoi(tree->label);
    (*iv_idx)++;
    //printf("%s, ", tree->label);
    //printf("Node: ");
    //printNode(tree);
  } else {
    // inner node
    assert(tree->next != NULL);

    // assign inner nodes to an array
    std::vector<pll_unode_t*> inner_nodes;
    pll_unode_t * temp_node = tree->next;
    while(temp_node->node_index != tree->node_index) {
      inner_nodes.push_back(temp_node);
      assert(temp_node->next != NULL);
      temp_node = temp_node->next;
    }
    assert(temp_node->node_index == tree->node_index);

    // sort the internal nodes
    std::sort (inner_nodes.begin(), inner_nodes.end(), innerNodeCompare);

    /*for (size_t i = 0; i < inner_nodes.size(); i++) {
      std::cout << (intptr_t) inner_nodes[i]->data << ", ";
    }
    std::cout << "\n";*/

    //printf("%i\t", *n);
    node_id_to_branch_id[tree->node_index] = *n;
    node_id_to_branch_id[tree->back->node_index] = *n;
    (*n)++;
    //printNode(tree);

    for (size_t i = 0; i < inner_nodes.size(); i++) {
      bp[*bp_idx] = 0;
      (*bp_idx)++;
      assignBranchNumbersRec(inner_nodes[i]->back, bp_idx, bp, iv_idx, iv, bl_idx, branch_lengths, n, node_id_to_branch_id);
      bp[*bp_idx] = 1;
      (*bp_idx)++;
    }


    /*int n1 = (intptr_t) tree->next->data;
    int n2 = (intptr_t) tree->next->next->data;
    if(n1 < n2) {
      //printNode(tree->next);
      bp[*bp_idx] = 0;
      (*bp_idx)++;
      assignBranchNumbersRec(tree->next->back, bp_idx, bp, iv_idx, iv, bl_idx, branch_lengths, n, node_id_to_branch_id);
      bp[*bp_idx] = 1;
      (*bp_idx)++;
      //printNode(tree->next->next);
      bp[*bp_idx] = 0;
      (*bp_idx)++;
      assignBranchNumbersRec(tree->next->next->back, bp_idx, bp, iv_idx, iv, bl_idx, branch_lengths, n, node_id_to_branch_id);
      bp[*bp_idx] = 1;
      (*bp_idx)++;
    } else {
      //printNode(tree->next->next);
      bp[*bp_idx] = 0;
      (*bp_idx)++;
      assignBranchNumbersRec(tree->next->next->back, bp_idx, bp, iv_idx, iv, bl_idx, branch_lengths, n, node_id_to_branch_id);
      bp[*bp_idx] = 1;
      (*bp_idx)++;
      //printNode(tree->next);
      bp[*bp_idx] = 0;
      (*bp_idx)++;
      assignBranchNumbersRec(tree->next->back, bp_idx, bp, iv_idx, iv, bl_idx, branch_lengths, n, node_id_to_branch_id);
      bp[*bp_idx] = 1;
      (*bp_idx)++;
    }*/

  }
}

void assignBranchNumbers(pll_unode_t * tree, sdsl::bit_vector &bp, sdsl::int_vector<> &iv,
                std::vector<double> &branch_lengths, unsigned int* node_id_to_branch_id) {
  assert(tree->next == NULL);
  assert(atoi(tree->label) == 1);
  bp[0] = 0;
  bp[1] = 0;
  bp[2] = 1;
  iv[0] = 1; // first node is always the root
  //printNode(tree);
  unsigned int n = 2;
  bp[3] = 0;
  unsigned int bp_idx = 4;
  unsigned int iv_idx = 1;
  unsigned int bl_idx = 1;
  assignBranchNumbersRec(tree->back, &bp_idx, bp, &iv_idx, iv, &bl_idx, branch_lengths, &n, node_id_to_branch_id);
  bp[bp_idx] = 1;
  bp_idx++;
  bp[bp_idx] = 1;
  bp_idx++;
  assert(bp_idx == bp.size());
  assert(iv_idx == iv.size());
}

/*
 * Searches the internal predecessor p of a node and returns a pointer to it.
 * The node p fulfills p->next = node.
 */
pll_unode_t * internalPredecessor(pll_unode_t * node) {
  pll_unode_t * predecessor = node;
  while(predecessor->next->node_index != node->node_index) {
    predecessor = predecessor->next;
  }
  assert(predecessor->next->node_index == node->node_index);
  return predecessor;
}

void contractEdge(pll_unode_t * node) {
  assert(node != NULL);
  assert(node->back != NULL);
  pll_unode_t * node_back_predecessor = internalPredecessor(node->back);
  pll_unode_t * node_back_successor = node->back->next;
  pll_unode_t * node_predecessor = internalPredecessor(node);
  pll_unode_t * node_successor = node->next;
  node_predecessor->next = node_back_successor;
  node_back_predecessor->next = node_successor;
}
