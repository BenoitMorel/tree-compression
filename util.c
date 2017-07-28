#include <assert.h>

#include "util.h"


void printNode(pll_unode_t * node) {
  assert(node != NULL);
  printf("Index: %i\t\tLabel: %s\t\tLength: %f\n", node->node_index, node->label, node->length);
}

pll_unode_t * searchRoot(pll_utree_t * tree) {
  //printf("Tip_count: %i\n", tree->tip_count);
  //printf("Inner_count: %i\n", tree->inner_count);
  int nodes = (tree->tip_count + tree->inner_count);
  /*for (size_t i = 0; i < nodes; i++) {
    printf("%i\t", i);
    printNode(tree->nodes[i]);
  }*/

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
    printf("Node: ");
    printNode(tree);
  } else {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next != NULL);
    assert(tree->next->back != NULL);
    assert(tree->next->next->back != NULL);
    printNode(tree);
    printNode(tree->next);
    printNode(tree->next->next);
    printTreeRec(tree->next->back);
    printTreeRec(tree->next->next->back);
  }
}

void printTree(pll_unode_t * tree) {
  assert(tree->next == NULL);

  printNode(tree);
  printTreeRec(tree->back);

}

int setTreeRec(pll_unode_t * tree) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
    int n = atoi(tree->label);
    tree->data = (void*) n;
    tree->back->data = (void*) n;
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
    tree->data = (void*) n;
    tree->back->data = (void*) n;
    return n;
  }
}

void setTree(pll_unode_t * tree) {
  assert(tree->next == NULL);
  setTreeRec(tree->back);
  //printf("\n\nSmallest: %i\n\n", setTreeRec(tree->back));

}

void assignBranchLengthsRec(pll_unode_t * tree, unsigned int * n, unsigned int* node_id_to_branch_id) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
    //printf("%i\t", *n);
    node_id_to_branch_id[tree->node_index] = *n;
    node_id_to_branch_id[tree->back->node_index] = *n;
    (*n)++;
    //printf("Node: ");
    //printNode(tree);
  } else {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next != NULL);
    assert(tree->next->back != NULL);
    assert(tree->next->next->back != NULL);
    //printf("%i\t", *n);
    node_id_to_branch_id[tree->node_index] = *n;
    node_id_to_branch_id[tree->back->node_index] = *n;
    (*n)++;
    //printNode(tree);
    int n1 = (int) tree->next->data;
    int n2 = (int) tree->next->next->data;
    if(n1 < n2) {
      //printNode(tree->next);
      assignBranchLengthsRec(tree->next->back, n, node_id_to_branch_id);
      //printNode(tree->next->next);
      assignBranchLengthsRec(tree->next->next->back, n, node_id_to_branch_id);
    } else {
      //printNode(tree->next->next);
      assignBranchLengthsRec(tree->next->next->back, n, node_id_to_branch_id);
      //printNode(tree->next);
      assignBranchLengthsRec(tree->next->back, n, node_id_to_branch_id);
    }

  }
}

void assignBranchLengths(pll_unode_t * tree, unsigned int* node_id_to_branch_id) {
  assert(tree->next == NULL);
  //printf("1\t");
  //printNode(tree);
  unsigned int n = 2;
  assignBranchLengthsRec(tree->back, &n, node_id_to_branch_id);
}
