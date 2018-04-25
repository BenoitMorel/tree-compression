#include <assert.h>

#include "util.h"

std::string toNewickRec(pll_unode_t * tree) {
  assert(tree != NULL);

  std::stringstream ss;
  if(tree->next == NULL) {
    //leaf
    ss << tree->label << ":" << tree->length;
  } else {
    assert(tree->next->next->next == tree); // tree is binary
    ss << "(" << toNewickRec(tree->next->back) << "," << toNewickRec(tree->next->next->back)
               << "):" << tree->length;
  }
  return ss.str();
}

std::string toNewick(pll_unode_t * tree) {
  assert(tree != NULL);

  if(tree->next == NULL) {
    assert(tree->back != NULL);
    assert(tree->back->next != NULL);
    assert(tree->back->next->next != NULL);

    std::stringstream ss;

    return "(" + toNewickRec(tree) + "," + toNewickRec(tree->back->next->back) + "," +
                      toNewickRec(tree->back->next->next->back) + ");";
  } else {
    assert(tree->next != NULL);
    assert(tree->next->next != NULL);

    std::stringstream ss;

    return "(" + toNewickRec(tree->back) + "," + toNewickRec(tree->next->back) + "," +
                      toNewickRec(tree->next->next->back) + ");";
  }
}

void printNode(pll_unode_t * node) {
  assert(node != NULL);
  printf("Node Index: %i, PMatrix Index: %i, Label: %s, Length: %f, Data: %i\n",
          node->node_index, node->pmatrix_index, node->label, node->length, (intptr_t) node->data);
}

pll_unode_t * searchRoot(pll_utree_t * tree) {

  for (size_t i = 0; i < tree->tip_count; i++) {
    if(strcmp(tree->nodes[i]->label, "1") == 0) {
      return tree->nodes[i];
    }
  }
  return NULL;
}

void printTreeRec(pll_unode_t * tree, std::string tabs) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
    std::cout << tabs;
    printNode(tree);
  } else {
    // inner node
    assert(tree->next != NULL);

    pll_unode_t * temp = tree->next;
    std::cout << tabs;
    printNode(tree);
    int ctr = 1;
    while(temp != tree) {
      std::cout << tabs;
      printNode(temp);
      printTreeRec(temp->back, (tabs + "\t"));
      temp = temp->next;
      assert(temp != NULL);
    }
  }
}

void printTree(pll_unode_t * tree) {

  if(tree->next == NULL) {
    printNode(tree);
    printTreeRec(tree->back, "\t");
  } else {
    printTreeRec(tree, "");
  }
}

void printTreeDescending(pll_unode_t * tree) {
    printNode(tree);
    if(tree->next != NULL) {
        printTreeRec(tree->back, "\t");
    }
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
}

void orderTreeRec(pll_unode_t * tree) {
  assert(tree != NULL);
  if(tree->next != NULL) {
    // inner node
    assert(tree->next != NULL);

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

    node_id_to_branch_id[tree->node_index] = *n;
    node_id_to_branch_id[tree->back->node_index] = *n;
    (*n)++;
    iv[*iv_idx] = atoi(tree->label);
    (*iv_idx)++;
  } else {
    // inner node

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

    node_id_to_branch_id[tree->node_index] = *n;
    node_id_to_branch_id[tree->back->node_index] = *n;
    (*n)++;

    for (size_t i = 0; i < inner_nodes.size(); i++) {
      bp[*bp_idx] = 0;
      (*bp_idx)++;
      assignBranchNumbersRec(inner_nodes[i]->back, bp_idx, bp, iv_idx, iv, bl_idx, branch_lengths, n, node_id_to_branch_id);
      bp[*bp_idx] = 1;
      (*bp_idx)++;
    }
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
  assert(bl_idx == branch_lengths.size());
}

/*
 * Searches the internal predecessor p of a node and returns a pointer to it.
 * The node p fulfills p->next = node.
 */
pll_unode_t * internalPredecessor(pll_unode_t * node) {
  pll_unode_t * predecessor = node;
  while(predecessor->next != node) {
    predecessor = predecessor->next;
  }
  assert(predecessor->next == node);
  return predecessor;
}

void contractEdge(pll_unode_t * node) {
  assert(node != NULL);
  assert(node->back != NULL);
  assert(node->next != NULL); // don't contract leafs

  pll_unode_t * node_back_predecessor = internalPredecessor(node->back);
  pll_unode_t * node_back_successor = node->back->next;
  pll_unode_t * node_predecessor = internalPredecessor(node);
  pll_unode_t * node_successor = node->next;

  node_predecessor->next = node_back_successor;
  node_back_predecessor->next = node_successor;
}

bool saveArray( const double* pdata, size_t length, const std::string& file_path ) {
    std::ofstream os(file_path.c_str(), std::ios::binary | std::ios::out);
    if ( !os.is_open() )
        return false;
    os.write(reinterpret_cast<const char*>(pdata), std::streamsize(length*sizeof(double)));
    os.close();
    return true;
}

bool loadArray( double* pdata, size_t length, const std::string& file_path) {
    std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
    if ( !is.is_open() )
        return false;
    is.read(reinterpret_cast<char*>(pdata), std::streamsize(length*sizeof(double)));
    is.close();
    return true;
}

void traverseTreeRec(pll_unode_t * tree, void (*leaf_func)(pll_unode_t *),
            void (*inner_node_func)(pll_unode_t *)) {
  assert(tree != NULL);
  assert(leaf_func != NULL);
  assert(inner_node_func != NULL);

  if(tree->next == NULL) {
    // leaf
    leaf_func(tree);
  } else {
    // inner node
    pll_unode_t * current_node = tree->next;
    while(current_node != tree) {
      pll_unode_t * temp_node = current_node;
      traverseTreeRec(current_node->back, leaf_func, inner_node_func);
      current_node = current_node->next;
      inner_node_func(temp_node);
    }
    inner_node_func(tree);
  }
}

void traverseTree(pll_unode_t * root, void (*leaf_func)(pll_unode_t *),
            void (*inner_node_func)(pll_unode_t *)) {
  assert(root->back != NULL);

  traverseTreeRec(root->back, leaf_func, inner_node_func);
}

void traverseConsensusRec(pll_unode_t * tree, std::vector<std::vector<int>> &perms) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
  } else {
    // inner node
    assert(tree->next != NULL);

    int ctr = 1;
    std::vector<int> perm;

    pll_unode_t * temp = tree->next;
    while(temp != tree) {
      traverseConsensusRec(temp->back, perms);

      perm.push_back((intptr_t) temp->back->data);

      temp = temp->next;

      assert(temp != NULL);
      ctr++;
    }

    if(ctr>3){
        perms.push_back(perm);
    }
  }
}

/**
 * Traverses the consensus tree, searches for nodes with outdegree > 2 and appends
 * the order of the children to a vector.
 * @param tree  root of the consensus tree
 * @param perms vector to store permutations of children of nodes with outdegree > 2
 */
void traverseConsensus(pll_unode_t * tree, std::vector<std::vector<int>> &perms) {
  assert(tree->next == NULL);
  assert(tree->back != NULL);

  traverseConsensusRec(tree->back, perms);
}

void traverseConsensusRec(pll_unode_t * tree, std::vector<pll_unode_t *> &subtree_roots,
                  std::vector<std::vector<pll_unode_t *>> &children) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
  } else {
    // inner node
    assert(tree->next != NULL);

    int ctr = 1;
    std::vector<pll_unode_t *> children_node;

    pll_unode_t * temp = tree->next;
    while(temp != tree) {
      traverseConsensusRec(temp->back, subtree_roots, children);

      children_node.push_back(temp);

      temp = temp->next;

      assert(temp != NULL);
      ctr++;
    }

    if(ctr>3){
        subtree_roots.push_back(tree);
        children.push_back(children_node);
    }
  }
}

void traverseConsensus(pll_unode_t * tree, std::vector<pll_unode_t *> &subtree_roots,
                      std::vector<std::vector<pll_unode_t *>> &perms) {
  assert(tree->next == NULL);
  assert(tree->back != NULL);

  traverseConsensusRec(tree->back, subtree_roots, perms);

  assert(subtree_roots.size() == perms.size());
}

void getNonBinaryNodesDFSRec(pll_unode_t * tree, std::vector<pll_unode_t *> &nodes) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf
  } else {
    // inner node
    assert(tree->next != NULL);

    int ctr = 1;

    pll_unode_t * temp = tree->next;
    while(temp != tree) {
      ctr++;
      if(ctr > 3) {
          nodes.push_back(tree);
          break;
      }
      temp = temp->next;
      assert(temp != NULL);
    }

    temp = tree->next;
    while(temp != tree) {
      getNonBinaryNodesDFSRec(temp->back, nodes);
      temp = temp->next;
      assert(temp != NULL);
    }
  }
}

/**
 * Traverses a given tree and returns a vector of all nodes with outdegree > 2,
 * i.e. all nodes that make the tree non-binary, in DFS order.
 * @param tree tree to traverse
 */
std::vector<pll_unode_t *> getNonBinaryNodesDFS(pll_unode_t * tree) {
  assert(tree->next == NULL);
  assert(tree->back != NULL);

  std::vector<pll_unode_t *> nodes;
  getNonBinaryNodesDFSRec(tree->back, nodes);
  return nodes;
}

/**
 * Check if vector a is a permutation of vector b, i.e. check if both vectors
 * contain the same elements but in a different order.
 * @param  a first vector
 * @param  b second vector
 * @return true iff a is permutation of b
 */
bool isPermutation(const std::vector<int> &a, const std::vector<int> &b) {
    if(a.size() != b.size()) {
        return false;
    }

    std::vector<int> a_copy(a.size());
    for (size_t i = 0; i < a.size(); i++) {
      a_copy[i] = a[i];
    }
    std::sort(a_copy.begin(), a_copy.end());

    std::vector<int> b_copy(b.size());
    for (size_t i = 0; i < b.size(); i++) {
      b_copy[i] = b[i];
    }
    std::sort(b_copy.begin(), b_copy.end());

    for (size_t i = 0; i < a_copy.size(); i++) {
      if(a_copy[i] != b_copy[i]) {
        return false;
      }
    }
    return true;
}

void printTreeEqualError(std::string message, pll_unode_t * node1, pll_unode_t * node2) {
  std::cout << "---------------------------------------------------\n" << message << "\n";
  std::cout << "Node 1:\n";
  printNode(node1);
  std::cout << "Node 2:\n";
  printNode(node2);
  std::cout << "---------------------------------------------------\n";
}

bool subnodesEqual(pll_unode_t * subnode1, pll_unode_t * subnode2) {
    assert(subnode1 != NULL);
    assert(subnode2 != NULL);

    if(subnode1->label == NULL) {
        if(subnode2->label != NULL) {
            printTreeEqualError("labels of subnodes are not identical", subnode1, subnode2);
            return false;
        }
    }
    if(subnode2->label == NULL) {
        if(subnode1->label != NULL) {
            printTreeEqualError("labels of subnodes are not identical", subnode1, subnode2);
            return false;
        }
    }
    if(subnode1->label != NULL && subnode2->label != NULL) {
        if(strcmp(subnode1->label, subnode2->label) != 0) {
          printTreeEqualError("labels of subnodes are not identical", subnode1, subnode2);
          return false;
        }
    }

    if(fabs(subnode1->length - subnode2->length) > 0.0000001) {
      printTreeEqualError("lengths of subnodes are not identical", subnode1, subnode2);
      return false;
    }
    return true;
}

bool treesEqualRec(pll_unode_t * node1, pll_unode_t * node2) {
    if(node1 == NULL) {
      if(node2 == NULL) {
        return true;
      } else {
        printTreeEqualError("node1 is NULL, node 2 not", node1, node2);
        return false;
      }
    }

    if(node2 == NULL) {
      if(node1 == NULL) {
        return true;
      } else {
        printTreeEqualError("node2 is NULL, node 1 not", node1, node2);
        return false;
      }
    }

    if(node1->next == NULL) {
      if (node2->next == NULL) {
        return subnodesEqual(node1, node2);
      } else {
        printTreeEqualError("node2 has more subnodes than node1", node1, node2);
        return false;
      }
    }

    if(node2->next == NULL) {
      if (node1->next == NULL) {
        return subnodesEqual(node1, node2);
      } else {
        printTreeEqualError("node1 has more subnodes than node1", node1, node2);
        return false;
      }
    }

    int subnodes_in_node1 = 1;
    pll_unode_t * temp = node1->next;
    while(temp != node1) {
      subnodes_in_node1++;
      temp = temp->next;
      assert(temp != NULL);
    }

    int subnodes_in_node2 = 1;
    temp = node2->next;
    while(temp != node2) {
      subnodes_in_node2++;
      temp = temp->next;
      assert(temp != NULL);
    }

    if(subnodes_in_node1 != subnodes_in_node2) {
      printTreeEqualError("different amount of subnodes", node1, node2);
      return false;
    }

    if(!subnodesEqual(node1, node2)) {
      return false;
    }
    pll_unode_t * temp1 = node1->next;
    pll_unode_t * temp2 = node2->next;
    while(temp1 != node1) {
      if(!subnodesEqual(temp1, temp2)) {
          return false;
      }
      if(!treesEqualRec(temp1->back, temp2->back)) {
          return false;
      }
      temp1 = temp1->next;
      temp2 = temp2->next;
      assert(temp != NULL);
    }
    return true;
}

bool treesEqual(pll_unode_t * node1, pll_unode_t * node2) {
    if(node1 == NULL) {
      if(node2 == NULL) {
        return true;
      } else {
        printTreeEqualError("node1 is NULL, node 2 not", node1, node2);
        return false;
      }
    }

    if(node2 == NULL) {
      if(node1 == NULL) {
        return true;
      } else {
        printTreeEqualError("node2 is NULL, node 1 not", node1, node2);
        return false;
      }
    }

    if(node1->next == NULL) {
      if (node2->next == NULL) {
        if(!subnodesEqual(node1, node2)) {
            return false;
        }
        if(!treesEqualRec(node1->back, node2->back)) {
            return false;
        }
      } else {
        printTreeEqualError("node2 has more subnodes than node1", node1, node2);
        return false;
      }
    }

    if(node2->next == NULL) {
      if (node1->next == NULL) {
        if(!subnodesEqual(node1, node2)) {
            return false;
        }
        if(!treesEqualRec(node1->back, node2->back)) {
            return false;
        }
      } else {
        printTreeEqualError("node1 has more subnodes than node1", node1, node2);
        return false;
      }
    }

    int subnodes_in_node1 = 1;
    pll_unode_t * temp = node1->next;
    while(temp != node1) {
      subnodes_in_node1++;
      temp = temp->next;
      assert(temp != NULL);
    }

    int subnodes_in_node2 = 1;
    temp = node2->next;
    while(temp != node2) {
      subnodes_in_node2++;
      temp = temp->next;
      assert(temp != NULL);
    }

    if(subnodes_in_node1 != subnodes_in_node2) {
      printTreeEqualError("different amount of subnodes", node1, node2);
      return false;
    }

    if(!subnodesEqual(node1, node2)) {
      return false;
    }
    pll_unode_t * temp1 = node1->next;
    pll_unode_t * temp2 = node2->next;
    while(temp1 != node1) {
      if(!subnodesEqual(temp1, temp2)) {
          return false;
      }
      if(!treesEqualRec(temp1->back, temp2->back)) {
          return false;
      }
      temp1 = temp1->next;
      temp2 = temp2->next;
      assert(temp != NULL);
    }
    return true;
}
