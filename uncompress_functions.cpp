 #include "uncompress_functions.h"

// pll_unode_t * createTreeRec(sdsl::bit_vector succinct_structure, unsigned int * succinct_idx,
//                         std::vector<pll_unode_t *> node_permutation, unsigned int * permutation_idx) {
//
//    assert(succinct_structure[*succinct_idx] == 0);
//    assert(*succinct_idx < succinct_structure.size() - 1);
//
//    if(succinct_structure[*succinct_idx + 1] == 0) {
//      // create a new node
//      pll_unode_t * new_innernode = pllmod_utree_create_node(0, 0, NULL, NULL);
//      (*succinct_idx)++;
//      new_innernode->next->back = createTreeRec(succinct_structure, succinct_idx, node_permutation, permutation_idx);
//      new_innernode->next->back->back = new_innernode->next;
//      assert(succinct_structure[*succinct_idx - 1] == 1);
//      assert(succinct_structure[*succinct_idx] == 0);
//      new_innernode->next->next->back = createTreeRec(succinct_structure, succinct_idx, node_permutation, permutation_idx);
//      new_innernode->next->next->back->back = new_innernode->next->next;
//      assert(succinct_structure[*succinct_idx - 1] == 1);
//      assert(succinct_structure[*succinct_idx] == 1);
//      (*succinct_idx)++;
//
//      return new_innernode;
//    } else {
//      assert(succinct_structure[*succinct_idx] == 0);
//      assert(succinct_structure[*succinct_idx + 1] == 1);
//
//      *succinct_idx = *succinct_idx + 2;
//
//      pll_unode_t * leaf = node_permutation[permutation_idx];
//      (*permutation_idx)++;
//      return leaf;
//    }
// }
//
// pll_unode_t * createTree(sdsl::bit_vector succinct_structure, unsigned int start,
//                         std::vector<pll_unode_t *> node_permutation) {
//
//     unsigned int succinct_idx = start;
//     unsigned int permutation_idx = 0;
//     pll_unode_t * tree = createTreeRec(succinct_structure, &succinct_idx, node_permutation, &permutation_idx);
//     assert(succinct_idx == succinct_structure.size());
//     assert(permutation_idx == node_permutation.size());
//     return tree;
// }

pll_unode_t * createLeaf(int index) {
    pll_unode_t * new_leaf = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
    if (!new_leaf) {
      // error
       assert(false);
       return NULL;
    }
    new_leaf->next = NULL;

    char * label = (char *) malloc(10);
    sprintf(label, "%d", index);

    new_leaf->label = label;
    new_leaf->data = NULL;
    new_leaf->length = 0;
    new_leaf->clv_index = 0;
    new_leaf->scaler_index = 0;

    return new_leaf;
}

pll_unode_t * createLeaf(double length, char * label) {
    pll_unode_t * new_leaf = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
    if (!new_leaf) {
      // error
       assert(false);
       return NULL;
    }
    new_leaf->next = NULL;

    new_leaf->label = label;
    int n = atoi(label);
    new_leaf->data = (void*)(intptr_t) n;
    new_leaf->length = length;
    new_leaf->clv_index = 0;
    new_leaf->scaler_index = 0;

    return new_leaf;
}

pll_unode_t * createTreeRec(sdsl::bit_vector succinct_structure, unsigned int * succinct_idx,
                        sdsl::int_vector<> node_permutation, unsigned int * node_idx) {

   assert(succinct_structure[*succinct_idx] == 0);
   assert(*succinct_idx < succinct_structure.size() - 1);

   if(succinct_structure[*succinct_idx + 1] == 0) {
     // create a new node
     pll_unode_t * new_innernode = pllmod_utree_create_node(0, 0, NULL, NULL);
     (*succinct_idx)++;
     new_innernode->next->back = createTreeRec(succinct_structure, succinct_idx, node_permutation, node_idx);
     new_innernode->next->back->back = new_innernode->next;
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 0);
     new_innernode->next->next->back = createTreeRec(succinct_structure, succinct_idx, node_permutation, node_idx);
     new_innernode->next->next->back->back = new_innernode->next->next;
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 1);
     (*succinct_idx)++;

     return new_innernode;
   } else {
     assert(succinct_structure[*succinct_idx] == 0);
     assert(succinct_structure[*succinct_idx + 1] == 1);

     // create new leaf
     pll_unode_t * new_leaf = createLeaf(node_permutation[*node_idx]);

     (*node_idx)++;
     *succinct_idx = *succinct_idx + 2;

     return new_leaf;
   }
}

pll_unode_t * createTreeRec(sdsl::bit_vector succinct_structure, unsigned int * succinct_idx,
                        std::vector<pll_unode_t *> leaves, unsigned int * node_idx) {

   assert(succinct_structure[*succinct_idx] == 0);
   assert(*succinct_idx < succinct_structure.size() - 1);

   if(succinct_structure[*succinct_idx + 1] == 0) {
     // create a new node
     pll_unode_t * new_innernode = pllmod_utree_create_node(0, 0, NULL, NULL);
     (*succinct_idx)++;
     new_innernode->next->back = createTreeRec(succinct_structure, succinct_idx, leaves, node_idx);
     new_innernode->next->back->back = new_innernode->next;
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 0);
     new_innernode->next->next->back = createTreeRec(succinct_structure, succinct_idx, leaves, node_idx);
     new_innernode->next->next->back->back = new_innernode->next->next;
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 1);
     (*succinct_idx)++;

     return new_innernode;
   } else {
     assert(succinct_structure[*succinct_idx] == 0);
     assert(succinct_structure[*succinct_idx + 1] == 1);

     // assign leaf
     pll_unode_t * new_leaf = leaves[*node_idx];

     (*node_idx)++;
     *succinct_idx = *succinct_idx + 2;

     return new_leaf;
   }
}

pll_unode_t * createTreeRec(std::vector<int> succinct_structure, unsigned int * succinct_idx,
                        std::vector<pll_unode_t *> leaves, unsigned int * node_idx) {

   assert(succinct_structure[*succinct_idx] == 0);
   assert(*succinct_idx < succinct_structure.size() - 1);

   if(succinct_structure[*succinct_idx + 1] == 0) {
     // create a new node
     pll_unode_t * new_innernode = pllmod_utree_create_node(0, 0, NULL, NULL);
     (*succinct_idx)++;
     new_innernode->next->back = createTreeRec(succinct_structure, succinct_idx, leaves, node_idx);
     new_innernode->next->back->back = new_innernode->next;
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 0);
     new_innernode->next->next->back = createTreeRec(succinct_structure, succinct_idx, leaves, node_idx);
     new_innernode->next->next->back->back = new_innernode->next->next;
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 1);
     (*succinct_idx)++;

     return new_innernode;
   } else {
     assert(succinct_structure[*succinct_idx] == 0);
     assert(succinct_structure[*succinct_idx + 1] == 1);

     // assign leaf
     pll_unode_t * new_leaf = leaves[*node_idx]->back;

     (*node_idx)++;
     *succinct_idx = *succinct_idx + 2;

     return new_leaf;
   }
}

/**
 * Creates a tree from the given structures.
 * @param  succinct_structure topology of the tree
 * @param  node_permutation   permutation of the nodes in the tree
 * @return                    root of the created tree
 */
pll_unode_t * createTree(sdsl::bit_vector succinct_structure, sdsl::int_vector<> node_permutation) {

      unsigned int succinct_idx = 0;
      unsigned int node_idx = 0;
      pll_unode_t * tree = createTreeRec(succinct_structure, &succinct_idx, node_permutation, &node_idx);
      assert(succinct_idx == succinct_structure.size());
      assert(node_idx == node_permutation.size());
      return tree;
}

/**
 * Creates a tree from the topology given by succinct_structure and assigns the
 * leaves given by leaves to the created tree.
 * @param  succinct_structure topology of the tree
 * @param  leaves             leaves of the tree, order is the order you would
 *                            see the leaves in depth-first search
 * @return                    root of the created tree
 */
pll_unode_t * createTree(sdsl::bit_vector succinct_structure, std::vector<pll_unode_t *> leaves) {

      unsigned int succinct_idx = 0;
      unsigned int node_idx = 0;
      pll_unode_t * tree = createTreeRec(succinct_structure, &succinct_idx, leaves, &node_idx);
      assert(succinct_idx == succinct_structure.size());
      assert(node_idx == leaves.size());
      return tree;
}

pll_unode_t * createTree(std::vector<int> succinct_structure, std::vector<pll_unode_t *> leaves) {
      unsigned int succinct_idx = 0;
      unsigned int node_idx = 0;
      pll_unode_t * tree = createTreeRec(succinct_structure, &succinct_idx, leaves, &node_idx);
      assert(succinct_idx == succinct_structure.size());

      assert(node_idx == leaves.size());
      return tree;
}

pll_unode_t * simple_uncompression(sdsl::bit_vector &succinct_structure, sdsl::int_vector<> &node_permutation) {

  pll_unode_t * tree = createTree(succinct_structure, node_permutation);

  assert(atoi(tree->next->back->label) == 1);

  pll_unode_t * root = tree->next->back;
  root->back = tree->next->next->back;
  tree->next->next->back->back = root;

  return root;
}

// void contractEdgesRec(pll_unode_t * node, unsigned int * edge_idx,
//         sdsl::int_vector<> &edges_to_contract, unsigned int * edge_to_contract_idx) {
//   assert(node->next = NULL || node->next->next = node);
//   if(*edge_idx == edges_to_contract[*edge_to_contract_idx]) {
//     (*edge_to_contract_idx)++;
//     contractEdgesRec(node->next, , edges_to_contract, edge_to_contract_idx);
//     contractEdgesRec(node->next->next, , edges_to_contract, edge_to_contract_idx);
//     contractEdge(node);
//   }
// }
//
// pll_unode_t * copyNode(pll_unode_t * node) {
//     char * copied_label;
//     strncpy(copied_label, node->label, sizeof(node->label));
//     return pllmod_utree_create_node(0, 0, copied_label, NULL);
// }

pll_unode_t * copyTreeRec(const pll_unode_t * original) {
    if(original->next == NULL) {
        // leaf
        return createLeaf(original->length, original->label);
    } else {
        assert(original->next != NULL);
        assert(original->next->next->next == original); // tree is binary

        pll_unode_t * copied_node = pllmod_utree_create_node(0, 0, NULL, NULL);

        copied_node->label = original->label;
        copied_node->data = original->data;
        copied_node->length = original->length;
        copied_node->next->label = original->next->label;
        copied_node->next->data = original->next->data;
        copied_node->next->length = original->next->length;
        copied_node->next->next->label = original->next->next->label;
        copied_node->next->next->data = original->next->next->data;
        copied_node->next->next->length = original->next->next->length;

        copied_node->next->back = copyTreeRec(original->next->back);
        copied_node->next->back->back = copied_node->next;

        copied_node->next->next->back = copyTreeRec(original->next->next->back);
        copied_node->next->next->back->back = copied_node->next->next;

        return copied_node;
    }
}

/**
 * Takes a binary tree and creates a copy of its topology,
 * including copies of the labels of each node
 * @param  tree the tree to copy
 * @return      creates copy of the trees topology
 */
pll_unode_t * copyTree(const pll_unode_t * tree) {
    if(tree->next == NULL && tree->back >= NULL) {
        // root of the tree is a leaf
        pll_unode_t * root = createLeaf(tree->length, tree->label);
        root->back = copyTreeRec(tree->back);
        root->back->back = root;
        return root;
    }
    return copyTreeRec(tree);
}

void traverseAndDeleteEdgesRec(pll_unode_t * tree, sdsl::int_vector<> &edges_to_contract,
                  unsigned int * edges_to_contract_idx, unsigned int * edges_idx,
                  std::vector<pll_unode_t *> &nodes_to_contract) {

  bool contract = false;

  if(*edges_to_contract_idx < edges_to_contract.size() &&
                *edges_idx == edges_to_contract[*edges_to_contract_idx]) {
      assert(tree->next != NULL);
      (*edges_to_contract_idx)++;
      nodes_to_contract.push_back(tree);
  }
  if(tree->next != NULL) {
      // inner node
      pll_unode_t * temp = tree->next;
      while(temp != tree) {
        (*edges_idx)++;
        traverseAndDeleteEdgesRec(temp->back, edges_to_contract, edges_to_contract_idx, edges_idx, nodes_to_contract);
        temp = temp->next;
        assert(temp != NULL);
      }
  }
  // if(contract) {
  //   contractEdge(tree);
  // }
}

/**
 * Traverses the given tree and contracts the edges given by the vector "edges_to_contract"
 * @param tree              tree
 * @param edges_to_contract edges to contract in the tree
 */
void traverseAndDeleteEdges(pll_unode_t * tree, sdsl::int_vector<> &edges_to_contract) {
    if(edges_to_contract.size() == 0) {
      return;
    }
    assert(edges_to_contract[0] > 2);
    unsigned int edges_to_contract_idx = 0;
    unsigned int edges_idx = 2;
    std::vector<pll_unode_t *> nodes_to_contract;
    traverseAndDeleteEdgesRec(tree->back, edges_to_contract, &edges_to_contract_idx, &edges_idx, nodes_to_contract);

    // std::cout << nodes_to_contract.size() << "\n";

    for(auto node: nodes_to_contract){
      contractEdge(node);
    }
}

pll_unode_t * buildSubtreeRec(sdsl::bit_vector &subtrees_succinct, unsigned int * idx) {
    assert(subtrees_succinct[*idx] != 1);
    if(subtrees_succinct[*idx] == 0) {
        if(subtrees_succinct[*idx+1] == 0) {
            pll_unode_t * new_innernode_0 = pllmod_utree_create_node(0, 0, NULL, NULL);
            pll_unode_t * new_innernode_1 = pllmod_utree_create_node(0, 0, NULL, NULL);
            pll_unode_t * new_innernode_2 = pllmod_utree_create_node(0, 0, NULL, NULL);

            new_innernode_0->next = new_innernode_1;
            new_innernode_1->next = new_innernode_2;
            new_innernode_2->next = new_innernode_0;

            (*idx)++;
            pll_unode_t * new_innernode_child_1 = buildSubtreeRec(subtrees_succinct, idx);
            pll_unode_t * new_innernode_child_2 = buildSubtreeRec(subtrees_succinct, idx);
            (*idx)++;
        } else {
            // sequence 01 found
            (*idx)++;
            (*idx)++;
            return pllmod_utree_create_node(0, 0, NULL, NULL);
            // TODO: add two childs
        }
    }
}

pll_unode_t * buildSubtree(std::vector<pll_unode_t *> children, sdsl::bit_vector &subtrees_succinct,
          sdsl::int_vector<> &succinct_permutations) {
    assert(subtrees_succinct[0] == 0);
    unsigned int idx = 0;
}

// void expandTreeRec(pll_unode_t * tree) {
//   if(tree->next == NULL) {
//       // leaf
//   } else {
//       assert(tree->next != NULL);
//
//       if(tree->next->next->next != tree) {
//           // tree is not binary
//
//           // children to array
//           std::vector<pll_unode_t *> children;
//           pll_unode_t * temp = tree->next;
//           while(temp != tree) {
//             children.push_back(temp);
//             traverseAndDeleteEdgesRec(temp->back, (*edges_idx)++);
//             temp = temp->next;
//             assert(temp != NULL);
//           }
//           assert(copy_temp->next == NULL);
//           copy_temp->next = copy;
//
//       } else {
//           // tree is binary
//           expandTreeRec(tree->next);
//           expandTreeRec(tree->next->next);
//       }
//
//       pll_unode_t * temp = tree->next;
//       while(temp != tree) {
//
//         traverseAndDeleteEdgesRec(temp->back, (*edges_idx)++);
//
//         temp = temp->next;
//         assert(temp != NULL);
//       }
//       assert(copy_temp->next == NULL);
//       copy_temp->next = copy;
//   }
// }

pll_unode_t * rf_distance_uncompression(const pll_unode_t * predecessor_tree, sdsl::int_vector<> &edges_to_contract,
          sdsl::bit_vector &subtrees_succinct, sdsl::int_vector<> &succinct_permutations) {

  // split subtrees_succinct
  std::vector<std::vector<int>> subtrees_split;

  if(!(subtrees_succinct.size() == 0)) {
    assert(subtrees_succinct[0] == 0);

    size_t temp_start_index = 0;
    unsigned int temp_ctr = 1;
    size_t temp_size = 1;

    while(temp_start_index + temp_size < subtrees_succinct.size()) {
      while(temp_ctr > 0) {
        if(subtrees_succinct[temp_start_index + temp_size] == 0) {
          temp_ctr++;
        } else {
          temp_ctr--;
        }
        temp_size++;
        assert(temp_ctr >= 0);
      }
      assert(temp_ctr == 0);
      subtrees_split.push_back(std::vector<int> (subtrees_succinct.begin() + temp_start_index, subtrees_succinct.begin() + temp_start_index + temp_size));

      temp_start_index = temp_start_index + temp_size;
      //assert(subtrees_succinct[temp_start_index] == 0);
      temp_size = 1;
      assert(temp_start_index <= subtrees_succinct.size());
      temp_ctr = 1;
    }
    assert(temp_ctr == 1);
  }

  // for(auto order: subtrees_split) {
  //   for(auto x: order) {
  //     std::cout << x;
  //   }
  //   std::cout << "  " << order.size();
  //   std::cout << "\n\n";
  // }

  // split succinct_permutations
  std::vector<std::vector<int>> permutations;
  size_t start_idx_succinct = 0;
  size_t end_idx_succinct;
  for (size_t i = 0; i < subtrees_split.size(); i++) {
      assert(subtrees_split[i].size() % 2 == 0);
      assert((subtrees_split[i].size() + 2) % 4 == 0);
      // std::cout << (subtrees_split[i].size() + 2) / 4 << "\n";
      end_idx_succinct = start_idx_succinct + (subtrees_split[i].size() + 2) / 4;
      assert(end_idx_succinct <= succinct_permutations.size());

      // for (size_t j = start_idx_succinct; j < end_idx_succinct; j++) {
      //     // TODO: Leafs to array
      //     createLeaf(subtrees_split[i][j]);
      // }


      permutations.push_back(std::vector<int> (succinct_permutations.begin() + start_idx_succinct,
                      succinct_permutations.begin() + end_idx_succinct));
      start_idx_succinct = end_idx_succinct;
  }
  assert(start_idx_succinct == succinct_permutations.size());


  // assert predecessor_tree ordered
  pll_unode_t * tree = copyTree(predecessor_tree);

  //printTree(tree);

  traverseAndDeleteEdges(tree, edges_to_contract);

  //printTree(tree);

  std::vector<pll_unode_t *> consensus_subtree_roots;
  std::vector<std::vector<pll_unode_t *>> consensus_orders;
  traverseConsensus(tree, consensus_subtree_roots, consensus_orders);

  // for(auto order: consensus_orders) {
  //   for(auto x: order) {
  //     printNode(x);
  //   }
  //   std::cout << "\n\n";
  // }
  //
  // for(auto perm: permutations) {
  //   for(auto x: perm) {
  //     std::cout << x << " ";
  //   }
  //   std::cout << ";    ";
  // }
  // std::cout << "\n";

  assert(consensus_orders.size() == permutations.size());
  for (size_t i = 0; i < consensus_orders.size(); i++) {
    // std::cout << "consensus_orders[" << i << "].size(): " << consensus_orders[i].size() << "\n";
    // std::cout << "permutations[" << i << "].size(): " << permutations[i].size() << "\n";

    assert(consensus_orders[i].size() == permutations[i].size());
  }

  //TODO: apply permutations

  std::vector<std::vector<pll_unode_t *>> new_orders;
  for (size_t i = 0; i < consensus_orders.size(); i++) {
      std::vector<pll_unode_t *> new_order(consensus_orders[i].size());
      for (size_t j = 0; j < consensus_orders[i].size(); j++) {
        new_order[j] = consensus_orders[i][permutations[i][j]];
      }

      for (size_t j = 0; j < consensus_orders[i].size(); j++) {
        assert(new_order[j] != NULL);
      }

      new_orders.push_back(new_order);
  }

  // for(auto order: new_orders) {
  //   for(auto x: order) {
  //     printNode(x);
  //   }
  //   std::cout << "\n\n";
  // }

  assert(subtrees_split.size() == new_orders.size());
  for (size_t i = 0; i < subtrees_split.size(); i++) {
    // for(auto x: subtrees_split[i]) {
    //   std::cout << x;
    // }
    // std::cout << " <---> ";
    // for(auto y: new_orders[i]) {
    //   std::cout << y<< ", ";
    // }
    // std::cout << "\n";
    pll_unode_t * subtree = createTree(subtrees_split[i], new_orders[i]);

    // printTree(subtree);
    // std::cout << "\n\n";

    // replace part in consensus tr
    assert(consensus_subtree_roots[i] != NULL);
    assert(consensus_subtree_roots[i]->back != NULL);
    assert(consensus_subtree_roots[i]->next != NULL); // don't contract leafs

    consensus_subtree_roots[i]->back->back = subtree;
    subtree->back = consensus_subtree_roots[i]->back;

    // TODO: free consensus_subtree_roots

    // pll_unode_t * old_root_node_predecessor = internalPredecessor(consensus_subtree_roots[i]);
    // pll_unode_t * old_root_node_successor = consensus_subtree_roots[i]->next;
    //
    // old_root_node_predecessor->next = subtree;
    // subtree->next = old_root_node_successor;

    // TODO: free consensus_subtree_roots[i]
  }

  //printTree(tree);

  // std::vector<pll_unode_t *> subtrees;
  //
  // subtrees: subtrees.push_back(createTree(subtree, tree2_order))

  // assert(subtrees.size() == nodes_to_replace.size());
  // for (size_t i = 0; i < subtrees.size(); i++) {
  //   // replace the nodes in consensus tree
  //   nodes_to_replace[i]->back->back = subtrees[i];
  //   subtrees[i]->back = nodes_to_replace[i]->back;
  //   free(nodes_to_replace[i]);
  // }

  return tree;
}
