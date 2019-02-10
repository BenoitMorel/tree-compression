 #include "uncompress_functions.h"

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
                        sdsl::int_vector<> node_permutation, unsigned int * node_idx,
                        std::vector<double> branch_lengths, unsigned int * branch_idx) {

   assert(succinct_structure[*succinct_idx] == 0);
   assert(*succinct_idx < succinct_structure.size() - 1);

   if(succinct_structure[*succinct_idx + 1] == 0) {
     // create a new node
     pll_unode_t * new_innernode = pllmod_utree_create_node(0, 0, NULL, NULL);
     double branch_length;
     (*succinct_idx)++;
     branch_length = branch_lengths[*branch_idx];
     (*branch_idx)++;
     new_innernode->next->back = createTreeRec(succinct_structure, succinct_idx, node_permutation, node_idx, branch_lengths, branch_idx);
     new_innernode->next->back->length = branch_length;
     new_innernode->next->back->back = new_innernode->next;
     new_innernode->next->back->back->length = branch_length;

     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 0);
     branch_length = branch_lengths[*branch_idx];
     (*branch_idx)++;
     new_innernode->next->next->back = createTreeRec(succinct_structure, succinct_idx, node_permutation, node_idx, branch_lengths, branch_idx);
     new_innernode->next->next->back->length = branch_length;
     new_innernode->next->next->back->back = new_innernode->next->next;
     new_innernode->next->next->back->back->length = branch_length;
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

pll_unode_t * createTreeRecSpecial(std::vector<int> succinct_structure, unsigned int * succinct_idx) {

   assert(succinct_structure[*succinct_idx] == 0);
   assert(*succinct_idx < succinct_structure.size() - 1);

   if(succinct_structure[*succinct_idx + 1] == 0) {
     // create a new node
     pll_unode_t * new_innernode = pllmod_utree_create_node(0, 0, NULL, NULL);
     (*succinct_idx)++;

     pll_unode_t * leaf = createTreeRecSpecial(succinct_structure, succinct_idx);
     new_innernode->next->back = leaf;
     if(leaf != NULL) {
        leaf->back = new_innernode->next;
     }
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 0);

     leaf = createTreeRecSpecial(succinct_structure, succinct_idx);
     new_innernode->next->next->back = leaf;
     if(leaf != NULL) {
        leaf->back = new_innernode->next->next;
     }
     assert(succinct_structure[*succinct_idx - 1] == 1);
     assert(succinct_structure[*succinct_idx] == 1);
     (*succinct_idx)++;

     assert(new_innernode != NULL);
     assert(new_innernode->next != NULL);
     assert(new_innernode->next->next != NULL);

     return new_innernode;
   } else {
     assert(succinct_structure[*succinct_idx] == 0);
     assert(succinct_structure[*succinct_idx + 1] == 1);

     // assign leaf
     *succinct_idx = *succinct_idx + 2;

     return NULL;
   }
}

/**
 * Creates a tree from the given structures.
 * @param  succinct_structure topology of the tree
 * @param  node_permutation   permutation of the nodes in the tree
 * @return                    root of the created tree
 */
pll_unode_t * createTree(sdsl::bit_vector succinct_structure, sdsl::int_vector<> node_permutation,
              std::vector<double> branch_lengths) {

      unsigned int succinct_idx = 0;
      unsigned int node_idx = 0;
      unsigned int branch_idx = 0;

      pll_unode_t * tree = createTreeRec(succinct_structure, &succinct_idx, node_permutation, &node_idx,
                  branch_lengths, &branch_idx);
      assert(succinct_idx == succinct_structure.size());
      assert(node_idx == node_permutation.size());
      assert(branch_idx == branch_lengths.size());

      return tree;
}

void assignBranchLengthsSubtreeRec(pll_unode_t * tree, std::vector<pll_unode_t *> &leaves, unsigned int * node_idx,
                  std::vector<double> branch_lengths, unsigned int * branch_idx){
  assert(tree != NULL);
  assert(tree->next != NULL);
  assert(tree->next->next->next == tree); // assert tree is binary

  if(tree->next->back != NULL) {
    double bl = branch_lengths[*branch_idx];
    tree->next->length = bl;
    tree->next->back->length = bl;
    (*branch_idx)++;

    assignBranchLengthsSubtreeRec(tree->next->back, leaves, node_idx, branch_lengths, branch_idx);
  } else {
    pll_unode_t * new_leaf = leaves[*node_idx]->back;
    (*node_idx)++;

    assert((intptr_t) new_leaf->back->data != 0);

    tree->next->back = new_leaf;
    new_leaf->back = tree->next;
    tree->next->length = new_leaf->length;
    tree->next->data = new_leaf->data;
  }

  if(tree->next->next->back != NULL) {
    double bl = branch_lengths[*branch_idx];
    tree->next->next->length = bl;
    tree->next->next->back->length = bl;
    (*branch_idx)++;

    assignBranchLengthsSubtreeRec(tree->next->next->back, leaves, node_idx, branch_lengths, branch_idx);
  } else {
    pll_unode_t * new_leaf = leaves[*node_idx]->back;
    (*node_idx)++;

    assert((intptr_t) new_leaf->back->data != 0);

    tree->next->next->back = new_leaf;
    new_leaf->back = tree->next->next;
    tree->next->next->length = new_leaf->length;
    tree->next->next->data = new_leaf->data;
  }
}

pll_unode_t * createTreeSpecial(std::vector<int> succinct_structure, std::vector<pll_unode_t *> &leaves,
                std::vector<double> branch_lengths) {

      unsigned int succinct_idx = 0;
      unsigned int node_idx = 0;

      pll_unode_t * tree = createTreeRecSpecial(succinct_structure, &succinct_idx);

      assert(succinct_idx == succinct_structure.size());

      unsigned int branch_idx = 0;
      assignBranchLengthsSubtreeRec(tree, leaves, &node_idx, branch_lengths, &branch_idx);
      assert(branch_idx == branch_lengths.size());



      return tree;
}

pll_unode_t * simple_uncompression(sdsl::bit_vector &succinct_structure, sdsl::int_vector<> &node_permutation, std::vector<double> branch_lengths) {

  pll_unode_t * tree = createTree(succinct_structure, node_permutation, branch_lengths);

  assert(atoi(tree->next->back->label) == 1);

  pll_unode_t * root = tree->next->back;
  root->back = tree->next->next->back;
  tree->next->next->back->back = root;
  root->length = root->back->length;

  return root;
}

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

    for(auto node: nodes_to_contract){
      contractEdge(node);
    }
}

void applyBranchLengthDiffsRec(pll_unode_t * tree, std::vector<double> consensus_branch_diffs,
                      unsigned int * branches_idx) {
    assert(tree != NULL);

    if(((intptr_t) tree->data != 0) || ((intptr_t) tree->back->data != 0)) {
        // apply branch diff
        double new_bl = tree->length + consensus_branch_diffs[*branches_idx];
        (*branches_idx)++;
        tree->length = new_bl;
        tree->back->length = new_bl;
    } else {
        assert((intptr_t) tree->data == 0 && (intptr_t) tree->back->data == 0);
    }

    if(tree->next == NULL) {
      return;
    }

    assert(tree->next->next->next == tree);

    applyBranchLengthDiffsRec(tree->next->back, consensus_branch_diffs, branches_idx);
    applyBranchLengthDiffsRec(tree->next->next->back, consensus_branch_diffs, branches_idx);
}

void applyBranchLengthDiffs(pll_unode_t * tree, std::vector<double> consensus_diffs) {
    unsigned int branches_idx = 1;
    applyBranchLengthDiffsRec(tree->back, consensus_diffs, &branches_idx);
    assert(branches_idx == consensus_diffs.size());
}

pll_unode_t * rf_distance_uncompression(const pll_unode_t * predecessor_tree, sdsl::int_vector<> &edges_to_contract,
          sdsl::bit_vector &subtrees_succinct, sdsl::int_vector<> &succinct_permutations,
          std::vector<double> consensus_branches, std::vector<double> non_consensus_branches) {

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
      temp_size = 1;
      assert(temp_start_index <= subtrees_succinct.size());
      temp_ctr = 1;
    }
    assert(temp_ctr == 1);
  }

  // split succinct_permutations
  std::vector<std::vector<int>> permutations;
  size_t start_idx_succinct = 0;
  size_t end_idx_succinct;
  for (size_t i = 0; i < subtrees_split.size(); i++) {
      assert(subtrees_split[i].size() % 2 == 0);
      assert((subtrees_split[i].size() + 2) % 4 == 0);
      end_idx_succinct = start_idx_succinct + (subtrees_split[i].size() + 2) / 4;
      assert(end_idx_succinct <= succinct_permutations.size());

      permutations.push_back(std::vector<int> (succinct_permutations.begin() + start_idx_succinct,
                      succinct_permutations.begin() + end_idx_succinct));
      start_idx_succinct = end_idx_succinct;
  }
  assert(start_idx_succinct == succinct_permutations.size());

  // split non_consensus_branch_lengths
  std::vector<std::vector<double>> branch_lengths;
  start_idx_succinct = 0;
  for (size_t i = 0; i < permutations.size(); i++) {
      end_idx_succinct = start_idx_succinct + (permutations[i].size() - 2);
      assert(end_idx_succinct <= succinct_permutations.size());

      branch_lengths.push_back(std::vector<double> (non_consensus_branches.begin() + start_idx_succinct,
                      non_consensus_branches.begin() + end_idx_succinct));
      start_idx_succinct = end_idx_succinct;
  }
  assert(start_idx_succinct == non_consensus_branches.size());

  // assert predecessor_tree ordered
  pll_unode_t * tree = copyTree(predecessor_tree);

  traverseAndDeleteEdges(tree, edges_to_contract);

  std::vector<pll_unode_t *> consensus_subtree_roots;
  std::vector<std::vector<pll_unode_t *>> consensus_orders;
  traverseConsensus(tree, consensus_subtree_roots, consensus_orders);

  assert(consensus_orders.size() == permutations.size());
  for (size_t i = 0; i < consensus_orders.size(); i++) {
    assert(consensus_orders[i].size() == permutations[i].size());
  }

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

  assert(subtrees_split.size() == new_orders.size());
  for (size_t i = 0; i < subtrees_split.size(); i++) {

    pll_unode_t * subtree = createTreeSpecial(subtrees_split[i], new_orders[i], branch_lengths[i]);

    assert(consensus_subtree_roots[i] != NULL);
    assert(consensus_subtree_roots[i]->back != NULL);
    assert(consensus_subtree_roots[i]->next != NULL); // don't contract leafs

    assert((intptr_t) consensus_subtree_roots[i]->data != 0);
    assert((intptr_t) subtree->data == 0);

    consensus_subtree_roots[i]->back->back = subtree;
    subtree->back = consensus_subtree_roots[i]->back;

    subtree->length = subtree->back->length;
    subtree->data = consensus_subtree_roots[i]->data;

    // TODO: free consensus_subtree_roots
    // TODO: free consensus_subtree_roots[i]
  }

  applyBranchLengthDiffs(tree, consensus_branches);

  return tree;
}
