#include "compress_functions.h"
#include "datastructure_compression_functions.h"

int simple_compression(const char * tree_file, const char * succinct_structure_file,
        const char * node_permutation_file, const char * branch_lengths_file, int flags) {

  /* tree properties */
  pll_utree_t * tree = NULL;
  unsigned int tip_count;

  /* parse the input tree */
  tree = pll_utree_parse_newick (tree_file);
  if(tree == NULL) {
      // ERROR: tree could not be parsed
      // --> syntax of newick file is not correct
      // --> tree has less than 3 leaves
      return -1;
  }
  tip_count = tree->tip_count;
  assert(tip_count >= 3);

  // search for the root of the tree
  pll_unode_t * root = searchRoot(tree);

  // set and order the tree
  setTree(root);
  orderTree(root);

  // succinct_structure stores the topology in balanced parantheses ("0=(, 1=)")
  sdsl::bit_vector succinct_structure(4 * tip_count - 2, 0);
  // succinct_structure stores the permutation of the taxa
  sdsl::int_vector<> node_permutation(tip_count, 0, 32);
  // branch_lengths stores all branch lengths of the tree
  std::vector<double> branch_lengths(2 * tip_count - 2);
  // create mapping from node_id to branch_id
  unsigned int* node_id_to_branch_id = (unsigned int*) malloc ((tree->inner_count * 3 + tree->tip_count) * sizeof(unsigned int));
  // fill the created structures with the given tree
  assignBranchNumbers(root, succinct_structure, node_permutation, branch_lengths, node_id_to_branch_id);

  auto size_topology = sdsl::size_in_bytes(succinct_structure);

  // TODO: compress permutation
  sdsl::util::bit_compress(node_permutation);
  auto size_node_permutation = sdsl::size_in_bytes(node_permutation);

  // TODO: compress branches
  auto size_branches = writeDoubleVector(branch_lengths, branch_lengths_file);

  if (flags & PRINT_COMPRESSION_STRUCTURES) {
    std::cout << "Succinct representation: " << succinct_structure << "\n";
    std::cout << "\tuncompressed size: " << size_topology << " bytes\n";

    std::cout << "Node permutation: " << node_permutation << "\n";
    std::cout << "\tcompressed size: " << size_node_permutation << " bytes\n";

    std::cout << "Branch Lengths: ";
    for(auto x: branch_lengths) {
      std::cout << x << " ";
    }
    std::cout << "\n\tcompressed size: " << size_branches << " bytes\n";
  }

  // write stuctures to file
  store_to_file(succinct_structure, succinct_structure_file);
  store_to_file(node_permutation, node_permutation_file);

  if(flags & PRINT_COMPRESSION) {
    std::cout << "Simple compression size: " << size_topology
    << " (topology) + " << size_node_permutation << " (leaves) + " << size_branches
    << " (branches) = " << size_topology + size_node_permutation + size_branches
    << " bytes\n";

    std::cout << "---------------------------------------------------------\n";
  }

  free(node_id_to_branch_id);

  return 0;
}

void commonBranchesOrderedRec(pll_unode_t * tree, const std::vector<bool> &edgeIncidentPresent2, std::vector<double> &branches) {
  assert(tree != NULL);

  if(!edgeIncidentPresent2[tree->node_index]) {
    branches.push_back(tree->length);
  }
  if(tree->next == NULL) {
    // leaf
  } else {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next->next == tree);

    // if(edgeIncidentPresent2[tree->next->node_index]) {
    //     branches.push_back(tree->next->length);
    // }
    commonBranchesOrderedRec(tree->next->back, edgeIncidentPresent2, branches);
    // if(edgeIncidentPresent2[tree->next->next->node_index]) {
    //     branches.push_back(tree->next->next->length);
    // }
    commonBranchesOrderedRec(tree->next->next->back, edgeIncidentPresent2, branches);

  }
}

std::vector<double> commonBranchesOrdered(pll_unode_t * tree, const std::vector<bool> &edgeIncidentPresent2) {
    std::vector<double> branches;
    branches.push_back(0);
    commonBranchesOrderedRec(tree, edgeIncidentPresent2, branches);
    return branches;
}



void commonBranchesOrderedCompareRec(pll_unode_t * tree, pll_unode_t * consensus, const std::vector<bool> &edgeIncidentPresent2, std::vector<double> &branches) {
  assert(tree != NULL);
  assert(consensus != NULL);

  if(tree->next == NULL) {
      assert(consensus->next == NULL);
      return;
  }

  assert(tree->next->next->next == tree); // tree is binary

  int i = 1;
  pll_unode_t * temp = consensus->next;
  while(temp != consensus) {
      temp = temp->next;
      i++;
  }

  if(i > 3) {
      int x = (intptr_t) tree->next->data;
      int y = (intptr_t) tree->next->next->data;
      pll_unode_t * subtree1;
      pll_unode_t * subtree2;

      pll_unode_t * temp = consensus->next;
      if((intptr_t) temp->data == x) {
          subtree1 = temp;
      }
      if((intptr_t) temp->data == y) {
          subtree2 = temp;
      }
      while(temp != consensus) {
          if((intptr_t) temp->data == x) {
              subtree1 = temp;
          }
          if((intptr_t) temp->data == y) {
              subtree2 = temp;
          }
          temp = temp->next;
      }
      assert(subtree1 != NULL);
      assert(subtree2 != NULL);

      if(!edgeIncidentPresent2[tree->next->node_index]) {
          assert((intptr_t) tree->next->data == (intptr_t) subtree1->data);
          branches.push_back(tree->next->length - subtree1->length);
          commonBranchesOrderedCompareRec(tree->next->back, subtree1->back, edgeIncidentPresent2, branches);
      } else {
          assert((intptr_t) tree->next->data == (intptr_t) subtree1->data);
          commonBranchesOrderedCompareRec(tree->next->back, consensus, edgeIncidentPresent2, branches);
      }

      if(!edgeIncidentPresent2[tree->next->next->node_index]) {
          assert((intptr_t) tree->next->next->data == (intptr_t) subtree2->data);
          branches.push_back(tree->next->next->length - subtree2->length);
          commonBranchesOrderedCompareRec(tree->next->next->back, subtree2->back, edgeIncidentPresent2, branches);
      } else {
          assert((intptr_t) tree->next->next->data == (intptr_t) subtree2->data);
          commonBranchesOrderedCompareRec(tree->next->next->back, consensus, edgeIncidentPresent2, branches);
      }

  } else {
      assert(i == 3);
      assert((intptr_t) tree->next->data == (intptr_t) consensus->next->data);
      assert((intptr_t) tree->next->next->data == (intptr_t) consensus->next->next->data);
      branches.push_back(tree->next->length - consensus->next->length);
      commonBranchesOrderedCompareRec(tree->next->back, consensus->next->back, edgeIncidentPresent2, branches);
      branches.push_back(tree->next->next->length - consensus->next->next->length);
      commonBranchesOrderedCompareRec(tree->next->next->back, consensus->next->next->back, edgeIncidentPresent2, branches);
  }
}

std::vector<double> commonBranchesOrderedCompare(pll_unode_t * tree, pll_unode_t * consensus, const std::vector<bool> &edgeIncidentPresent2) {
    std::vector<double> branches;
    branches.push_back(0);
    branches.push_back(tree->length - consensus->length);
    commonBranchesOrderedCompareRec(tree, consensus, edgeIncidentPresent2, branches);
    return branches;
}


std::tuple<std::vector<int>, std::vector<int>, std::vector<double>> findRFSubtreesRecL(pll_unode_t * tree, const std::vector<bool> &edgeIncidentPresent2,
            std::stack<pll_unode_t *> &tasks) {
  assert(tree != NULL);
  std::vector<int> return_vector_topology;
  std::vector<int> return_vector_order;
  std::vector<double> return_branches;

  if(tree->next == NULL) {
    // leaf; edge incident is always present in both trees
    //std::cout << "Leaf: " << tree->node_index << "\t    Edge incident present: " << edgeIncidentPresent2[tree->node_index] << "\n";
  } else {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next != NULL);
    assert(tree->next->back != NULL);
    assert(tree->next->next->back != NULL);
    assert(tree->next->next->next = tree);

    //std::cout << "Inner node: " << tree->node_index << "\t    Edge incident present: " << edgeIncidentPresent2[tree->node_index] << "\n";

    //std::cout << "1. Inner node: " << tree->next->node_index <<
    //   "\t    Edge incident present: " << edgeIncidentPresent2[tree->next->node_index] << "\n";
    if(edgeIncidentPresent2[tree->next->node_index] && edgeIncidentPresent2[tree->next->next->node_index]) {
          std::vector<int> temp_topology_right;
          std::vector<int> temp_order_right;
          std::vector<double> temp_branches_right;
          std::tie (temp_topology_right, temp_order_right, temp_branches_right) = findRFSubtreesRecL(tree->next->next->back, edgeIncidentPresent2, tasks);

          std::vector<int> temp_topology_left;
          std::vector<int> temp_order_left;
          std::vector<double> temp_branches_left;
          std::tie (temp_topology_left, temp_order_left, temp_branches_left) = findRFSubtreesRecL(tree->next->back, edgeIncidentPresent2, tasks);

          return_vector_topology.push_back(0);
          return_vector_topology.insert(return_vector_topology.end(), temp_topology_left.begin(), temp_topology_left.end());
          return_vector_order.insert(return_vector_order.end(), temp_order_left.begin(), temp_order_left.end());
          return_vector_topology.push_back(1);
          return_vector_topology.push_back(0);
          return_vector_topology.insert(return_vector_topology.end(), temp_topology_right.begin(), temp_topology_right.end());
          return_vector_order.insert(return_vector_order.end(), temp_order_right.begin(), temp_order_right.end());
          return_vector_topology.push_back(1);
          return_branches.push_back(tree->next->back->length);
          return_branches.insert(return_branches.end(), temp_branches_left.begin(), temp_branches_left.end());
          return_branches.push_back(tree->next->next->back->length);
          return_branches.insert(return_branches.end(), temp_branches_right.begin(), temp_branches_right.end());
    } else if(edgeIncidentPresent2[tree->next->node_index]) {
          assert(!edgeIncidentPresent2[tree->next->next->node_index]);

          tasks.push(tree->next->next->back);

          return_vector_topology.push_back(0);
          std::vector<int> temp_topology;
          std::vector<int> temp_order;
          std::vector<double> temp_branches;
          std::tie (temp_topology, temp_order, temp_branches) = findRFSubtreesRecL(tree->next->back, edgeIncidentPresent2, tasks);
          return_vector_topology.insert(return_vector_topology.end(), temp_topology.begin(), temp_topology.end());
          return_vector_order.insert(return_vector_order.end(), temp_order.begin(), temp_order.end());
          return_vector_topology.push_back(1);

          return_vector_topology.push_back(0);
          return_vector_order.push_back((intptr_t) tree->next->next->back->data);
          return_vector_topology.push_back(1);
          return_branches.push_back(tree->next->back->length);
          return_branches.insert(return_branches.end(), temp_branches.begin(), temp_branches.end());
    } else if(edgeIncidentPresent2[tree->next->next->node_index]) {
          assert(!edgeIncidentPresent2[tree->next->node_index]);

          return_vector_topology.push_back(0);
          return_vector_order.push_back((intptr_t) tree->next->back->data);
          return_vector_topology.push_back(1);

          return_vector_topology.push_back(0);
          std::vector<int> temp_topology;
          std::vector<int> temp_order;
          std::vector<double> temp_branches;
          std::tie (temp_topology, temp_order, temp_branches) = findRFSubtreesRecL(tree->next->next->back, edgeIncidentPresent2, tasks);
          return_vector_topology.insert(return_vector_topology.end(), temp_topology.begin(), temp_topology.end());
          return_vector_order.insert(return_vector_order.end(), temp_order.begin(), temp_order.end());
          return_vector_topology.push_back(1);
          return_branches.push_back(tree->next->next->back->length);
          return_branches.insert(return_branches.end(), temp_branches.begin(), temp_branches.end());

          tasks.push(tree->next->back);
    } else {
          assert(!edgeIncidentPresent2[tree->next->node_index] && !edgeIncidentPresent2[tree->next->next->node_index]);

          tasks.push(tree->next->next->back);
          tasks.push(tree->next->back);
          return_vector_topology.push_back(0);
          return_vector_order.push_back((intptr_t) tree->next->back->data);
          return_vector_topology.push_back(1);
          return_vector_topology.push_back(0);
          return_vector_order.push_back((intptr_t) tree->next->next->back->data);
          return_vector_topology.push_back(1);
    }
  }
  return std::make_tuple(return_vector_topology, return_vector_order, return_branches);
}

/**
 * If the queue of tasks is not empty, the method pops a subtree of the queue
 * and searches for a node with degree > 3. The method returns the rf-subtree
 * topology and permutation of the nodes.
 * @param tasks                queue of subtrees to check
 * @param edgeIncidentPresent2 vector indicating which edges are are present in
 * both trees (edges of the consensus tree)
 */
std::tuple<std::vector<int>, std::vector<int>, std::vector<double>> findRFSubtreesL(std::stack<pll_unode_t *> &tasks,
                      const std::vector<bool> &edgeIncidentPresent2) {
  std::vector<int> return_topology;
  std::vector<int> return_order;
  std::vector<double> return_branches;
  if(tasks.empty()) {
    return std::make_tuple(return_topology, return_order, return_branches);
  }
  pll_unode_t * tree = tasks.top();
  tasks.pop();

  if(tree == NULL) {
    return std::make_tuple(return_topology, return_order, return_branches);
  } else {
      return findRFSubtreesRecL(tree, edgeIncidentPresent2, tasks);
  }
}

void free_leaf(pll_unode_t * leaf) {
    if (leaf->label)
      free(leaf->label);
    free(leaf);
}

void free_inner_node(pll_unode_t * node) {
  if (node->label)
    free(node->label);
  free(node);
}

bool innerNodeCompare_(pll_unode_t * node1, pll_unode_t * node2) {
  return ((intptr_t) node1->data) < ((intptr_t) node2->data);
}

void consensusDiffRec(pll_utree_t * tree1, pll_utree_t * tree2, pll_unode_t * node,
              unsigned int * bl_idx, std::vector<double> &branch_lengths) {
  assert(node != NULL);

  double diff = tree2->nodes[node->pmatrix_index]->length - node->length;

  branch_lengths[*bl_idx] = diff;
  (*bl_idx)++;
  if(node->next == NULL) {
    // leaf

  } else {
    // inner node

    // assign inner nodes to an array
    std::vector<pll_unode_t*> inner_nodes;
    pll_unode_t * temp_node = node->next;
    while(temp_node->node_index != node->node_index) {
      inner_nodes.push_back(temp_node);
      assert(temp_node->next != NULL);
      temp_node = temp_node->next;
    }
    assert(temp_node->node_index == node->node_index);

    // sort the internal nodes
    std::sort (inner_nodes.begin(), inner_nodes.end(), innerNodeCompare_);


    for (size_t i = 0; i < inner_nodes.size(); i++) {
        consensusDiffRec(tree1, tree2, inner_nodes[i]->back, bl_idx, branch_lengths);
    }
  }
}

void consensusDiff(pll_utree_t * tree1, pll_utree_t * tree2, pll_unode_t * node,
            std::vector<double> &branch_lengths) {
  assert(node->next == NULL);
  assert(atoi(node->label) == 1);
  unsigned int bl_idx = 1;
  consensusDiffRec(tree1, tree2, node->back, &bl_idx, branch_lengths);
  assert(bl_idx == branch_lengths.size());
}

void nonConsensusBranchLengthsRec(pll_unode_t * tree, const std::vector<bool> node_incident, std::vector<double> &branch_lengths) {
  assert(tree != NULL);

  if(tree->next == NULL) {
    // leaf
  } else {
    // inner node
    if(node_incident[tree->node_index] == false) {
      branch_lengths.push_back(tree->length);
    }
    pll_unode_t * current_node = tree->next;
    while(current_node != tree) {
      pll_unode_t * temp_node = current_node;
      nonConsensusBranchLengthsRec(current_node->back, node_incident, branch_lengths);
      current_node = current_node->next;
    }
  }
}

/**
 * Stores all branch lengths of the branches that are not part of the consensus
 * tree in an array
 * @param root           root of the tree
 * @param node_incident  vector telling which node is incident to an edge in the consensus tree
 * @param branch_lengths return vector to write branch lengths in
 */
void nonConsensusBranchLengths(pll_unode_t * root, const std::vector<bool> node_incident, std::vector<double> &branch_lengths) {
  assert(root->back != NULL);
  assert(node_incident[root->node_index] == true);
  assert(node_incident[root->back->node_index] == true);

  nonConsensusBranchLengthsRec(root->back, node_incident, branch_lengths);
}

/**
 * Compares the sum of the element of two arrays
 * @param  a first array
 * @param  b second array
 * @return   true if sum of first array is greater than the sum of the second array
 */
bool arraySumComp(const std::vector<int> &a, const std::vector<int> &b) {
  int sum_a = 0;
  for (size_t i = 0; i < a.size(); i++) {
    sum_a += a[i];
  }
  int sum_b = 0;
  for (size_t i = 0; i < b.size(); i++) {
    sum_b += b[i];
  }
  return sum_a > sum_b ;
}

/**
 * Find the permuation to get from vector a to vector. Element i in the resulting
 * vector tells on which position in array b element a[i] is found.
 *
 * @param  a source vector
 * @param  b destination vector
 * @return   permuations to permutate vector a into vector b
 */
// Precondition: a and b must have the same length and contain the same, unique elements.
std::vector<int> findPermutation(const std::vector<int> &a, const std::vector<int> &b) {
    std::vector<int> permutation(a.size());
    for (size_t i = 0; i < a.size(); i++) {
      for (size_t j = 0; j < b.size(); j++) {
        if(a[i] == b[j]) {
          permutation[i] = j;
          break;
        }
      }
    }
    return permutation;
}

int rf_distance_compression(const char * tree1_file, const char * tree2_file,
        const char * edges_to_contract_file, const char * subtrees_succinct_file,
        const char * node_permutations_file, const char * branch_lengths_consensus_file,
        const char * branch_lengths_non_consensus_file, int flags) {

  /* tree properties */
  pll_utree_t * tree1 = NULL,
              * tree2 = NULL;
  unsigned int tip_count;

  /* parse the input trees */
  tree1 = pll_utree_parse_newick (tree1_file);
  if(tree1 == NULL) {
      // ERROR: tree could not be parsed
      // --> syntax of newick file is not correct
      // --> tree has less than 3 leaves
      return -1;
  }
  tree2 = pll_utree_parse_newick (tree2_file);
  if(tree2 == NULL) {
      // ERROR: tree could not be parsed
      // --> syntax of newick file is not correct
      // --> tree has less than 3 leaves
      return -1;
  }
  tip_count = tree1->tip_count;

  if (tip_count != tree2->tip_count) {
    // ERROR: Trees have different number of tips!
    return -1;
  }

  if (!pllmod_utree_consistency_set(tree1, tree2)) {
    // ERROR: Cannot set trees consistent!
    return -1;
  }

  if (!pllmod_utree_consistency_check(tree1, tree2)) {
    // ERROR: Tip node IDs are not consistent!
    return -1;
  }

  // search for the root of the trees
  pll_unode_t * root1 = searchRoot(tree1);
  pll_unode_t * root2 = searchRoot(tree2);

  // set and order both trees
  setTree(root1);
  setTree(root2);
  orderTree(root1);
  orderTree(root2);

  // succinct_structure stores the topology in balanced parantheses ("0=(, 1=)")
  sdsl::bit_vector succinct_structure1(4 * tip_count - 2, 0);
  // succinct_structure stores the permutation of the taxa
  sdsl::int_vector<> node_permutation1(tip_count, 0, 32);
  // branch_lengths stores all branch lengths of the tree
  // TODO: lots of potential here!
  std::vector<double> branch_lengths1(2 * tip_count - 2);
  // create mapping from node_id to branch_id
  unsigned int* node_id_to_branch_id1 = (unsigned int*) malloc ((tree1->inner_count * 3 + tree1->tip_count) * sizeof(unsigned int));
  // fill the created structures with the given tree
  assignBranchNumbers(root1, succinct_structure1, node_permutation1, branch_lengths1, node_id_to_branch_id1);

  sdsl::bit_vector succinct_structure2(4 * tip_count - 2, 0);
  sdsl::int_vector<> node_permutation2(tip_count, 0, 32);
  std::vector<double> branch_lengths2(2 * tip_count - 2);
  unsigned int* node_id_to_branch_id2 = (unsigned int*) malloc ((tree2->inner_count * 3 + tree2->tip_count) * sizeof(unsigned int));
  assignBranchNumbers(root2, succinct_structure2, node_permutation2, branch_lengths2, node_id_to_branch_id2);

  if(flags & PRINT_COMPRESSION_STRUCTURES) {
    std::cout << "Succinct representation tree 1: " << succinct_structure1 << "\n";
    std::cout << "Succinct representation tree 2: " << succinct_structure2 << "\n";
  }

  /* 1. creating the split sets manually */
  unsigned int n_splits = tip_count - 3;
  pll_unode_t ** splits_to_node1 = (pll_unode_t **) malloc(n_splits * sizeof(pll_unode_t *));
  pll_split_t * splits1 = pllmod_utree_split_create(tree1->nodes[tip_count],
                                                    tip_count,
                                                    splits_to_node1);

  if (flags & PRINT_SPLITS) {
    unsigned int i;
    for (i=0; i<n_splits; ++i)
    {
      pllmod_utree_split_show(splits1[i], tip_count);
      printf("\n");
    }
    printf("\n");
  }

  /* compute the splits, but also the nodes corresponding to each split */
  pll_unode_t ** splits_to_node2 = (pll_unode_t **) malloc(n_splits * sizeof(pll_unode_t *));
  pll_split_t * splits2 = pllmod_utree_split_create(tree2->nodes[tip_count],
                                                    tip_count,
                                                    splits_to_node2);

  if(flags & PRINT_SPLITS) {
    unsigned int i;
    for (i=0; i<n_splits; ++i)
    {
      pllmod_utree_split_show(splits2[i], tip_count);
      printf(" node: Pmatrix:%d Nodes:%d<->%d Length:%lf\n",
      splits_to_node2[i]->pmatrix_index,
      splits_to_node2[i]->node_index,
      splits_to_node2[i]->back->node_index,
      splits_to_node2[i]->length);
    }
  }

  // create arrays indicating whether a split is common in both trees or not
  int * s1_present = (int*) calloc(n_splits, sizeof(int));
  int * s2_present = (int*) calloc(n_splits, sizeof(int));

  // fill the arrays s1_present and s2_present
  int rf_distance = pllmod_utree_split_rf_distance_extended(splits1, splits2, s1_present, s2_present, tip_count);

  // vector storing if node is incident to edge in consensus tree
  // true -> is incident
  // false -> is not incident
  std::vector<bool> node_incident (3 * (tip_count - 2) + tip_count, true);

  for (unsigned int i=0; i<n_splits; ++i) {
    if (s2_present[i] != 1) {
      node_incident[splits_to_node2[i]->node_index] = false;
      node_incident[splits_to_node2[i]->back->node_index] = false;
    }
  }

  // std::vector<double> non_consensus_branch_lengths;
  // nonConsensusBranchLengths(root2, node_incident, non_consensus_branch_lengths);
  //
  // std::cout << "non_consensus_branch_lengths: ";
  // for (auto nbl: non_consensus_branch_lengths) {
  //   std::cout << nbl << " ";
  // }
  // std::cout << '\n';
  //
  // auto size_non_consensus_branch_lengths = compressBranchLengthsAndStore(non_consensus_branch_lengths, branch_lengths_non_consensus_file);

  // take minimum
  //auto size_non_consensus_branch_lengths = (size_compressed_non_consensus_branch_lengths < (non_consensus_branch_lengths.size()*8) ? size_compressed_non_consensus_branch_lengths : (non_consensus_branch_lengths.size()*8));
  //auto size_non_consensus_branch_lengths = size_compressed_non_consensus_branch_lengths;

  if(flags & PRINT_COMPRESSION_STRUCTURES) {
    std::cout << "RF-distance: " << rf_distance << "\n";
  }

  assert(rf_distance % 2 == 0);
  // create array containing all edges to contract in tree1
  sdsl::int_vector<> edges_to_contract(rf_distance / 2, 0);
  size_t idx = 0;

  // create consensus tree
  for (size_t i = 0; i < n_splits; i++) {
      if(s1_present[i] == 0) {
          // contract edge in tree 1 to get to the conseneus tree
          contractEdge(splits_to_node1[i]);
          edges_to_contract[idx] = node_id_to_branch_id1[splits_to_node1[i]->node_index];
          idx++;
      }
  }

  // sort all edges to contract
  sort(edges_to_contract.begin(), edges_to_contract.end());

  // TODO: elias-fano encoding for the ordered edges?
  sdsl::util::bit_compress(edges_to_contract);
  auto size_edges_to_contract = sdsl::size_in_bytes(edges_to_contract);

  if(flags & PRINT_COMPRESSION_STRUCTURES) {
    std::cout << "Edges to contract in tree 1: " << edges_to_contract << "\n";
    std::cout << "\tcompressed size: " << size_edges_to_contract << " bytes\n\n";
  }

  // succinct_structure stores the topology for the consensus tree in balanced parantheses ("0=(, 1=)")
  sdsl::bit_vector consensus_succinct_structure(4 * tip_count - rf_distance - 2, 0);
  sdsl::int_vector<> consensus_node_permutation(tip_count, 0, 32);
  std::vector<double> consensus_branch_lengths(2 * tip_count - 2 - (rf_distance / 2));

  unsigned int* consensus_node_id_to_branch_id = (unsigned int*) malloc ((tree1->inner_count * 3 + tree1->tip_count) * sizeof(unsigned int));
  assignBranchNumbers(root1, consensus_succinct_structure, consensus_node_permutation, consensus_branch_lengths, consensus_node_id_to_branch_id);

  std::vector<double> consensus_branch_diff_lengths(2 * tip_count - 2 - (rf_distance / 2));
  consensusDiff(tree1, tree2, root1, consensus_branch_diff_lengths);

  // TODO: later, new method
  //auto size_consensus_branch_lengths = compressBranchLengthsAndStore(consensus_branch_diff_lengths, branch_lengths_consensus_file);
  //auto size_consensus_branch_lengths = writeDoubleVector(consensus_branch_diff_lengths, branch_lengths_consensus_file);


  /*std::cout << "size consensus= " << consensus_branch_diff_lengths.size() * 8 << std::endl;
  std::cout << "complete size = " << sdsl::size_in_bytes(wt) << std::endl;
  std::cout << "bits per number = " << (sdsl::size_in_bytes(wt)*8.0)/seq.size() << std::endl;*/

  if(flags & PRINT_COMPRESSION_STRUCTURES) {
    std::cout << "Consensus tree after edge contraction: " << consensus_succinct_structure << "\n";
  }

  // create an array that indicates whether the edge incident to the node_id in tree 2 is common in both trees or not
  std::vector<bool> edgeIncidentPresent2 (tree2->inner_count * 3 + tree2->tip_count);
  for (unsigned int i=0; i<n_splits; ++i)
  {

    if(s2_present[i] == 0) {
      // edge is not common in both trees
      edgeIncidentPresent2[splits_to_node2[i]->node_index] = true;
      edgeIncidentPresent2[splits_to_node2[i]->back->node_index] = true;
    } else {
      // edge is common in both trees
      edgeIncidentPresent2[splits_to_node2[i]->node_index] = false;
      edgeIncidentPresent2[splits_to_node2[i]->back->node_index] = false;
    }
  }

  std::vector<double> branches_common_tree2 = commonBranchesOrdered(root2->back, edgeIncidentPresent2);

  std::vector<double> branches_tree2_compare = commonBranchesOrderedCompare(root2->back, root1->back, edgeIncidentPresent2);

  auto size_consensus_branch_lengths = writeDoubleVector(branches_tree2_compare, branch_lengths_consensus_file);

  std::stack<pll_unode_t *> tasks;
  tasks.push(root2->back);

  std::vector<std::vector<int>> subtrees;
  std::vector<std::vector<int>> permutations;
  std::vector<std::vector<double>> branch_lengths;

  // recursively find all subtrees that need to be inserted into the consensus tree
  while(!tasks.empty()) {
      std::vector<int> subtree;
      std::vector<int> leaf_order;
      std::vector<double> branches;

      std::tie(subtree, leaf_order, branches) = findRFSubtreesL(tasks, edgeIncidentPresent2);

      if(!subtree.empty() && leaf_order.size() > 2) {
          std::vector<int> subtree_extended;
          subtree_extended.push_back(0);
          subtree_extended.insert(subtree_extended.end(), subtree.begin(), subtree.end());
          subtree_extended.push_back(1);
          subtrees.push_back(subtree_extended);
          permutations.push_back(leaf_order);
          branch_lengths.push_back(branches);
      }
  }

    std::vector<int> permutations_index(permutations.size(), 0);
    for (int i = 0 ; i != permutations_index.size() ; i++) {
        permutations_index[i] = i;
    }

    std::sort(permutations_index.begin(), permutations_index.end(),
                                        [&](const int& a, const int& b) {
                                              return (arraySumComp(permutations[a], permutations[b]));
                                          });

    std::vector<int> permutations_index2;


    std::sort(permutations.begin(), permutations.end(), arraySumComp);

    std::vector<std::vector<int>> tree2_perms;
    traverseConsensus(root1, tree2_perms);

    if(flags & PRINT_COMPRESSION_STRUCTURES) {
      std::cout << "\nPermutations:\n";
      std::cout << "tree 2 <---> consensus tree\n";
    }

    // find corresponding permutations of nodes
    std::vector<std::vector<int>> normalized_permutations;
    int permutation_elements = 0;
    for(auto perm: tree2_perms) {
        permutation_elements += perm.size();

        size_t index = std::lower_bound(permutations.begin(), permutations.end(), perm, arraySumComp) - permutations.begin();

        assert(!arraySumComp(permutations[index], perm)); // permutation must be present

        size_t temp_index = index;
        while(!arraySumComp(permutations[temp_index], perm) && temp_index < permutations.size()) {
          if(isPermutation(permutations[temp_index], perm)) {
            if(flags & PRINT_COMPRESSION_STRUCTURES) {
              for(auto x: permutations[temp_index]) {
                  std::cout << x << " ";
              }
              std::cout << "<---> ";
              for(auto x: perm) {
                  std::cout << x << " ";
              }

              std::cout << "\t\tpermutation: ";
              for(auto x: findPermutation(permutations[temp_index], perm)) {
                  std::cout << x << " ";
              }
              std::cout <<  '\n';
            }
            normalized_permutations.push_back(findPermutation(permutations[temp_index], perm));

            permutations_index2.push_back(temp_index);

            break;
          }
          temp_index++;
        }
        assert(temp_index < permutations.size()); // no matching permutation found
    }

    size_t permutation_index = 0;
    // succinct_permutations stores all permutations according to the subtrees, split by a "0"
    // TODO: prefix sum?
    sdsl::int_vector<> succinct_permutations(permutation_elements, 0);
    for (std::vector<int> i: normalized_permutations) {
      for (auto j : i) {
        succinct_permutations[permutation_index] = j;
        permutation_index++;
      }
    }
    assert(permutation_index = succinct_permutations.size());



    if(!subtrees.empty()) {
      if(flags & PRINT_COMPRESSION_STRUCTURES) {
        std::cout << "\nSubtrees: \n";
      }
      int subtree_elements = 0;
      for (size_t i = 0; i < subtrees.size(); i++) {
        subtree_elements += subtrees[i].size();
        if(flags & PRINT_COMPRESSION_STRUCTURES) {
        }
      }

      // reorder subtrees
      std::vector<std::vector<int>> subtrees_perms_1 (subtrees.size());
      for (size_t i = 0; i < subtrees_perms_1.size(); i++) {
        subtrees_perms_1[i] = subtrees[permutations_index[i]];
      }

      std::vector<std::vector<int>> subtrees_perms_2 (subtrees.size());
      for (size_t i = 0; i < subtrees_perms_2.size(); i++) {
        subtrees_perms_2[i] = subtrees_perms_1[permutations_index2[i]];
      }
      for (size_t i = 0; i < subtrees.size(); i++) {
        if(flags & PRINT_COMPRESSION_STRUCTURES) {
          for (auto j : subtrees_perms_2[i]) {
            std::cout << j;
          }
          std::cout << "\n";
        }
      }

      // reorder branches
      std::vector<std::vector<double>> branches_perms_1 (branch_lengths.size());
      for (size_t i = 0; i < branches_perms_1.size(); i++) {
        branches_perms_1[i] = branch_lengths[permutations_index[i]];
      }

      std::vector<std::vector<double>> branches_perms_2 (branch_lengths.size());
      for (size_t i = 0; i < branches_perms_2.size(); i++) {
        branches_perms_2[i] = branches_perms_1[permutations_index2[i]];
      }



      std::vector<double> non_consensus_branch_lengths;

      for (size_t i = 0; i < branches_perms_2.size(); i++) {
        non_consensus_branch_lengths.insert(non_consensus_branch_lengths.end(),
        branches_perms_2[i].begin(), branches_perms_2[i].end());
      }
      
      // auto size_non_consensus_branch_lengths = compressBranchLengthsAndStore(non_consensus_branch_lengths, branch_lengths_non_consensus_file);
      auto size_non_consensus_branch_lengths = writeDoubleVector(non_consensus_branch_lengths, branch_lengths_non_consensus_file);


      subtrees = subtrees_perms_2;

      size_t subtrees_index = 0;
      // subtrees_succinct stores all subtrees to insert into the consensus tree
      sdsl::bit_vector subtrees_succinct(subtree_elements, 1);
      for (std::vector<int> i: subtrees) {
        for (auto j : i) {
          subtrees_succinct[subtrees_index] = j;
          subtrees_index++;
        }
      }
      assert(subtrees_index = subtrees_succinct.size());

      auto size_subtrees = sdsl::size_in_bytes(subtrees_succinct);

      if(flags & PRINT_COMPRESSION_STRUCTURES) {
        std::cout << "\nSuccinct subtree representation: " << subtrees_succinct << "\n";
        std::cout << "\tcompressed size: " << size_subtrees << " bytes\n";
      }

    // TODO: better compression?
    sdsl::util::bit_compress(succinct_permutations);
    auto size_permutations = sdsl::size_in_bytes(succinct_permutations);

    if(flags & PRINT_COMPRESSION_STRUCTURES) {
      std::cout << "\nSuccinct permutation representation: " << succinct_permutations << "\n";
      std::cout << "\tcompressed size: " << size_permutations << " bytes\n";
    }

    if(flags & PRINT_COMPRESSION) {
      std::cout << "\nRF compression size: " << size_edges_to_contract
      << " (edges to contract) + " << size_subtrees << " (subtrees) + "
      << size_permutations << " (permutations) + " << branch_lengths1.size() * 8 << " (branches) = "
      << size_edges_to_contract + size_subtrees + size_permutations + branch_lengths1.size() * 8
      << " bytes\n";

      std::cout << "\nRF compression size: " << size_edges_to_contract
      << " (edges to contract) + " << size_subtrees << " (subtrees) + "
      << size_permutations << " (permutations) + " << size_consensus_branch_lengths
      << " (consensus branches) + " << size_non_consensus_branch_lengths << " (non-consensus branches) = "
      << size_edges_to_contract + size_subtrees + size_permutations + size_consensus_branch_lengths + size_non_consensus_branch_lengths
      << " bytes\n";

      std::cout << "---------------------------------------------------------\n";
      /*std::cout << size_edges_to_contract << ";" << size_subtrees << ";"
      << size_permutations << ";" << size_consensus_branch_lengths
      << ";" << size_non_consensus_branch_lengths << ";"
      << size_edges_to_contract + size_subtrees + size_permutations
      + size_consensus_branch_lengths + size_non_consensus_branch_lengths
      << ";\n";*/

    }

    // write stuctures to file
    store_to_file(edges_to_contract, edges_to_contract_file);
    store_to_file(subtrees_succinct, subtrees_succinct_file);
    store_to_file(succinct_permutations, node_permutations_file);
  }

  // TODO: free procs segmentation fault

  //printf("RF [manual]\n");
  //printf("distance = %d\n", rf_dist);
  //printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  //printf("Amount of branchs with same lengths = %d\n", same_branchs);

  pllmod_utree_split_destroy(splits1);
  pllmod_utree_split_destroy(splits2);
  free(splits_to_node1);
  free(splits_to_node2);

  /* clean */
  //pll_utree_destroy_consensus (tree1);

  traverseTree(root1, free_leaf, free_inner_node);
  free(tree1->nodes);
  free(tree1);
  pll_utree_destroy (tree2, NULL);

  free(node_id_to_branch_id1);
  free(node_id_to_branch_id2);
  free(consensus_node_id_to_branch_id);
  free(s1_present);
  free(s2_present);

  return 0;
}
