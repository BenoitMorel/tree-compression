#include "compress_functions.h"

/* set to 1 for printing splits */
#define PRINT_SPLITS 1

/* static functions */
static void fatal (const char * format, ...);

void traverseTreeRec(pll_unode_t * tree, const std::vector<bool> &edgeIncidentPresent2, std::queue<pll_unode_t *> &tasks) {
  assert(tree != NULL);
  if(tree->next == NULL) {
    // leaf; edge incident is always present in both trees
    //std::cout << "Leaf: " << tree->node_index << "\t    Edge incident present: " << edgeIncidentPresent2[tree->node_index] << "\n";
  } else {
    // inner node
    assert(tree->next != NULL);
    assert(tree->next->next != NULL);
    assert(tree->next->back != NULL);
    assert(tree->next->next->back != NULL);

    //std::cout << "Inner node: " << tree->node_index << "\t    Edge incident present: " << edgeIncidentPresent2[tree->node_index] << "\n";

    //std::cout << "\n1. Inner node: " << tree->next->node_index <<
    //    "\t    Edge incident present: " << edgeIncidentPresent2[tree->next->node_index] << "\n";
    if(edgeIncidentPresent2[tree->next->node_index]) {
          std::cout << "0";
          traverseTreeRec(tree->next->back, edgeIncidentPresent2, tasks);
          std::cout << "1";
    } else {
          tasks.push(tree->next->back);
    }
    //std::cout << "\n2. Inner node: " << tree->next->next->node_index <<
    //    "\t    Edge incident present: " << edgeIncidentPresent2[tree->next->next->node_index] << "\n";
    if(edgeIncidentPresent2[tree->next->next->node_index]) {
          std::cout << "0";
          traverseTreeRec(tree->next->next->back, edgeIncidentPresent2, tasks);
          std::cout << "1";
    } else {
          tasks.push(tree->next->next->back);
    }
  }
}

void traverseTree(std::queue<pll_unode_t *> &tasks, const std::vector<bool> &edgeIncidentPresent2) {
  if(tasks.empty()) {
    return;
  }
  pll_unode_t * tree = tasks.front();
  tasks.pop();

  if(tree == NULL) {
    return;
  } else {
      traverseTreeRec(tree, edgeIncidentPresent2, tasks);
  }
}

void simple_compression(char * tree_file) {
  /* tree properties */
  pll_utree_t * tree = NULL;
  unsigned int tip_count;

  /* parse the input tree */
  tree = pll_utree_parse_newick (tree_file);
  tip_count = tree->tip_count;

  // search for the root of the tree
  pll_unode_t * root = searchRoot(tree);
  // set the tree
  setTree(root);
  // create a mapping from node_ids in tree to branch numbers
  sdsl::bit_vector succinct_structure(4 * tip_count - 2, 0);
  sdsl::int_vector<> node_permutation(tip_count, 0, 32);
  std::vector<double> branch_lengths(2 * tip_count - 2);
  unsigned int* node_id_to_branch_id = (unsigned int*) malloc ((tree->inner_count * 3 + tree->tip_count) * sizeof(unsigned int));
  assignBranchNumbers(root, succinct_structure, node_permutation, branch_lengths, node_id_to_branch_id);

  std::cout << "Succinct representation: " << succinct_structure << "\n";
  std::cout << "\tuncompressed size: " << sdsl::size_in_bytes(succinct_structure) << " bytes\n";
  sdsl::util::bit_compress(succinct_structure);
  std::cout << "\tcompressed size: " << sdsl::size_in_bytes(succinct_structure) << " bytes\n";

  std::cout << "Node permutation: " << node_permutation << "\n";
  std::cout << "\tuncompressed size: " << sdsl::size_in_bytes(node_permutation) << " bytes\n";
  sdsl::util::bit_compress(node_permutation);
  std::cout << "\tcompressed size: " << sdsl::size_in_bytes(node_permutation) << " bytes\n";

  std::cout << "Branch Lengths: ";
  for(auto x: branch_lengths) {
    std::cout << x << " ";
  }
  std::cout << "\n\n";
}


void rf_distance_compression(char * tree1_file, char * tree2_file) {
  /* tree properties */
  pll_utree_t * tree1 = NULL,
              * tree2 = NULL;
  unsigned int tip_count;

  /* parse the input trees */
  tree1 = pll_utree_parse_newick (tree1_file);
  tree2 = pll_utree_parse_newick (tree2_file);
  tip_count = tree1->tip_count;

  if (tip_count != tree2->tip_count)
    fatal("Trees have different number of tips!");

  if (!pllmod_utree_consistency_set(tree1, tree2))
    fatal("Cannot set trees consistent!");

  if (!pllmod_utree_consistency_check(tree1, tree2))
    fatal("Tip node IDs are not consistent!");


  // search for the root of the trees
  pll_unode_t * root1 = searchRoot(tree1);
  pll_unode_t * root2 = searchRoot(tree2);

  // order both trees
  setTree(root1);
  setTree(root2);

  // create a mapping from node_ids in tree1 to branch numbers
  sdsl::bit_vector succinct_structure1(4 * tip_count - 2, 0);
  size_t zeros = sdsl::rank_support_v<0>(&succinct_structure1)(succinct_structure1.size());
  sdsl::bit_vector::select_0_type b_sel(&succinct_structure1);
  //for (size_t i=1; i <= 2 * tip_count - 1; ++i)
  // std::cout << i << ": " << b_sel(i) << "\n";

  sdsl::int_vector<> node_permutation1(tip_count, 0, 32);
  std::vector<double> branch_lengths1(2 * tip_count - 2);
  unsigned int* node_id_to_branch_id1 = (unsigned int*) malloc ((tree1->inner_count * 3 + tree1->tip_count) * sizeof(unsigned int));
  assignBranchNumbers(root1, succinct_structure1, node_permutation1, branch_lengths1, node_id_to_branch_id1);

  sdsl::bit_vector succinct_structure2(4 * tip_count - 2, 0);
  sdsl::int_vector<> node_permutation2(tip_count, 0, 32);
  std::vector<double> branch_lengths2(2 * tip_count - 2);
  unsigned int* node_id_to_branch_id2 = (unsigned int*) malloc ((tree2->inner_count * 3 + tree2->tip_count) * sizeof(unsigned int));
  assignBranchNumbers(root2, succinct_structure2, node_permutation2, branch_lengths2, node_id_to_branch_id2);


  std::cout << "Succinct representation tree 1: " << succinct_structure1 << "\n";

  std::cout << "Succinct representation tree 2: " << succinct_structure2 << "\n";

  //std::cout << "Node permutation tree 1: " << node_permutation1 << "\n";
  /*std::cout << "Branch Lengths tree 1: ";
  for(auto x: branch_lengths1) {
    std::cout << x << " ";
  }
  std::cout << "\n";*/

  /*for (size_t i = 0; i < (tree1->inner_count * 3 + tree1->tip_count); i++) {
      printf("Index: %i, Edge: %u\n", i, node_id_to_branch_id[i]);
  }*/


  /* uncomment lines below for displaying the trees in ASCII format */
  // pll_utree_show_ascii(tree1, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);
  // pll_utree_show_ascii(tree2, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

  /* next, we compute the RF distance in 2 different ways: */

  /* 1. creating the split sets manually */
  unsigned int n_splits = tip_count - 3;
  pll_unode_t ** splits_to_node1 = (pll_unode_t **) malloc(n_splits * sizeof(pll_unode_t *));
  pll_split_t * splits1 = pllmod_utree_split_create(tree1->nodes[tip_count],
                                                    tip_count,
                                                    splits_to_node1);

#if(PRINT_SPLITS)
  {
    unsigned int i;
    for (i=0; i<n_splits; ++i)
    {
      pllmod_utree_split_show(splits1[i], tip_count);
      printf("\n");
    }
    printf("\n");
  }
#endif

  /* now we compute the splits, but also the nodes corresponding to each split */
  pll_unode_t ** splits_to_node2 = (pll_unode_t **) malloc(n_splits * sizeof(pll_unode_t *));
  pll_split_t * splits2 = pllmod_utree_split_create(tree2->nodes[tip_count],
                                                    tip_count,
                                                    splits_to_node2);


#if(PRINT_SPLITS)
  {
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
#endif

  // create arrays indicating whether a split is common in both trees or not
  int * s1_present = (int*) calloc(n_splits, sizeof(int));
  int * s2_present = (int*) calloc(n_splits, sizeof(int));

  // fill the arrays s1_present and s2_present
  int rf_distance = pllmod_utree_split_rf_distance_extended(splits1, splits2, s1_present, s2_present, tip_count);
  std::cout << "RF-distance: " << rf_distance << "\n";

  // determine number of branches with the same length
  unsigned int same_branchs = same_branch_lengths(splits1, splits2, splits_to_node1, splits_to_node2, tip_count);

  assert(rf_distance % 2 == 0);
  sdsl::int_vector<> edges_to_contract(rf_distance / 2, 0);
  size_t idx = 0;

  std::cout << "\ns1_present\n";
  for (size_t i = 0; i < n_splits; i++) {
      std::cout << i << ":\t" << s1_present[i] << "\t" << node_id_to_branch_id1[splits_to_node1[i]->node_index] << "\n";
  }
  std::cout << "\ns2_present\n";
  for (size_t i = 0; i < n_splits; i++) {
      std::cout << i << ":\t" << s2_present[i] << "\t" << node_id_to_branch_id2[splits_to_node2[i]->node_index] << "\n";
  }


  for (size_t i = 0; i < n_splits; i++) {
      if(s1_present[i] == 0) {
          contractEdge(splits_to_node1[i]);
          edges_to_contract[idx] = node_id_to_branch_id1[splits_to_node1[i]->node_index];
          idx++;
      }
  }

  sort(edges_to_contract.begin(), edges_to_contract.end());
  std::cout << "Edges to contract in tree 1: " << edges_to_contract;
  printf("\n\n");

  // create a mapping from node_ids in tree1 to branch numbers
  sdsl::bit_vector consensus_succinct_structure1(4 * tip_count - rf_distance - 2, 0);
  size_t consensus_zeros = sdsl::rank_support_v<0>(&consensus_succinct_structure1)(consensus_succinct_structure1.size());
  sdsl::bit_vector::select_0_type consensus_b_sel(&consensus_succinct_structure1);
  //for (size_t i=1; i <= 2 * tip_count - 1; ++i)
  // std::cout << i << ": " << b_sel(i) << "\n";

  sdsl::int_vector<> consensus_node_permutation1(tip_count, 0, 32);
  std::vector<double> consensus_branch_lengths1(2 * tip_count - 2 - (rf_distance / 2));
  unsigned int* consensus_node_id_to_branch_id1 = (unsigned int*) malloc ((tree1->inner_count * 3 + tree1->tip_count) * sizeof(unsigned int));
  assignBranchNumbers(root1, consensus_succinct_structure1, consensus_node_permutation1, consensus_branch_lengths1, consensus_node_id_to_branch_id1);

  std::cout << "Consensus tree after edge contraction: " << consensus_succinct_structure1 << "\n";




  // create an array that indicates whether the edge incident to the node_id in tree 2 is common
  // in both trees or not
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

  /*std::cout << "edgeIncidentPresent2:\n";
  for(int i = 0; i < edgeIncidentPresent2.size(); i++) {
    std::cout << i << ": " << edgeIncidentPresent2[i] << "\n";
  }*/

  std::queue<pll_unode_t *> tasks;
  tasks.push(root2->back);

  //printTree(root2);
  while(!tasks.empty()) {
      std::cout << " ";
      traverseTree(tasks, edgeIncidentPresent2);
  }


  //printf("RF [manual]\n");
  //printf("distance = %d\n", rf_dist);
  //printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  //printf("Amount of branchs with same lengths = %d\n", same_branchs);

  pllmod_utree_split_destroy(splits1);
  pllmod_utree_split_destroy(splits2);
  free(splits_to_node1);
  free(splits_to_node2);

  /* clean */
  pll_utree_destroy (tree1, NULL);
  pll_utree_destroy (tree2, NULL);

  free(node_id_to_branch_id1);
  free(s1_present);
  free(s2_present);

}



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void fatal (const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf (stderr, format, argptr);
  va_end(argptr);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}
