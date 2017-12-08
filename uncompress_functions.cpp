#include "uncompress_functions.h"

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

     // leaf is reached
     pll_unode_t * new_leaf = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
     if (!new_leaf) {
      // error
        assert(false);
        return NULL;
     }

     new_leaf->next = NULL;

     char * label = (char *) malloc(10);
     int label_int = node_permutation[*node_idx];
     sprintf(label, "%d", label_int);
     new_leaf->label = label;
     (*node_idx)++;
     new_leaf->data = NULL;
     new_leaf->length = 0;
     new_leaf->clv_index = 0;
     new_leaf->scaler_index = 0;

     *succinct_idx = *succinct_idx + 2;

     return new_leaf;
   }
}

void simple_uncompression(sdsl::bit_vector &succinct_structure, sdsl::int_vector<> &node_permutation) {
  unsigned int succinct_idx = 0;
  unsigned int node_idx = 0;
  pll_unode_t * tree = createTreeRec(succinct_structure, &succinct_idx, node_permutation, &node_idx);
  assert(succinct_idx == succinct_structure.size());
  assert(node_idx == node_permutation.size());

  assert(atoi(tree->next->back->label) == 1);

  pll_unode_t * root = tree->next->back;
  root->back = tree->next->next->back;
  tree->next->next->back->back = root;

  printTree(root);
}
