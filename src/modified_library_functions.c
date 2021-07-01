#include "modified_library_functions.h"

// copied because "compare_splits" uses it internally
/*
unsigned int bitv_length(unsigned int bit_count)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = bit_count % split_size;

  return bit_count / split_size + (split_offset>0);
}
*/

/*
 * The order of the splits is not really significant, as long as the two
 * following agree.
 *
 * _cmp_splits is used for sorting.
 * compare_splits is used for comparing splits from different trees
 */
int compare_splits (pll_split_t s1,
                           pll_split_t s2,
                           unsigned int split_len)
{
  unsigned int i;

  for (i=0; i<split_len; ++i)
  {
    if (s1[i] != s2[i])
      return (int) (s1[i] > s2[i]?1:-1);
  }
  return 0;
}

/*
 * Precondition: splits must be normalized and sorted!
 */
unsigned int pllmod_utree_split_rf_distance_extended(pll_split_t * s1,
                                                       pll_split_t * s2,
                                                       int * s1_present,
                                                       int * s2_present,
                                                       unsigned int tip_count)
{
  unsigned int split_count = tip_count - 3;
  unsigned int split_len   = bitv_length(tip_count);
  unsigned int equal = 0;
  unsigned int s1_idx = 0,
               s2_idx = 0;

  for (s1_idx=0; s1_idx<split_count && s2_idx<split_count; ++s1_idx)
  {
    int cmp = compare_splits(s1[s1_idx], s2[s2_idx], split_len);
    if (!cmp)
    {
      equal++;
      s1_present[s1_idx] = 1;
      s2_present[s2_idx] = 1;
      //printf("equal: %i, s1_idx: %i, s2_idx: %i\n", equal, s1_idx, s2_idx);
      s2_idx++;
    }
    else
    {
      if (cmp > 0)
      {
        while(++s2_idx < split_count &&
              (cmp = compare_splits(s1[s1_idx], s2[s2_idx], split_len)) > 0);
        if (!cmp)
        {
           equal++;
           s1_present[s1_idx] = 1;
           s2_present[s2_idx] = 1;
           //printf("equal: %i, s1_idx: %i, s2_idx: %i\n", equal, s1_idx, s2_idx);
           //s2_idx++;
        }
      }
    }
  }

  assert(equal <= (tip_count-3));
  //printf("tip_count: %i, equal: %i\n", tip_count, equal);
  return 2*(tip_count - 3 - equal);
}
