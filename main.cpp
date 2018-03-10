#include <assert.h>
#include <stdarg.h>

#include "compress_functions.h"

/* static functions */
static void fatal (const char * format, ...);

int main (int argc, const char * argv[])
{
  if (argc != 3)
    fatal (" syntax: %s [newick] [newick]", argv[0]);

  std::cout << "Simple compression of tree 1: \n";
  simple_compression(argv[1], "output_files/succinct_tree.sdsl", "output_files/node_permutation.sdsl");
  std::cout << "\n";

  // std::cout << "Simple compression of tree 2: \n";
  // simple_compression(argv[2]);
  // std::cout << "\n";
  //
  // std::cout << "RF-compression: \n";
  // rf_distance_compression(argv[1], argv[2]);
  // std::cout << "\n";

  /*for (int i = 1; i < 1001; i++) {
      //std::cout << "RF-compression: \n";

      std::stringstream ss1;
      ss1 << "500\ \(500000\ gen\)/tree_" << i << ".nwk";

      std::stringstream ss2;
      ss2 << "500\ \(500000\ gen\)/tree_" << (i+1) << ".nwk";

      //std::cout << ss1.str() << " & " << ss2.str() << "\n";

      //simple_compression(strdup(ss1.str().c_str()));
      rf_distance_compression(strdup(ss1.str().c_str()), strdup(ss2.str().c_str()));
  }*/

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
