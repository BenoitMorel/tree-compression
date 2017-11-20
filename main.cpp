#include <assert.h>
#include <stdarg.h>

#include "compress_functions.h"

/* static functions */
static void fatal (const char * format, ...);

int main (int argc, char * argv[])
{
  if (argc != 3)
    fatal (" syntax: %s [newick] [newick]", argv[0]);

  std::cout << "Simple compression of tree 1: \n";
  simple_compression(argv[1]);
  std::cout << "\n";

  std::cout << "Simple compression of tree 2: \n";
  simple_compression(argv[2]);
  std::cout << "\n";

  std::cout << "RF-compression: \n";
  rf_distance_compression(argv[1], argv[2]);
  std::cout << "\n";

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
