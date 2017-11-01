/*
 Copyright (C) 2016 Diego Darriba

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */
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

  std::cout << "RF-compression with branch lengths: \n";
  rf_distance_compression_with_branches(argv[1], argv[2]);
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
