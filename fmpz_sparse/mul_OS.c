/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
     
    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Authored 2015 by A. Whitman Groves; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"
#include "fmpz_vec.h"
#include "math.h"

void
fmpz_sparse_mul_OS(fmpz_sparse_t res, flint_rand_t state, const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
  fmpz * test;
  slong length;
  fmpz_t f_height, g_height, C, p;

  fmpz_init(f_height);
  fmpz_init(g_height);
  fmpz_init(f_deg);
  fmpz_init(g_deg);
  fmpz_init(C);
  fmpz_init(p);

  test = NULL;

  flint_printf("\nthere is no god\n");
  length = fmpz_sparse_sumset(&test, state, poly1, poly2);

  /*SPARSE_MUL_COEFFS*/
  
  /*C = F-height*g_height*MAX_DEGREE(f, g)*/
  fmpz_sparse_get_height(f_height, poly1);
  fmpz_sparse_get_height(g_height, poly2);

  fmpz_sparse_get_degree(f_deg, poly1);
  fmpz_sparse_get_degree(g_deg, poly2);

  if(fmpz_cmp(f_deg, g_deg) > 0)
    fmpz_set(C, f_deg);
  else
    fmpz_set(C, g_deg);

  fmpz_mul(C, f_height, g_height);

  fmpz_van_prime(p, state, length, fmpz_bits(C), .25);  
    /*GET PRIM ROOTS*/

    /*CHINESE REMAINDER THEOREM*/

  fmpz_sparse_zero(res);

  fmpz_clear(f_height);
  fmpz_clear(g_height);
  fmpz_clear(f_deg);
  fmpz_clear(g_deg);
  fmpz_clear(C);
}
