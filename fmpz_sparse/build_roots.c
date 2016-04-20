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

    Authored 2015 by A. Whitman Groves and Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"

void 
_fmpz_mod_poly_build_roots(fmpz_mod_poly_t res, const fmpz * roots, slong len)
{
  fmpz_poly_struct ** tree;
  slong i, j;

  tree = _fmpz_mod_poly_tree_alloc(len);

  _fmpz_mod_poly_tree_build(tree, roots, len, &(res->p));

  for(i = 0, j = WORD(2); j < len; i++, j*=2);

  j = j/2 + 1;

  fmpz_mod_poly_fit_length(res, len + 1);

  _fmpz_mod_poly_mul(res->coeffs, tree[i][0].coeffs, j, tree[i][1].coeffs, len + 2 - j, &(res->p));

  _fmpz_mod_poly_set_length(res, len+1);
}
