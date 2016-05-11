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

        Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_spoly.h"
#include "fmpz_mod_poly.h"

void _fmpz_spoly_transp_vandermonde_precomp(fmpz* bb,
        fmpz_poly_struct * const * tree, const fmpz* tree_root,
        const fmpz* xx, slong len, const fmpz_t p)
{
    /* Evaluates the Vandermonde system V(vv)^T * xx mod p
     * len is the size of xx, vv, and bb.
     * tree and root are built from entries in vv (not given explicitly).
     */
    /* fmpz_poly_struct *A; 
    fmpz_poly_t t1, t2;
    slong i, j, k; */

    /* calculate sum_i( -bb[i]/vv[i] prod_{i != j}( x - vv[j] */

    /* TODO FIXME */
}

void _fmpz_spoly_transp_vandermonde(fmpz* bb,
        const fmpz* vv, const fmpz* xx, slong len, const fmpz_t p)
{
    fmpz_poly_struct ** tree;
    fmpz_poly_t root;
    slong ind;

    if (len == 0) return;
    else if (len == 1) 
    {
        fmpz_set(bb + 0, vv + 0);
        return;
    }

    tree = _fmpz_mod_poly_tree_alloc(len);
    _fmpz_mod_poly_tree_build(tree, vv, len, p);

    fmpz_poly_init2(root, len + 1);
    ind = FLINT_CLOG2(len);
    FLINT_ASSERT(len ==
            fmpz_poly_degree(tree[ind] + 0) + fmpz_poly_degree(tree[ind] + 1));
    fmpz_poly_mul(root, tree[ind] + 0, tree[ind] + 1);

    _fmpz_spoly_transp_vandermonde_precomp(bb, tree, root->coeffs, xx, len, p);

    _fmpz_mod_poly_tree_free(tree, len);
}
