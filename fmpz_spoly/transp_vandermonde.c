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
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_spoly_transp_vandermonde_precomp(fmpz* bb, const fmpz* vv_inv, 
        fmpz_poly_struct * const * tree, const fmpz* tree_root,
        const fmpz* xx, slong len, const fmpz_t p)
{
    /* Evaluates the Vandermonde system V(vv)^T * xx mod p
     * len is the size of xx, vv, and bb.
     * tree and root are built from entries in vv inverted
     * (which are also provided explicitly).
     */
    fmpz_poly_struct *A; 
    fmpz_poly_t t1, t2;
    slong i, j, lenA;

    fmpz_poly_init(t1);
    fmpz_poly_init(t2);

    /* calculate sum_i( xx[i]/(1 - vv[i]*x) = -xx[i]*vv_inv[i]/(x - vv_inv[i]) ) */

    /* A is used to calculate the numerator, by following up the tree. */
    lenA = len;
    A = flint_malloc(sizeof(fmpz_poly_struct) * lenA);
    for (i = 0; i < len; ++i)
    {
        fmpz* coeff;
        fmpz_poly_init(A + i);
        fmpz_poly_set_ui(A + i, UWORD(1));
        coeff = fmpz_poly_get_coeff_ptr(A + i, 0);
        FLINT_ASSERT(coeff != NULL);

        fmpz_mul(coeff, vv_inv + i, xx + i);
        fmpz_neg(coeff, coeff);
        fmpz_mod(coeff, coeff, p);
    }

    for (i = 0; lenA > 1; ++i)
    {
        /* invariant: length of tree[i] == lenA */
        slong next_lenA = 0;

        for (j = 0; j + 1 < lenA; j += 2)
        {
            /* cross product */
            fmpz_poly_mul(t1, A + j, tree[i] + (j + 1));
            fmpz_poly_mul(t2, A + (j + 1), tree[i] + j);
            fmpz_poly_add(A + next_lenA, t1, t2);
            fmpz_poly_scalar_mod_fmpz(A + next_lenA, A + next_lenA, p);
            next_lenA++;
        }

        if (j < lenA)
        {
            /* odd case */
            fmpz_poly_swap(A + next_lenA, A + j);
            next_lenA++;
        }

        lenA = next_lenA;
    }

    /* now divide computed numerator by the product tree root */
    _fmpz_mod_poly_div_series(bb, A[0].coeffs, len, tree_root, len + 1, p, len);

    /* clean-up */
    for (i = 0; i < len; ++i)
    {
        fmpz_poly_clear(A + i);
    }
    flint_free(A);

    fmpz_poly_clear(t1);
    fmpz_poly_clear(t2);
}

void _fmpz_spoly_transp_vandermonde(fmpz* bb,
        const fmpz* vv, const fmpz* xx, slong len, const fmpz_t p)
{
    fmpz_poly_struct ** tree;
    fmpz* vv_inv;
    fmpz_poly_t root;
    slong ind;

    if (len == 0) return;
    else if (len == 1) 
    {
        fmpz_mod(bb + 0, xx + 0, p);
        return;
    }

    vv_inv = _fmpz_vec_init(len);
    for (ind = 0; ind < len; ++ind)
    {
        fmpz_invmod(vv_inv + ind, vv + ind, p);
    }

    tree = _fmpz_mod_poly_tree_alloc(len);
    _fmpz_mod_poly_tree_build(tree, vv_inv, len, p);

    fmpz_poly_init2(root, len + 1);
    ind = FLINT_CLOG2(len) - 1;
    FLINT_ASSERT(len ==
            fmpz_poly_degree(tree[ind] + 0) + fmpz_poly_degree(tree[ind] + 1));
    fmpz_poly_mul(root, tree[ind] + 0, tree[ind] + 1);
    fmpz_poly_scalar_mod_fmpz(root, root, p);

    _fmpz_spoly_transp_vandermonde_precomp(bb, vv_inv, tree, root->coeffs, xx, len, p);

    _fmpz_mod_poly_tree_free(tree, len);
    _fmpz_vec_clear(vv_inv, len);
    fmpz_poly_clear(root);
}
