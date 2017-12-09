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

#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_transp_vandermonde_inv_precomp(fmpz* xx,
        fmpz_poly_struct * const * tree, const fmpz* tree_root,
        const fmpz* bb, slong len, const fmpz_t p)
{
    /* Solves the Vandermonde system V(vv)^T * xx = bb mod p
     * for the vector xx.
     * len is the size of xx, vv, and bb.
     * tree and root are built from entries in vv (not given explicitly).
     * Entries in xx are reduced in the symmetric range [-p/2 .. p/2].
     */
    fmpz *D, *Q;
    slong i;

    if (len <= 0) return;

    D = _fmpz_vec_init(3*len);
    Q = D + len;

    /* D = reversal of bb */
    for (i=0; i<len; ++i) fmpz_set(D + i, bb + (len-i-1));

    /* Q = root*D / x^len. TODO: make this faster using mullow? */
    _fmpz_mod_poly_mul(Q, tree_root, len+1, D, len, p);
    Q = Q + len;

    /* xx = evals of Q at points in vv */
    _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(xx,
            Q, len, tree, len, p);

    /* overwrite D with evals of (d poly/dx) at points in vv */
    _fmpz_mod_poly_derivative(Q, tree_root, len+1, p);
    _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(D,
            Q, len, tree, len, p);

    /* pairwise divide the two sets of evaluations */
    for (i=0; i<len; ++i)
    {
        fmpz_invmod(D+i, D+i, p);
        fmpz_mul(xx+i, xx+i, D+i);
    }
    _fmpz_vec_scalar_smod_fmpz(xx, xx, len, p);

    _fmpz_vec_clear(D, 3*len);
}

void _fmpz_mod_poly_transp_vandermonde_inv(fmpz* xx,
        const fmpz* vv, const fmpz* bb, slong len, const fmpz_t p)
{
    fmpz_poly_struct ** tree;
    fmpz_poly_t root;
    slong ind;

    if (len == 0)
        return;
    else if (len == 1)
    {
        fmpz_mods(xx + 0, bb + 0, p);
        return;
    }

    tree = _fmpz_mod_poly_tree_alloc(len);
    _fmpz_mod_poly_tree_build(tree, vv, len, p);

    fmpz_poly_init2(root, len + 1);
    ind = FLINT_CLOG2(len) - 1;
    FLINT_ASSERT(len ==
            fmpz_poly_degree(tree[ind] + 0) + fmpz_poly_degree(tree[ind] + 1));
    fmpz_poly_mul(root, tree[ind] + 0, tree[ind] + 1);
    fmpz_poly_scalar_mod_fmpz(root, root, p);

    _fmpz_mod_poly_transp_vandermonde_inv_precomp(xx, tree, root->coeffs, bb, len, p);

    _fmpz_mod_poly_tree_free(tree, len);
    fmpz_poly_clear(root);
}
