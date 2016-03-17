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

#include "fmpz_sparse.h"
#include "fq.h"
#include "fq_mat.h"
#include "long_extras.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_powmod_x_2exp(fmpz* res, 
        const fmpz* poly, slong len, ulong k, const fmpz_t p)
{
    /* Computes x^(2^k) mod poly.
     * Assumes res has size len-1 and poly is monic.
     */
    fmpz *T, *Q;
    const fmpz *one = poly + (len-1);
    ulong i;
    slong lenT, lenQ;

    FLINT_ASSERT(len >= 1);
    FLINT_ASSERT(fmpz_is_one(one));

    if (len == 1) return;
    else if (len == 2)
    {
        fmpz_t e;
        fmpz_init(e);
        fmpz_setbit(e, k);

        fmpz_neg(res+0, poly+0);
        fmpz_powm(res+0, res+0, e, p);

        fmpz_clear(e);
        return;
    }

    lenT = 2 * len - 3;
    lenQ = len - 2;

    T = _fmpz_vec_init(lenT + lenQ);
    Q = T + lenT;

    i = FLINT_MIN(k, z_sizeinbase(lenQ, 2));
    FLINT_ASSERT((lenT-1) >> i >= WORD(1)); /* 2^i <= lenT-1 */ 
    fmpz_one(T + (UWORD(1) << i));
    _fmpz_mod_poly_divrem(Q, res, T, lenT, poly, len, one, p);
    /* invariant: res = x^(2^i) mod f */

    for (; i < k; ++i)
    {
        _fmpz_mod_poly_sqr(T, res, len-1, p);
        _fmpz_mod_poly_divrem(Q, res, T, lenT, poly, len, one, p);
    }

    _fmpz_vec_clear(T, lenT + lenQ);
}


slong _fmpz_mod_poly_binary_roots(fmpz* roots, fmpz* expons, 
        const fmpz* poly, slong len, const fmpz_t theta, slong k, const fmpz_t p)
{
    /* Computes the roots of the given polynomial mod p, under certain
     * assumptions:
     * (1) poly is monic
     * (1) poly splits completely into distinct degree-1 factors modulo p
     * (2) every root is a power of the given field element theta
     * (3) theta has order 2^k in Z[p]
     *
     * The returned value is the number of distinct roots found.
     * Each entry in roots is one of the roots of the polynomial,
     * and the corresponding entry in expons gives the power of theta
     * that root is equal to.
     *
     * roots and expons should already be allocated to size at least
     * len-1 each.
     */
    slong i, nroots;
    const fmpz *one = poly + (len-1);
    
    if (len <= 1) return 0;
    FLINT_ASSERT(fmpz_is_one(one));
    FLINT_ASSERT(k >= 0);
    
    if (k < FLINT_BITS && (len-1) >> k > 0)
    {
        /* degree exceeds 2^k, so every root must be present. */
        nroots = WORD(1) << k;
        FLINT_ASSERT(len == nroots + 1);

        fmpz_one(roots+0);
        fmpz_zero(expons+0);
        for (i=1; i<nroots; ++i)
        {
            fmpz_mul(roots+i, roots+(i-1), theta);
            fmpz_mod(roots+i, roots+i, p);
            fmpz_set_si(expons+i, i);
        }
    }
    else
    {
        fmpz *g, *h;
        fmpz_t t2, temp;
        slong leng = len-1, lenh = len, oddroots;

        g = _fmpz_vec_init(leng + lenh);
        h = g + leng;

        fmpz_init(temp);
        fmpz_init(t2);
        /* t2 = theta^2 mod p */
        fmpz_mul(t2, theta, theta);
        fmpz_mod(t2, t2, p);

        /* compute h = x^(2^(k-1)) - 1 mod poly */
        _fmpz_mod_poly_powmod_x_2exp(h, poly, len, k-1, p);
        if (fmpz_is_zero(h+0)) fmpz_set(h+0, p);
        fmpz_sub_ui(h+0, h+0, UWORD(1));

        FMPZ_VEC_NORM(h, lenh);
        if (lenh == 0)
        {
            /* all root powers are even; only need one recursive call. */
            nroots = _fmpz_mod_poly_binary_roots(roots, expons, poly, len, t2, k-1, p);
            _fmpz_vec_scalar_mul_2exp(expons, expons, nroots, UWORD(1));
        }
        else
        {
            /* split off even roots into g */
            fmpz_invmod(temp, h + (lenh-1), p);
            leng = _fmpz_mod_poly_gcd(g, poly, len, h, lenh, temp, p);
            FLINT_ASSERT(leng >= 1);

            if (leng == 1)
            {
                /* all root powers are odd; only need one recursive call */
                nroots = 0;

                /* set h = poly */
                h = g + leng;
                _fmpz_vec_set(h, poly, len);
                lenh = len;
            }
            else
            {
                /* make g monic */
                fmpz_invmod(temp, g + (leng-1), p);
                _fmpz_mod_poly_scalar_mul_fmpz(g, g, leng, temp, p);

                /* recursive call on even roots */
                nroots = _fmpz_mod_poly_binary_roots(roots, expons, g, leng, t2, k-1, p);
                _fmpz_vec_scalar_mul_2exp(expons, expons, nroots, UWORD(1));

                /* set h = poly/g */
                h = g + leng;
                lenh = len - leng + 1;
                _fmpz_mod_poly_divrem(h, h+lenh, poly, len, g, leng, one, p);
            }

            /* split off odd roots into h.
             * set h = h(theta*x) / theta^(deg h)
             */

            fmpz_set(temp, theta);
            i = 1;

            for (i=1; i < lenh-1; ++i)
            {
                fmpz_mul(h+i, h+i, temp);
                fmpz_mod(h+i, h+i, p);
                fmpz_mul(temp, temp, theta);
            }

            fmpz_invmod(temp, temp, p);
            _fmpz_mod_poly_scalar_mul_fmpz(h, h, lenh-1, temp, p);

            /* recursive call on odd roots */
            oddroots = _fmpz_mod_poly_binary_roots(roots+nroots, expons+nroots, 
                    h, lenh, t2, k-1, p);
            _fmpz_mod_poly_scalar_mul_fmpz(roots+nroots, roots+nroots, oddroots, theta, p);

            for (i=nroots; i < nroots+oddroots; ++i)
            {
                fmpz_mul_2exp(expons+i, expons+i, 1);
                fmpz_add_ui(expons+i, expons+i, UWORD(1));
            }
            nroots += oddroots;
        }

        _fmpz_vec_clear(g, 2*len-1);
        fmpz_clear(t2);
        fmpz_clear(temp);
    }

    FLINT_ASSERT(nroots == len-1);
    return nroots;
}

void _fmpz_mod_poly_transposed_vandermonde(fmpz* xx,
        const fmpz* vv, const fmpz* bb, slong len, const fmpz* poly, const fmpz_t p)
{
    /* Solves the Vandermonde system V(vv)^T * xx = bb mod p
     * for the vector xx.
     * len is the size of xx, vv, and bb.
     * (poly, len+1) is a polynomial whose roots are the vv entries.
     * Entries in xx are reduced in the symmetric range [-p/2 .. p/2].
     */
    fmpz *D, *Q;
    slong i;

    if (len <= 0) return;

    D = _fmpz_vec_init(3*len);
    Q = D + len;

    /* D = reversal of bb */
    for (i=0; i<len; ++i) fmpz_set(D + i, bb + (len-i-1));

    /* Q = poly*D / x^len. TODO: make this faster using mullow? */
    _fmpz_mod_poly_mul(Q, poly, len+1, D, len, p);
    Q = Q + len;

    /* xx = evals of Q at points in vv */
    _fmpz_mod_poly_evaluate_fmpz_vec(xx, Q, len, vv, len, p);

    /* overwrite D with evals of (d poly/dx) at points in vv */
    /* TODO: re-use subproduct tree to make this faster? */
    _fmpz_mod_poly_derivative(Q, poly, len+1, p);
    _fmpz_mod_poly_evaluate_fmpz_vec(D, Q, len, vv, len, p);

    /* pairwise divide the two sets of evaluations */
    for (i=0; i<len; ++i)
    {
        fmpz_invmod(D+i, D+i, p);
        fmpz_mul(xx+i, xx+i, D+i);
    }
    _fmpz_vec_scalar_smod_fmpz(xx, xx, len, p);

    _fmpz_vec_clear(D, 3*len);
}

int fmpz_sparse_bp_interp(fmpz_sparse_t res,
        const fmpz_sparse_bp_interp_t evals)
{
    /* Berlekamp-Massey to find the Prony polynomial
      * Compute roots of Prony poly, along with discrete logs
      * Compute coefficients from transposed Vandermode
      */

    fmpz_mod_poly_t G;
    fmpz * roots;
    slong i, t;
    const fmpz * w = evals->sample_points + 1;

    fmpz_sparse_zero(res);

    if (evals->length == 0) 
    {
        return 1;
    }

    /* Berlekamp-Massey to discover Prony polynomial */

    fmpz_mod_poly_init(G, evals->q);
    fmpz_mod_poly_minpoly(G, evals->evaluations, evals->length);

    t = fmpz_mod_poly_degree(G);
    if (t > evals->length / 2)
    {
        /* sparsity estimate was too low */
        fmpz_mod_poly_clear(G);
        return 0;
    }

    _fmpz_sparse_reserve(res, t);
    for (i=res->length; i<t; ++i)
    {
        fmpz_init(res->coeffs + i);
        fmpz_init(res->expons + i);
    }
    _fmpz_sparse_set_length(res, t);

    /* find roots of Prony polynomial, and their orders */

    roots = _fmpz_vec_init(t);
    _fmpz_mod_poly_binary_roots(roots, res->expons, 
            G->coeffs, G->length, w, evals->log2_order, evals->q);

    /* solve transposed Vandermode to get coeffs */
    /* Varndermonde(roots)^T * x = evals, truncated to length t */

    _fmpz_mod_poly_transposed_vandermonde(res->coeffs,
            roots, evals->evaluations, t, G->coeffs, evals->q);

    /* sort terms and remove zero coeffs */
    _fmpz_sparse_normalise(res);

    /* clean-up */
    _fmpz_vec_clear(roots, t);
    fmpz_mod_poly_clear(G);

    return 1;
}
