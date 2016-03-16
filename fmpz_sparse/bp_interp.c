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

int OHNO = 0; /* FIXME */

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
    
if (OHNO && len <= 1)
{ /* FIXME */
flint_printf("hit base case k=%wd\n", k);
fflush(stdout);
}
    if (len <= 1) return 0;
    FLINT_ASSERT(fmpz_is_one(one));
    FLINT_ASSERT(k >= 0);
    
if (OHNO) {
flint_printf("\nk=%wd, theta=", k);
fmpz_print(theta); flint_printf("\n");
_fmpz_mod_poly_print(poly, len, p); flint_printf("\n");
fflush(stdout); /* FIXME */
}
    if ((len-1) >> k > 0)
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
if (OHNO) { 
flint_printf("**even roots only k=%wd**\n", k);
fflush(stdout); /* FIXME */
}
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

if (OHNO) { 
flint_printf("**odd roots only k=%wd**\n", k);
fflush(stdout); /* FIXME */
}
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

if (OHNO) { 
flint_printf("**even roots REC CALL k=%wd**\n", k);
flint_printf("leng=%wd\n", leng);
fflush(stdout); /* FIXME */
}
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

if (OHNO) { 
flint_printf("**odd roots REC CALL k=%wd**\n", k);
fflush(stdout); /* FIXME */
}
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
#ifdef WANT_ASSERT
    { /* FIXME very expensive assertion here */
        fmpz *evals = _fmpz_vec_init(nroots);
        _fmpz_mod_poly_evaluate_fmpz_vec(evals, poly, len, roots, nroots, p);
if(OHNO)
{
flint_printf("\nchecking k=%wd p=",k); fmpz_print(p); flint_printf("\n"); 
_fmpz_vec_print(roots, nroots); flint_printf("\n");
_fmpz_vec_print(expons, nroots); flint_printf("\n");
fflush(stdout);
}
        FLINT_ASSERT(_fmpz_vec_is_zero(evals, nroots));
        _fmpz_vec_clear(evals, nroots);
    }
#endif
    return nroots;
}

int fmpz_sparse_bp_interp(fmpz_sparse_t res,
        const fmpz_sparse_bp_interp_t evals)
{
    /* Berlekamp-Massey to find the Prony polynomial
      * Compute roots of Prony poly, along with discrete logs
      * Compute coefficients from transposed Vandermode
      */

    fq_ctx_t Fq;
    fq_mat_t M;
    fq_mat_t b;
    fmpz_mod_poly_t G;
    fmpz_t temp_fmpz;
    fmpz * roots;
    fmpz * expons;
    slong T = evals->length / 2;
    slong t;
    slong i, j;
    slong temp_slong;
    slong * perm;
    const fmpz * w = evals->sample_points + 1;

/* FIXME
flint_printf("\nstart with T=%wd q=",T); fmpz_print(evals->q); flint_printf("\n"); fflush(stdout);
*/

    fmpz_sparse_zero(res);

    if (T == 0) return 1;

    /* Berlekamp-Massey to discover Prony polynomial */

    fmpz_mod_poly_init(G, evals->q);
    fmpz_mod_poly_minpoly(G, evals->evaluations, evals->length);

    t = fmpz_mod_poly_degree(G);
    if (t > T)
    {
        /* sparsity estimate was too low */
        fmpz_mod_poly_clear(G);
        return 0;
    }

    /* find roots of Prony polynomial, and their orders */

    roots = _fmpz_vec_init(t);
    expons = _fmpz_vec_init(t);
{
fmpz_t checkp;
fmpz_init(checkp);
fmpz_set_str(checkp, "0", 10);
if (fmpz_equal(checkp, evals->q)) OHNO=1;
if (OHNO) {
flint_printf("------------binary roots---------------------------------------------\n\n"); fflush(stdout); /* FIXME */
}
fmpz_clear(checkp);
}
    _fmpz_mod_poly_binary_roots(roots, expons, 
            G->coeffs, G->length, w, evals->log2_order, evals->q);

    /* solve transposed Vandermode to get coeffs */
    fmpz_init(temp_fmpz);

    /* build finite field */
    fq_ctx_init(Fq, evals->q, WORD(1), "x");

    fq_mat_init(M, t, t, Fq);
    perm = flint_calloc(t, sizeof(slong));

    for (i=0; i<t; ++i)
    {
        fq_set_ui(fq_mat_entry(M,0,i), UWORD(1), Fq);
    }

    for (i=1; i<t; ++i)
    {
        for (j=0; j<t; ++j) 
        {
            fmpz_mul(temp_fmpz, 
                    fmpz_poly_get_coeff_ptr(fq_mat_entry(M,i-1,j), 0),
                    roots+j);
            fq_set_fmpz(fq_mat_entry(M,i,j), temp_fmpz, Fq);
        }
    }

    temp_slong = fq_mat_lu(perm, M, 1, Fq);
    FLINT_ASSERT (temp_slong);

#ifdef WANT_ASSERT
    for (i=0; i<t; ++i) FLINT_ASSERT(perm[i] == i);
#endif

    fq_mat_init(b, t, 1, Fq);
    for (i=0; i<t; ++i)
    {
        fq_set_fmpz(fq_mat_entry(b,i,0), evals->evaluations+i, Fq);
    }

    fq_mat_solve_tril(b, M, b, 1, Fq);
    fq_mat_solve_triu(b, M, b, 0, Fq);

    /* set coeffs and expons of the actual polynomial */
    
    if (evals->laurent)
    {
        fmpz_set_ui(temp_fmpz, UWORD(1));
        fmpz_mul_2exp(temp_fmpz, temp_fmpz, evals->log2_order-1);

        for (i=0; i<t; ++i)
        {
            if (fmpz_cmp(expons+i, temp_fmpz) > 0)
            {
                fmpz_submul_ui(expons+i, temp_fmpz, UWORD(2));
            }
        }
    }

    for (i=0; i<t; ++i)
    {
        fq_struct * coeff = fq_mat_entry(b,i,0);

        if (!fq_is_zero(coeff, Fq))
        {
            FLINT_ASSERT(fmpz_poly_length(coeff) == UWORD(1));
            fmpz_mods(temp_fmpz, fmpz_poly_get_coeff_ptr(coeff,0), evals->q);
            fmpz_sparse_set_coeff(res, temp_fmpz, expons+i);
        }
    }

    fq_mat_clear(b, Fq);
    fq_mat_clear(M, Fq);

    /* clean-up */
    _fmpz_vec_clear(roots, t);
    _fmpz_vec_clear(expons, t);
    fmpz_mod_poly_clear(G);
    flint_free(perm);
    fq_ctx_clear(Fq);
    fmpz_clear(temp_fmpz);

    return 1;
}
