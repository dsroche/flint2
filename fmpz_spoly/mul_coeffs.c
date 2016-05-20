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

        Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_spoly.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

/* bit vector operations */
#define FMPZ_SPOLY_SETBIT(vec, ind) \
    ( (vec)[(ind) / FLINT_BITS] |= (UWORD(1) << ((ind) % FLINT_BITS)) )

#define FMPZ_SPOLY_GETBIT(vec, ind) \
    ( (vec)[(ind) / FLINT_BITS] & (UWORD(1) << ((ind) % FLINT_BITS)) )

void _fmpz_spoly_mul_coeffs(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
    /* NOTE: no aliasing between inputs and output */
    ulong* known;
    slong num_remain = res->length;
    slong last_success = 1;
    fmpz_t allp;
    slong* frompos;
    fmpz *poly1_p, *poly2_p, *res_p, *known_p;
    slong vecs_len = 3 * num_remain;
    slong will_get;
    slong i, j; 
    slong p = 2;

    if (num_remain < 10)
    {
        fmpz_spoly_mul_heaps(res, poly1, poly2);
        return;
    }

    /* bit field to store which coeffs are known */
    known = flint_calloc((num_remain - 1) / FLINT_BITS + 1, sizeof(ulong));

    /* stores products of all p's used */
    fmpz_init_set_ui(allp, UWORD(1));

    /* indicates from which original exponent index this one came */
    frompos = flint_malloc(vecs_len * sizeof(slong));

    poly1_p = _fmpz_vec_init(vecs_len);
    poly2_p = _fmpz_vec_init(vecs_len);
    res_p = _fmpz_vec_init(vecs_len * 2);
    known_p = _fmpz_vec_init(vecs_len);

    while (num_remain)
    {
        /* choose prime p for exponents mod */
        if (last_success)
        {
            /* choose new p based on num_remain */
            p = (slong) FLINT_MAX(n_nextprime(((ulong)num_remain) * 2, 0), 31);
        }
        else {
            p = (slong) n_nextprime((ulong) p, 0);
        }

        while (fmpz_divisible_si(allp, p))
        {
            /* skip over primes already tried */
            p = (slong) n_nextprime((ulong) p, 0);
        }

        if (p > vecs_len) 
        {
            _fmpz_vec_clear(poly1_p, vecs_len);
            _fmpz_vec_clear(poly2_p, vecs_len);
            _fmpz_vec_clear(res_p, vecs_len * 2);
            _fmpz_vec_clear(known_p, vecs_len);
            vecs_len = p * 2;
            frompos = (slong*) flint_realloc(frompos, vecs_len * sizeof(slong));
            poly1_p = _fmpz_vec_init(vecs_len);
            poly2_p = _fmpz_vec_init(vecs_len);
            res_p = _fmpz_vec_init(vecs_len * 2);
            known_p = _fmpz_vec_init(vecs_len);
        }

        /* find out where each exponent maps mod p */
        for (i = 0; i < p; ++i) frompos[i] = WORD(-1);

        _fmpz_vec_zero(known_p, p);
        will_get = 0;

        for (i = 0; i < res->length; ++i)
        {
            j = (slong) fmpz_fdiv_ui(res->expons + i, (ulong) p);
            if (FMPZ_SPOLY_GETBIT(known, i))
            {
                /* add to known result */
                fmpz_add(known_p + j, known_p + j, res->coeffs + i);
            }
            else
            {
                /* keep track of index */
                if (frompos[j] < 0)
                {
                    /* nothing there yet; put this there */
                    frompos[j] = i;
                    will_get += 1;
                }
                else if (frompos[j] < res->length)
                {
                    /* collision; set to invalid index res->length */
                    frompos[j] = res->length;
                    will_get -= 1;
                }
            }
        }

        if (will_get <= num_remain / 2) 
        {
            /* not worth it */
            last_success = 0;
            continue;
        }

        /* reduce inputs mod x^p - 1 */
        _fmpz_vec_zero(poly1_p, p);
        for (i = 0; i < poly1->length; ++i)
        {
            j = (slong) fmpz_fdiv_ui(poly1->expons + i, (ulong) p);
            fmpz_add(poly1_p + j, poly1_p + j, poly1->coeffs + i);
        }

        _fmpz_vec_zero(poly2_p, p);
        for (i = 0; i < poly2->length; ++i)
        {
            j = (slong) fmpz_fdiv_ui(poly2->expons + i, (ulong) p);
            fmpz_add(poly2_p + j, poly2_p + j, poly2->coeffs + i);
        }

        /* actually multiply */
        _fmpz_poly_mul(res_p, poly1_p, p, poly2_p, p);
        fmpz_zero(res_p + (p * 2 - 1));

        /* recover coefficients */
        for (i = 0; i < p; ++i)
        {
            j = frompos[i];
            if (j >= 0 && j < res->length)
            {
                fmpz_add(res->coeffs + j, res_p + i, res_p + (i + p));
                fmpz_sub(res->coeffs + j, res->coeffs + j, known_p + i);
                FMPZ_SPOLY_SETBIT(known, j);
                --num_remain;
            }
        }

        /* update product of all primes used */
        fmpz_mul_si(allp, allp, p);

        last_success = 1;
    }

    /* clean up */
    flint_free(known);
    flint_free(frompos);
    fmpz_clear(allp);
    _fmpz_vec_clear(poly1_p, vecs_len);
    _fmpz_vec_clear(poly2_p, vecs_len);
    _fmpz_vec_clear(res_p, vecs_len * 2);
    _fmpz_vec_clear(known_p, vecs_len);
}

#undef FMPZ_SPOLY_SETBIT
#undef FMPZ_SPOLY_GETBIT
