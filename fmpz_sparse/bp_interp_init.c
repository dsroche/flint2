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

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_sparse.h"

void fmpz_sparse_bp_interp_init(fmpz_sparse_bp_interp_t res,
        slong terms, const fmpz_t height, const fmpz_t degree, flint_rand_t rs)
{
    fmpz_t k; /* we will have q = order*k + 1 */
    fmpz_t w;
    fmpz_t test;
    fmpz_t order;
    slong i;

    FLINT_ASSERT(terms >= 0);

    if (terms == 0) {
        fmpz_init_set_ui(res->q, UWORD(2));
        res->log2_order = UWORD(0);
        res->sample_points = NULL;
        res->evaluations = NULL;
        res->length = WORD(0);
        return;
    }

    fmpz_init_set_ui(res->q, UWORD(1));
    fmpz_init_set_ui(order, UWORD(1));
    fmpz_init(k);

    /* order = next higher power of 2 from degree */
    res->log2_order = FLINT_MAX(1, fmpz_bits(degree));

    if (fmpz_cmp_si(degree, WORD(0)) >= 0)
    {
        res->laurent = 0;
    }
    else
    {
        /* negative degree means Laurent polynomial */
        res->log2_order += 1;
        res->laurent = 1;
    }

    fmpz_mul_2exp(order, order, res->log2_order);

    /* Starting point for q = ceil(2*height/order)*order + 1 */
    fmpz_mul_2exp(k, height, UWORD(1));
    fmpz_cdiv_q(k, k, order);
    fmpz_addmul(res->q, k, order);

    /* Search for next prime in the arithmetic progression. */
    while (! fmpz_is_probabprime(res->q)) 
    {
        fmpz_add(res->q, res->q, order);
        fmpz_add_ui(k, k, UWORD(1));
    }

    /* Now find a PRU with the given order modulo q. */
    fmpz_init2(w, fmpz_size(res->q));
    fmpz_init2(test, fmpz_size(res->q));
    fmpz_fdiv_q_2exp(order, order, UWORD(1)); /* actual order / 2 */
    do {
        fmpz_randm(w, rs, res->q);
        fmpz_powm(w, w, k, res->q);
        fmpz_powm(test, w, order, res->q);
        fmpz_add_ui(test, test, UWORD(1));
    } while (!fmpz_equal(test, res->q));

    fmpz_clear(k);
    fmpz_clear(test);
    fmpz_clear(order);

    /* Need 2*terms evaluation points. */
    res->length = 2*terms;
    res->sample_points = _fmpz_vec_init(res->length);
    res->evaluations = _fmpz_vec_init(res->length);

    /* Sample point i = w^i mod q */
    fmpz_set_ui (res->sample_points+0, UWORD(1));
    for (i = 1; i < res->length; ++i) {
        fmpz_mul(res->sample_points+i, res->sample_points+(i-1), w);
        fmpz_mod(res->sample_points+i, res->sample_points+i, res->q);
    }

    fmpz_clear(w);
}
