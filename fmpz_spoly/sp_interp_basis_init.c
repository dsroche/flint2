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

#include "flint.h"
#include "nmod_vec.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_spoly.h"

const slong COEFF_PBITS = FLINT_BITS - 1;

/* table lookup for number of primes with given bitlength */
const slong _n_prime_bits_table[35] = {0, 0, 2, 2, 2, 5, 7, 13,
    23, 43, 75, 137, 255, 464, 872, 1612,
    3030, 5709, 10749, 20390, 38635, 73586, 140336, 268216,
    513708, 985818, 1894120, 3645744, 7027290, 13561907, 26207278, 50697537,
    98182656, 190335585, 369323305};
const slong _n_prime_bits_tlen = 35;

/* computes lower bound for number of primes with given bitlength */
static ulong _n_prime_bits_lb(slong bits)
{
    if (bits < 0) 
    {
        return 0;
    }
    else if (bits < _n_prime_bits_tlen) 
    {
        return _n_prime_bits_table[bits];
    }
    else
    {
        ulong lo1, hi1, lo2, hi2;
        ulong n = (UWORD(1) << (bits - 1)) - 1;

        n_prime_pi_bounds(&lo1, &hi1, n);

        n = n * 2 + 1;
        n_prime_pi_bounds(&lo2, &hi2, n);

        return lo2 - hi1;
    }
}

void fmpz_spoly_sp_interp_basis_init(fmpz_spoly_sp_interp_basis_t res, flint_rand_t state,
        slong terms, mp_bitcnt_t d, mp_bitcnt_t h)
{
    slong num_rounds, groups_per, coeffs_per, pbits;
    slong i, j, round_start, groups_remain = 0, coeffs_remain = 0;
    fmpz_t cmprod, emprod;

    if (terms == 0 || h == 0)
    {
        /* edge case: zero polynomial */
        res->length = 0;
        res->cimg_per_round = res->eimg_per_round = 0;
        res->cimg_needed = res->eimg_needed = 0;
        return;
    }
    else if (d <= FLINT_BIT_COUNT(terms) + 2)
    {
        /* edge case: dense polynomial */
        num_rounds = 1;
        pbits = FLINT_MAX(2, d + 1);
        groups_per = 1;
        coeffs_per = 1 + h / (COEFF_PBITS - 1);
    }
    else
    {
        /* default case: sparse polynomial */
        num_rounds = FLINT_BIT_COUNT(FLINT_BIT_COUNT(terms) + 11);
        pbits = FLINT_BIT_COUNT(terms) + 2;
        groups_per = 1 + (d * 2 - 1) / (pbits - 1);
        
        /* ensure sufficiently many primes available */
        while (_n_prime_bits_lb(pbits) < UWORD(2) * num_rounds * groups_per)
        {
            ++pbits;
            groups_per = 1 + (d * 2 - 1) / (pbits - 1);
        }

        coeffs_per = 1 + (h * 2 + 1) / (COEFF_PBITS - 1);
    }

    res->length = num_rounds * (groups_per + coeffs_per);
    res->cimg_per_round = terms * coeffs_per;
    res->cimg_needed = 1 + h / (COEFF_PBITS - 1);
    res->eimg_per_round = terms * groups_per;
    res->eimg_needed = FLINT_MAX(1, 1 + (d - 1) / (pbits - 1));

    res->cmods = flint_malloc(res->length * sizeof *res->cmods);
    res->shifts = flint_malloc(res->length * sizeof *res->shifts);
    res->emods = flint_malloc(res->length * sizeof *res->emods);

    fmpz_init(cmprod);
    fmpz_init(emprod);

    i = 0;
    while (i < res->length)
    {
        if (groups_remain == 0)
        {
            /* start new round */
            FLINT_ASSERT(coeffs_remain == 0);
            --num_rounds;
            groups_remain = groups_per;
            coeffs_remain = coeffs_per;
            round_start = i;

            fmpz_set_si(cmprod, WORD(1));
            fmpz_set_si(emprod, WORD(1));

            nmod_init(res->cmods + i, n_randprime(state, COEFF_PBITS, 0));
            do
            {
                res->shifts[i] = n_randint(state, res->cmods[i].n);
            }
            while (res->shifts[i] <= 1);
        }
        else
        {
            res->cmods[i] = res->cmods[round_start];
            res->shifts[i] = res->shifts[round_start];
        }

        /* start new group */
        do
        {
            nmod_init(res->emods + i, n_randprime(state, pbits, 0));
        } 
        while (fmpz_fdiv_ui(emprod, res->emods[i].n) == 0);
        fmpz_mul_ui(emprod, emprod, res->emods[i].n);

        for (j = 1; j <= (coeffs_remain + groups_remain - 1) / groups_remain; ++j)
        {
            do
            {
                nmod_init(res->cmods + (i + j), n_randprime(state, COEFF_PBITS, 0));
            }
            while (fmpz_fdiv_ui(cmprod, res->cmods[i + j].n) == 0);
            fmpz_mul_ui(cmprod, cmprod, res->cmods[i + j].n);

            res->shifts[i + j] = 1;
            res->emods[i + j] = res->emods[i];
        }

        --groups_remain;
        coeffs_remain -= (j - 1);
        i += j;
    }

    FLINT_ASSERT(num_rounds == 0);
    FLINT_ASSERT(groups_remain == 0);
    FLINT_ASSERT(coeffs_remain == 0);
    FLINT_ASSERT(i == res->length);

    fmpz_clear(emprod);
    fmpz_clear(cmprod);
}
