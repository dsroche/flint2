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

    Authored 2016 by A. Whitman Groves; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_spoly.h"
#include "fmpz_vec.h"

void fmpz_spoly_randtest_kron(fmpz_spoly_t res, flint_rand_t state, 
        slong terms, const fmpz_t degree, mp_bitcnt_t bits, ulong limit,
        ulong vars)
{
    slong i, j;
    fmpz_t temp;

    if (fmpz_cmp_si(degree, WORD(0)) < 0 || fmpz_bits(degree) > limit)
    {
        flint_printf("Exception (fmpz_spoly_randtest_kron). Degree negative or too large.\n");
        abort();
    }

    if (vars < 1)
    {
        flint_printf("Exception (fmpz_spoly_randtest_kron). Need at least one variable.\n");
        abort();
    }

    fmpz_spoly_zero(res);

    /*res is zero if no terms are desired*/
    if(terms == 0)
    {
        return;
    }

    _fmpz_spoly_reserve(res, terms);

    fmpz_init(temp);

    /*The height is guaranteed by the first term and then bitsize becomes random*/
    fmpz_randbits(res->coeffs, state, bits);
   
    /*Make the first term and include the degree*/
    fmpz_set(res->expons, degree);
    for(i = 1; (ulong)i < vars; i++)
    {
        fmpz_mul_2exp(res->expons, res->expons, limit);
        fmpz_add(res->expons, res->expons, degree);
    }

    if(terms == 1 || fmpz_is_zero(degree))
    {
        _fmpz_spoly_set_length(res, 1);
        return;
    }

    i = 1;

    fmpz_add_ui(temp, degree, 1);
    fmpz_pow_ui(temp, temp, vars);
    if (fmpz_cmp_ui(temp, ((ulong) terms) * 2) < 0)
    {
        /* dense polynomial */
        fmpz_spoly_t tpoly;
        ulong maxterms = fmpz_get_ui(temp);
        slong deg;
        slong* multideg;
        
        FLINT_ASSERT(fmpz_fits_si(degree));
        deg = fmpz_get_si(degree);

        multideg = flint_malloc(vars * sizeof *multideg);
        for (j = 0; (ulong)j < vars; ++j) multideg[j] = deg;

        fmpz_spoly_init(tpoly);
        fmpz_spoly_randtest_kron(tpoly, state, FLINT_MAX(0, (slong)(maxterms - terms)) + 1, degree, 1, limit, vars);

        while (i < terms)
        {
            fmpz_set_si(res->expons + i, multideg[0]);
            for (j = 1; (ulong)j < vars; ++j)
            {
                fmpz_mul_2exp(res->expons + i, res->expons + i, limit);
                fmpz_add_ui(res->expons + i, res->expons + i, (ulong)multideg[j]);
            }

            if (_fmpz_spoly_index(tpoly, res->expons + i) < 0)
            {
                fmpz_randbits(res->coeffs + i, state, bits);
                ++i;
            }

            for (j = (slong)vars - 1; j >= 0; --j)
            {
                if (multideg[j] > 0) break;
            }

            if (j < 0) break;

            multideg[j]--;
            for (++j; (ulong)j < vars; ++j) multideg[j] = deg;
        }

        FLINT_ASSERT((ulong)i == FLINT_MIN((ulong)terms, maxterms));
        
        fmpz_spoly_clear(tpoly);
        flint_free(multideg);
    }
    else
    {
        while (i < terms)
        {
            fmpz_randbits(res->coeffs + i, state, bits);
            
            fmpz_randm(res->expons + i, state, degree);
            for(j = 1; (ulong)j < vars; j++)
            {
                fmpz_mul_2exp(res->expons + i, res->expons + i, limit);
                fmpz_randm(temp, state, degree);
                fmpz_add(res->expons + i, res->expons + i, temp);
            }

            for (j = 0; j < i; ++j)
            {
                if (fmpz_equal(res->expons + j, res->expons + i)) break;
            }

            if (j == i) ++i;
        }
    }

    _fmpz_spoly_set_length(res, terms);
    _fmpz_spoly_normalise(res);

    fmpz_clear(temp);
}
