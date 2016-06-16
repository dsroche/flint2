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

    Authored 2015 by Daniel S. Roche and A. Whitman Groves;
    US Government work in the public domain. 

******************************************************************************/

#include "fmpz_spoly.h"
#include "fmpz_vec.h"

void fmpz_spoly_randtest(fmpz_spoly_t res, flint_rand_t state, 
        slong terms, const fmpz_t degree, mp_bitcnt_t bits)
{
    slong i, j;
    fmpz_t temp;

    if (fmpz_cmp_si(degree, WORD(0)) < 0)
    {
        flint_printf("Exception (fmpz_spoly_randtest). Degree negative\n");
        abort();
    }

    fmpz_spoly_zero(res);

    /*res is zero if no terms are desired*/
    if(terms == 0)
    {
        return;
    }

    if(fmpz_cmp_si(degree, terms - 1) < 0)
        terms = fmpz_get_ui(degree) + 1;

    _fmpz_spoly_reserve(res,terms);

    _fmpz_spoly_set_length(res, terms);
    
    /*Make the first term and include the degree*/
    fmpz_randbits(fmpz_spoly_get_term_coeff_ptr(res, 0), state, bits);
    fmpz_set(fmpz_spoly_get_term_expon_ptr(res, 0), degree);

    /*res is one term of the desired degree if one terms are desired
      * or degree is zero*/
    if(terms == 1 || fmpz_is_zero(degree))
    {
        _fmpz_spoly_set_length(res, 1);
        return;
    }


    i = 1;

    if(fmpz_cmp_ui(degree, ((ulong) terms) * 2) < 0)
    {
        fmpz_spoly_t tpoly;
        slong deg;

        FLINT_ASSERT(fmpz_fits_si(degree));
        deg = fmpz_get_si(degree);

        fmpz_spoly_init(tpoly);
        fmpz_spoly_randtest(tpoly, state, FLINT_MAX(0, (slong)(deg - terms)), degree, 1);

        while(i < terms)
        {
            fmpz_set_si(fmpz_spoly_get_term_expon_ptr(res, i), deg);

            if(_fmpz_spoly_index(tpoly, fmpz_spoly_get_term_expon_ptr(res, i)) < 0)
            {
                fmpz_randbits(fmpz_spoly_get_term_coeff_ptr(res, i), state, bits);
                i++;
            }

            deg--;
        }

        fmpz_spoly_clear(tpoly);
    }
    else
    {
        while(i < terms)
        {
            fmpz_randbits(fmpz_spoly_get_term_coeff_ptr(res, i), state, bits);

            fmpz_randm(fmpz_spoly_get_term_expon_ptr(res, i), state, degree);

            for(j = 0; j < i; j++)
            {
                if(fmpz_equal(fmpz_spoly_get_term_expon_ptr(res,j), fmpz_spoly_get_term_expon_ptr(res,i))) break;
            }

            if(i == j) i++;
        }
    }

    _fmpz_spoly_normalise(res);

    fmpz_clear(temp);
}
