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
    slong i, j, k, length, unique;
    fmpz * rands;
    fmpz_t half, temp;

    unique = 0;

    fmpz_spoly_zero(res);

    /*res is zero if no terms are desired*/
    if(terms == 0)
        return;

    if(fmpz_cmp_si(degree, terms - 1) < 0)
    {
        terms = fmpz_get_ui(degree) + 1;
    }

    _fmpz_spoly_reserve(res, terms);

    fmpz_randbits(res->coeffs, state, bits);
    fmpz_set(res->expons, degree);

    /*res is one term of the desired degree if one terms are desired
      * or degree is zero*/
    if(terms == 1 || fmpz_is_zero(degree))
    {
        _fmpz_spoly_set_length(res, 1);
        return;
    }

    fmpz_init(half);
    fmpz_init(temp);
    fmpz_fdiv_q_ui(half, degree, 2);
    fmpz_set_ui(temp, terms);

    /*if terms > 1/2 * degree choose D + 1 - T unique exponents*/
    if(fmpz_cmp(temp, half) > 0)
    {
        length = fmpz_get_ui(degree) + 1 - terms;
    }
    /*if terms <= 1/2 * degree choose T - 1 unique exponents*/
    else
    {
        length = terms - 1;
    }

    /*if a completely dense polynomial is desired*/
    if(length == 0)
    {
        for(i = 1; i < terms; i++)
        {
            fmpz_randbits(res->coeffs + i, state, bits);
            fmpz_set_ui(res->expons + i, terms - i - 1);
        }
    }
    else
    {
        rands = _fmpz_vec_init(length);

        for(i = 0; i < length; ++i)
        {
            fmpz_randm(rands + i, state, degree);
        }

        /*checks to make sure that all random exponents are unique*/
        while(unique == 0)
        {
            unique = 1;

            /*fmpz_vec sorts in ascending order*/
            _fmpz_vec_sort(rands, length);

            for(i = 0; i < length - 1; ++i)
            {
                if(fmpz_equal(rands + i, rands + i + 1))
                {
                    fmpz_randm(rands + i, state, degree);
                    unique = 0;
                }
            }
        }

        /*if terms > .5 * degree*/
        if(fmpz_cmp(temp, half) > 0)
        {
            k = 1;
            j = length - 1;
            
            /*fmpz_vec rands is sorted in ascending order*/
            for(i = fmpz_get_si(degree) - 1; i >= 0; --i)
            {
                if(j >= 0 && fmpz_equal_si(rands + j, i))
                    j--;
                else     
                {
                    fmpz_set_si(res->expons + k, i);
                    if(fmpz_sgn(res->expons + k) == -1)
                        fmpz_mul_si(res->expons + k, res->expons + k, -1);
                    k++;
                }
            }
        }
        /*if terms <= .5 * degree*/
        else
        {
            for(i = 1; i < terms; ++i)
            {
                fmpz_set(res->expons + i, rands + i - 1);
            }
        }

        for (i = 1; i < terms; ++i)
        {
            fmpz_randbits(res->coeffs + i, state, bits);
        }

        _fmpz_vec_clear(rands, length);
    }
    
    _fmpz_spoly_set_length(res, terms);
    _fmpz_spoly_normalise(res);
    
    fmpz_clear(half);
    fmpz_clear(temp);
}

void fmpz_spoly_randtest_kron(fmpz_spoly_t res, flint_rand_t state, slong terms, 
    const fmpz_t degree, const fmpz_t limit, mp_bitcnt_t bits, slong vars)
{
    /*slong i, j, k, length, unique;
    fmpz * rands;
    fmpz_t abs, half, temp;

    unique = 0;

    fmpz_spoly_zero(res);

    if(terms == 0)
        return;

    if(fmpz_cmp_si(degree, terms - 1) < 0)
    {
        terms = fmpz_get_ui(degree) + 1;
    }

    _fmpz_spoly_reserve(res, terms);

    */

    return;
}
