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
        slong terms, const fmpz_t degree, slong limit, mp_bitcnt_t bits,
        slong vars)
{
    slong i, j;
    fmpz_t temp;

    fmpz_init(temp);

    fmpz_spoly_zero(res);

    /*res is zero if no terms are desired*/
    if(terms == 0)
    {
        return;
    }

    _fmpz_spoly_reserve(res, terms);

    /*The height is guaranteed by the first term and then bitsize becomes random*/
    fmpz_randbits(res->coeffs, state, bits);
   
    /*Make the first term and include the degree*/
    fmpz_set(res->expons, degree);
    for(i = 1; i < vars; i++)
    {
        fmpz_mul_2exp(res->expons, res->expons, limit);
        fmpz_randm(temp, state, degree);
        fmpz_add(res->expons, res->expons, temp);
    }

    if(terms == 1 || fmpz_is_zero(degree))
    {
        _fmpz_spoly_set_length(res, 1);
        return;
    }

    for(i = 1; i < terms; i++)
    {
        fmpz_randbits(res->coeffs + i, state, bits);
        
        fmpz_randm(res->expons + i, state, degree);
        for(j = 0; j < vars; j++)
        {
            fmpz_mul_2exp(res->expons + i, res->expons + i, limit);
            fmpz_randm(temp, state, degree);
            fmpz_add(res->expons + i, res->expons + i, temp);
        }
    }

    _fmpz_spoly_set_length(res, terms);
    _fmpz_spoly_normalise(res);

    fmpz_clear(temp);
}
