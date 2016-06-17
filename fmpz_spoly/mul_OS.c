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

    Authored 2015 by A. Whitman Groves; US Government work in the public domain. 
    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_spoly.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

void
fmpz_spoly_mul_OS(fmpz_spoly_t res, flint_rand_t state, const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
    fmpz * test;
    slong length;
    ulong p, pinv, x, eval_actual, eval_check;

    if (fmpz_spoly_is_zero(poly1) || fmpz_spoly_is_zero(poly2))
    {
        fmpz_spoly_zero(res);
        return;
    }

    if (poly1->length == 1 && poly2->length == 1)
    {
        fmpz_add(res->expons, poly1->expons + 0, poly2->expons + 0);
        fmpz_mul(res->coeffs, poly1->coeffs + 0, poly2->coeffs + 0);
        res->length = 1;
        return;
    }

    if (res == poly1 || res == poly2)
    {
        /* No aliasing between output and inputs */
        fmpz_spoly_t temp;
        fmpz_spoly_init(temp);
        fmpz_spoly_mul_OS(temp, state, poly1, poly2);
        fmpz_spoly_swap(res, temp);
        fmpz_spoly_clear(temp);
        return;
    }

    do
    {
        test = NULL;
        length = fmpz_spoly_sumset(&test, state, poly1, poly2);
        FLINT_ASSERT(length > 1);

        _fmpz_spoly_reserve(res, length);
        _fmpz_poly_reverse(res->expons, test, length, length);
        _fmpz_spoly_set_length(res, length);
        _fmpz_vec_clear(test, length);
        _fmpz_spoly_mul_coeffs(res, poly1, poly2);
        
        /* check the result mod a random prime */
        p = n_randprime(state, FLINT_BITS, 0);
        pinv = n_preinvert_limb(p);
        x = n_randint(state, p);
        eval_actual = n_mulmod2_preinv(
                fmpz_spoly_evaluate_mod_ui(poly1, x, p),
                fmpz_spoly_evaluate_mod_ui(poly2, x, p),
                p, pinv);
        eval_check = fmpz_spoly_evaluate_mod_ui(res, x, p);

        if (eval_actual != eval_check) {
            flint_printf("DETECTED mul_OS failure %wu %wu %wu\n", poly1->length, poly2->length, length);
        }
    }
    while (eval_actual != eval_check);
}
