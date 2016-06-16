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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_mod....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c, n;
        fmpz_init(n);
        fmpz_randtest_unsigned_not_zero(n, state, 200);

        fmpz_init(c);
        fmpz_randtest_unsigned(c, state, 200);
        
        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), c, 200);

        fmpz_spoly_scalar_mod(b, a, n);
        fmpz_spoly_scalar_mod(a, a, n);

        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL (1):\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(n);
        fmpz_clear(c);
        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
    }

    /* Compare with fmpz_spoly_scalar_mul_si */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c, n1;
        slong n;
        fmpz_init(n1);
        n = (slong) n_randbits(state, FLINT_BITS - 1);
        if (n_randint(state, 2))
            n = -n;
        fmpz_set_si(n1, n);

        fmpz_init(c);
        fmpz_randtest_unsigned(c, state, 100);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), c, 200);

        fmpz_spoly_scalar_mod(b, a, n1);
        _fmpz_vec_scalar_mod_fmpz(a->coeffs, a->coeffs, a->length, n1);
        _fmpz_spoly_normalise(a);

        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL (2):\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(n1);
        fmpz_clear(c);
        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
