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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_spoly.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("sparse_mul_si....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c;
        slong n = z_randtest(state);

        fmpz_init(c);
        fmpz_randtest_unsigned(c, state, 100);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), c, 200);

        fmpz_spoly_scalar_mul_si(b, a, n);
        fmpz_spoly_scalar_mul_si(a, a, n);

        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL MUL:\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_clear(c);
    }

    /* Compare with fmpz_spoly_sparse_mul_ui */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c;
        ulong n = n_randbits(state, FLINT_BITS - 1);

        fmpz_init(c);
        fmpz_randtest_unsigned(c, state, 100);
        
        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), c, 200);

        fmpz_spoly_scalar_mul_ui(b, a, n);
        fmpz_spoly_scalar_mul_si(a, a, n);

        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL COMPARE WITH UI:\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_clear(c);
    }

    /* Check (a*n1)*n2 = a*(n1*n2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c;
        fmpz_t d;
        slong n1 = (slong) n_randbits(state, (FLINT_BITS - 2) / 2);
        slong n2 = (slong) n_randbits(state, (FLINT_BITS - 2) / 2);
        if (n_randint(state, 2))
            n1 = -n1;
        if (n_randint(state, 2))
            n2 = -n2;

        fmpz_init(d);
        fmpz_randtest_unsigned(d, state, 100);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), d, 200);

        fmpz_spoly_scalar_mul_si(b, a, n1);
        fmpz_spoly_scalar_mul_si(c, b, n2);
        fmpz_spoly_scalar_mul_si(b, a, n1 * n2);

        result = (fmpz_spoly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL n1 = %wd, n2 = %wd:\n", n1, n2);
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            fmpz_spoly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
        fmpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
