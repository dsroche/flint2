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

    flint_printf("scalar_addmul....");
    fflush(stdout);

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c, n;
        fmpz_init(n);
        fmpz_randtest(n, state, 200);

        fmpz_init(c);
        fmpz_randtest(c, state, 200);
        
        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), c, 200);
        fmpz_spoly_set(b, a);

        fmpz_spoly_scalar_addmul(b, a, n);
        fmpz_spoly_scalar_addmul(a, a, n);

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

    /* Check that b += x*a equals c = b + x*a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c;
        fmpz_t d, e, x;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_init(x);
        fmpz_randtest(d, state, 100);
        fmpz_randtest(e, state, 100);

        fmpz_randtest(x, state, n_randint(state, 100));

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), d, 200);
        fmpz_spoly_randtest(b, state, n_randint(state, 100), e, 200);

        fmpz_spoly_scalar_mul(c, a, x);
        fmpz_spoly_add(c, b, c);

        fmpz_spoly_scalar_addmul(b, a, x);

        result = (fmpz_spoly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (2):\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            fmpz_spoly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(e);
        fmpz_clear(x);
        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
