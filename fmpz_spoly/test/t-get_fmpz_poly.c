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
#include "fmpz_poly.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("get_fmpz_poly....");
    fflush(stdout);

    

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a;
        fmpz_poly_t b, c;
        fmpz_t n;

        fmpz_spoly_init(a);
        
        fmpz_poly_init(b);
        fmpz_poly_init(c);


        fmpz_init(n);

        fmpz_randtest_unsigned(n, state, 15);

        fmpz_spoly_randtest(a, state, n_randint(state, 50), n, 200);
        fmpz_spoly_get_fmpz_poly(b, a);
        fmpz_poly_set(c, b);

        result = (fmpz_poly_equal(b, c));

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_spoly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_clear(n);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        fmpz_spoly_t m, n;
        fmpz_t c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        fmpz_spoly_init(m);
        fmpz_spoly_init(n);

        fmpz_init(c);

        fmpz_randtest_unsigned(c, state, 15);

        fmpz_spoly_randtest(m, state, n_randint(state, 50), c, 200);
        fmpz_spoly_randtest(n, state, n_randint(state, 50), c, 200);
        while (fmpz_spoly_equal(m, n))
        {
            if(fmpz_is_zero(c))
              fmpz_randtest_unsigned(c, state, 15);
            fmpz_spoly_randtest(n, state, n_randint(state, 50), c, 200);
        }
        fmpz_spoly_get_fmpz_poly(a, m);
        fmpz_spoly_get_fmpz_poly(b, n);

        result = (!fmpz_poly_equal(a, b));
        if (!result && (fmpz_sgn(m->expons) > 0 && fmpz_sgn(n->expons) > 0))
        {
            flint_printf("FAIL:\n");
            flint_printf("m = "), fmpz_spoly_print(m), flint_printf("\n\n");
            flint_printf("n = "), fmpz_spoly_print(n), flint_printf("\n\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(m);
        fmpz_spoly_clear(n);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
