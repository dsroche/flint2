/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPinterpE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

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

    flint_printf("mul_interp....");
    fflush(stdout);

    /* Compare against mul_heaps */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c1, c2;
        fmpz_t d, e;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest_unsigned(d, state, 20);
        fmpz_randtest_unsigned(e, state, 20);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c1);
        fmpz_spoly_init(c2);
        fmpz_spoly_randtest(a, state, n_randint(state, 30), d, 50);
        fmpz_spoly_randtest(b, state, n_randint(state, 30), e, 50);

        fmpz_spoly_mul_heaps(c1, a, b);
        fmpz_spoly_mul_interp(c2, state, a, b, fmpz_spoly_terms(c1));

        result = fmpz_spoly_equal(c1, c2);
        if (!result)
        {
            flint_printf("FAIL PHASE 1:\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            fmpz_spoly_print(c1), flint_printf("\n\n");
            fmpz_spoly_print(c2), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c1);
        fmpz_spoly_clear(c2);
        fmpz_clear(d);
        fmpz_clear(e);
    }
    
    /* Check aliasing of a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c;
        fmpz_t d, e;
        slong outsize;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest_unsigned(d, state, 20);
        fmpz_randtest_unsigned(e, state, 20);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_randtest(b, state, n_randint(state, 30), d, 50);
        fmpz_spoly_randtest(c, state, n_randint(state, 30), e, 50);
        outsize = fmpz_spoly_terms(b) * fmpz_spoly_terms(c);

        fmpz_spoly_mul_interp(a, state, b, c, outsize);
        fmpz_spoly_mul_interp(b, state, b, c, outsize);

        result = fmpz_spoly_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL PHASE 2:\n");
            fmpz_spoly_print(a); flint_printf("\n\n");
            fmpz_spoly_print(b); flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c;
        fmpz_t d, e;
        slong outsize;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest_unsigned(d, state, 20);
        fmpz_randtest_unsigned(e, state, 20);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_randtest(b, state, n_randint(state, 30), d, 50);
        fmpz_spoly_randtest(c, state, n_randint(state, 30), e, 50);
        outsize = fmpz_spoly_terms(b) * fmpz_spoly_terms(c);

        fmpz_spoly_mul_interp(a, state, b, c, outsize);
        fmpz_spoly_mul_interp(c, state, b, c, outsize);

        result = (fmpz_spoly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL PHASE 3:\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
